#include "platform/backend/spcc/VDBSDFField.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <openvdb/io/File.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/MeshToVolume.h>

namespace platform {
namespace backend {
namespace spcc {

namespace {

bool IsFiniteScalar(double v) {
    return std::isfinite(v);
}

bool IsFiniteVec(const chrono::ChVector3d& v) {
    return IsFiniteScalar(v.x()) && IsFiniteScalar(v.y()) && IsFiniteScalar(v.z());
}

chrono::ChVector3d SafeNormalize(const chrono::ChVector3d& v,
                                 const chrono::ChVector3d& fallback) {
    const double n = v.Length();
    if (!std::isfinite(n) || n <= 1.0e-12) {
        return fallback;
    }
    return v * (1.0 / n);
}

float ToFiniteFloat(double value, const char* error_context) {
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string(error_context) + ": non-finite scalar");
    }
    if (value < -static_cast<double>(std::numeric_limits<float>::max()) ||
        value > static_cast<double>(std::numeric_limits<float>::max())) {
        throw std::runtime_error(std::string(error_context) + ": scalar out of float range");
    }
    return static_cast<float>(value);
}

void EnsureOpenVDBInitialized() {
    static const bool initialized = []() {
        openvdb::initialize();
        return true;
    }();
    (void)initialized;
}

using FloatGridSampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>;

double SamplePhiWorld(const openvdb::FloatGrid& grid, const chrono::ChVector3d& x_M) {
    const FloatGridSampler sampler(grid);
    return static_cast<double>(sampler.wsSample(openvdb::Vec3d(x_M.x(), x_M.y(), x_M.z())));
}

chrono::ChVector3d GetVoxelStep(const openvdb::FloatGrid& grid) {
    const auto voxel = grid.voxelSize();
    return chrono::ChVector3d(std::max(1.0e-6, static_cast<double>(voxel[0])),
                              std::max(1.0e-6, static_cast<double>(voxel[1])),
                              std::max(1.0e-6, static_cast<double>(voxel[2])));
}

}  // namespace

struct VDBSDFField::Impl {
    openvdb::FloatGrid::Ptr open_grid;
    chrono::ChVector3d bounding_center_M{0.0, 0.0, 0.0};
    double bounding_radius = 0.0;
    bool has_bounding_sphere = false;
    bool ready = false;
    std::string last_error;
};

VDBSDFField::VDBSDFField() : impl_(new Impl()) {}

VDBSDFField::~VDBSDFField() = default;

bool VDBSDFField::BuildFromTriangleMesh(const chrono::ChTriangleMeshConnected& mesh_M) {
    return BuildFromTriangleMesh(mesh_M, BuildOptions{});
}

bool VDBSDFField::BuildFromTriangleMesh(const chrono::ChTriangleMeshConnected& mesh_M,
                                        const BuildOptions& options) {
    EnsureOpenVDBInitialized();

    impl_->last_error.clear();
    impl_->open_grid.reset();
    impl_->bounding_center_M = chrono::ChVector3d(0.0, 0.0, 0.0);
    impl_->bounding_radius = 0.0;
    impl_->has_bounding_sphere = false;
    impl_->ready = false;

    if (!IsFiniteScalar(options.voxel_size) || options.voxel_size <= 0.0) {
        impl_->last_error = "VDBSDFField::BuildFromTriangleMesh invalid voxel_size";
        return false;
    }
    if (!IsFiniteScalar(options.half_band_width_voxels) || options.half_band_width_voxels <= 0.0) {
        impl_->last_error = "VDBSDFField::BuildFromTriangleMesh invalid half_band_width_voxels";
        return false;
    }

    const auto& vertices = mesh_M.GetCoordsVertices();
    const auto& triangles = mesh_M.GetIndicesVertexes();
    if (vertices.empty() || triangles.empty()) {
        impl_->last_error = "VDBSDFField::BuildFromTriangleMesh mesh has no vertices or triangles";
        return false;
    }

    try {
        std::vector<openvdb::Vec3s> points;
        points.reserve(vertices.size());
        chrono::ChVector3d bbox_min(std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::infinity());
        chrono::ChVector3d bbox_max(-std::numeric_limits<double>::infinity(),
                                    -std::numeric_limits<double>::infinity(),
                                    -std::numeric_limits<double>::infinity());
        for (const auto& v : vertices) {
            if (!IsFiniteVec(v)) {
                impl_->last_error = "VDBSDFField::BuildFromTriangleMesh mesh has non-finite vertex";
                return false;
            }
            bbox_min.x() = std::min(bbox_min.x(), v.x());
            bbox_min.y() = std::min(bbox_min.y(), v.y());
            bbox_min.z() = std::min(bbox_min.z(), v.z());
            bbox_max.x() = std::max(bbox_max.x(), v.x());
            bbox_max.y() = std::max(bbox_max.y(), v.y());
            bbox_max.z() = std::max(bbox_max.z(), v.z());
            points.emplace_back(ToFiniteFloat(v.x(), "VDBSDFField::BuildFromTriangleMesh vertex x"),
                                ToFiniteFloat(v.y(), "VDBSDFField::BuildFromTriangleMesh vertex y"),
                                ToFiniteFloat(v.z(), "VDBSDFField::BuildFromTriangleMesh vertex z"));
        }

        impl_->bounding_center_M = 0.5 * (bbox_min + bbox_max);
        impl_->bounding_radius = 0.0;
        for (const auto& v : vertices) {
            impl_->bounding_radius =
                std::max(impl_->bounding_radius, (v - impl_->bounding_center_M).Length());
        }
        impl_->has_bounding_sphere = std::isfinite(impl_->bounding_radius);

        std::vector<openvdb::Vec3I> tri_faces;
        tri_faces.reserve(triangles.size());
        for (const auto& tri : triangles) {
            const int i0 = tri.x();
            const int i1 = tri.y();
            const int i2 = tri.z();
            if (i0 < 0 || i1 < 0 || i2 < 0) {
                impl_->last_error = "VDBSDFField::BuildFromTriangleMesh negative triangle index";
                return false;
            }
            if (static_cast<std::size_t>(i0) >= vertices.size() ||
                static_cast<std::size_t>(i1) >= vertices.size() ||
                static_cast<std::size_t>(i2) >= vertices.size()) {
                impl_->last_error = "VDBSDFField::BuildFromTriangleMesh triangle index out of range";
                return false;
            }
            tri_faces.emplace_back(i0, i1, i2);
        }

        const auto transform = openvdb::math::Transform::createLinearTransform(options.voxel_size);
        std::vector<openvdb::Vec4I> quad_faces;
        const float half_band = ToFiniteFloat(options.half_band_width_voxels,
                                              "VDBSDFField::BuildFromTriangleMesh half_band_width_voxels");

        impl_->open_grid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(
            *transform, points, tri_faces, quad_faces, half_band, half_band);

        if (!impl_->open_grid) {
            impl_->last_error = "VDBSDFField::BuildFromTriangleMesh failed to create OpenVDB float grid";
            return false;
        }

        if (!options.grid_name.empty()) {
            impl_->open_grid->setName(options.grid_name);
        }

        return FinalizeOpenGrid();
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] OpenVDB Build Exception: " << e.what() << std::endl;
        impl_->last_error = std::string("VDBSDFField::BuildFromTriangleMesh exception: ") + e.what();
        impl_->ready = false;
        impl_->open_grid.reset();
        return false;
    }
}

bool VDBSDFField::LoadFromVDBFile(const std::string& vdb_path,
                                  const std::string& grid_name) {
    EnsureOpenVDBInitialized();

    impl_->last_error.clear();
    impl_->ready = false;
    impl_->open_grid.reset();
    impl_->bounding_center_M = chrono::ChVector3d(0.0, 0.0, 0.0);
    impl_->bounding_radius = 0.0;
    impl_->has_bounding_sphere = false;

    if (vdb_path.empty()) {
        impl_->last_error = "VDBSDFField::LoadFromVDBFile empty vdb_path";
        return false;
    }

    try {
        openvdb::io::File file(vdb_path);
        file.open();

        std::string selected_grid_name = grid_name;
        if (selected_grid_name.empty()) {
            const auto metadata = file.readAllGridMetadata();
            if (!metadata || metadata->empty()) {
                file.close();
                impl_->last_error = "VDBSDFField::LoadFromVDBFile no grids in file";
                return false;
            }
            selected_grid_name = metadata->front()->getName();
        }

        openvdb::GridBase::Ptr base_grid = file.readGrid(selected_grid_name);
        file.close();

        if (!base_grid) {
            impl_->last_error = "VDBSDFField::LoadFromVDBFile grid not found: " + selected_grid_name;
            return false;
        }

        auto float_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
        if (!float_grid) {
            impl_->last_error = "VDBSDFField::LoadFromVDBFile grid is not FloatGrid: " + selected_grid_name;
            return false;
        }

        impl_->open_grid = float_grid;
        return FinalizeOpenGrid();
    } catch (const std::exception& e) {
        impl_->last_error = std::string("VDBSDFField::LoadFromVDBFile exception: ") + e.what();
        return false;
    }
}

bool VDBSDFField::IsReady() const {
    return impl_->ready;
}

const std::string& VDBSDFField::LastError() const {
    return impl_->last_error;
}

bool VDBSDFField::QueryPhiGradM(const chrono::ChVector3d& x_M,
                                double& phi,
                                chrono::ChVector3d& grad_M) const {
    if (!impl_->ready || !impl_->open_grid || !IsFiniteVec(x_M)) {
        return false;
    }

    phi = SamplePhiWorld(*impl_->open_grid, x_M);
    const chrono::ChVector3d h = GetVoxelStep(*impl_->open_grid);

    const double phi_xp = SamplePhiWorld(*impl_->open_grid, x_M + chrono::ChVector3d(h.x(), 0.0, 0.0));
    const double phi_xm = SamplePhiWorld(*impl_->open_grid, x_M - chrono::ChVector3d(h.x(), 0.0, 0.0));
    const double phi_yp = SamplePhiWorld(*impl_->open_grid, x_M + chrono::ChVector3d(0.0, h.y(), 0.0));
    const double phi_ym = SamplePhiWorld(*impl_->open_grid, x_M - chrono::ChVector3d(0.0, h.y(), 0.0));
    const double phi_zp = SamplePhiWorld(*impl_->open_grid, x_M + chrono::ChVector3d(0.0, 0.0, h.z()));
    const double phi_zm = SamplePhiWorld(*impl_->open_grid, x_M - chrono::ChVector3d(0.0, 0.0, h.z()));

    const chrono::ChVector3d grad_raw((phi_xp - phi_xm) / (2.0 * h.x()),
                                      (phi_yp - phi_ym) / (2.0 * h.y()),
                                      (phi_zp - phi_zm) / (2.0 * h.z()));
    grad_M = SafeNormalize(grad_raw, chrono::ChVector3d(1.0, 0.0, 0.0));

    return IsFiniteScalar(phi) && IsFiniteVec(grad_M);
}

bool VDBSDFField::QueryPhiM(const chrono::ChVector3d& x_M, double& phi) const {
    if (!impl_->ready || !impl_->open_grid || !IsFiniteVec(x_M)) {
        return false;
    }

    phi = SamplePhiWorld(*impl_->open_grid, x_M);
    return IsFiniteScalar(phi);
}

bool VDBSDFField::GetBoundingSphereM(chrono::ChVector3d& center_M, double& radius) const {
    if (!impl_->has_bounding_sphere || !std::isfinite(impl_->bounding_radius)) {
        center_M = chrono::ChVector3d(0.0, 0.0, 0.0);
        radius = 0.0;
        return false;
    }

    center_M = impl_->bounding_center_M;
    radius = impl_->bounding_radius;
    return true;
}

bool VDBSDFField::QueryPhiGradHessianM(const chrono::ChVector3d& x_M,
                                       double& phi,
                                       chrono::ChVector3d& grad_M,
                                       chrono::ChMatrix33<>& hessian_M) const {
    if (!QueryPhiGradM(x_M, phi, grad_M) || !impl_->open_grid) {
        return false;
    }

    const chrono::ChVector3d h = GetVoxelStep(*impl_->open_grid);
    const chrono::ChVector3d dx(h.x(), 0.0, 0.0);
    const chrono::ChVector3d dy(0.0, h.y(), 0.0);
    const chrono::ChVector3d dz(0.0, 0.0, h.z());

    auto phi_func = [&](const chrono::ChVector3d& point) {
        return SamplePhiWorld(*impl_->open_grid, point);
    };

    const double phi0 = phi;
    const double phi_xp = phi_func(x_M + dx);
    const double phi_xm = phi_func(x_M - dx);
    const double phi_yp = phi_func(x_M + dy);
    const double phi_ym = phi_func(x_M - dy);
    const double phi_zp = phi_func(x_M + dz);
    const double phi_zm = phi_func(x_M - dz);

    hessian_M(0, 0) = (phi_xp - 2.0 * phi0 + phi_xm) / (h.x() * h.x());
    hessian_M(1, 1) = (phi_yp - 2.0 * phi0 + phi_ym) / (h.y() * h.y());
    hessian_M(2, 2) = (phi_zp - 2.0 * phi0 + phi_zm) / (h.z() * h.z());

    hessian_M(0, 1) = hessian_M(1, 0) =
        (phi_func(x_M + dx + dy) - phi_func(x_M + dx - dy) -
         phi_func(x_M - dx + dy) + phi_func(x_M - dx - dy)) /
        (4.0 * h.x() * h.y());
    hessian_M(0, 2) = hessian_M(2, 0) =
        (phi_func(x_M + dx + dz) - phi_func(x_M + dx - dz) -
         phi_func(x_M - dx + dz) + phi_func(x_M - dx - dz)) /
        (4.0 * h.x() * h.z());
    hessian_M(1, 2) = hessian_M(2, 1) =
        (phi_func(x_M + dy + dz) - phi_func(x_M + dy - dz) -
         phi_func(x_M - dy + dz) + phi_func(x_M - dy - dz)) /
        (4.0 * h.y() * h.z());

    return true;
}

bool VDBSDFField::FinalizeOpenGrid() {
    if (!impl_->open_grid) {
        impl_->last_error = "VDBSDFField::FinalizeOpenGrid missing OpenVDB float grid";
        impl_->ready = false;
        return false;
    }

    impl_->ready = true;
    return true;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
