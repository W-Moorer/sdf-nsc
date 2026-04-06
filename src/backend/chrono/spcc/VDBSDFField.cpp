#include "platform/backend/spcc/VDBSDFField.h"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <openvdb/io/File.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/math/SampleFromVoxels.h>
#include <nanovdb/tools/CreateNanoGrid.h>

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

double GetEnvDouble(const char* name, double fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || !std::isfinite(parsed)) {
        return fallback;
    }
    return parsed;
}

chrono::ChVector3d SafeNormalize(const chrono::ChVector3d& v,
                                 const chrono::ChVector3d& fallback) {
    const double n = v.Length();
    if (!std::isfinite(n) || n <= 1e-12) {
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

}  // namespace

struct VDBSDFField::Impl {
    openvdb::FloatGrid::Ptr open_grid;
    nanovdb::GridHandle<nanovdb::HostBuffer> nano_handle;
    const nanovdb::NanoGrid<float>* nano_grid = nullptr;
    bool direct_phi_hessian = false;
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
    impl_->nano_handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
    impl_->nano_grid = nullptr;
    impl_->direct_phi_hessian = options.direct_phi_hessian;
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
        for (const auto& v : vertices) {
            if (!IsFiniteVec(v)) {
                impl_->last_error = "VDBSDFField::BuildFromTriangleMesh mesh has non-finite vertex";
                return false;
            }
            points.emplace_back(
                ToFiniteFloat(v.x(), "VDBSDFField::BuildFromTriangleMesh vertex x"),
                ToFiniteFloat(v.y(), "VDBSDFField::BuildFromTriangleMesh vertex y"),
                ToFiniteFloat(v.z(), "VDBSDFField::BuildFromTriangleMesh vertex z"));
        }

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
            *transform,
            points,
            tri_faces,
            quad_faces,
            half_band,
            half_band);

        if (!impl_->open_grid) {
            impl_->last_error = "VDBSDFField::BuildFromTriangleMesh failed to create OpenVDB float grid";
            return false;
        }

        if (!options.grid_name.empty()) {
            impl_->open_grid->setName(options.grid_name);
        }

        return BuildNanoFromOpenGrid();
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] OpenVDB Build Exception: " << e.what() << std::endl;
        impl_->last_error = std::string("VDBSDFField::BuildFromTriangleMesh exception: ") + e.what();
        impl_->ready = false;
        impl_->open_grid.reset();
        impl_->nano_handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
        impl_->nano_grid = nullptr;
        return false;
    }
}

bool VDBSDFField::LoadFromVDBFile(const std::string& vdb_path,
                                  const std::string& grid_name) {
    EnsureOpenVDBInitialized();

    impl_->last_error.clear();
    impl_->ready = false;
    impl_->open_grid.reset();
    impl_->nano_handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
    impl_->nano_grid = nullptr;

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
        return BuildNanoFromOpenGrid();
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
    if (!impl_->ready || !impl_->nano_grid || !IsFiniteVec(x_M)) {
        return false;
    }

    const nanovdb::Vec3d x_world(x_M.x(), x_M.y(), x_M.z());
    const auto x_index = impl_->nano_grid->worldToIndex(x_world);

    const auto accessor = impl_->nano_grid->getAccessor();
    const auto sampler = nanovdb::math::createSampler<1>(accessor);

    const float phi_value = sampler(x_index);
    
    // SampleFromVoxels<3> doesn't have .gradient(), use finite difference in index space
    double hd = 0.5; // Half voxel step
    nanovdb::Vec3d grad_index(
        (sampler(nanovdb::Vec3d(x_index[0] + hd, x_index[1], x_index[2])) - sampler(nanovdb::Vec3d(x_index[0] - hd, x_index[1], x_index[2]))) / (2.0 * hd),
        (sampler(nanovdb::Vec3d(x_index[0], x_index[1] + hd, x_index[2])) - sampler(nanovdb::Vec3d(x_index[0], x_index[1] - hd, x_index[2]))) / (2.0 * hd),
        (sampler(nanovdb::Vec3d(x_index[0], x_index[1], x_index[2] + hd)) - sampler(nanovdb::Vec3d(x_index[0], x_index[1], x_index[2] - hd))) / (2.0 * hd)
    );
    const auto grad_world = impl_->nano_grid->indexToWorldGrad(grad_index);

    phi = static_cast<double>(phi_value);
    const chrono::ChVector3d grad_raw(
        static_cast<double>(grad_world[0]),
        static_cast<double>(grad_world[1]),
        static_cast<double>(grad_world[2]));
    grad_M = SafeNormalize(grad_raw, chrono::ChVector3d(1.0, 0.0, 0.0));

    if (!IsFiniteScalar(phi) || !IsFiniteVec(grad_M)) {
        return false;
    }
    return true;
}

bool VDBSDFField::QueryPhiGradHessianM(const chrono::ChVector3d& x_M,
                                       double& phi,
                                       chrono::ChVector3d& grad_M,
                                       chrono::ChMatrix33<>& hessian_M) const {
    if (!QueryPhiGradM(x_M, phi, grad_M)) {
        return false;
    }

    // Finite difference step (2*voxel_size) for Hessian
    double h = impl_->nano_grid->voxelSize()[0] * 2.0;

    if (impl_->direct_phi_hessian) {
        auto phi_func = [&](const chrono::ChVector3d& point) -> double {
            double d = 0.0;
            chrono::ChVector3d g(0, 1, 0);
            QueryPhiGradM(point, d, g);
            return d;
        };

        const chrono::ChVector3d dx(h, 0, 0);
        const chrono::ChVector3d dy(0, h, 0);
        const chrono::ChVector3d dz(0, 0, h);
        const double phi0 = phi;

        const double phi_xp = phi_func(x_M + dx);
        const double phi_xm = phi_func(x_M - dx);
        const double phi_yp = phi_func(x_M + dy);
        const double phi_ym = phi_func(x_M - dy);
        const double phi_zp = phi_func(x_M + dz);
        const double phi_zm = phi_func(x_M - dz);

        hessian_M(0, 0) = (phi_xp - 2.0 * phi0 + phi_xm) / (h * h);
        hessian_M(1, 1) = (phi_yp - 2.0 * phi0 + phi_ym) / (h * h);
        hessian_M(2, 2) = (phi_zp - 2.0 * phi0 + phi_zm) / (h * h);

        hessian_M(0, 1) = hessian_M(1, 0) =
            (phi_func(x_M + dx + dy) - phi_func(x_M + dx - dy) -
             phi_func(x_M - dx + dy) + phi_func(x_M - dx - dy)) / (4.0 * h * h);
        hessian_M(0, 2) = hessian_M(2, 0) =
            (phi_func(x_M + dx + dz) - phi_func(x_M + dx - dz) -
             phi_func(x_M - dx + dz) + phi_func(x_M - dx - dz)) / (4.0 * h * h);
        hessian_M(1, 2) = hessian_M(2, 1) =
            (phi_func(x_M + dy + dz) - phi_func(x_M + dy - dz) -
             phi_func(x_M - dy + dz) + phi_func(x_M - dy - dz)) / (4.0 * h * h);
        return true;
    }

    auto grad_func = [&](const chrono::ChVector3d& point) -> chrono::ChVector3d {
        double d;
        chrono::ChVector3d g(0, 1, 0); // default fallback
        QueryPhiGradM(point, d, g);
        return g;
    };

    chrono::ChVector3d dx(h, 0, 0);
    chrono::ChVector3d dy(0, h, 0);
    chrono::ChVector3d dz(0, 0, h);

    chrono::ChVector3d gx_p = grad_func(x_M + dx);
    chrono::ChVector3d gx_m = grad_func(x_M - dx);
    chrono::ChVector3d gy_p = grad_func(x_M + dy);
    chrono::ChVector3d gy_m = grad_func(x_M - dy);
    chrono::ChVector3d gz_p = grad_func(x_M + dz);
    chrono::ChVector3d gz_m = grad_func(x_M - dz);

    chrono::ChVector3d h_x = (gx_p - gx_m) / (2.0 * h);
    chrono::ChVector3d h_y = (gy_p - gy_m) / (2.0 * h);
    chrono::ChVector3d h_z = (gz_p - gz_m) / (2.0 * h);

    hessian_M(0, 0) = h_x.x(); hessian_M(0, 1) = h_x.y(); hessian_M(0, 2) = h_x.z();
    hessian_M(1, 0) = h_y.x(); hessian_M(1, 1) = h_y.y(); hessian_M(1, 2) = h_y.z();
    hessian_M(2, 0) = h_z.x(); hessian_M(2, 1) = h_z.y(); hessian_M(2, 2) = h_z.z();

    // symmetrize
    hessian_M = 0.5 * (hessian_M + hessian_M.transpose());
    return true;
}

bool VDBSDFField::BuildNanoFromOpenGrid() {
    if (!impl_->open_grid) {
        impl_->last_error = "VDBSDFField::BuildNanoFromOpenGrid missing OpenVDB float grid";
        impl_->ready = false;
        impl_->nano_grid = nullptr;
        impl_->nano_handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
        return false;
    }

    try {
        impl_->nano_handle = nanovdb::tools::createNanoGrid(*impl_->open_grid);
        impl_->nano_grid = impl_->nano_handle.grid<float>();
        impl_->ready = (impl_->nano_grid != nullptr);
        if (!impl_->ready) {
            impl_->last_error = "VDBSDFField::BuildNanoFromOpenGrid produced no float NanoVDB grid";
        }
        return impl_->ready;
    } catch (const std::exception& e) {
        impl_->last_error = std::string("VDBSDFField::BuildNanoFromOpenGrid exception: ") + e.what();
        impl_->ready = false;
        impl_->nano_grid = nullptr;
        impl_->nano_handle = nanovdb::GridHandle<nanovdb::HostBuffer>();
        return false;
    }
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform

