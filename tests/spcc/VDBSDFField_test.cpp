#include "platform/backend/spcc/VDBSDFField.h"

#include <cmath>

#include <gtest/gtest.h>

#if defined(SPCC_ENABLE_VDB)

namespace {

chrono::ChTriangleMeshConnected BuildUnitCubeMesh() {
    chrono::ChTriangleMeshConnected mesh;

    auto& v = mesh.GetCoordsVertices();
    v = {
        chrono::ChVector3d(-0.5, -0.5, -0.5),
        chrono::ChVector3d(0.5, -0.5, -0.5),
        chrono::ChVector3d(0.5, 0.5, -0.5),
        chrono::ChVector3d(-0.5, 0.5, -0.5),
        chrono::ChVector3d(-0.5, -0.5, 0.5),
        chrono::ChVector3d(0.5, -0.5, 0.5),
        chrono::ChVector3d(0.5, 0.5, 0.5),
        chrono::ChVector3d(-0.5, 0.5, 0.5),
    };

    auto& f = mesh.GetIndicesVertexes();
    f = {
        chrono::ChVector3i(0, 1, 2), chrono::ChVector3i(0, 2, 3),
        chrono::ChVector3i(4, 6, 5), chrono::ChVector3i(4, 7, 6),
        chrono::ChVector3i(0, 4, 5), chrono::ChVector3i(0, 5, 1),
        chrono::ChVector3i(1, 5, 6), chrono::ChVector3i(1, 6, 2),
        chrono::ChVector3i(2, 6, 7), chrono::ChVector3i(2, 7, 3),
        chrono::ChVector3i(3, 7, 4), chrono::ChVector3i(3, 4, 0),
    };

    return mesh;
}

TEST(VDBSDFFieldTest, BuildFromMeshAndQueryPhiGrad) {
    platform::backend::spcc::VDBSDFField sdf;

    platform::backend::spcc::VDBSDFField::BuildOptions options;
    options.voxel_size = 0.05;
    options.half_band_width_voxels = 4.0;

    const auto mesh = BuildUnitCubeMesh();
    ASSERT_TRUE(sdf.BuildFromTriangleMesh(mesh, options)) << sdf.LastError();
    ASSERT_TRUE(sdf.IsReady());

    double phi_center = 0.0;
    chrono::ChVector3d grad_center;
    ASSERT_TRUE(sdf.QueryPhiGradM(chrono::ChVector3d(0.0, 0.0, 0.0), phi_center, grad_center));
    EXPECT_LT(phi_center, 0.0);

    double phi_outside = 0.0;
    chrono::ChVector3d grad_outside;
    ASSERT_TRUE(sdf.QueryPhiGradM(chrono::ChVector3d(0.8, 0.0, 0.0), phi_outside, grad_outside));
    EXPECT_GT(phi_outside, 0.0);

    double phi_surface = 0.0;
    chrono::ChVector3d grad_surface;
    ASSERT_TRUE(sdf.QueryPhiGradM(chrono::ChVector3d(0.5, 0.0, 0.0), phi_surface, grad_surface));
    EXPECT_NEAR(phi_surface, 0.0, 0.08);
    EXPECT_TRUE(std::isfinite(grad_surface.x()));
    EXPECT_TRUE(std::isfinite(grad_surface.y()));
    EXPECT_TRUE(std::isfinite(grad_surface.z()));
    EXPECT_GT(grad_surface.Length(), 0.1);
}

TEST(VDBSDFFieldTest, RejectsInvalidOptions) {
    platform::backend::spcc::VDBSDFField sdf;
    const auto mesh = BuildUnitCubeMesh();

    platform::backend::spcc::VDBSDFField::BuildOptions bad_options;
    bad_options.voxel_size = 0.0;
    EXPECT_FALSE(sdf.BuildFromTriangleMesh(mesh, bad_options));
}

}  // namespace

#else

TEST(VDBSDFFieldTest, BuildDisabled) {
    SUCCEED();
}

#endif
