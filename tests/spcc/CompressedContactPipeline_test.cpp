#include "platform/backend/spcc/CompressedContactPipeline.h"

#include <gtest/gtest.h>

namespace {

class PlaneSDF final : public platform::backend::spcc::FirstOrderSDF {
  public:
    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        phi = x_M.y();
        grad_M = chrono::ChVector3d(0.0, 1.0, 0.0);
        return true;
    }
};

platform::backend::spcc::RigidBodyStateW MakeIdentityState() {
    platform::backend::spcc::RigidBodyStateW state;
    state.x_com_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    state.x_ref_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    state.v_com_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    state.w_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    return state;
}

}  // namespace

TEST(CompressedContactPipelineTest, ReducesDensePatchIntoLimitedSupportPoints) {
    platform::backend::spcc::CompressedContactPipeline pipeline;
    platform::backend::spcc::CompressedContactConfig cfg;
    cfg.delta_on = 0.02;
    cfg.delta_off = 0.03;
    cfg.max_active_dense = 64;
    cfg.patch_radius = 0.2;
    cfg.normal_cos_min = 0.95;
    cfg.max_patch_diameter = 0.4;
    cfg.max_reduced_points_per_patch = 4;
    pipeline.Configure(cfg);

    std::vector<platform::backend::spcc::DenseSurfaceSample> samples;
    for (int ix = -2; ix <= 2; ++ix) {
        for (int iz = -2; iz <= 2; ++iz) {
            platform::backend::spcc::DenseSurfaceSample sample;
            sample.xi_slave_S = chrono::ChVector3d(0.02 * ix, -0.01, 0.02 * iz);
            sample.normal_slave_S = chrono::ChVector3d(0.0, 1.0, 0.0);
            sample.area_weight = 1.0;
            samples.push_back(sample);
        }
    }
    pipeline.SetSlaveSurfaceSamples(samples);

    PlaneSDF sdf;
    const auto master_state = MakeIdentityState();
    const auto slave_state = MakeIdentityState();

    std::vector<platform::backend::spcc::ReducedContactPoint> reduced;
    platform::backend::spcc::CompressionStats stats;
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, 0.2, 1.0e-3, reduced, &stats);

    ASSERT_FALSE(reduced.empty());
    EXPECT_LE(reduced.size(), 4u);
    EXPECT_EQ(stats.patch_count, 1u);
    for (const auto& contact : reduced) {
        EXPECT_LT(contact.phi, 0.0);
        EXPECT_NEAR(contact.n_W.y(), 1.0, 1.0e-9);
        EXPECT_NEAR(contact.x_master_surface_W.y(), 0.0, 1.0e-9);
    }
}
