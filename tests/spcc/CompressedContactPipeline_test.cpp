#include "platform/backend/spcc/CompressedContactPipeline.h"

#include <algorithm>

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

class SphereSDF final : public platform::backend::spcc::FirstOrderSDF {
  public:
    explicit SphereSDF(double radius) : radius_(radius) {}

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        const double dist = x_M.Length();
        phi = dist - radius_;
        if (dist > 1.0e-12) {
            grad_M = x_M * (1.0 / dist);
        } else {
            grad_M = chrono::ChVector3d(1.0, 0.0, 0.0);
        }
        return true;
    }

    bool GetBoundingSphereM(chrono::ChVector3d& center_M, double& radius) const override {
        center_M = chrono::ChVector3d(0.0, 0.0, 0.0);
        radius = radius_;
        return true;
    }

  private:
    double radius_ = 0.0;
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
    cfg.max_subpatch_diameter = 0.0;
    cfg.max_plane_error = 0.0;
    cfg.sentinel_spacing = 0.0;
    cfg.sentinel_margin = 0.0;
    cfg.max_subpatch_depth = 0;
    cfg.min_dense_points_per_subpatch = 0;
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
    EXPECT_EQ(stats.subpatch_count, 1u);
    for (const auto& contact : reduced) {
        EXPECT_LT(contact.phi, 0.0);
        EXPECT_NEAR(contact.n_W.y(), 1.0, 1.0e-9);
        EXPECT_NEAR(contact.x_master_surface_W.y(), 0.0, 1.0e-9);
    }
}

TEST(CompressedContactPipelineTest, UsesBVHToReduceExactDenseQueries) {
    platform::backend::spcc::CompressedContactPipeline pipeline;
    platform::backend::spcc::CompressedContactConfig cfg;
    cfg.delta_on = 0.02;
    cfg.delta_off = 0.03;
    cfg.max_active_dense = 0;
    cfg.bvh_leaf_size = 12;
    cfg.bvh_query_margin = 0.0;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 0.2;
    cfg.normal_cos_min = 0.8;
    cfg.max_patch_diameter = 0.4;
    cfg.max_subpatch_diameter = 0.0;
    cfg.max_plane_error = 0.0;
    cfg.sentinel_spacing = 0.0;
    cfg.sentinel_margin = 0.0;
    cfg.max_subpatch_depth = 0;
    cfg.min_dense_points_per_subpatch = 0;
    cfg.max_reduced_points_per_patch = 4;
    pipeline.Configure(cfg);

    std::vector<platform::backend::spcc::DenseSurfaceSample> samples;
    for (int iy = -30; iy <= 30; ++iy) {
        for (int iz = -30; iz <= 30; ++iz) {
            platform::backend::spcc::DenseSurfaceSample sample;
            sample.xi_slave_S = chrono::ChVector3d(0.075, 0.008 * iy, 0.008 * iz);
            sample.normal_slave_S = chrono::ChVector3d(-1.0, 0.0, 0.0);
            sample.area_weight = 1.0;
            samples.push_back(sample);
        }
    }
    pipeline.SetSlaveSurfaceSamples(samples);

    SphereSDF sdf(0.08);
    const auto master_state = MakeIdentityState();
    const auto slave_state = MakeIdentityState();

    std::vector<platform::backend::spcc::ReducedContactPoint> reduced;
    platform::backend::spcc::CompressionStats stats;
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, 0.2, 1.0e-3, reduced, &stats);

    ASSERT_FALSE(reduced.empty());
    EXPECT_EQ(stats.total_samples, samples.size());
    EXPECT_LT(stats.candidate_count, samples.size());
    EXPECT_GT(stats.bvh_nodes_pruned_obb + stats.bvh_nodes_pruned_sdf, 0u);
}

TEST(CompressedContactPipelineTest, SplitsLongStripIntoMultipleSubpatches) {
    platform::backend::spcc::CompressedContactPipeline pipeline;
    platform::backend::spcc::CompressedContactConfig cfg;
    cfg.delta_on = 0.02;
    cfg.delta_off = 0.03;
    cfg.max_active_dense = 0;
    cfg.bvh_leaf_size = 12;
    cfg.bvh_enable_sdf_node_bound = false;
    cfg.patch_radius = 0.5;
    cfg.normal_cos_min = 0.95;
    cfg.max_patch_diameter = 0.5;
    cfg.max_subpatch_diameter = 0.08;
    cfg.max_plane_error = 1.0e-4;
    cfg.sentinel_spacing = 0.02;
    cfg.sentinel_margin = 0.002;
    cfg.max_subpatch_depth = 4;
    cfg.min_dense_points_per_subpatch = 12;
    cfg.max_reduced_points_per_patch = 4;
    pipeline.Configure(cfg);

    std::vector<platform::backend::spcc::DenseSurfaceSample> samples;
    for (int ix = -20; ix <= 20; ++ix) {
        for (int iz = -2; iz <= 2; ++iz) {
            platform::backend::spcc::DenseSurfaceSample sample;
            sample.xi_slave_S = chrono::ChVector3d(0.01 * ix, -0.01, 0.01 * iz);
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
    EXPECT_EQ(stats.patch_count, 1u);
    EXPECT_GT(stats.subpatch_count, 1u);
    EXPECT_GT(stats.reduced_count, 4u);
}

TEST(CompressedContactPipelineTest, KeepsChainConnectedContactRegionAsSinglePatch) {
    platform::backend::spcc::CompressedContactPipeline pipeline;
    platform::backend::spcc::CompressedContactConfig cfg;
    cfg.delta_on = 0.02;
    cfg.delta_off = 0.03;
    cfg.max_active_dense = 0;
    cfg.patch_radius = 0.06;
    cfg.normal_cos_min = 0.95;
    cfg.max_patch_diameter = 0.0;
    cfg.max_subpatch_diameter = 0.0;
    cfg.max_plane_error = 0.0;
    cfg.sentinel_spacing = 0.0;
    cfg.sentinel_margin = 0.0;
    cfg.max_subpatch_depth = 0;
    cfg.min_dense_points_per_subpatch = 0;
    cfg.max_reduced_points_per_patch = 4;
    pipeline.Configure(cfg);

    std::vector<platform::backend::spcc::DenseSurfaceSample> samples;
    for (int ix = 0; ix < 4; ++ix) {
        platform::backend::spcc::DenseSurfaceSample sample;
        sample.xi_slave_S = chrono::ChVector3d(0.05 * ix, -0.01, 0.0);
        sample.normal_slave_S = chrono::ChVector3d(0.0, 1.0, 0.0);
        sample.area_weight = 1.0;
        samples.push_back(sample);
    }
    pipeline.SetSlaveSurfaceSamples(samples);

    PlaneSDF sdf;
    const auto master_state = MakeIdentityState();
    const auto slave_state = MakeIdentityState();

    std::vector<platform::backend::spcc::ReducedContactPoint> reduced;
    platform::backend::spcc::CompressionStats stats;
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, 0.2, 1.0e-3, reduced, &stats);

    ASSERT_FALSE(reduced.empty());
    EXPECT_EQ(stats.patch_count, 1u);
    EXPECT_EQ(stats.subpatch_count, 1u);
}

TEST(CompressedContactPipelineTest, PreservesPersistentIdAcrossSmallMotion) {
    platform::backend::spcc::CompressedContactPipeline pipeline;
    platform::backend::spcc::CompressedContactConfig cfg;
    cfg.delta_on = 0.02;
    cfg.delta_off = 0.03;
    cfg.max_active_dense = 0;
    cfg.patch_radius = 0.2;
    cfg.normal_cos_min = 0.95;
    cfg.max_patch_diameter = 0.4;
    cfg.max_subpatch_diameter = 0.0;
    cfg.max_plane_error = 0.0;
    cfg.sentinel_spacing = 0.0;
    cfg.sentinel_margin = 0.0;
    cfg.max_subpatch_depth = 0;
    cfg.min_dense_points_per_subpatch = 0;
    cfg.max_reduced_points_per_patch = 4;
    cfg.warm_start_match_radius = 0.02;
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
    auto slave_state = MakeIdentityState();

    std::vector<platform::backend::spcc::ReducedContactPoint> reduced_a;
    platform::backend::spcc::CompressionStats stats_a;
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, 0.2, 1.0e-3, reduced_a, &stats_a);

    slave_state.x_ref_W = chrono::ChVector3d(5.0e-4, 0.0, 0.0);
    std::vector<platform::backend::spcc::ReducedContactPoint> reduced_b;
    platform::backend::spcc::CompressionStats stats_b;
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, 0.2, 1.0e-3, reduced_b, &stats_b);

    ASSERT_FALSE(reduced_a.empty());
    ASSERT_FALSE(reduced_b.empty());
    EXPECT_EQ(stats_a.patch_count, 1u);
    EXPECT_EQ(stats_b.patch_count, 1u);
    ASSERT_EQ(reduced_a.size(), reduced_b.size());
    for (const auto& contact : reduced_a) {
        EXPECT_NE(contact.persistent_id, 0u);
    }

    std::vector<std::size_t> support_ids_a;
    std::vector<std::size_t> support_ids_b;
    for (const auto& contact : reduced_a) {
        support_ids_a.push_back(contact.support_id);
    }
    for (const auto& contact : reduced_b) {
        EXPECT_EQ(contact.persistent_id, reduced_a.front().persistent_id);
        support_ids_b.push_back(contact.support_id);
    }
    std::sort(support_ids_a.begin(), support_ids_a.end());
    std::sort(support_ids_b.begin(), support_ids_b.end());
    EXPECT_EQ(support_ids_b, support_ids_a);
}
