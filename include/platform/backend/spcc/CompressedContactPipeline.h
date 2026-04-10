#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/DenseContactCloud.h"
#include "platform/backend/spcc/DenseSampleBVH.h"
#include "platform/backend/spcc/DenseSurfaceSampler.h"
#include "platform/backend/spcc/FirstOrderSDF.h"

namespace platform {
namespace backend {
namespace spcc {

struct ReducedContactPoint {
    std::size_t persistent_id = 0;
    std::size_t patch_id = 0;
    std::size_t subpatch_id = 0;
    std::size_t support_id = 0;
    std::size_t dense_members = 0;
    int emission_count = 1;

    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d n_W;
    chrono::ChVector3d v_rel_W;
    chrono::ChVector3d stencil_axis_W;
    chrono::ChVector3d stencil_axis_secondary_W;

    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
    double support_weight = 0.0;
    double allocated_load = 0.0;
    chrono::ChVector3d allocated_force_W;
    double stencil_half_extent = 0.0;
    double stencil_half_extent_secondary = 0.0;
    std::array<double, 5> stencil_gap_offsets{};
    std::array<chrono::ChVector3d, 5> slot_x_W{};
    std::array<chrono::ChVector3d, 5> slot_x_master_surface_W{};
    std::array<chrono::ChVector3d, 5> slot_n_W{};
    std::array<double, 5> slot_weight{};
    double mu = 0.0;
    // Fixed cache slots: center, primary negative/positive, secondary negative/positive.
    std::array<float, 6> reaction_cache_primary{};
    std::array<float, 6> reaction_cache_secondary{};
    std::array<float, 6> reaction_cache_tertiary{};
    std::array<float, 6> reaction_cache_quaternary{};
    std::array<float, 6> reaction_cache_quinary{};
};

struct CompressionStats {
    std::size_t total_samples = 0;
    std::size_t candidate_count = 0;
    std::size_t phi_prefilter_count = 0;
    std::size_t exact_count = 0;
    std::size_t dense_count = 0;
    std::size_t reduced_count = 0;
    std::size_t patch_count = 0;
    std::size_t subpatch_count = 0;
    std::size_t bvh_nodes_visited = 0;
    std::size_t bvh_nodes_pruned_obb = 0;
    std::size_t bvh_nodes_pruned_sdf = 0;
    std::size_t bvh_leaf_samples_tested = 0;
    double epsilon_F = 0.0;
    double epsilon_M = 0.0;
    double epsilon_CoP = 0.0;
    double epsilon_gap = 0.0;
    double dense_worst_gap = 0.0;
    double reduced_worst_gap = 0.0;
    double max_subpatch_plane_error = 0.0;
    double max_subpatch_second_moment_error = 0.0;
    double max_subpatch_cone_error = 0.0;
    double max_subpatch_gap_error = 0.0;
    double max_subpatch_force_residual = 0.0;
    double max_subpatch_moment_residual = 0.0;
    double max_subpatch_reference_wrench_error = 0.0;
    double max_subpatch_reference_cop_error = 0.0;
    double max_dense_micro_force_residual = 0.0;
    double max_dense_micro_moment_residual = 0.0;
    double dense_cloud_ms = 0.0;
    double patch_ms = 0.0;
    double subpatch_ms = 0.0;
    double dense_micro_ms = 0.0;
    double support_build_ms = 0.0;
    double allocation_ms = 0.0;
    double reinjection_ms = 0.0;
    double total_pipeline_ms = 0.0;
};

struct TemporalSubpatchState {
    std::size_t persistent_id = 0;
    chrono::ChVector3d centroid_W;
    chrono::ChVector3d avg_normal_W;
    double diameter = 0.0;
    chrono::ChVector3d reference_origin_W;
    chrono::ChVector3d reference_force_W;
    chrono::ChVector3d reference_moment_W;
    double reference_total_load = 0.0;
    chrono::ChVector3d impulse_origin_W;
    chrono::ChVector3d impulse_force_W;
    chrono::ChVector3d impulse_moment_W;
    double impulse_total_load = 0.0;
    bool has_impulse_wrench = false;
    std::vector<ReducedContactPoint> contacts;
};

class CompressedContactPipeline {
  public:
    void Configure(const CompressedContactConfig& cfg);
    void SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples);
    void SyncTemporalWarmStart(const std::vector<ReducedContactPoint>& emitted_contacts) const;

    void BuildReducedContacts(const RigidBodyStateW& master_state,
                              const RigidBodyStateW& slave_state,
                              const FirstOrderSDF& sdf,
                              double mu_default,
                              double step_size,
                              std::vector<ReducedContactPoint>& out_contacts,
                              CompressionStats* out_stats = nullptr) const;

  private:
    CompressedContactConfig cfg_;
    std::vector<DenseSurfaceSample> slave_surface_samples_;
    DenseSampleBVH dense_sample_bvh_;
    mutable std::vector<ReducedContactPoint> previous_contacts_;
    mutable std::vector<TemporalSubpatchState> previous_subpatches_;
    mutable std::size_t next_persistent_id_ = 1;
    mutable std::size_t build_step_counter_ = 0;
    mutable double previous_step_size_ = 0.0;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
