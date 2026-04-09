#pragma once

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

    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
    double support_weight = 0.0;
    double allocated_load = 0.0;
    double mu = 0.0;
};

struct CompressionStats {
    std::size_t total_samples = 0;
    std::size_t candidate_count = 0;
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
    double max_subpatch_plane_error = 0.0;
    double max_subpatch_gap_error = 0.0;
    double max_subpatch_force_residual = 0.0;
    double max_subpatch_moment_residual = 0.0;
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
    std::vector<ReducedContactPoint> contacts;
};

class CompressedContactPipeline {
  public:
    void Configure(const CompressedContactConfig& cfg);
    void SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples);

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
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
