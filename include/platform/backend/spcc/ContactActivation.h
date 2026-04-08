#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/SPCCProblem.h"
#include "platform/backend/spcc/ContactTuning.h"

namespace platform {
namespace backend {
namespace spcc {

// Minimal SDF query interface for ContactActivation.
// Real SDF implementations can inherit from this interface.
class SDFField {
  public:
    virtual ~SDFField() = default;

    // Input : x_M (master local point)
    // Output: phi and grad_M
    virtual bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                               double& phi,
                               chrono::ChVector3d& grad_M) const = 0;

    // Output: phi, grad_M, and hessian_M
    virtual bool QueryPhiGradHessianM(const chrono::ChVector3d& x_M,
                                      double& phi,
                                      chrono::ChVector3d& grad_M,
                                      chrono::ChMatrix33<>& hessian_M) const {
        // default fallback
        hessian_M.setZero();
        return QueryPhiGradM(x_M, phi, grad_M);
    }
};

class ContactActivation {
  public:
    struct Stats {
        std::size_t queried = 0;
        std::size_t accepted_before_cap = 0;
        std::size_t accepted_after_cap = 0;
        std::size_t rejected_invalid = 0;
        std::size_t local_fit_attempted = 0;
        std::size_t local_fit_applied = 0;
        std::size_t local_fit_rejected_positive_gap = 0;
    };

    void SetEnvPrefix(const std::string& env_prefix);
    void SetPolicy(const ContactRegimeConfig& policy);
    void Configure(double delta_on, double delta_off, int hold_steps, std::size_t max_active_keep);

    void Reset(std::size_t sample_count);

    // first-pass assumption: local_samples_S.size() == world_samples_W.size()
    void BuildActiveSet(const RigidBodyStateW& master_pred,
                        const RigidBodyStateW& slave_pred,
                        const SDFField& sdf,
                        const std::vector<chrono::ChVector3d>& local_samples_S,
                        const std::vector<chrono::ChVector3d>& world_samples_W,
                        double mu_default,
                        double step_size,
                        bool need_hessian,
                        std::vector<ActiveContactSample>& out_active);

    const Stats& GetStats() const { return stats_; }

  private:
    struct PersistentPatchState {
        std::size_t manifold_id = 0;
        chrono::ChVector3d x_W;
        chrono::ChVector3d x_master_M;
        chrono::ChVector3d n_W;
        double phi = 0.0;
        int age = 0;
    };

    void EnsureLocalNeighborGraph(const std::vector<chrono::ChVector3d>& local_samples_S, double radius);
    void ApplySlidingPersistentManifold(const RigidBodyStateW& master_pred,
                                        const RigidBodyStateW& slave_pred,
                                        const SDFField& sdf,
                                        bool need_hessian,
                                        std::vector<ActiveContactSample>& patches);
    void RefreshPatchKinematics(const RigidBodyStateW& master_pred,
                                const RigidBodyStateW& slave_pred,
                                ActiveContactSample& patch) const;
    bool QueryPatchGeometry(const RigidBodyStateW& master_pred,
                            const SDFField& sdf,
                            bool need_hessian,
                            const chrono::ChVector3d& query_master_M,
                            const chrono::ChVector3d& reference_normal_W,
                            ActiveContactSample& patch) const;
    bool QuerySweptPatchGeometry(const RigidBodyStateW& master_pred,
                                 const SDFField& sdf,
                                 bool need_hessian,
                                 const chrono::ChVector3d& start_master_M,
                                 const chrono::ChVector3d& end_master_M,
                                 const chrono::ChVector3d& reference_normal_W,
                                 int path_samples,
                                 ActiveContactSample& patch) const;
    void SelectSlidingCoveragePatches(std::vector<ActiveContactSample>& patches) const;

    std::string env_prefix_;
    ContactRegimeType regime_ = ContactRegimeType::CompactDiscrete;
    PatchGeometryMode patch_geometry_mode_ = PatchGeometryMode::DeepestPoint;
    double delta_on_ = 0.0;
    double delta_off_ = 0.0;
    int hold_steps_ = 0;
    std::size_t max_active_keep_ = 0;
    int full_scan_period_default_ = 1;
    double local_scan_radius_default_ = 0.0;
    double cluster_angle_deg_default_ = 30.0;
    double separating_cutoff_default_ = 1.0e-4;
    double cluster_radius_default_ = 0.0;
    bool avg_point_default_ = false;
    double persistent_match_radius_default_ = 0.0;
    double persistent_normal_cos_min_default_ = -1.0;
    double persistent_blend_alpha_default_ = 1.0;
    int persistent_path_samples_default_ = 1;
    double coverage_spacing_radius_default_ = 0.0;
    int onset_refine_steps_default_ = 0;
    int onset_refine_path_samples_default_ = 1;
    double onset_refine_backtrack_scale_default_ = 1.0;
    int local_fit_onset_steps_default_ = 0;
    int local_fit_min_cluster_size_default_ = 3;
    int local_fit_path_samples_default_ = 1;
    double local_fit_max_shift_ratio_default_ = 0.5;
    double local_fit_blend_default_ = 1.0;
    double local_fit_reject_positive_phi_default_ = -1.0;
    double onset_gate_current_phi_max_default_ = -1.0;
    int single_point_local_fit_path_samples_default_ = 0;
    double single_point_local_fit_backtrack_scale_default_ = 1.0;
    bool curvature_gate_enabled_ = false;
    bool curvature_tangential_only_default_ = false;
    double normal_alignment_cos_min_default_ = -1.0;
    double max_hessian_frobenius_default_ = 0.0;
    double max_curvature_term_abs_default_ = 0.0;
    double max_curvature_term_ratio_default_ = 0.0;
    double curvature_gap_floor_default_ = 0.0;
    int curvature_ramp_steps_default_ = 0;

    std::vector<uint8_t> prev_active_;
    std::vector<int> hold_counter_;
    std::vector<int> active_age_;
    std::vector<std::vector<std::uint32_t>> local_neighbors_;
    std::vector<PersistentPatchState> persistent_patches_;
    double local_neighbor_radius_ = -1.0;
    std::size_t build_step_ = 0;
    std::size_t next_manifold_id_ = 1;
    Stats stats_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
