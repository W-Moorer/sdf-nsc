#include "platform/backend/spcc/ContactActivation.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <numeric>

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

bool NormalizeOrReject(chrono::ChVector3d& v) {
    const double n = v.Length();
    if (!std::isfinite(n) || n <= 1e-12) {
        return false;
    }
    v *= (1.0 / n);
    return true;
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

chrono::ChMatrix33<> BuildProjectorFromNormal(const chrono::ChVector3d& n_W) {
    chrono::ChMatrix33<> P(1);
    P(0, 0) -= n_W.x() * n_W.x();
    P(0, 1) -= n_W.x() * n_W.y();
    P(0, 2) -= n_W.x() * n_W.z();
    P(1, 0) -= n_W.y() * n_W.x();
    P(1, 1) -= n_W.y() * n_W.y();
    P(1, 2) -= n_W.y() * n_W.z();
    P(2, 0) -= n_W.z() * n_W.x();
    P(2, 1) -= n_W.z() * n_W.y();
    P(2, 2) -= n_W.z() * n_W.z();
    return P;
}

}  // namespace

void ContactActivation::Configure(double delta_on,
                                  double delta_off,
                                  int hold_steps,
                                  std::size_t max_active_keep) {
    delta_on_ = delta_on;
    delta_off_ = delta_off;
    hold_steps_ = std::max(0, hold_steps);
    max_active_keep_ = max_active_keep;

    if (delta_off_ < delta_on_) {
        delta_off_ = delta_on_;
    }
}

void ContactActivation::Reset(std::size_t sample_count) {
    prev_active_.assign(sample_count, 0);
    hold_counter_.assign(sample_count, 0);
    active_age_.assign(sample_count, 0);
    stats_ = Stats{};
}

void ContactActivation::BuildActiveSet(const RigidBodyStateW& master_pred,
                                       const RigidBodyStateW& slave_pred,
                                       const SDFField& sdf,
                                       const std::vector<chrono::ChVector3d>& local_samples_S,
                                       const std::vector<chrono::ChVector3d>& world_samples_W,
                                       double mu_default,
                                       std::vector<ActiveContactSample>& out_active) {
    out_active.clear();
    stats_ = Stats{};

    // first-pass assumption: local and world samples are pre-aligned by index.
    if (local_samples_S.size() != world_samples_W.size()) {
        return;
    }

    const std::size_t n_samples = local_samples_S.size();
    if (prev_active_.size() != n_samples) {
        Reset(n_samples);
    }

    const double mu_use = (std::isfinite(mu_default) && mu_default >= 0.0) ? mu_default : 0.0;

    std::vector<ActiveContactSample> candidates;
    candidates.reserve(n_samples);

    for (std::size_t i = 0; i < n_samples; ++i) {
        stats_.queried += 1;

        const chrono::ChVector3d& xi_S = local_samples_S[i];
        const chrono::ChVector3d& x_W = world_samples_W[i];

        if (!IsFiniteVec(xi_S) || !IsFiniteVec(x_W)) {
            stats_.rejected_invalid += 1;
            prev_active_[i] = 0;
            hold_counter_[i] = 0;
            active_age_[i] = 0;
            continue;
        }

        const chrono::ChVector3d x_master_M = master_pred.R_WRef.transpose() * (x_W - master_pred.x_ref_W);

        double phi = 0.0;
        chrono::ChVector3d grad_M;
        chrono::ChMatrix33<> hessian_M(0);
        if (!sdf.QueryPhiGradHessianM(x_master_M, phi, grad_M, hessian_M)) {
            stats_.rejected_invalid += 1;
            prev_active_[i] = 0;
            hold_counter_[i] = 0;
            active_age_[i] = 0;
            continue;
        }

        if (!IsFiniteScalar(phi) || !IsFiniteVec(grad_M)) {
            stats_.rejected_invalid += 1;
            prev_active_[i] = 0;
            hold_counter_[i] = 0;
            active_age_[i] = 0;
            continue;
        }

        chrono::ChVector3d n_W = master_pred.R_WRef * grad_M;
        if (!NormalizeOrReject(n_W)) {
            stats_.rejected_invalid += 1;
            prev_active_[i] = 0;
            hold_counter_[i] = 0;
            active_age_[i] = 0;
            continue;
        }

        const chrono::ChMatrix33<> P_W = BuildProjectorFromNormal(n_W);

        const bool was_active = (prev_active_[i] != 0);
        bool is_active = false;

        if (!was_active) {
            is_active = (phi <= delta_on_);
        } else {
            is_active = (phi <= delta_off_) || (hold_counter_[i] > 0);
        }

        if (!is_active) {
            prev_active_[i] = 0;
            hold_counter_[i] = 0;
            active_age_[i] = 0;
            continue;
        }

        if (!was_active) {
            hold_counter_[i] = hold_steps_;
            active_age_[i] = 1;
        } else {
            if (hold_counter_[i] > 0) {
                hold_counter_[i] -= 1;
            }
            active_age_[i] += 1;
        }
        prev_active_[i] = 1;

        ActiveContactSample s;
        s.sample_id = i;
        s.xi_slave_S = xi_S;
        s.x_W = x_W;
        s.x_master_M = x_master_M;
        s.phi = phi;
        s.grad_M = grad_M;
        s.hessian_M = hessian_M;
        s.n_W = n_W;
        s.hessian_W = master_pred.R_WRef * hessian_M * master_pred.R_WRef.transpose();
        s.P_W = P_W;

        s.rA_W = x_W - master_pred.x_com_W;
        s.rB_W = x_W - slave_pred.x_com_W;

        s.mu = mu_use;

        const chrono::ChVector3d vA_point_W = master_pred.v_com_W + chrono::Vcross(master_pred.w_W, s.rA_W);
        const chrono::ChVector3d vB_point_W = slave_pred.v_com_W + chrono::Vcross(slave_pred.w_W, s.rB_W);

        s.u_pred_W = vB_point_W - vA_point_W;
        s.u_n_pred = chrono::Vdot(s.n_W, s.u_pred_W);
        s.u_tau_pred_W = s.P_W * s.u_pred_W;

        s.active_age = active_age_[i];
        s.cluster_size = 1;

        candidates.push_back(s);
    }

    stats_.accepted_before_cap = candidates.size();

    // ==========================================
    // TOPOLOGY-AWARE CONTACT MANIFOLD SUBSUMPTION 
    // (Multi-Patch Manifold Extraction with Normal Consistency Filtering)
    // ==========================================
    if (candidates.size() > 0) {
        // Sort all candidates by penetration depth (deepest first)
        std::sort(candidates.begin(), candidates.end(), [](const ActiveContactSample& a, const ActiveContactSample& b) {
            return a.phi < b.phi;
        });

        std::vector<ActiveContactSample> active_patches;
        const double angle_threshold_deg = GetEnvDouble("SPCC_GEAR_CLUSTER_ANGLE_DEG", 30.0);
        const double angle_threshold_cos = std::cos(angle_threshold_deg * chrono::CH_PI / 180.0);

        std::vector<ActiveContactSample> filtered_candidates;
        filtered_candidates.reserve(candidates.size());
        const double separating_cutoff = GetEnvDouble("SPCC_GEAR_SEPARATING_CUTOFF", 1.0e-4);
        for (const auto& candidate : candidates) {
            if (candidate.u_n_pred <= separating_cutoff) {
                filtered_candidates.push_back(candidate);
            }
        }

        const double distance_threshold = GetEnvDouble("SPCC_GEAR_CLUSTER_RADIUS", 3.0e-4);
        for (const auto& candidate : filtered_candidates) {
            bool is_new_patch = true;

            for (auto& patch : active_patches) {
                double dist = (candidate.x_W - patch.x_W).Length();
                double normal_dot = chrono::Vdot(candidate.n_W, patch.n_W);
                if (dist < distance_threshold && normal_dot > angle_threshold_cos) {
                    is_new_patch = false;
                    const int prev_cluster_size = patch.cluster_size;
                    patch.cluster_size += 1;

                    if (GetEnvDouble("SPCC_GEAR_AVG_POINT", 0.0) > 0.5) {
                        const double inv = 1.0 / static_cast<double>(patch.cluster_size);
                        patch.x_W = (patch.x_W * static_cast<double>(prev_cluster_size) + candidate.x_W) * inv;
                        patch.x_master_M =
                            (patch.x_master_M * static_cast<double>(prev_cluster_size) + candidate.x_master_M) * inv;
                    }

                    patch.n_W = patch.n_W * static_cast<double>(prev_cluster_size) + candidate.n_W;
                    patch.n_W.Normalize();

                    // Keep the deepest point/Hessian in the patch while only averaging the normal.
                    if (candidate.phi < patch.phi) {
                        patch.sample_id = candidate.sample_id;
                        patch.xi_slave_S = candidate.xi_slave_S;
                        patch.x_W = candidate.x_W;
                        patch.x_master_M = candidate.x_master_M;
                        patch.phi = candidate.phi;
                        patch.grad_M = candidate.grad_M;
                        patch.hessian_M = candidate.hessian_M;
                        patch.hessian_W = candidate.hessian_W;
                        patch.rA_W = candidate.rA_W;
                        patch.rB_W = candidate.rB_W;
                        patch.mu = candidate.mu;
                        patch.u_pred_W = candidate.u_pred_W;
                        patch.active_age = candidate.active_age;
                    }
                    break;
                }
            }

            if (is_new_patch) {
                active_patches.push_back(candidate);
            }
        }

        for (auto& patch : active_patches) {
            if (GetEnvDouble("SPCC_GEAR_AVG_POINT", 0.0) > 0.5) {
                patch.rA_W = patch.x_W - master_pred.x_com_W;
                patch.rB_W = patch.x_W - slave_pred.x_com_W;
                const chrono::ChVector3d vA_point_W =
                    master_pred.v_com_W + chrono::Vcross(master_pred.w_W, patch.rA_W);
                const chrono::ChVector3d vB_point_W =
                    slave_pred.v_com_W + chrono::Vcross(slave_pred.w_W, patch.rB_W);
                patch.u_pred_W = vB_point_W - vA_point_W;
            }
            patch.P_W = BuildProjectorFromNormal(patch.n_W);
            patch.u_n_pred = chrono::Vdot(patch.n_W, patch.u_pred_W);
            patch.u_tau_pred_W = patch.P_W * patch.u_pred_W;
        }

        candidates = active_patches;
    }

    out_active.insert(out_active.end(), candidates.begin(), candidates.end());
    stats_.accepted_after_cap = out_active.size();
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform

