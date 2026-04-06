#include "platform/backend/spcc/ContactActivation.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>

namespace platform {
namespace backend {
namespace spcc {

namespace {

struct CellKey {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const CellKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& key) const {
        const std::size_t h1 = static_cast<std::size_t>(key.x) * 73856093u;
        const std::size_t h2 = static_cast<std::size_t>(key.y) * 19349663u;
        const std::size_t h3 = static_cast<std::size_t>(key.z) * 83492791u;
        return h1 ^ h2 ^ h3;
    }
};

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

int GetEnvInt(const char* name, int fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value) {
        return fallback;
    }
    return static_cast<int>(parsed);
}

std::string MakeScopedEnvName(const std::string& prefix, const char* suffix) {
    if (prefix.empty()) {
        return std::string(suffix);
    }
    return prefix + "_" + suffix;
}

double GetScopedEnvDouble(const std::string& prefix, const char* suffix, double fallback) {
    return GetEnvDouble(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

int GetScopedEnvInt(const std::string& prefix, const char* suffix, int fallback) {
    return GetEnvInt(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

CellKey MakeCellKey(const chrono::ChVector3d& p, double cell_size) {
    return CellKey{
        static_cast<int>(std::floor(p.x() / cell_size)),
        static_cast<int>(std::floor(p.y() / cell_size)),
        static_cast<int>(std::floor(p.z() / cell_size)),
    };
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

double MatrixFrobeniusNorm(const chrono::ChMatrix33<>& M) {
    double sum_sq = 0.0;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            const double v = M(r, c);
            sum_sq += v * v;
        }
    }
    return std::sqrt(sum_sq);
}

struct PatchLocalFitPoint {
    chrono::ChVector3d x_master_M;
    double phi = 0.0;
};

bool BuildOrthonormalTangentBasis(const chrono::ChVector3d& n_M,
                                  chrono::ChVector3d& t1_M,
                                  chrono::ChVector3d& t2_M) {
    chrono::ChVector3d seed = (std::abs(n_M.z()) < 0.9) ? chrono::ChVector3d(0, 0, 1) : chrono::ChVector3d(0, 1, 0);
    t1_M = chrono::Vcross(seed, n_M);
    if (!NormalizeOrReject(t1_M)) {
        seed = chrono::ChVector3d(1, 0, 0);
        t1_M = chrono::Vcross(seed, n_M);
        if (!NormalizeOrReject(t1_M)) {
            return false;
        }
    }
    t2_M = chrono::Vcross(n_M, t1_M);
    return NormalizeOrReject(t2_M);
}

bool SolveLinear3x3(double A[3][4], double x[3]) {
    for (int pivot = 0; pivot < 3; ++pivot) {
        int best = pivot;
        double best_abs = std::abs(A[pivot][pivot]);
        for (int row = pivot + 1; row < 3; ++row) {
            const double cand = std::abs(A[row][pivot]);
            if (cand > best_abs) {
                best = row;
                best_abs = cand;
            }
        }
        if (!(best_abs > 1.0e-14) || !std::isfinite(best_abs)) {
            return false;
        }
        if (best != pivot) {
            for (int col = pivot; col < 4; ++col) {
                std::swap(A[pivot][col], A[best][col]);
            }
        }
        const double inv_pivot = 1.0 / A[pivot][pivot];
        for (int col = pivot; col < 4; ++col) {
            A[pivot][col] *= inv_pivot;
        }
        for (int row = 0; row < 3; ++row) {
            if (row == pivot) {
                continue;
            }
            const double factor = A[row][pivot];
            if (std::abs(factor) <= 1.0e-16) {
                continue;
            }
            for (int col = pivot; col < 4; ++col) {
                A[row][col] -= factor * A[pivot][col];
            }
        }
    }
    x[0] = A[0][3];
    x[1] = A[1][3];
    x[2] = A[2][3];
    return IsFiniteScalar(x[0]) && IsFiniteScalar(x[1]) && IsFiniteScalar(x[2]);
}

bool EstimatePatchLocalFitTarget(const std::vector<PatchLocalFitPoint>& points,
                                 const chrono::ChVector3d& center_master_M,
                                 const chrono::ChVector3d& normal_master_M,
                                 double max_shift_ratio,
                                 double shift_blend,
                                 double* out_predicted_min_phi,
                                 chrono::ChVector3d& out_target_master_M) {
    if (points.size() < 3) {
        return false;
    }

    chrono::ChVector3d n_M = normal_master_M;
    if (!NormalizeOrReject(n_M)) {
        return false;
    }

    chrono::ChVector3d t1_M;
    chrono::ChVector3d t2_M;
    if (!BuildOrthonormalTangentBasis(n_M, t1_M, t2_M)) {
        return false;
    }

    double c00 = 0.0;
    double c01 = 0.0;
    double c11 = 0.0;
    for (const auto& p : points) {
        const chrono::ChVector3d d_M = p.x_master_M - center_master_M;
        const double u = chrono::Vdot(d_M, t1_M);
        const double v = chrono::Vdot(d_M, t2_M);
        c00 += u * u;
        c01 += u * v;
        c11 += v * v;
    }

    chrono::ChVector3d tangent_master_M = t1_M;
    if (std::abs(c01) > 1.0e-16 || std::abs(c00 - c11) > 1.0e-16) {
        const double trace = c00 + c11;
        const double disc = std::sqrt(std::max(0.0, (c00 - c11) * (c00 - c11) + 4.0 * c01 * c01));
        const double lambda = 0.5 * (trace + disc);
        double ex = c01;
        double ey = lambda - c00;
        if (std::abs(ex) <= 1.0e-16 && std::abs(ey) <= 1.0e-16) {
            ex = (c00 >= c11) ? 1.0 : 0.0;
            ey = (c00 >= c11) ? 0.0 : 1.0;
        }
        tangent_master_M = t1_M * ex + t2_M * ey;
        if (!NormalizeOrReject(tangent_master_M)) {
            tangent_master_M = t1_M;
        }
    }

    double s_abs_max = 0.0;
    double s2_sum = 0.0;
    double s3_sum = 0.0;
    double s4_sum = 0.0;
    double sy_sum = 0.0;
    double s2y_sum = 0.0;
    double s_sum = 0.0;
    const double count = static_cast<double>(points.size());

    for (const auto& p : points) {
        const double s = chrono::Vdot(p.x_master_M - center_master_M, tangent_master_M);
        const double s2 = s * s;
        s_abs_max = std::max(s_abs_max, std::abs(s));
        s_sum += s;
        s2_sum += s2;
        s3_sum += s2 * s;
        s4_sum += s2 * s2;
        sy_sum += s * p.phi;
        s2y_sum += s2 * p.phi;
    }

    if (!(s_abs_max > 1.0e-6) || !std::isfinite(s_abs_max)) {
        return false;
    }

    double A[3][4] = {
        {count, s_sum, s2_sum, 0.0},
        {s_sum, s2_sum, s3_sum, sy_sum},
        {s2_sum, s3_sum, s4_sum, s2y_sum},
    };

    double rhs0 = 0.0;
    for (const auto& p : points) {
        rhs0 += p.phi;
    }
    A[0][3] = rhs0;

    double coeffs[3] = {0.0, 0.0, 0.0};
    if (!SolveLinear3x3(A, coeffs)) {
        return false;
    }

    const double c2 = coeffs[2];
    if (!(c2 > 1.0e-6) || !std::isfinite(c2)) {
        return false;
    }

    const double s_star_raw = -coeffs[1] / (2.0 * c2);
    if (!std::isfinite(s_star_raw)) {
        return false;
    }

    const double ratio = std::clamp(max_shift_ratio, 0.0, 1.0);
    const double max_shift = std::max(1.0e-6, ratio * s_abs_max);
    const double blend = std::clamp(shift_blend, 0.0, 1.0);
    const double s_star = std::clamp(s_star_raw, -max_shift, max_shift) * blend;
    if (!(std::abs(s_star) > 1.0e-6) || !std::isfinite(s_star)) {
        return false;
    }

    const double phi_min = coeffs[0] + coeffs[1] * s_star + coeffs[2] * s_star * s_star;
    if (out_predicted_min_phi != nullptr) {
        *out_predicted_min_phi = phi_min;
    }
    out_target_master_M = center_master_M + tangent_master_M * s_star;
    return IsFiniteVec(out_target_master_M);
}

}  // namespace

void ContactActivation::SetEnvPrefix(const std::string& env_prefix) {
    env_prefix_ = env_prefix;
    local_neighbors_.clear();
    local_neighbor_radius_ = -1.0;
    build_step_ = 0;
}

void ContactActivation::SetPolicy(const ContactRegimeConfig& policy) {
    regime_ = policy.regime;
    patch_geometry_mode_ = policy.patch_geometry_mode;
    full_scan_period_default_ = std::max(1, policy.activation.full_scan_period);
    local_scan_radius_default_ = policy.activation.local_scan_radius;
    cluster_angle_deg_default_ = policy.activation.cluster_angle_deg;
    separating_cutoff_default_ = policy.activation.separating_cutoff;
    cluster_radius_default_ = policy.activation.cluster_radius;
    avg_point_default_ = policy.activation.avg_point;
    persistent_match_radius_default_ = policy.activation.persistent_match_radius;
    persistent_normal_cos_min_default_ = policy.activation.persistent_normal_cos_min;
    persistent_blend_alpha_default_ = policy.activation.persistent_blend_alpha;
    persistent_path_samples_default_ = std::max(1, policy.activation.persistent_path_samples);
    coverage_spacing_radius_default_ = policy.activation.coverage_spacing_radius;
    onset_refine_steps_default_ = std::max(0, policy.activation.onset_refine_steps);
    onset_refine_path_samples_default_ = std::max(1, policy.activation.onset_refine_path_samples);
    onset_refine_backtrack_scale_default_ = policy.activation.onset_refine_backtrack_scale;
    local_fit_onset_steps_default_ = std::max(0, policy.activation.local_fit_onset_steps);
    local_fit_min_cluster_size_default_ = std::max(3, policy.activation.local_fit_min_cluster_size);
    local_fit_path_samples_default_ = std::max(1, policy.activation.local_fit_path_samples);
    local_fit_max_shift_ratio_default_ = policy.activation.local_fit_max_shift_ratio;
    local_fit_blend_default_ = policy.activation.local_fit_blend;
    local_fit_reject_positive_phi_default_ = policy.activation.local_fit_reject_positive_phi;
    curvature_gate_enabled_ = policy.curvature.enabled;
    curvature_tangential_only_default_ = policy.curvature.tangential_only;
    normal_alignment_cos_min_default_ = policy.curvature.normal_alignment_cos_min;
    max_hessian_frobenius_default_ = policy.curvature.max_hessian_frobenius;
    max_curvature_term_abs_default_ = policy.curvature.max_curvature_term_abs;
    max_curvature_term_ratio_default_ = policy.curvature.max_curvature_term_ratio;
    curvature_gap_floor_default_ = policy.curvature.gap_floor;
    curvature_ramp_steps_default_ = std::max(0, policy.curvature.ramp_steps);
    local_neighbors_.clear();
    local_neighbor_radius_ = -1.0;
    build_step_ = 0;
    persistent_patches_.clear();
    next_manifold_id_ = 1;
    Configure(policy.activation.delta_on,
              policy.activation.delta_off,
              policy.activation.hold_steps,
              static_cast<std::size_t>(std::max(0, policy.activation.max_active_keep)));
}

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
    persistent_patches_.clear();
    build_step_ = 0;
    next_manifold_id_ = 1;
    stats_ = Stats{};
}

void ContactActivation::EnsureLocalNeighborGraph(const std::vector<chrono::ChVector3d>& local_samples_S, double radius) {
    if (!(radius > 0.0) || !std::isfinite(radius)) {
        local_neighbors_.clear();
        local_neighbor_radius_ = -1.0;
        return;
    }

    if (local_neighbors_.size() == local_samples_S.size() &&
        std::abs(local_neighbor_radius_ - radius) <= 1.0e-12) {
        return;
    }

    local_neighbors_.assign(local_samples_S.size(), {});
    local_neighbor_radius_ = radius;

    std::unordered_map<CellKey, std::vector<std::uint32_t>, CellKeyHash> grid;
    grid.reserve(local_samples_S.size());

    for (std::uint32_t i = 0; i < local_samples_S.size(); ++i) {
        const auto& p = local_samples_S[i];
        if (!IsFiniteVec(p)) {
            continue;
        }
        grid[MakeCellKey(p, radius)].push_back(i);
    }

    const double radius_sq = radius * radius;
    for (std::uint32_t i = 0; i < local_samples_S.size(); ++i) {
        const auto& p = local_samples_S[i];
        if (!IsFiniteVec(p)) {
            continue;
        }

        const CellKey key = MakeCellKey(p, radius);
        auto& neighbors = local_neighbors_[i];
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const CellKey probe{key.x + dx, key.y + dy, key.z + dz};
                    const auto it = grid.find(probe);
                    if (it == grid.end()) {
                        continue;
                    }
                    for (std::uint32_t j : it->second) {
                        const auto d = local_samples_S[j] - p;
                        if (chrono::Vdot(d, d) <= radius_sq) {
                            neighbors.push_back(j);
                        }
                    }
                }
            }
        }
    }
}

void ContactActivation::RefreshPatchKinematics(const RigidBodyStateW& master_pred,
                                               const RigidBodyStateW& slave_pred,
                                               ActiveContactSample& patch) const {
    patch.rA_W = patch.x_W - master_pred.x_com_W;
    patch.rB_W = patch.x_W - slave_pred.x_com_W;

    const chrono::ChVector3d vA_point_W = master_pred.v_com_W + chrono::Vcross(master_pred.w_W, patch.rA_W);
    const chrono::ChVector3d vB_point_W = slave_pred.v_com_W + chrono::Vcross(slave_pred.w_W, patch.rB_W);
    patch.u_pred_W = vB_point_W - vA_point_W;

    patch.P_W = BuildProjectorFromNormal(patch.n_W);
    patch.u_n_pred = chrono::Vdot(patch.n_W, patch.u_pred_W);
    patch.u_tau_pred_W = patch.P_W * patch.u_pred_W;
}

bool ContactActivation::QueryPatchGeometry(const RigidBodyStateW& master_pred,
                                           const SDFField& sdf,
                                           bool need_hessian,
                                           const chrono::ChVector3d& query_master_M,
                                           const chrono::ChVector3d& reference_normal_W,
                                           ActiveContactSample& patch) const {
    patch.curvature_gate = 1.0;
    patch.curvature_tangential_only = curvature_gate_enabled_ && curvature_tangential_only_default_;
    patch.curvature_term_abs_max = curvature_gate_enabled_ ? max_curvature_term_abs_default_ : 0.0;
    patch.curvature_term_ratio_max = curvature_gate_enabled_ ? max_curvature_term_ratio_default_ : 0.0;
    patch.curvature_gap_floor = curvature_gate_enabled_ ? curvature_gap_floor_default_ : 0.0;
    patch.queried_normal_alignment = 1.0;
    patch.hessian_frobenius = 0.0;
    patch.hessian_M.setZero();
    patch.hessian_W.setZero();

    double phi_q = patch.phi;
    chrono::ChVector3d grad_q = patch.grad_M;
    if (!sdf.QueryPhiGradM(query_master_M, phi_q, grad_q)) {
        patch.curvature_gate = 0.0;
        return false;
    }

    chrono::ChVector3d n_q_W = master_pred.R_WRef * grad_q;
    if (!NormalizeOrReject(n_q_W)) {
        patch.curvature_gate = 0.0;
        return false;
    }

    patch.x_master_M = query_master_M;
    patch.phi = phi_q;
    patch.grad_M = grad_q;
    patch.n_W = n_q_W;

    if (curvature_gate_enabled_) {
        patch.queried_normal_alignment = chrono::Vdot(reference_normal_W, n_q_W);
        if (!std::isfinite(patch.queried_normal_alignment)) {
            patch.queried_normal_alignment = -1.0;
            patch.curvature_gate = 0.0;
        } else if (normal_alignment_cos_min_default_ > -1.0 &&
                   patch.queried_normal_alignment < normal_alignment_cos_min_default_) {
            patch.curvature_gate = 0.0;
        }
    }

    if (need_hessian) {
        if (patch.curvature_gate > 0.0 || !curvature_gate_enabled_) {
            double phi_h = patch.phi;
            chrono::ChVector3d grad_h = patch.grad_M;
            chrono::ChMatrix33<> hessian_h(0);
            if (sdf.QueryPhiGradHessianM(query_master_M, phi_h, grad_h, hessian_h)) {
                chrono::ChVector3d n_h_W = master_pred.R_WRef * grad_h;
                if (NormalizeOrReject(n_h_W)) {
                    patch.phi = phi_h;
                    patch.grad_M = grad_h;
                    patch.n_W = n_h_W;
                    if (curvature_gate_enabled_) {
                        patch.queried_normal_alignment = chrono::Vdot(reference_normal_W, n_h_W);
                        if (!std::isfinite(patch.queried_normal_alignment)) {
                            patch.queried_normal_alignment = -1.0;
                            patch.curvature_gate = 0.0;
                        } else if (normal_alignment_cos_min_default_ > -1.0 &&
                                   patch.queried_normal_alignment < normal_alignment_cos_min_default_) {
                            patch.curvature_gate = 0.0;
                        }
                    }
                    if (patch.curvature_gate > 0.0 || !curvature_gate_enabled_) {
                        patch.hessian_M = hessian_h;
                        patch.hessian_W = master_pred.R_WRef * hessian_h * master_pred.R_WRef.transpose();
                    }
                } else {
                    patch.curvature_gate = 0.0;
                }
            } else {
                patch.curvature_gate = 0.0;
            }
        }

        patch.hessian_frobenius = MatrixFrobeniusNorm(patch.hessian_W);
        if (!std::isfinite(patch.hessian_frobenius)) {
            patch.curvature_gate = 0.0;
            patch.hessian_frobenius = 0.0;
            patch.hessian_M.setZero();
            patch.hessian_W.setZero();
        } else if (curvature_gate_enabled_ && max_hessian_frobenius_default_ > 0.0 &&
                   patch.hessian_frobenius > max_hessian_frobenius_default_) {
            patch.curvature_gate *= (max_hessian_frobenius_default_ / patch.hessian_frobenius);
        }
    }

    if (!std::isfinite(patch.curvature_gate)) {
        patch.curvature_gate = 0.0;
    }
    patch.curvature_gate = std::clamp(patch.curvature_gate, 0.0, 1.0);
    patch.P_W = BuildProjectorFromNormal(patch.n_W);
    return true;
}

bool ContactActivation::QuerySweptPatchGeometry(const RigidBodyStateW& master_pred,
                                                const SDFField& sdf,
                                                bool need_hessian,
                                                const chrono::ChVector3d& start_master_M,
                                                const chrono::ChVector3d& end_master_M,
                                                const chrono::ChVector3d& reference_normal_W,
                                                int path_samples,
                                                ActiveContactSample& patch) const {
    const int sample_count = std::max(1, path_samples);
    double phi_sum = 0.0;
    chrono::ChVector3d grad_sum(0, 0, 0);
    int valid = 0;

    for (int i = 0; i < sample_count; ++i) {
        const double alpha = (sample_count <= 1) ? 1.0 : static_cast<double>(i) / static_cast<double>(sample_count - 1);
        const chrono::ChVector3d query_master_M = start_master_M * (1.0 - alpha) + end_master_M * alpha;
        double phi_q = 0.0;
        chrono::ChVector3d grad_q;
        if (!sdf.QueryPhiGradM(query_master_M, phi_q, grad_q)) {
            continue;
        }
        if (!IsFiniteScalar(phi_q) || !IsFiniteVec(grad_q)) {
            continue;
        }
        phi_sum += phi_q;
        grad_sum += grad_q;
        ++valid;
    }

    if (valid <= 0) {
        patch.curvature_gate = 0.0;
        return false;
    }

    const chrono::ChVector3d query_master_M = start_master_M * 0.5 + end_master_M * 0.5;
    patch.curvature_gate = 1.0;
    patch.curvature_tangential_only = curvature_gate_enabled_ && curvature_tangential_only_default_;
    patch.curvature_term_abs_max = curvature_gate_enabled_ ? max_curvature_term_abs_default_ : 0.0;
    patch.curvature_term_ratio_max = curvature_gate_enabled_ ? max_curvature_term_ratio_default_ : 0.0;
    patch.curvature_gap_floor = curvature_gate_enabled_ ? curvature_gap_floor_default_ : 0.0;
    patch.queried_normal_alignment = 1.0;
    patch.hessian_frobenius = 0.0;
    patch.hessian_M.setZero();
    patch.hessian_W.setZero();

    patch.x_master_M = query_master_M;
    patch.phi = phi_sum / static_cast<double>(valid);
    patch.grad_M = grad_sum * (1.0 / static_cast<double>(valid));

    chrono::ChVector3d n_avg_W = master_pred.R_WRef * patch.grad_M;
    if (!NormalizeOrReject(n_avg_W)) {
        patch.curvature_gate = 0.0;
        return false;
    }
    patch.n_W = n_avg_W;

    if (curvature_gate_enabled_) {
        patch.queried_normal_alignment = chrono::Vdot(reference_normal_W, patch.n_W);
        if (!std::isfinite(patch.queried_normal_alignment)) {
            patch.queried_normal_alignment = -1.0;
            patch.curvature_gate = 0.0;
        } else if (normal_alignment_cos_min_default_ > -1.0 &&
                   patch.queried_normal_alignment < normal_alignment_cos_min_default_) {
            patch.curvature_gate = 0.0;
        }
    }

    if (need_hessian) {
        if (patch.curvature_gate > 0.0 || !curvature_gate_enabled_) {
            double phi_h = patch.phi;
            chrono::ChVector3d grad_h = patch.grad_M;
            chrono::ChMatrix33<> hessian_h(0);
            if (sdf.QueryPhiGradHessianM(query_master_M, phi_h, grad_h, hessian_h)) {
                chrono::ChVector3d n_h_W = master_pred.R_WRef * grad_h;
                if (NormalizeOrReject(n_h_W)) {
                    if (curvature_gate_enabled_) {
                        const double alignment = chrono::Vdot(reference_normal_W, n_h_W);
                        if (!std::isfinite(alignment)) {
                            patch.curvature_gate = 0.0;
                        } else {
                            patch.queried_normal_alignment = alignment;
                            if (normal_alignment_cos_min_default_ > -1.0 && alignment < normal_alignment_cos_min_default_) {
                                patch.curvature_gate = 0.0;
                            }
                        }
                    }
                    patch.hessian_M = hessian_h;
                    patch.hessian_W = master_pred.R_WRef * hessian_h * master_pred.R_WRef.transpose();
                } else {
                    patch.curvature_gate = 0.0;
                }
            } else {
                patch.curvature_gate = 0.0;
            }
        }

        patch.hessian_frobenius = MatrixFrobeniusNorm(patch.hessian_W);
        if (!std::isfinite(patch.hessian_frobenius)) {
            patch.curvature_gate = 0.0;
            patch.hessian_frobenius = 0.0;
            patch.hessian_M.setZero();
            patch.hessian_W.setZero();
        } else if (curvature_gate_enabled_ && max_hessian_frobenius_default_ > 0.0 &&
                   patch.hessian_frobenius > max_hessian_frobenius_default_) {
            patch.curvature_gate *= (max_hessian_frobenius_default_ / patch.hessian_frobenius);
        }
    }

    if (!std::isfinite(patch.curvature_gate)) {
        patch.curvature_gate = 0.0;
    }
    patch.curvature_gate = std::clamp(patch.curvature_gate, 0.0, 1.0);
    patch.P_W = BuildProjectorFromNormal(patch.n_W);
    return true;
}

void ContactActivation::ApplySlidingPersistentManifold(const RigidBodyStateW& master_pred,
                                                       const RigidBodyStateW& slave_pred,
                                                       const SDFField& sdf,
                                                       bool need_hessian,
                                                       std::vector<ActiveContactSample>& patches) {
    if (patches.empty()) {
        return;
    }

    const double match_radius =
        GetScopedEnvDouble(env_prefix_, "PERSISTENT_MATCH_RADIUS", persistent_match_radius_default_);
    const double normal_cos_min =
        GetScopedEnvDouble(env_prefix_, "PERSISTENT_NORMAL_COS_MIN", persistent_normal_cos_min_default_);
    const double blend_alpha =
        std::clamp(GetScopedEnvDouble(env_prefix_, "PERSISTENT_BLEND_ALPHA", persistent_blend_alpha_default_), 0.0,
                   1.0);
    const int path_samples =
        std::max(1, GetScopedEnvInt(env_prefix_, "PERSISTENT_PATH_SAMPLES", persistent_path_samples_default_));

    std::vector<uint8_t> prev_used(persistent_patches_.size(), 0);
    for (auto& patch : patches) {
        patch.manifold_id = 0;
        patch.manifold_matched = false;

        int best_prev = -1;
        double best_score = std::numeric_limits<double>::infinity();
        if (match_radius > 0.0 && std::isfinite(match_radius)) {
            for (std::size_t i = 0; i < persistent_patches_.size(); ++i) {
                if (prev_used[i] != 0) {
                    continue;
                }
                const auto& prev = persistent_patches_[i];
                const double normal_dot = chrono::Vdot(prev.n_W, patch.n_W);
                if (normal_cos_min > -1.0 && normal_dot < normal_cos_min) {
                    continue;
                }
                const double dist = (prev.x_W - patch.x_W).Length();
                if (!std::isfinite(dist) || dist > match_radius) {
                    continue;
                }
                const double score = dist + 0.25 * match_radius * (1.0 - std::clamp(normal_dot, -1.0, 1.0));
                if (score < best_score) {
                    best_score = score;
                    best_prev = static_cast<int>(i);
                }
            }
        }

        if (best_prev >= 0) {
            prev_used[static_cast<std::size_t>(best_prev)] = 1;
            const auto& prev = persistent_patches_[static_cast<std::size_t>(best_prev)];
            patch.manifold_id = prev.manifold_id;
            patch.manifold_matched = true;
            patch.active_age = std::max(patch.active_age, prev.age + 1);

            if (blend_alpha < 1.0) {
                const double keep = 1.0 - blend_alpha;
                const chrono::ChVector3d blended_x_W = prev.x_W * keep + patch.x_W * blend_alpha;
                const chrono::ChVector3d current_x_master_M = patch.x_master_M;
                const chrono::ChVector3d blended_x_master_M =
                    prev.x_master_M * keep + current_x_master_M * blend_alpha;
                chrono::ChVector3d blended_normal_W = prev.n_W * keep + patch.n_W * blend_alpha;
                if (!NormalizeOrReject(blended_normal_W)) {
                    blended_normal_W = patch.n_W;
                }

                patch.x_W = blended_x_W;
                const bool query_ok =
                    (path_samples > 1)
                        ? QuerySweptPatchGeometry(master_pred, sdf, need_hessian, prev.x_master_M, current_x_master_M,
                                                  blended_normal_W, path_samples, patch)
                        : QueryPatchGeometry(master_pred, sdf, need_hessian, blended_x_master_M, blended_normal_W,
                                             patch);
                if (!query_ok) {
                    patch.x_master_M = blended_x_master_M;
                    patch.n_W = blended_normal_W;
                    patch.P_W = BuildProjectorFromNormal(patch.n_W);
                }
                RefreshPatchKinematics(master_pred, slave_pred, patch);
            }
        } else {
            patch.manifold_id = next_manifold_id_++;
        }
    }
}

void ContactActivation::SelectSlidingCoveragePatches(std::vector<ActiveContactSample>& patches) const {
    if (max_active_keep_ == 0 || patches.size() <= max_active_keep_) {
        return;
    }

    const double spacing_radius =
        GetScopedEnvDouble(env_prefix_, "COVERAGE_SPACING_RADIUS", coverage_spacing_radius_default_);

    std::vector<std::size_t> selected;
    selected.reserve(max_active_keep_);
    std::vector<uint8_t> used(patches.size(), 0);

    auto deepest_it = std::min_element(patches.begin(), patches.end(), [](const auto& a, const auto& b) {
        return a.phi < b.phi;
    });
    const std::size_t deepest_idx = static_cast<std::size_t>(std::distance(patches.begin(), deepest_it));
    selected.push_back(deepest_idx);
    used[deepest_idx] = 1;

    while (selected.size() < max_active_keep_) {
        int best_idx = -1;
        double best_score = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < patches.size(); ++i) {
            if (used[i] != 0) {
                continue;
            }

            double min_dist = 0.0;
            bool has_selected = false;
            for (std::size_t selected_idx : selected) {
                const double d = (patches[i].x_W - patches[selected_idx].x_W).Length();
                if (!has_selected || d < min_dist) {
                    min_dist = d;
                    has_selected = true;
                }
            }

            const double depth_score = -patches[i].phi;
            const double spacing_score =
                (spacing_radius > 0.0 && std::isfinite(spacing_radius)) ? std::min(1.0, min_dist / spacing_radius) : 0.0;
            const double persistence_bonus = patches[i].manifold_matched ? 0.2 : 0.0;
            const double age_bonus = 0.02 * std::min(patches[i].active_age, 10);
            const double score = depth_score + 5.0e-4 * spacing_score + 2.5e-4 * persistence_bonus + 1.0e-5 * age_bonus;

            if (score > best_score) {
                best_score = score;
                best_idx = static_cast<int>(i);
            }
        }

        if (best_idx < 0) {
            break;
        }
        used[static_cast<std::size_t>(best_idx)] = 1;
        selected.push_back(static_cast<std::size_t>(best_idx));
    }

    std::sort(selected.begin(), selected.end());
    std::vector<ActiveContactSample> kept;
    kept.reserve(selected.size());
    for (std::size_t idx : selected) {
        kept.push_back(patches[idx]);
    }
    std::sort(kept.begin(), kept.end(), [](const auto& a, const auto& b) {
        return a.phi < b.phi;
    });
    patches.swap(kept);
}

void ContactActivation::BuildActiveSet(const RigidBodyStateW& master_pred,
                                       const RigidBodyStateW& slave_pred,
                                       const SDFField& sdf,
                                       const std::vector<chrono::ChVector3d>& local_samples_S,
                                       const std::vector<chrono::ChVector3d>& world_samples_W,
                                       double mu_default,
                                       double step_size,
                                       bool need_hessian,
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
    const int full_scan_period = std::max(1, GetScopedEnvInt(env_prefix_, "FULL_SCAN_PERIOD", full_scan_period_default_));
    const double local_scan_radius =
        GetScopedEnvDouble(env_prefix_, "LOCAL_SCAN_RADIUS", local_scan_radius_default_);

    std::vector<std::size_t> query_indices;
    bool full_scan = true;
    if (n_samples > 0 && full_scan_period > 1 && local_scan_radius > 0.0 && std::isfinite(local_scan_radius)) {
        EnsureLocalNeighborGraph(local_samples_S, local_scan_radius);
        full_scan = ((build_step_ % static_cast<std::size_t>(full_scan_period)) == 0);
        if (!full_scan) {
            std::vector<uint8_t> query_mask(n_samples, 0);
            std::size_t seed_count = 0;
            for (std::size_t i = 0; i < n_samples; ++i) {
                if (prev_active_[i] == 0 && hold_counter_[i] <= 0) {
                    continue;
                }
                seed_count += 1;
                for (std::uint32_t j : local_neighbors_[i]) {
                    query_mask[j] = 1;
                }
            }

            if (seed_count == 0) {
                full_scan = true;
            } else {
                query_indices.reserve(n_samples);
                for (std::size_t i = 0; i < n_samples; ++i) {
                    if (query_mask[i] != 0) {
                        query_indices.push_back(i);
                    }
                }
                if (query_indices.empty()) {
                    full_scan = true;
                }
            }
        }
    }

    if (full_scan) {
        query_indices.resize(n_samples);
        std::iota(query_indices.begin(), query_indices.end(), std::size_t{0});
    }
    build_step_ += 1;

    candidates.reserve(query_indices.size());

    for (std::size_t i : query_indices) {
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
        if (!sdf.QueryPhiGradM(x_master_M, phi, grad_M)) {
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
        s.x_deepest_master_M = x_master_M;
        s.phi = phi;
        s.grad_M = grad_M;
        s.hessian_M.setZero();
        s.n_W = n_W;
        s.hessian_W.setZero();
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
        std::vector<std::vector<PatchLocalFitPoint>> patch_member_points;
        const double angle_threshold_deg =
            GetScopedEnvDouble(env_prefix_, "CLUSTER_ANGLE_DEG", cluster_angle_deg_default_);
        const double angle_threshold_cos = std::cos(angle_threshold_deg * chrono::CH_PI / 180.0);

        std::vector<ActiveContactSample> filtered_candidates;
        filtered_candidates.reserve(candidates.size());
        const double separating_cutoff =
            GetScopedEnvDouble(env_prefix_, "SEPARATING_CUTOFF", separating_cutoff_default_);
        for (const auto& candidate : candidates) {
                    if (candidate.u_n_pred <= separating_cutoff) {
                        filtered_candidates.push_back(candidate);
                    }
        }

        const double distance_threshold =
            GetScopedEnvDouble(env_prefix_, "CLUSTER_RADIUS", cluster_radius_default_);
        const bool avg_point = (GetScopedEnvDouble(env_prefix_, "AVG_POINT", avg_point_default_ ? 1.0 : 0.0) > 0.5);
        for (const auto& candidate : filtered_candidates) {
            bool is_new_patch = true;

            for (auto& patch : active_patches) {
                const std::size_t patch_index = static_cast<std::size_t>(&patch - active_patches.data());
                double dist = (candidate.x_W - patch.x_W).Length();
                double normal_dot = chrono::Vdot(candidate.n_W, patch.n_W);
                if (dist < distance_threshold && normal_dot > angle_threshold_cos) {
                    is_new_patch = false;
                    const int prev_cluster_size = patch.cluster_size;
                    patch.cluster_size += 1;

                    if (avg_point) {
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
                        if (!(avg_point && patch_geometry_mode_ == PatchGeometryMode::RepresentativeQuery)) {
                            patch.x_W = candidate.x_W;
                            patch.x_master_M = candidate.x_master_M;
                        }
                        patch.x_deepest_master_M = candidate.x_deepest_master_M;
                        patch.phi = candidate.phi;
                        patch.grad_M = candidate.grad_M;
                        patch.rA_W = candidate.rA_W;
                        patch.rB_W = candidate.rB_W;
                        patch.mu = candidate.mu;
                        patch.u_pred_W = candidate.u_pred_W;
                        patch.active_age = candidate.active_age;
                    }
                    patch_member_points[patch_index].push_back(PatchLocalFitPoint{candidate.x_master_M, candidate.phi});
                    break;
                }
            }

            if (is_new_patch) {
                active_patches.push_back(candidate);
                patch_member_points.push_back({PatchLocalFitPoint{candidate.x_master_M, candidate.phi}});
            }
        }

        for (std::size_t patch_index = 0; patch_index < active_patches.size(); ++patch_index) {
            auto& patch = active_patches[patch_index];
            if (avg_point) {
                patch.rA_W = patch.x_W - master_pred.x_com_W;
                patch.rB_W = patch.x_W - slave_pred.x_com_W;
                const chrono::ChVector3d vA_point_W =
                    master_pred.v_com_W + chrono::Vcross(master_pred.w_W, patch.rA_W);
                const chrono::ChVector3d vB_point_W =
                    slave_pred.v_com_W + chrono::Vcross(slave_pred.w_W, patch.rB_W);
                patch.u_pred_W = vB_point_W - vA_point_W;
            }
            patch.curvature_gate = 1.0;
            patch.curvature_tangential_only = curvature_gate_enabled_ && curvature_tangential_only_default_;
            patch.curvature_term_abs_max = curvature_gate_enabled_ ? max_curvature_term_abs_default_ : 0.0;
            patch.curvature_term_ratio_max = curvature_gate_enabled_ ? max_curvature_term_ratio_default_ : 0.0;
            patch.curvature_gap_floor = curvature_gate_enabled_ ? curvature_gap_floor_default_ : 0.0;
            patch.queried_normal_alignment = 1.0;
            patch.hessian_frobenius = 0.0;
            const bool representative_query = (patch_geometry_mode_ == PatchGeometryMode::RepresentativeQuery);
            if (representative_query) {
                const chrono::ChVector3d averaged_normal_W = patch.n_W;
                const chrono::ChVector3d query_master_M = patch.x_master_M;
                double phi_q = patch.phi;
                chrono::ChVector3d grad_q = patch.grad_M;
                const bool query_ok = sdf.QueryPhiGradM(query_master_M, phi_q, grad_q);
                if (query_ok) {
                    chrono::ChVector3d n_q_W = master_pred.R_WRef * grad_q;
                    if (NormalizeOrReject(n_q_W)) {
                        patch.phi = phi_q;
                        patch.grad_M = grad_q;
                        patch.n_W = n_q_W;
                        if (curvature_gate_enabled_) {
                            patch.queried_normal_alignment = chrono::Vdot(averaged_normal_W, n_q_W);
                            if (!std::isfinite(patch.queried_normal_alignment)) {
                                patch.queried_normal_alignment = -1.0;
                                patch.curvature_gate = 0.0;
                            } else if (normal_alignment_cos_min_default_ > -1.0 &&
                                       patch.queried_normal_alignment < normal_alignment_cos_min_default_) {
                                patch.curvature_gate = 0.0;
                            }
                        }
                        if (need_hessian) {
                            if (patch.curvature_gate > 0.0) {
                                double phi_h = patch.phi;
                                chrono::ChVector3d grad_h = patch.grad_M;
                                chrono::ChMatrix33<> hessian_h(0);
                                if (sdf.QueryPhiGradHessianM(query_master_M, phi_h, grad_h, hessian_h)) {
                                    chrono::ChVector3d n_h_W = master_pred.R_WRef * grad_h;
                                    if (NormalizeOrReject(n_h_W)) {
                                        patch.phi = phi_h;
                                        patch.grad_M = grad_h;
                                        patch.n_W = n_h_W;
                                        if (curvature_gate_enabled_) {
                                            patch.queried_normal_alignment = chrono::Vdot(averaged_normal_W, n_h_W);
                                            if (!std::isfinite(patch.queried_normal_alignment)) {
                                                patch.queried_normal_alignment = -1.0;
                                                patch.curvature_gate = 0.0;
                                            } else if (normal_alignment_cos_min_default_ > -1.0 &&
                                                       patch.queried_normal_alignment < normal_alignment_cos_min_default_) {
                                                patch.curvature_gate = 0.0;
                                            }
                                        }
                                        if (patch.curvature_gate > 0.0) {
                                            patch.hessian_M = hessian_h;
                                            patch.hessian_W = master_pred.R_WRef * hessian_h * master_pred.R_WRef.transpose();
                                        } else {
                                            patch.hessian_M.setZero();
                                            patch.hessian_W.setZero();
                                        }
                                    } else {
                                        patch.curvature_gate = 0.0;
                                        patch.hessian_M.setZero();
                                        patch.hessian_W.setZero();
                                    }
                                } else {
                                    patch.curvature_gate = 0.0;
                                    patch.hessian_M.setZero();
                                    patch.hessian_W.setZero();
                                }
                            } else {
                                patch.hessian_M.setZero();
                                patch.hessian_W.setZero();
                            }
                        } else {
                            patch.hessian_M.setZero();
                            patch.hessian_W.setZero();
                        }
                    } else {
                        patch.curvature_gate = 0.0;
                        patch.hessian_M.setZero();
                        patch.hessian_W.setZero();
                    }
                } else {
                    patch.curvature_gate = 0.0;
                    if (!need_hessian) {
                        patch.hessian_M.setZero();
                        patch.hessian_W.setZero();
                    } else {
                        patch.hessian_M.setZero();
                        patch.hessian_W.setZero();
                    }
                }
            } else if (need_hessian) {
                double phi_h = patch.phi;
                chrono::ChVector3d grad_h = patch.grad_M;
                chrono::ChMatrix33<> hessian_h(0);
                if (sdf.QueryPhiGradHessianM(patch.x_deepest_master_M, phi_h, grad_h, hessian_h)) {
                    patch.phi = phi_h;
                    patch.grad_M = grad_h;
                    patch.hessian_M = hessian_h;
                    patch.hessian_W = master_pred.R_WRef * hessian_h * master_pred.R_WRef.transpose();
                } else {
                    patch.hessian_M.setZero();
                    patch.hessian_W.setZero();
                }
            } else {
                patch.hessian_M.setZero();
                patch.hessian_W.setZero();
            }

            if (need_hessian) {
                patch.hessian_frobenius = MatrixFrobeniusNorm(patch.hessian_W);
                if (!std::isfinite(patch.hessian_frobenius)) {
                    patch.curvature_gate = 0.0;
                    patch.hessian_frobenius = 0.0;
                    patch.hessian_M.setZero();
                    patch.hessian_W.setZero();
                } else if (curvature_gate_enabled_ && max_hessian_frobenius_default_ > 0.0 &&
                           patch.hessian_frobenius > max_hessian_frobenius_default_) {
                    patch.curvature_gate *= (max_hessian_frobenius_default_ / patch.hessian_frobenius);
                }
            }
            if (!std::isfinite(patch.curvature_gate)) {
                patch.curvature_gate = 0.0;
            }
            if (regime_ == ContactRegimeType::SlidingPatch) {
                const int onset_steps =
                    std::max(0, GetScopedEnvInt(env_prefix_, "ONSET_REFINE_STEPS", onset_refine_steps_default_));
                const int onset_path_samples = std::max(
                    1, GetScopedEnvInt(env_prefix_, "ONSET_REFINE_PATH_SAMPLES", onset_refine_path_samples_default_));
                const double onset_backtrack_scale =
                    GetScopedEnvDouble(env_prefix_, "ONSET_REFINE_BACKTRACK_SCALE", onset_refine_backtrack_scale_default_);
                if (onset_steps > 0 && patch.active_age <= onset_steps && onset_path_samples > 1 &&
                    std::isfinite(onset_backtrack_scale) && onset_backtrack_scale > 0.0 &&
                    std::isfinite(step_size) && step_size > 0.0) {
                    const chrono::ChVector3d dx_master_M =
                        master_pred.R_WRef.transpose() * (patch.u_pred_W * (step_size * onset_backtrack_scale));
                    const chrono::ChVector3d start_master_M = patch.x_master_M - dx_master_M;
                    const chrono::ChVector3d end_master_M = patch.x_master_M;
                    QuerySweptPatchGeometry(master_pred, sdf, need_hessian, start_master_M, end_master_M, patch.n_W,
                                            onset_path_samples, patch);
                }

                const int local_fit_onset_steps = std::max(
                    0, GetScopedEnvInt(env_prefix_, "LOCAL_FIT_ONSET_STEPS", local_fit_onset_steps_default_));
                const int local_fit_min_cluster_size = std::max(
                    3, GetScopedEnvInt(env_prefix_, "LOCAL_FIT_MIN_CLUSTER_SIZE", local_fit_min_cluster_size_default_));
                const int local_fit_path_samples = std::max(
                    1, GetScopedEnvInt(env_prefix_, "LOCAL_FIT_PATH_SAMPLES", local_fit_path_samples_default_));
                const double local_fit_max_shift_ratio = GetScopedEnvDouble(
                    env_prefix_, "LOCAL_FIT_MAX_SHIFT_RATIO", local_fit_max_shift_ratio_default_);
                const double local_fit_blend =
                    GetScopedEnvDouble(env_prefix_, "LOCAL_FIT_BLEND", local_fit_blend_default_);
                const double local_fit_reject_positive_phi = GetScopedEnvDouble(
                    env_prefix_, "LOCAL_FIT_REJECT_POSITIVE_PHI", local_fit_reject_positive_phi_default_);

                if (need_hessian && local_fit_onset_steps > 0 && patch.active_age <= local_fit_onset_steps &&
                    patch.cluster_size >= local_fit_min_cluster_size && patch_index < patch_member_points.size()) {
                    chrono::ChVector3d normal_master_M = patch.grad_M;
                    if (NormalizeOrReject(normal_master_M)) {
                        chrono::ChVector3d fitted_master_M;
                        double predicted_min_phi = std::numeric_limits<double>::quiet_NaN();
                        if (EstimatePatchLocalFitTarget(patch_member_points[patch_index], patch.x_master_M,
                                                        normal_master_M, local_fit_max_shift_ratio, local_fit_blend,
                                                        &predicted_min_phi, fitted_master_M)) {
                            if (std::isfinite(local_fit_reject_positive_phi) && local_fit_reject_positive_phi >= 0.0 &&
                                std::isfinite(predicted_min_phi) && predicted_min_phi > local_fit_reject_positive_phi) {
                                patch.cluster_size = 0;
                                continue;
                            }
                            const chrono::ChVector3d current_master_M = patch.x_master_M;
                            if (local_fit_path_samples > 1) {
                                QuerySweptPatchGeometry(master_pred, sdf, need_hessian, current_master_M,
                                                        fitted_master_M, patch.n_W, local_fit_path_samples, patch);
                            } else {
                                QueryPatchGeometry(master_pred, sdf, need_hessian, fitted_master_M, patch.n_W, patch);
                            }
                        }
                    }
                }
            }
            if (curvature_gate_enabled_ && curvature_ramp_steps_default_ > 0) {
                const double ramp =
                    std::clamp(static_cast<double>(std::max(0, patch.active_age)) /
                                   static_cast<double>(curvature_ramp_steps_default_),
                               0.0, 1.0);
                patch.curvature_gate *= ramp;
            }
            patch.curvature_gate = std::clamp(patch.curvature_gate, 0.0, 1.0);

            patch.P_W = BuildProjectorFromNormal(patch.n_W);
            patch.u_n_pred = chrono::Vdot(patch.n_W, patch.u_pred_W);
            patch.u_tau_pred_W = patch.P_W * patch.u_pred_W;
        }

        active_patches.erase(std::remove_if(active_patches.begin(), active_patches.end(),
                                            [](const ActiveContactSample& patch) { return patch.cluster_size <= 0; }),
                             active_patches.end());
        candidates = active_patches;
    }

    if (regime_ == ContactRegimeType::SlidingPatch) {
        ApplySlidingPersistentManifold(master_pred, slave_pred, sdf, need_hessian, candidates);
    }

    if (max_active_keep_ > 0 && candidates.size() > max_active_keep_) {
        if (regime_ == ContactRegimeType::SlidingPatch) {
            SelectSlidingCoveragePatches(candidates);
        } else {
            std::nth_element(candidates.begin(),
                             candidates.begin() + static_cast<std::ptrdiff_t>(max_active_keep_),
                             candidates.end(),
                             [](const ActiveContactSample& a, const ActiveContactSample& b) {
                                 return a.phi < b.phi;
                             });
            candidates.resize(max_active_keep_);
            std::sort(candidates.begin(), candidates.end(),
                      [](const ActiveContactSample& a, const ActiveContactSample& b) { return a.phi < b.phi; });
        }
    }

    out_active.insert(out_active.end(), candidates.begin(), candidates.end());
    if (regime_ == ContactRegimeType::SlidingPatch) {
        persistent_patches_.clear();
        persistent_patches_.reserve(out_active.size());
        for (const auto& patch : out_active) {
            PersistentPatchState state;
            state.manifold_id = patch.manifold_id;
            state.x_W = patch.x_W;
            state.x_master_M = patch.x_master_M;
            state.n_W = patch.n_W;
            state.phi = patch.phi;
            state.age = std::max(1, patch.active_age);
            persistent_patches_.push_back(state);
        }
    } else {
        persistent_patches_.clear();
    }
    stats_.accepted_after_cap = out_active.size();
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform

