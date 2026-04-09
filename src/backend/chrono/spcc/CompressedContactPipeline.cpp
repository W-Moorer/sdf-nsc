#include "platform/backend/spcc/CompressedContactPipeline.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_set>

namespace platform {
namespace backend {
namespace spcc {

namespace {

struct DenseContactPoint {
    std::size_t sample_id = 0;
    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d n_W;
    chrono::ChVector3d v_rel_W;
    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
};

struct PatchAccumulator {
    std::vector<std::size_t> members;
    chrono::ChVector3d centroid_W;
    chrono::ChVector3d avg_normal_W;
};

chrono::ChVector3d SafeNormalized(const chrono::ChVector3d& v, const chrono::ChVector3d& fallback) {
    const double len = v.Length();
    if (!(len > 1.0e-12) || !std::isfinite(len)) {
        return fallback;
    }
    return v * (1.0 / len);
}

void BuildOrthonormalBasis(const chrono::ChVector3d& n_W,
                           chrono::ChVector3d& t1_W,
                           chrono::ChVector3d& t2_W) {
    const chrono::ChVector3d seed = (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0)
                                                               : chrono::ChVector3d(1.0, 0.0, 0.0);
    t1_W = SafeNormalized(chrono::Vcross(seed, n_W), chrono::ChVector3d(1.0, 0.0, 0.0));
    t2_W = SafeNormalized(chrono::Vcross(n_W, t1_W), chrono::ChVector3d(0.0, 1.0, 0.0));
}

double DistanceToLine(const chrono::ChVector3d& point_W,
                      const chrono::ChVector3d& origin_W,
                      const chrono::ChVector3d& dir_W) {
    const chrono::ChVector3d rel = point_W - origin_W;
    return chrono::Vcross(rel, dir_W).Length();
}

double ProxyLoad(const DenseContactPoint& point) {
    return std::max(0.0, -point.phi_eff) * std::max(0.0, point.area_weight);
}

chrono::ChVector3d ProxyForce(const std::vector<DenseContactPoint>& points) {
    chrono::ChVector3d force(0.0, 0.0, 0.0);
    for (const auto& point : points) {
        force += ProxyLoad(point) * point.n_W;
    }
    return force;
}

chrono::ChVector3d ProxyMoment(const std::vector<DenseContactPoint>& points, const chrono::ChVector3d& origin_W) {
    chrono::ChVector3d moment(0.0, 0.0, 0.0);
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        moment += chrono::Vcross(point.x_W - origin_W, load * point.n_W);
    }
    return moment;
}

chrono::ChVector3d ProxyCoP(const std::vector<DenseContactPoint>& points) {
    double weight_sum = 0.0;
    chrono::ChVector3d cop(0.0, 0.0, 0.0);
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        cop += load * point.x_W;
        weight_sum += load;
    }
    if (!(weight_sum > 1.0e-12)) {
        return cop;
    }
    return cop * (1.0 / weight_sum);
}

std::vector<std::size_t> SelectSupportIndices(const std::vector<DenseContactPoint>& dense_points,
                                              const PatchAccumulator& patch,
                                              int target_points) {
    std::vector<std::size_t> selected;
    if (patch.members.empty() || target_points <= 0) {
        return selected;
    }

    target_points = std::max(1, target_points);

    std::size_t deepest_index = patch.members.front();
    double deepest_phi = dense_points[deepest_index].phi_eff;
    for (const auto point_index : patch.members) {
        if (dense_points[point_index].phi_eff < deepest_phi) {
            deepest_phi = dense_points[point_index].phi_eff;
            deepest_index = point_index;
        }
    }
    selected.push_back(deepest_index);
    if (target_points == 1 || patch.members.size() == 1) {
        return selected;
    }

    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildOrthonormalBasis(patch.avg_normal_W, t1_W, t2_W);

    auto pick_extreme = [&](const chrono::ChVector3d& axis_W, bool maximize) {
        std::size_t best_index = patch.members.front();
        double best_value = maximize ? -std::numeric_limits<double>::infinity()
                                     : std::numeric_limits<double>::infinity();
        for (const auto point_index : patch.members) {
            const double value = chrono::Vdot(dense_points[point_index].x_W - patch.centroid_W, axis_W);
            if ((maximize && value > best_value) || (!maximize && value < best_value)) {
                best_value = value;
                best_index = point_index;
            }
        }
        return best_index;
    };

    std::vector<std::size_t> candidates = {
        pick_extreme(t1_W, false),
        pick_extreme(t1_W, true),
        pick_extreme(t2_W, false),
        pick_extreme(t2_W, true),
    };

    for (const auto candidate : candidates) {
        if (std::find(selected.begin(), selected.end(), candidate) == selected.end()) {
            selected.push_back(candidate);
            if (static_cast<int>(selected.size()) >= target_points) {
                return selected;
            }
        }
    }

    while (static_cast<int>(selected.size()) < target_points &&
           selected.size() < patch.members.size()) {
        std::size_t best_index = patch.members.front();
        double best_score = -1.0;
        for (const auto point_index : patch.members) {
            if (std::find(selected.begin(), selected.end(), point_index) != selected.end()) {
                continue;
            }
            double score = 0.0;
            for (const auto chosen_index : selected) {
                score += (dense_points[point_index].x_W - dense_points[chosen_index].x_W).Length();
            }
            score += DistanceToLine(dense_points[point_index].x_W, patch.centroid_W, patch.avg_normal_W);
            if (score > best_score) {
                best_score = score;
                best_index = point_index;
            }
        }
        if (std::find(selected.begin(), selected.end(), best_index) == selected.end()) {
            selected.push_back(best_index);
        } else {
            break;
        }
    }

    return selected;
}

}  // namespace

void CompressedContactPipeline::Configure(const CompressedContactConfig& cfg) {
    cfg_ = cfg;
}

void CompressedContactPipeline::SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples) {
    slave_surface_samples_ = std::move(samples);
}

void CompressedContactPipeline::BuildReducedContacts(const RigidBodyStateW& master_state,
                                                     const RigidBodyStateW& slave_state,
                                                     const FirstOrderSDF& sdf,
                                                     double mu_default,
                                                     double step_size,
                                                     std::vector<ReducedContactPoint>& out_contacts,
                                                     CompressionStats* out_stats) const {
    out_contacts.clear();

    std::vector<DenseContactPoint> dense_points;
    dense_points.reserve(slave_surface_samples_.size());

    for (std::size_t sample_id = 0; sample_id < slave_surface_samples_.size(); ++sample_id) {
        const auto& sample = slave_surface_samples_[sample_id];
        const chrono::ChVector3d x_W = slave_state.x_ref_W + slave_state.R_WRef * sample.xi_slave_S;
        const chrono::ChVector3d x_master_M = master_state.R_WRef.transpose() * (x_W - master_state.x_ref_W);

        double phi = 0.0;
        chrono::ChVector3d grad_M;
        if (!sdf.QueryPhiGradM(x_master_M, phi, grad_M)) {
            continue;
        }

        chrono::ChVector3d n_W = master_state.R_WRef * grad_M;
        n_W = SafeNormalized(n_W, chrono::ChVector3d(0.0, 1.0, 0.0));

        const chrono::ChVector3d rA_W = x_W - master_state.x_com_W;
        const chrono::ChVector3d rB_W = x_W - slave_state.x_com_W;
        const chrono::ChVector3d v_master_W = master_state.v_com_W + chrono::Vcross(master_state.w_W, rA_W);
        const chrono::ChVector3d v_slave_W = slave_state.v_com_W + chrono::Vcross(slave_state.w_W, rB_W);
        const chrono::ChVector3d v_rel_W = v_slave_W - v_master_W;

        double phi_eff = phi;
        if (cfg_.predictive_gap) {
            phi_eff += step_size * std::min(0.0, chrono::Vdot(n_W, v_rel_W));
        }

        if (!(phi_eff <= cfg_.delta_on || phi <= cfg_.delta_on)) {
            continue;
        }

        DenseContactPoint point;
        point.sample_id = sample_id;
        point.x_W = x_W;
        point.x_master_M = x_master_M;
        point.x_master_surface_W = x_W - phi * n_W;
        point.n_W = n_W;
        point.v_rel_W = v_rel_W;
        point.phi = phi;
        point.phi_eff = phi_eff;
        point.area_weight = sample.area_weight;
        dense_points.push_back(point);
    }

    if (cfg_.max_active_dense > 0 && static_cast<int>(dense_points.size()) > cfg_.max_active_dense) {
        std::stable_sort(dense_points.begin(), dense_points.end(), [](const DenseContactPoint& a,
                                                                      const DenseContactPoint& b) {
            return a.phi_eff < b.phi_eff;
        });
        dense_points.resize(static_cast<std::size_t>(cfg_.max_active_dense));
    }

    std::vector<PatchAccumulator> patches;
    for (std::size_t point_index = 0; point_index < dense_points.size(); ++point_index) {
        const auto& point = dense_points[point_index];
        bool assigned = false;
        for (auto& patch : patches) {
            const chrono::ChVector3d avg_normal = SafeNormalized(patch.avg_normal_W, point.n_W);
            if (chrono::Vdot(point.n_W, avg_normal) < cfg_.normal_cos_min) {
                continue;
            }

            const double distance = (point.x_W - patch.centroid_W).Length();
            if (cfg_.patch_radius > 0.0 && distance > cfg_.patch_radius) {
                continue;
            }

            if (cfg_.max_patch_diameter > 0.0) {
                bool too_far = false;
                for (const auto member_index : patch.members) {
                    if ((dense_points[member_index].x_W - point.x_W).Length() > cfg_.max_patch_diameter) {
                        too_far = true;
                        break;
                    }
                }
                if (too_far) {
                    continue;
                }
            }

            patch.members.push_back(point_index);
            const double inv_count = 1.0 / static_cast<double>(patch.members.size());
            patch.centroid_W += (point.x_W - patch.centroid_W) * inv_count;
            patch.avg_normal_W = SafeNormalized(patch.avg_normal_W + point.n_W, point.n_W);
            assigned = true;
            break;
        }

        if (!assigned) {
            PatchAccumulator patch;
            patch.members.push_back(point_index);
            patch.centroid_W = point.x_W;
            patch.avg_normal_W = point.n_W;
            patches.push_back(patch);
        }
    }

    out_contacts.reserve(dense_points.size());
    for (std::size_t patch_id = 0; patch_id < patches.size(); ++patch_id) {
        const auto& patch = patches[patch_id];
        const int target_points = std::max(1, std::min(cfg_.max_reduced_points_per_patch,
                                                       static_cast<int>(patch.members.size())));
        const auto support_indices = SelectSupportIndices(dense_points, patch, target_points);
        double patch_area = 0.0;
        double worst_gap = 0.0;
        for (const auto member_index : patch.members) {
            patch_area += dense_points[member_index].area_weight;
            worst_gap = std::min(worst_gap, dense_points[member_index].phi_eff);
        }
        const double per_support_weight =
            support_indices.empty() ? 0.0 : patch_area / static_cast<double>(support_indices.size());

        for (std::size_t support_id = 0; support_id < support_indices.size(); ++support_id) {
            const auto dense_index = support_indices[support_id];
            const auto& dense = dense_points[dense_index];

            ReducedContactPoint reduced;
            reduced.patch_id = patch_id;
            reduced.support_id = support_id;
            reduced.dense_members = patch.members.size();
            reduced.x_W = dense.x_W;
            reduced.x_master_M = dense.x_master_M;
            reduced.x_master_surface_W = dense.x_master_surface_W;
            reduced.n_W = dense.n_W;
            reduced.v_rel_W = dense.v_rel_W;
            reduced.phi = dense.phi;
            reduced.phi_eff = dense.phi_eff;
            reduced.area_weight = dense.area_weight;
            reduced.support_weight = per_support_weight;
            reduced.mu = mu_default;
            out_contacts.push_back(reduced);
        }
    }

    if (!out_stats) {
        return;
    }

    out_stats->dense_count = dense_points.size();
    out_stats->reduced_count = out_contacts.size();
    out_stats->patch_count = patches.size();

    if (dense_points.empty() || out_contacts.empty()) {
        out_stats->epsilon_F = 0.0;
        out_stats->epsilon_M = 0.0;
        out_stats->epsilon_CoP = 0.0;
        out_stats->epsilon_gap = 0.0;
        return;
    }

    std::vector<DenseContactPoint> reduced_proxy;
    reduced_proxy.reserve(out_contacts.size());
    for (const auto& reduced : out_contacts) {
        DenseContactPoint proxy;
        proxy.x_W = reduced.x_W;
        proxy.x_master_M = reduced.x_master_M;
        proxy.x_master_surface_W = reduced.x_master_surface_W;
        proxy.n_W = reduced.n_W;
        proxy.v_rel_W = reduced.v_rel_W;
        proxy.phi = reduced.phi;
        proxy.phi_eff = reduced.phi_eff;
        proxy.area_weight = reduced.support_weight;
        reduced_proxy.push_back(proxy);
    }

    const chrono::ChVector3d dense_force = ProxyForce(dense_points);
    const chrono::ChVector3d reduced_force = ProxyForce(reduced_proxy);
    const chrono::ChVector3d dense_cop = ProxyCoP(dense_points);
    const chrono::ChVector3d reduced_cop = ProxyCoP(reduced_proxy);
    const chrono::ChVector3d dense_moment = ProxyMoment(dense_points, dense_cop);
    const chrono::ChVector3d reduced_moment = ProxyMoment(reduced_proxy, reduced_cop);

    const double dense_force_norm = std::max(1.0e-12, dense_force.Length());
    const double dense_moment_norm = std::max(1.0e-12, dense_moment.Length());
    const double dense_gap = std::max(1.0e-12, std::abs(std::min(0.0, dense_points.front().phi_eff)));

    out_stats->epsilon_F = (reduced_force - dense_force).Length() / dense_force_norm;
    out_stats->epsilon_M = (reduced_moment - dense_moment).Length() / dense_moment_norm;
    out_stats->epsilon_CoP = (reduced_cop - dense_cop).Length();

    double reduced_worst_gap = 0.0;
    for (const auto& reduced : out_contacts) {
        reduced_worst_gap = std::min(reduced_worst_gap, reduced.phi_eff);
    }
    double dense_worst_gap = 0.0;
    for (const auto& dense : dense_points) {
        dense_worst_gap = std::min(dense_worst_gap, dense.phi_eff);
    }
    out_stats->epsilon_gap = std::abs(reduced_worst_gap - dense_worst_gap) / dense_gap;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
