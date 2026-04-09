#include "platform/backend/spcc/CompressedContactPipeline.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace platform {
namespace backend {
namespace spcc {

namespace {

struct PatchAccumulator {
    std::vector<std::size_t> members;
    chrono::ChVector3d centroid_W;
    chrono::ChVector3d avg_normal_W;
};

struct ReducedSupportAggregate {
    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d n_W;
    chrono::ChVector3d v_rel_W;
    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
    double support_weight = 0.0;
    std::size_t dense_members = 0;
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

double ReductionWeight(const DenseContactPoint& point) {
    return std::max(ProxyLoad(point), std::max(1.0e-12, point.area_weight));
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

chrono::ChVector3d PatchReferenceCenter(const std::vector<DenseContactPoint>& dense_points,
                                        const PatchAccumulator& patch) {
    chrono::ChVector3d center_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto point_index : patch.members) {
        const auto& point = dense_points[point_index];
        const double weight = ProxyLoad(point);
        center_W += weight * point.x_W;
        weight_sum += weight;
    }
    if (weight_sum > 1.0e-12) {
        return center_W * (1.0 / weight_sum);
    }
    return patch.centroid_W;
}

double TangentialRadiusSquared(const chrono::ChVector3d& x_W,
                               const chrono::ChVector3d& origin_W,
                               const chrono::ChVector3d& n_W) {
    const chrono::ChVector3d rel = x_W - origin_W;
    const chrono::ChVector3d tangential = rel - chrono::Vdot(rel, n_W) * n_W;
    return tangential.Length2();
}

std::vector<std::size_t> SelectSupportIndices(const std::vector<DenseContactPoint>& dense_points,
                                              const PatchAccumulator& patch,
                                              const std::vector<ReducedContactPoint>& previous_contacts,
                                              double match_radius,
                                              int target_points) {
    std::vector<std::size_t> selected;
    if (patch.members.empty() || target_points <= 0) {
        return selected;
    }

    target_points = std::max(1, target_points);

    std::size_t max_load_index = patch.members.front();
    double max_load = ProxyLoad(dense_points[max_load_index]);
    for (const auto point_index : patch.members) {
        const double load = ProxyLoad(dense_points[point_index]);
        if (load > max_load) {
            max_load = load;
            max_load_index = point_index;
        }
    }
    selected.push_back(max_load_index);
    if (target_points == 1 || patch.members.size() == 1) {
        return selected;
    }

    const chrono::ChVector3d center_W = PatchReferenceCenter(dense_points, patch);
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildOrthonormalBasis(patch.avg_normal_W, t1_W, t2_W);

    if (match_radius > 0.0) {
        for (const auto& previous : previous_contacts) {
            if (chrono::Vdot(previous.n_W, patch.avg_normal_W) < 0.85) {
                continue;
            }

            std::size_t best_index = patch.members.front();
            double best_distance = std::numeric_limits<double>::infinity();
            for (const auto point_index : patch.members) {
                const double distance = (dense_points[point_index].x_W - previous.x_W).Length();
                if (distance < best_distance) {
                    best_distance = distance;
                    best_index = point_index;
                }
            }

            if (best_distance <= match_radius &&
                std::find(selected.begin(), selected.end(), best_index) == selected.end()) {
                selected.push_back(best_index);
                if (static_cast<int>(selected.size()) >= target_points) {
                    return selected;
                }
            }
        }
    }

    const int direction_count = std::max(8, target_points * 2);
    for (int direction_index = 0; direction_index < direction_count; ++direction_index) {
        const double theta =
            (2.0 * std::acos(-1.0) * static_cast<double>(direction_index)) / static_cast<double>(direction_count);
        const chrono::ChVector3d axis_W = std::cos(theta) * t1_W + std::sin(theta) * t2_W;

        std::size_t best_index = patch.members.front();
        double best_score = -std::numeric_limits<double>::infinity();
        for (const auto point_index : patch.members) {
            const auto& point = dense_points[point_index];
            const chrono::ChVector3d rel = point.x_W - center_W;
            const double projected = chrono::Vdot(rel, axis_W);
            const double radial = std::sqrt(TangentialRadiusSquared(point.x_W, center_W, patch.avg_normal_W));
            const double score = projected + 0.1 * radial + 1.0e-3 * ProxyLoad(point);
            if (score > best_score) {
                best_score = score;
                best_index = point_index;
            }
        }
        if (std::find(selected.begin(), selected.end(), best_index) == selected.end()) {
            selected.push_back(best_index);
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
                score += std::sqrt(
                    TangentialRadiusSquared(dense_points[point_index].x_W, dense_points[chosen_index].x_W,
                                            patch.avg_normal_W));
            }
            score += std::sqrt(TangentialRadiusSquared(dense_points[point_index].x_W, center_W, patch.avg_normal_W));
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

std::vector<ReducedSupportAggregate> BuildReducedSupportsForPatch(const std::vector<DenseContactPoint>& dense_points,
                                                                  const PatchAccumulator& patch,
                                                                  const std::vector<ReducedContactPoint>& previous_contacts,
                                                                  double match_radius,
                                                                  int target_points) {
    target_points = std::max(1, std::min(target_points, static_cast<int>(patch.members.size())));
    const auto support_indices = SelectSupportIndices(dense_points, patch, previous_contacts, match_radius,
                                                      target_points);
    if (support_indices.empty()) {
        return {};
    }

    std::vector<chrono::ChVector3d> centers_W;
    centers_W.reserve(support_indices.size());
    for (const auto dense_index : support_indices) {
        centers_W.push_back(dense_points[dense_index].x_W);
    }

    std::vector<std::size_t> assignments(patch.members.size(), 0);
    for (std::size_t local_index = 0; local_index < patch.members.size(); ++local_index) {
        const auto member_index = patch.members[local_index];
        const auto& dense = dense_points[member_index];

        std::size_t best_cluster = 0;
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t cluster = 0; cluster < centers_W.size(); ++cluster) {
            const double distance = (dense.x_W - centers_W[cluster]).Length2();
            if (distance < best_distance) {
                best_distance = distance;
                best_cluster = cluster;
            }
        }
        assignments[local_index] = best_cluster;
    }

    const chrono::ChVector3d patch_center_W = PatchReferenceCenter(dense_points, patch);
    std::vector<ReducedSupportAggregate> supports;
    supports.reserve(centers_W.size());
    for (std::size_t cluster = 0; cluster < centers_W.size(); ++cluster) {
        double sum_area = 0.0;
        double sum_load = 0.0;
        double sum_weight = 0.0;
        std::size_t count = 0;
        std::size_t representative_member_index = patch.members.front();
        double best_rep_score = -std::numeric_limits<double>::infinity();

        chrono::ChVector3d anchor_W = centers_W[cluster];
        bool anchored_to_history = false;
        if (match_radius > 0.0) {
            double best_anchor_distance = std::numeric_limits<double>::infinity();
            for (const auto& previous : previous_contacts) {
                const double distance = (previous.x_W - centers_W[cluster]).Length();
                if (distance < best_anchor_distance) {
                    best_anchor_distance = distance;
                    anchor_W = previous.x_W;
                }
            }
            anchored_to_history = best_anchor_distance <= match_radius;
        }

        chrono::ChVector3d rep_dir_W = centers_W[cluster] - patch_center_W;
        rep_dir_W -= chrono::Vdot(rep_dir_W, patch.avg_normal_W) * patch.avg_normal_W;
        rep_dir_W = SafeNormalized(rep_dir_W, chrono::ChVector3d(0.0, 0.0, 0.0));

        for (std::size_t local_index = 0; local_index < patch.members.size(); ++local_index) {
            if (assignments[local_index] != cluster) {
                continue;
            }

            const auto member_index = patch.members[local_index];
            const auto& dense = dense_points[member_index];
            const double weight = ReductionWeight(dense);

            sum_area += dense.area_weight;
            sum_load += ProxyLoad(dense);
            sum_weight += weight;
            ++count;

            const chrono::ChVector3d rel = dense.x_W - patch_center_W;
            chrono::ChVector3d tangential_rel = rel - chrono::Vdot(rel, patch.avg_normal_W) * patch.avg_normal_W;
            const double radial = tangential_rel.Length();
            double rep_score = 0.0;
            if (anchored_to_history) {
                rep_score = -(dense.x_W - anchor_W).Length() + 0.02 * radial + 1.0e-3 * ProxyLoad(dense);
            } else {
                const double aligned =
                    rep_dir_W.Length2() > 0.0 ? chrono::Vdot(tangential_rel, rep_dir_W) : radial;
                rep_score = aligned + 0.1 * radial + 1.0e-3 * ProxyLoad(dense);
            }
            if (rep_score > best_rep_score) {
                best_rep_score = rep_score;
                representative_member_index = member_index;
            }
        }

        if (!(sum_weight > 1.0e-12) || count == 0) {
            continue;
        }

        ReducedSupportAggregate support;
        const auto& representative = dense_points[representative_member_index];
        support.x_W = representative.x_W;
        support.x_master_M = representative.x_master_M;
        support.x_master_surface_W = representative.x_master_surface_W;
        support.n_W = representative.n_W;
        support.v_rel_W = representative.v_rel_W;
        support.phi = representative.phi;
        support.phi_eff = representative.phi_eff;
        support.area_weight = sum_area;
        support.support_weight =
            (support.phi_eff < -1.0e-12) ? (sum_load / std::max(1.0e-12, -support.phi_eff)) : sum_area;
        support.dense_members = count;
        supports.push_back(support);
    }

    return supports;
}

}  // namespace

void CompressedContactPipeline::Configure(const CompressedContactConfig& cfg) {
    cfg_ = cfg;
    if (!slave_surface_samples_.empty()) {
        dense_sample_bvh_.Build(slave_surface_samples_, cfg_.bvh_leaf_size);
    }
    previous_contacts_.clear();
}

void CompressedContactPipeline::SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples) {
    slave_surface_samples_ = std::move(samples);
    dense_sample_bvh_.Build(slave_surface_samples_, cfg_.bvh_leaf_size);
    previous_contacts_.clear();
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
    DenseContactCloudStats cloud_stats;
    DenseContactCloudBuilder::Build(cfg_, slave_surface_samples_, dense_sample_bvh_, master_state, slave_state, sdf,
                                    step_size, dense_points, &cloud_stats);

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
        std::vector<ReducedContactPoint> previous_patch_contacts;
        if (cfg_.warm_start_match_radius > 0.0) {
            const double patch_gate =
                std::max({cfg_.warm_start_match_radius, cfg_.patch_radius, cfg_.max_patch_diameter, 1.0e-6});
            for (const auto& previous : previous_contacts_) {
                if (chrono::Vdot(previous.n_W, patch.avg_normal_W) < cfg_.normal_cos_min) {
                    continue;
                }
                if ((previous.x_W - patch.centroid_W).Length() <= patch_gate) {
                    previous_patch_contacts.push_back(previous);
                }
            }
        }

        auto supports = BuildReducedSupportsForPatch(dense_points, patch, previous_patch_contacts,
                                                     cfg_.warm_start_match_radius,
                                                     cfg_.max_reduced_points_per_patch);

        for (std::size_t support_id = 0; support_id < supports.size(); ++support_id) {
            const auto& support = supports[support_id];
            ReducedContactPoint reduced;
            reduced.patch_id = patch_id;
            reduced.support_id = support_id;
            reduced.dense_members = support.dense_members;
            reduced.x_W = support.x_W;
            reduced.x_master_M = support.x_master_M;
            reduced.x_master_surface_W = support.x_master_surface_W;
            reduced.n_W = support.n_W;
            reduced.v_rel_W = support.v_rel_W;
            reduced.phi = support.phi;
            reduced.phi_eff = support.phi_eff;
            reduced.area_weight = support.area_weight;
            reduced.support_weight = support.support_weight;
            reduced.mu = mu_default;
            out_contacts.push_back(reduced);
        }
    }

    previous_contacts_ = out_contacts;

    if (!out_stats) {
        return;
    }

    out_stats->total_samples = cloud_stats.total_samples;
    out_stats->candidate_count = cloud_stats.candidate_samples;
    out_stats->dense_count = dense_points.size();
    out_stats->reduced_count = out_contacts.size();
    out_stats->patch_count = patches.size();
    out_stats->bvh_nodes_visited = cloud_stats.bvh.nodes_visited;
    out_stats->bvh_nodes_pruned_obb = cloud_stats.bvh.nodes_pruned_obb;
    out_stats->bvh_nodes_pruned_sdf = cloud_stats.bvh.nodes_pruned_sdf;
    out_stats->bvh_leaf_samples_tested = cloud_stats.bvh.leaf_samples_tested;

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
    const chrono::ChVector3d reduced_moment = ProxyMoment(reduced_proxy, dense_cop);

    const double dense_force_norm = std::max(1.0e-12, dense_force.Length());
    const double dense_gap = std::max(1.0e-12, std::abs(std::min(0.0, dense_points.front().phi_eff)));
    double characteristic_radius = 0.0;
    for (const auto& dense : dense_points) {
        characteristic_radius = std::max(characteristic_radius, (dense.x_W - dense_cop).Length());
    }
    const double dense_moment_scale =
        std::max(1.0e-12, dense_force_norm * std::max(1.0e-6, characteristic_radius));

    out_stats->epsilon_F = (reduced_force - dense_force).Length() / dense_force_norm;
    out_stats->epsilon_M = (reduced_moment - dense_moment).Length() / dense_moment_scale;
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
