#include "platform/backend/spcc/CompressedContactPipeline.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_set>

#include "platform/backend/spcc/LocalWrenchAllocator.h"
#include "platform/backend/spcc/SubpatchRefiner.h"

namespace platform {
namespace backend {
namespace spcc {

namespace {

struct ReducedSupportAggregate {
    std::size_t patch_id = 0;
    std::size_t subpatch_id = 0;
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
    double coverage_radius = 0.0;
    std::size_t dense_members = 0;
};

struct SubpatchMatchCandidate {
    std::size_t current_index = 0;
    std::size_t previous_index = 0;
    double score = 0.0;
};

struct SupportMatchCandidate {
    std::size_t current_index = 0;
    std::size_t previous_index = 0;
    double score = 0.0;
};

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

chrono::ChVector3d SubpatchReferenceCenter(const std::vector<DenseContactPoint>& dense_points,
                                           const DenseSubpatch& subpatch) {
    chrono::ChVector3d center_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = ProxyLoad(point);
        center_W += weight * point.x_W;
        weight_sum += weight;
    }
    if (weight_sum > 1.0e-12) {
        return center_W * (1.0 / weight_sum);
    }
    return subpatch.centroid_W;
}

double TangentialRadiusSquared(const chrono::ChVector3d& x_W,
                               const chrono::ChVector3d& origin_W,
                               const chrono::ChVector3d& n_W) {
    const chrono::ChVector3d rel = x_W - origin_W;
    const chrono::ChVector3d tangential = rel - chrono::Vdot(rel, n_W) * n_W;
    return tangential.Length2();
}

chrono::ChVector3d ProjectToTangentUnit(const chrono::ChVector3d& v_W, const chrono::ChVector3d& n_W) {
    const chrono::ChVector3d tangent_W = v_W - chrono::Vdot(v_W, n_W) * n_W;
    const double tangent_len = tangent_W.Length();
    if (!(tangent_len > 1.0e-12)) {
        return chrono::ChVector3d(0.0, 0.0, 0.0);
    }
    return tangent_W * (1.0 / tangent_len);
}

chrono::ChVector3d ChooseStencilAxis(const ReducedSupportAggregate& support, const DenseSubpatch& subpatch) {
    const chrono::ChVector3d tangential_velocity_W = support.v_rel_W -
                                                     chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W) *
                                                         subpatch.avg_normal_W;
    const chrono::ChVector3d from_velocity =
        ProjectToTangentUnit(chrono::Vcross(subpatch.avg_normal_W, tangential_velocity_W), subpatch.avg_normal_W);
    if (from_velocity.Length2() > 0.0) {
        return from_velocity;
    }

    const chrono::ChVector3d radial_W = support.x_W - subpatch.centroid_W;
    const chrono::ChVector3d from_radial = ProjectToTangentUnit(radial_W, subpatch.avg_normal_W);
    if (from_radial.Length2() > 0.0) {
        return from_radial;
    }

    return subpatch.t1_W;
}

double SubpatchMatchGate(const CompressedContactConfig& cfg,
                         const DenseSubpatch& current,
                         const TemporalSubpatchState& previous) {
    const double base =
        std::max({cfg.warm_start_match_radius, cfg.max_subpatch_diameter, cfg.max_patch_diameter,
                  0.5 * current.diameter, 0.5 * previous.diameter, 1.0e-6});
    return 1.5 * base;
}

ReferenceWrench BlendTemporalReference(const ReferenceWrench& current,
                                       const TemporalSubpatchState* previous,
                                       double blend) {
    if (!previous) {
        return current;
    }

    const double alpha = std::clamp(blend, 0.0, 1.0);
    if (!(alpha > 0.0)) {
        return current;
    }

    ReferenceWrench blended = current;
    const chrono::ChVector3d previous_force_W = previous->reference_force_W;
    const chrono::ChVector3d previous_moment_W =
        previous->reference_moment_W +
        chrono::Vcross(previous->reference_origin_W - current.origin_W, previous_force_W);

    blended.force_W = (1.0 - alpha) * current.force_W + alpha * previous_force_W;
    blended.moment_W = (1.0 - alpha) * current.moment_W + alpha * previous_moment_W;
    blended.total_load = (1.0 - alpha) * current.total_load + alpha * previous->reference_total_load;
    return blended;
}

std::vector<std::ptrdiff_t> MatchSubpatches(const std::vector<DenseSubpatch>& current_subpatches,
                                            const std::vector<TemporalSubpatchState>& previous_subpatches,
                                            const CompressedContactConfig& cfg) {
    std::vector<std::ptrdiff_t> matched_previous(current_subpatches.size(), -1);
    if (current_subpatches.empty() || previous_subpatches.empty() || !(cfg.warm_start_match_radius > 0.0)) {
        return matched_previous;
    }

    std::vector<SubpatchMatchCandidate> candidates;
    for (std::size_t current_index = 0; current_index < current_subpatches.size(); ++current_index) {
        const auto& current = current_subpatches[current_index];
        for (std::size_t previous_index = 0; previous_index < previous_subpatches.size(); ++previous_index) {
            const auto& previous = previous_subpatches[previous_index];
            const double normal_cos = chrono::Vdot(current.avg_normal_W, previous.avg_normal_W);
            if (normal_cos < cfg.normal_cos_min) {
                continue;
            }

            const double gate = SubpatchMatchGate(cfg, current, previous);
            const double distance = (current.centroid_W - previous.centroid_W).Length();
            if (distance > gate) {
                continue;
            }

            SubpatchMatchCandidate candidate;
            candidate.current_index = current_index;
            candidate.previous_index = previous_index;
            candidate.score = distance / gate + 0.5 * (1.0 - normal_cos);
            candidates.push_back(candidate);
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const SubpatchMatchCandidate& a, const SubpatchMatchCandidate& b) {
                  return a.score < b.score;
              });

    std::vector<char> previous_taken(previous_subpatches.size(), 0);
    std::vector<char> current_taken(current_subpatches.size(), 0);
    for (const auto& candidate : candidates) {
        if (current_taken[candidate.current_index] || previous_taken[candidate.previous_index]) {
            continue;
        }
        current_taken[candidate.current_index] = 1;
        previous_taken[candidate.previous_index] = 1;
        matched_previous[candidate.current_index] = static_cast<std::ptrdiff_t>(candidate.previous_index);
    }

    return matched_previous;
}

double LookupPreviousAllocatedLoad(const ReducedSupportAggregate& support,
                                   const std::vector<ReducedContactPoint>& previous_contacts,
                                   double match_radius) {
    if (previous_contacts.empty()) {
        return std::max(0.0, support.allocated_load);
    }

    double best_distance = std::numeric_limits<double>::infinity();
    double best_load = std::max(0.0, support.allocated_load);
    for (const auto& previous : previous_contacts) {
        const double distance = (previous.x_W - support.x_W).Length();
        if (distance < best_distance) {
            best_distance = distance;
            best_load = std::max(0.0, previous.allocated_load);
        }
    }

    if (match_radius > 0.0 && best_distance > 2.0 * match_radius) {
        return std::max(0.0, support.allocated_load);
    }
    return best_load;
}

std::vector<std::ptrdiff_t> MatchSupportsToPrevious(const std::vector<ReducedSupportAggregate>& supports,
                                                    const std::vector<ReducedContactPoint>& previous_contacts,
                                                    double match_radius) {
    std::vector<std::ptrdiff_t> matched_previous(supports.size(), -1);
    if (supports.empty() || previous_contacts.empty() || !(match_radius > 0.0)) {
        return matched_previous;
    }

    std::vector<SupportMatchCandidate> candidates;
    for (std::size_t current_index = 0; current_index < supports.size(); ++current_index) {
        const auto& support = supports[current_index];
        for (std::size_t previous_index = 0; previous_index < previous_contacts.size(); ++previous_index) {
            const auto& previous = previous_contacts[previous_index];
            const double normal_cos = chrono::Vdot(support.n_W, previous.n_W);
            if (normal_cos < 0.85) {
                continue;
            }

            const double distance = (support.x_W - previous.x_W).Length();
            if (distance > 2.0 * match_radius) {
                continue;
            }

            SupportMatchCandidate candidate;
            candidate.current_index = current_index;
            candidate.previous_index = previous_index;
            candidate.score = distance / std::max(2.0 * match_radius, 1.0e-6) + 0.25 * (1.0 - normal_cos);
            candidates.push_back(candidate);
        }
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const SupportMatchCandidate& a, const SupportMatchCandidate& b) {
                  return a.score < b.score;
              });

    std::vector<char> previous_taken(previous_contacts.size(), 0);
    std::vector<char> current_taken(supports.size(), 0);
    for (const auto& candidate : candidates) {
        if (current_taken[candidate.current_index] || previous_taken[candidate.previous_index]) {
            continue;
        }
        current_taken[candidate.current_index] = 1;
        previous_taken[candidate.previous_index] = 1;
        matched_previous[candidate.current_index] = static_cast<std::ptrdiff_t>(candidate.previous_index);
    }
    return matched_previous;
}

std::vector<std::size_t> AssignSupportIds(const std::vector<std::ptrdiff_t>& matched_previous,
                                          const std::vector<ReducedContactPoint>& previous_contacts) {
    std::vector<std::size_t> support_ids(matched_previous.size(), 0);
    std::unordered_set<std::size_t> used_ids;
    for (std::size_t i = 0; i < matched_previous.size(); ++i) {
        if (matched_previous[i] >= 0) {
            const std::size_t support_id =
                previous_contacts[static_cast<std::size_t>(matched_previous[i])].support_id;
            support_ids[i] = support_id;
            used_ids.insert(support_id);
        }
    }

    std::size_t next_support_id = 0;
    for (std::size_t i = 0; i < matched_previous.size(); ++i) {
        if (matched_previous[i] >= 0) {
            continue;
        }
        while (used_ids.find(next_support_id) != used_ids.end()) {
            ++next_support_id;
        }
        support_ids[i] = next_support_id;
        used_ids.insert(next_support_id);
        ++next_support_id;
    }

    return support_ids;
}

std::vector<std::size_t> SelectSupportIndices(const std::vector<DenseContactPoint>& dense_points,
                                              const DenseSubpatch& subpatch,
                                              const std::vector<ReducedContactPoint>& previous_contacts,
                                              double match_radius,
                                              int target_points) {
    std::vector<std::size_t> selected;
    if (subpatch.members.empty() || target_points <= 0) {
        return selected;
    }

    target_points = std::max(1, target_points);

    std::size_t max_load_index = subpatch.members.front();
    double max_load = ProxyLoad(dense_points[max_load_index]);
    for (const auto point_index : subpatch.members) {
        const double load = ProxyLoad(dense_points[point_index]);
        if (load > max_load) {
            max_load = load;
            max_load_index = point_index;
        }
    }
    selected.push_back(max_load_index);
    if (target_points == 1 || subpatch.members.size() == 1) {
        return selected;
    }

    const chrono::ChVector3d center_W = SubpatchReferenceCenter(dense_points, subpatch);

    if (match_radius > 0.0) {
        for (const auto& previous : previous_contacts) {
            if (chrono::Vdot(previous.n_W, subpatch.avg_normal_W) < 0.85) {
                continue;
            }

            std::size_t best_index = subpatch.members.front();
            double best_distance = std::numeric_limits<double>::infinity();
            for (const auto point_index : subpatch.members) {
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
        const chrono::ChVector3d axis_W =
            std::cos(theta) * subpatch.t1_W + std::sin(theta) * subpatch.t2_W;

        std::size_t best_index = subpatch.members.front();
        double best_score = -std::numeric_limits<double>::infinity();
        for (const auto point_index : subpatch.members) {
            const auto& point = dense_points[point_index];
            const chrono::ChVector3d rel = point.x_W - center_W;
            const double projected = chrono::Vdot(rel, axis_W);
            const double radial = std::sqrt(TangentialRadiusSquared(point.x_W, center_W, subpatch.avg_normal_W));
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
           selected.size() < subpatch.members.size()) {
        std::size_t best_index = subpatch.members.front();
        double best_score = -1.0;
        for (const auto point_index : subpatch.members) {
            if (std::find(selected.begin(), selected.end(), point_index) != selected.end()) {
                continue;
            }
            double score = 0.0;
            for (const auto chosen_index : selected) {
                score += std::sqrt(
                    TangentialRadiusSquared(dense_points[point_index].x_W, dense_points[chosen_index].x_W,
                                            subpatch.avg_normal_W));
            }
            score +=
                std::sqrt(TangentialRadiusSquared(dense_points[point_index].x_W, center_W, subpatch.avg_normal_W));
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

std::vector<ReducedSupportAggregate> BuildReducedSupportsForSubpatch(const std::vector<DenseContactPoint>& dense_points,
                                                                     const DenseSubpatch& subpatch,
                                                                     const std::vector<ReducedContactPoint>& previous_contacts,
                                                                     double match_radius,
                                                                     int target_points) {
    target_points = std::max(1, std::min(target_points, static_cast<int>(subpatch.members.size())));
    const auto support_indices = SelectSupportIndices(dense_points, subpatch, previous_contacts, match_radius,
                                                      target_points);
    if (support_indices.empty()) {
        return {};
    }

    std::vector<chrono::ChVector3d> centers_W;
    centers_W.reserve(support_indices.size());
    for (const auto dense_index : support_indices) {
        centers_W.push_back(dense_points[dense_index].x_W);
    }

    std::vector<std::size_t> assignments(subpatch.members.size(), 0);
    for (std::size_t local_index = 0; local_index < subpatch.members.size(); ++local_index) {
        const auto member_index = subpatch.members[local_index];
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

    std::vector<ReducedSupportAggregate> supports;
    supports.reserve(centers_W.size());
    for (std::size_t cluster = 0; cluster < centers_W.size(); ++cluster) {
        double sum_area = 0.0;
        double sum_load = 0.0;
        double sum_geom_weight = 0.0;
        double sum_load_weight = 0.0;
        double max_coverage_radius = 0.0;
        std::size_t count = 0;
        chrono::ChVector3d avg_x_W(0.0, 0.0, 0.0);
        chrono::ChVector3d avg_x_master_M(0.0, 0.0, 0.0);
        chrono::ChVector3d avg_n_W(0.0, 0.0, 0.0);
        chrono::ChVector3d avg_v_rel_W(0.0, 0.0, 0.0);
        double avg_phi = 0.0;
        double avg_phi_eff = 0.0;
        double avg_phi_eff_load = 0.0;

        chrono::ChVector3d anchor_W = centers_W[cluster];
        bool anchored_to_history = false;
        const ReducedContactPoint* anchor_contact = nullptr;
        if (match_radius > 0.0) {
            double best_anchor_distance = std::numeric_limits<double>::infinity();
            for (const auto& previous : previous_contacts) {
                const double distance = (previous.x_W - centers_W[cluster]).Length();
                if (distance < best_anchor_distance) {
                    best_anchor_distance = distance;
                    anchor_W = previous.x_W;
                    anchor_contact = &previous;
                }
            }
            anchored_to_history = best_anchor_distance <= match_radius;
        }

        for (std::size_t local_index = 0; local_index < subpatch.members.size(); ++local_index) {
            if (assignments[local_index] != cluster) {
                continue;
            }

            const auto member_index = subpatch.members[local_index];
            const auto& dense = dense_points[member_index];
            const double load = ProxyLoad(dense);
            const double geom_weight = (load > 1.0e-12) ? load : std::max(1.0e-12, dense.area_weight);

            sum_area += dense.area_weight;
            sum_load += load;
            sum_geom_weight += geom_weight;
            if (load > 1.0e-12) {
                sum_load_weight += load;
                avg_phi_eff_load += load * dense.phi_eff;
            }
            ++count;

            avg_x_W += geom_weight * dense.x_W;
            avg_x_master_M += geom_weight * dense.x_master_M;
            avg_n_W += geom_weight * dense.n_W;
            avg_v_rel_W += geom_weight * dense.v_rel_W;
            avg_phi += geom_weight * dense.phi;
            avg_phi_eff += geom_weight * dense.phi_eff;
        }

        if (!(sum_geom_weight > 1.0e-12) || count == 0) {
            continue;
        }

        avg_x_W *= (1.0 / sum_geom_weight);
        avg_x_master_M *= (1.0 / sum_geom_weight);
        avg_n_W *= (1.0 / sum_geom_weight);
        avg_v_rel_W *= (1.0 / sum_geom_weight);
        avg_phi *= (1.0 / sum_geom_weight);
        avg_phi_eff *= (1.0 / sum_geom_weight);
        if (sum_load_weight > 1.0e-12) {
            avg_phi_eff = avg_phi_eff_load * (1.0 / sum_load_weight);
        }

        if (anchored_to_history && anchor_contact) {
            constexpr double kAnchorBlend = 0.5;
            avg_x_W = (1.0 - kAnchorBlend) * avg_x_W + kAnchorBlend * anchor_contact->x_W;
            avg_x_master_M =
                (1.0 - kAnchorBlend) * avg_x_master_M + kAnchorBlend * anchor_contact->x_master_M;
            avg_n_W = (1.0 - kAnchorBlend) * avg_n_W + kAnchorBlend * anchor_contact->n_W;
            avg_v_rel_W = (1.0 - kAnchorBlend) * avg_v_rel_W + kAnchorBlend * anchor_contact->v_rel_W;
            avg_phi = (1.0 - kAnchorBlend) * avg_phi + kAnchorBlend * anchor_contact->phi;
            avg_phi_eff = (1.0 - kAnchorBlend) * avg_phi_eff + kAnchorBlend * anchor_contact->phi_eff;
        }

        for (std::size_t local_index = 0; local_index < subpatch.members.size(); ++local_index) {
            if (assignments[local_index] != cluster) {
                continue;
            }
            const auto member_index = subpatch.members[local_index];
            max_coverage_radius =
                std::max(max_coverage_radius,
                         std::sqrt(TangentialRadiusSquared(dense_points[member_index].x_W, avg_x_W,
                                                           subpatch.avg_normal_W)));
        }

        ReducedSupportAggregate support;
        support.patch_id = subpatch.patch_id;
        support.subpatch_id = subpatch.subpatch_id;
        support.x_W = avg_x_W;
        support.x_master_M = avg_x_master_M;
        support.n_W = avg_n_W;
        const double normal_len = support.n_W.Length();
        if (normal_len > 1.0e-12) {
            support.n_W *= (1.0 / normal_len);
        } else {
            support.n_W = subpatch.avg_normal_W;
        }
        support.v_rel_W = avg_v_rel_W;
        support.phi = avg_phi;
        support.phi_eff = avg_phi_eff;
        support.x_master_surface_W = support.x_W - support.phi * support.n_W;
        support.area_weight = sum_area;
        support.support_weight =
            (support.phi_eff < -1.0e-12) ? (sum_load / std::max(1.0e-12, -support.phi_eff)) : sum_area;
        support.allocated_load = sum_load;
        support.coverage_radius = max_coverage_radius;
        support.dense_members = count;
        supports.push_back(support);
    }

    return supports;
}

double ComputeSubpatchGapError(const DenseSubpatch& subpatch,
                               const std::vector<ReducedSupportAggregate>& supports,
                               const RigidBodyStateW& master_state,
                               const FirstOrderSDF& sdf,
                               const CompressedContactConfig& cfg) {
    if (supports.empty()) {
        return 0.0;
    }

    double max_gap_error = 0.0;
    for (const auto& sentinel_W : subpatch.sentinel_W) {
        const chrono::ChVector3d sentinel_M = master_state.R_WRef.transpose() * (sentinel_W - master_state.x_ref_W);
        double phi = 0.0;
        if (!sdf.QueryPhiM(sentinel_M, phi) || !(phi < 0.0)) {
            continue;
        }

        double min_uncovered_distance = std::numeric_limits<double>::infinity();
        for (const auto& support : supports) {
            const double tangential =
                std::sqrt(TangentialRadiusSquared(sentinel_W, support.x_W, subpatch.avg_normal_W));
            min_uncovered_distance =
                std::min(min_uncovered_distance, tangential - support.coverage_radius - cfg.sentinel_margin);
        }

        if (min_uncovered_distance > 0.0) {
            max_gap_error = std::max(max_gap_error, -phi);
        }
    }

    return max_gap_error;
}

}  // namespace

void CompressedContactPipeline::Configure(const CompressedContactConfig& cfg) {
    cfg_ = cfg;
    if (!slave_surface_samples_.empty()) {
        dense_sample_bvh_.Build(slave_surface_samples_, cfg_.bvh_leaf_size);
    }
    previous_contacts_.clear();
    previous_subpatches_.clear();
    next_persistent_id_ = 1;
}

void CompressedContactPipeline::SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples) {
    slave_surface_samples_ = std::move(samples);
    dense_sample_bvh_.Build(slave_surface_samples_, cfg_.bvh_leaf_size);
    previous_contacts_.clear();
    previous_subpatches_.clear();
    next_persistent_id_ = 1;
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

    std::vector<DensePatch> patches;
    SubpatchRefiner::BuildPatches(dense_points, cfg_, patches);

    std::vector<DenseSubpatch> subpatches;
    SubpatchRefiner::BuildSubpatches(dense_points, patches, cfg_, subpatches);
    const auto matched_previous = MatchSubpatches(subpatches, previous_subpatches_, cfg_);

    double max_subpatch_plane_error = 0.0;
    double max_subpatch_gap_error = 0.0;
    double max_subpatch_force_residual = 0.0;
    double max_subpatch_moment_residual = 0.0;

    out_contacts.reserve(dense_points.size());
    std::vector<TemporalSubpatchState> current_subpatches;
    current_subpatches.reserve(subpatches.size());
    for (std::size_t subpatch_index = 0; subpatch_index < subpatches.size(); ++subpatch_index) {
        const auto& subpatch = subpatches[subpatch_index];
        max_subpatch_plane_error = std::max(max_subpatch_plane_error, subpatch.plane_error);

        std::vector<ReducedContactPoint> previous_local_contacts;
        std::size_t persistent_id = 0;
        const TemporalSubpatchState* matched_previous_state = nullptr;
        if (matched_previous[subpatch_index] >= 0) {
            matched_previous_state =
                &previous_subpatches_[static_cast<std::size_t>(matched_previous[subpatch_index])];
            persistent_id = matched_previous_state->persistent_id;
            previous_local_contacts = matched_previous_state->contacts;
        } else if (cfg_.warm_start_match_radius > 0.0) {
            const double gate = std::max(
                {cfg_.warm_start_match_radius, cfg_.max_subpatch_diameter, cfg_.max_patch_diameter, 1.0e-6});
            for (const auto& previous : previous_contacts_) {
                if (chrono::Vdot(previous.n_W, subpatch.avg_normal_W) < cfg_.normal_cos_min) {
                    continue;
                }
                if ((previous.x_W - subpatch.centroid_W).Length() <= gate) {
                    previous_local_contacts.push_back(previous);
                }
            }
        }
        std::sort(previous_local_contacts.begin(), previous_local_contacts.end(),
                  [](const ReducedContactPoint& a, const ReducedContactPoint& b) {
                      return a.support_id < b.support_id;
                  });
        if (persistent_id == 0) {
            persistent_id = next_persistent_id_++;
        }
        auto selection_previous_contacts = previous_local_contacts;
        std::sort(selection_previous_contacts.begin(), selection_previous_contacts.end(),
                  [](const ReducedContactPoint& a, const ReducedContactPoint& b) {
                      if (a.allocated_load != b.allocated_load) {
                          return a.allocated_load > b.allocated_load;
                      }
                      return a.support_id < b.support_id;
                  });

        int target_points = cfg_.max_reduced_points_per_patch;
        if (cfg_.max_subpatch_diameter > 0.0 && subpatch.diameter <= 0.75 * cfg_.max_subpatch_diameter) {
            target_points = std::min(target_points, 3);
        }

        auto supports = BuildReducedSupportsForSubpatch(dense_points, subpatch, selection_previous_contacts,
                                                        cfg_.warm_start_match_radius, target_points);
        const auto matched_supports =
            MatchSupportsToPrevious(supports, previous_local_contacts, cfg_.warm_start_match_radius);
        const auto support_ids = AssignSupportIds(matched_supports, previous_local_contacts);
        std::vector<std::size_t> support_order(supports.size());
        std::iota(support_order.begin(), support_order.end(), 0);
        std::sort(support_order.begin(), support_order.end(),
                  [&support_ids](std::size_t a, std::size_t b) { return support_ids[a] < support_ids[b]; });
        const ReferenceWrench dense_reference =
            LocalWrenchAllocator::BuildDenseReference(dense_points, subpatch.members, subpatch.centroid_W);
        const ReferenceWrench reference =
            BlendTemporalReference(dense_reference, matched_previous_state, cfg_.temporal_reference_blend);

        if (!supports.empty()) {
            std::vector<SupportWrenchPoint> support_points;
            support_points.reserve(supports.size());
            for (std::size_t support_index = 0; support_index < supports.size(); ++support_index) {
                const auto& support = supports[support_index];
                SupportWrenchPoint point;
                point.x_W = support.x_W;
                point.n_W = support.n_W;
                if (matched_supports[support_index] >= 0) {
                    point.initial_load = std::max(
                        0.0, previous_local_contacts[static_cast<std::size_t>(matched_supports[support_index])]
                                 .allocated_load);
                } else {
                    point.initial_load =
                        LookupPreviousAllocatedLoad(support, previous_local_contacts, cfg_.warm_start_match_radius);
                }
                support_points.push_back(point);
            }

            WrenchAllocationResult allocation;
            LocalWrenchAllocator::Allocate(support_points, reference, cfg_.temporal_load_regularization, allocation);
            max_subpatch_force_residual = std::max(max_subpatch_force_residual, allocation.force_residual);
            max_subpatch_moment_residual = std::max(max_subpatch_moment_residual, allocation.moment_residual);

            if (allocation.feasible && allocation.loads.size() == supports.size()) {
                for (std::size_t i = 0; i < supports.size(); ++i) {
                    supports[i].allocated_load = allocation.loads[i];
                    if (supports[i].phi_eff < -1.0e-12) {
                        supports[i].support_weight = allocation.loads[i] / std::max(1.0e-12, -supports[i].phi_eff);
                    } else {
                        supports[i].support_weight = supports[i].area_weight;
                    }
                }
            }
        }

        double mean_allocated_load = 0.0;
        for (const auto& support : supports) {
            mean_allocated_load += std::max(0.0, support.allocated_load);
        }
        if (!supports.empty()) {
            mean_allocated_load /= static_cast<double>(supports.size());
        }

        max_subpatch_gap_error = std::max(max_subpatch_gap_error,
                                          ComputeSubpatchGapError(subpatch, supports, master_state, sdf, cfg_));

        const std::size_t contacts_begin = out_contacts.size();
        for (const auto support_index : support_order) {
            const auto& support = supports[support_index];
            ReducedContactPoint reduced;
            reduced.persistent_id = persistent_id;
            reduced.patch_id = support.patch_id;
            reduced.subpatch_id = support.subpatch_id;
            reduced.support_id = support_ids[support_index];
            reduced.dense_members = support.dense_members;
            if (mean_allocated_load > 1.0e-12) {
                reduced.emission_count = std::clamp(
                    static_cast<int>(std::lround(std::max(0.0, support.allocated_load) / mean_allocated_load)), 1, 2);
            } else {
                reduced.emission_count = 1;
            }
            reduced.x_W = support.x_W;
            reduced.x_master_M = support.x_master_M;
            reduced.x_master_surface_W = support.x_master_surface_W;
            reduced.n_W = support.n_W;
            reduced.v_rel_W = support.v_rel_W;
            reduced.stencil_axis_W = ChooseStencilAxis(support, subpatch);
            reduced.phi = support.phi;
            reduced.phi_eff = support.phi_eff;
            reduced.area_weight = support.area_weight;
            reduced.support_weight = support.support_weight;
            reduced.allocated_load = support.allocated_load;
            reduced.stencil_half_extent =
                (reduced.emission_count > 1) ? std::min(0.35 * support.coverage_radius, 0.15 * subpatch.diameter)
                                             : 0.0;
            reduced.mu = mu_default;
            out_contacts.push_back(reduced);
        }

        TemporalSubpatchState temporal_state;
        temporal_state.persistent_id = persistent_id;
        temporal_state.centroid_W = subpatch.centroid_W;
        temporal_state.avg_normal_W = subpatch.avg_normal_W;
        temporal_state.diameter = subpatch.diameter;
        temporal_state.reference_origin_W = reference.origin_W;
        temporal_state.reference_force_W = reference.force_W;
        temporal_state.reference_moment_W = reference.moment_W;
        temporal_state.reference_total_load = reference.total_load;
        for (std::size_t contact_index = contacts_begin; contact_index < out_contacts.size(); ++contact_index) {
            temporal_state.contacts.push_back(out_contacts[contact_index]);
        }
        current_subpatches.push_back(std::move(temporal_state));
    }

    previous_contacts_ = out_contacts;
    previous_subpatches_ = std::move(current_subpatches);

    if (!out_stats) {
        return;
    }

    out_stats->total_samples = cloud_stats.total_samples;
    out_stats->candidate_count = cloud_stats.candidate_samples;
    out_stats->dense_count = dense_points.size();
    out_stats->reduced_count = out_contacts.size();
    out_stats->patch_count = patches.size();
    out_stats->subpatch_count = subpatches.size();
    out_stats->bvh_nodes_visited = cloud_stats.bvh.nodes_visited;
    out_stats->bvh_nodes_pruned_obb = cloud_stats.bvh.nodes_pruned_obb;
    out_stats->bvh_nodes_pruned_sdf = cloud_stats.bvh.nodes_pruned_sdf;
    out_stats->bvh_leaf_samples_tested = cloud_stats.bvh.leaf_samples_tested;
    out_stats->max_subpatch_plane_error = max_subpatch_plane_error;
    out_stats->max_subpatch_gap_error = max_subpatch_gap_error;
    out_stats->max_subpatch_force_residual = max_subpatch_force_residual;
    out_stats->max_subpatch_moment_residual = max_subpatch_moment_residual;

    if (dense_points.empty() || out_contacts.empty()) {
        out_stats->epsilon_F = 0.0;
        out_stats->epsilon_M = 0.0;
        out_stats->epsilon_CoP = 0.0;
        out_stats->epsilon_gap = max_subpatch_gap_error;
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
    double dense_worst_gap = 0.0;
    for (const auto& dense : dense_points) {
        dense_worst_gap = std::min(dense_worst_gap, dense.phi_eff);
    }
    const double dense_gap = std::max(1.0e-12, std::abs(dense_worst_gap));
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
    const double worst_gap_mismatch = std::abs(reduced_worst_gap - dense_worst_gap) / dense_gap;
    out_stats->epsilon_gap = std::max(worst_gap_mismatch, max_subpatch_gap_error / dense_gap);
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
