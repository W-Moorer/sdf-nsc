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
    chrono::ChVector3d allocated_force_W;
    double coverage_radius = 0.0;
    std::size_t dense_members = 0;
    std::vector<std::size_t> member_indices;
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

struct SupportTransportSeed {
    double load = 0.0;
    chrono::ChVector3d force_W;
    double confidence = 0.0;
};

struct TransportBlendWeights {
    double scalar = 0.0;
    double normal = 0.0;
    double tangential = 0.0;
};

struct SecondMoment2D {
    double uu = 0.0;
    double uv = 0.0;
    double vv = 0.0;
};

double Clamp01(double value) {
    return std::clamp(value, 0.0, 1.0);
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

chrono::ChMatrix33<> BuildContactPlane(const chrono::ChVector3d& n_W) {
    chrono::ChMatrix33<> contact_plane;
    contact_plane.SetFromAxisX(n_W, chrono::VECT_Y);
    return contact_plane;
}

void ClearReactionCache(std::array<float, 6>& reaction_cache) {
    reaction_cache.fill(0.0f);
}

void EncodeReactionCacheWorldImpulse(std::array<float, 6>& reaction_cache,
                                     const chrono::ChVector3d& impulse_W,
                                     const chrono::ChVector3d& n_W) {
    ClearReactionCache(reaction_cache);
    const chrono::ChMatrix33<> contact_plane = BuildContactPlane(n_W);
    const chrono::ChVector3d local_impulse = contact_plane.transpose() * impulse_W;
    reaction_cache[0] = static_cast<float>(local_impulse.x());
    reaction_cache[1] = static_cast<float>(local_impulse.y());
    reaction_cache[2] = static_cast<float>(local_impulse.z());
}

chrono::ChVector3d DecodeReactionCacheWorldForce(const std::array<float, 6>& reaction_cache,
                                                 const chrono::ChVector3d& n_W) {
    const chrono::ChMatrix33<> contact_plane = BuildContactPlane(n_W);
    const chrono::ChVector3d local_force(reaction_cache[0], reaction_cache[1], reaction_cache[2]);
    return contact_plane * local_force;
}

chrono::ChVector3d DecodeTotalReactionForce(const ReducedContactPoint& contact) {
    const chrono::ChVector3d center_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_primary, contact.n_W);
    const chrono::ChVector3d negative_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_secondary, contact.n_W);
    const chrono::ChVector3d positive_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_tertiary, contact.n_W);
    const chrono::ChVector3d secondary_negative_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_quaternary, contact.n_W);
    const chrono::ChVector3d secondary_positive_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_quinary, contact.n_W);

    if (contact.emission_count <= 1) {
        return center_force_W;
    }
    if (contact.emission_count == 2) {
        return negative_force_W + positive_force_W;
    }
    if (contact.emission_count >= 5) {
        return center_force_W + negative_force_W + positive_force_W + secondary_negative_force_W +
               secondary_positive_force_W;
    }
    return center_force_W + negative_force_W + positive_force_W;
}

chrono::ChVector3d DecodeTotalReactionEquivalentForce(const ReducedContactPoint& contact, double step_size) {
    const chrono::ChVector3d cached_impulse_W = DecodeTotalReactionForce(contact);
    if (!(step_size > 1.0e-12)) {
        return cached_impulse_W;
    }
    return cached_impulse_W * (1.0 / step_size);
}

double DenseMetricWeight(const DenseContactPoint& point) {
    const double proxy = ProxyLoad(point);
    return (proxy > 1.0e-12) ? proxy : std::max(1.0e-12, point.area_weight);
}

double SupportMetricWeight(const ReducedSupportAggregate& support) {
    if (support.allocated_load > 1.0e-12) {
        return support.allocated_load;
    }
    if (support.support_weight > 1.0e-12) {
        return support.support_weight;
    }
    return std::max(1.0e-12, support.area_weight);
}

double ComputeSupportAxisSpread(const std::vector<DenseContactPoint>& dense_points,
                                const ReducedSupportAggregate& support,
                                const chrono::ChVector3d& center_W,
                                const chrono::ChVector3d& axis_W) {
    if (support.member_indices.empty()) {
        return 0.0;
    }

    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    for (const auto dense_index : support.member_indices) {
        const auto& point = dense_points[dense_index];
        const double weight = DenseMetricWeight(point);
        const double offset = chrono::Vdot(point.x_W - center_W, axis_W);
        weighted_sum += weight * offset * offset;
        weight_sum += weight;
    }
    if (!(weight_sum > 1.0e-12)) {
        return 0.0;
    }
    return std::sqrt(weighted_sum / weight_sum);
}

SecondMoment2D ComputeDenseSecondMoment(const std::vector<DenseContactPoint>& dense_points,
                                        const DenseSubpatch& subpatch,
                                        const chrono::ChVector3d& center_W) {
    SecondMoment2D moment;
    double weight_sum = 0.0;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = DenseMetricWeight(point);
        const chrono::ChVector3d rel = point.x_W - center_W;
        const double u = chrono::Vdot(rel, subpatch.t1_W);
        const double v = chrono::Vdot(rel, subpatch.t2_W);
        moment.uu += weight * u * u;
        moment.uv += weight * u * v;
        moment.vv += weight * v * v;
        weight_sum += weight;
    }
    if (weight_sum > 1.0e-12) {
        moment.uu /= weight_sum;
        moment.uv /= weight_sum;
        moment.vv /= weight_sum;
    }
    return moment;
}

SecondMoment2D ComputeReducedSecondMoment(const std::vector<ReducedSupportAggregate>& supports,
                                          const DenseSubpatch& subpatch,
                                          const chrono::ChVector3d& center_W) {
    SecondMoment2D moment;
    double weight_sum = 0.0;
    for (const auto& support : supports) {
        const double weight = SupportMetricWeight(support);
        const chrono::ChVector3d rel = support.x_W - center_W;
        const double u = chrono::Vdot(rel, subpatch.t1_W);
        const double v = chrono::Vdot(rel, subpatch.t2_W);
        moment.uu += weight * u * u;
        moment.uv += weight * u * v;
        moment.vv += weight * v * v;
        weight_sum += weight;
    }
    if (weight_sum > 1.0e-12) {
        moment.uu /= weight_sum;
        moment.uv /= weight_sum;
        moment.vv /= weight_sum;
    }
    return moment;
}

double ComputeSecondMomentError(const std::vector<DenseContactPoint>& dense_points,
                                const DenseSubpatch& subpatch,
                                const std::vector<ReducedSupportAggregate>& supports) {
    if (subpatch.members.empty() || supports.empty()) {
        return 0.0;
    }

    const chrono::ChVector3d center_W = SubpatchReferenceCenter(dense_points, subpatch);
    const auto dense_moment = ComputeDenseSecondMoment(dense_points, subpatch, center_W);
    const auto reduced_moment = ComputeReducedSecondMoment(supports, subpatch, center_W);
    const double dense_norm = std::sqrt(dense_moment.uu * dense_moment.uu + 2.0 * dense_moment.uv * dense_moment.uv +
                                        dense_moment.vv * dense_moment.vv);
    const double diff_norm = std::sqrt((reduced_moment.uu - dense_moment.uu) * (reduced_moment.uu - dense_moment.uu) +
                                       2.0 * (reduced_moment.uv - dense_moment.uv) *
                                           (reduced_moment.uv - dense_moment.uv) +
                                       (reduced_moment.vv - dense_moment.vv) *
                                           (reduced_moment.vv - dense_moment.vv));
    return diff_norm / std::max(dense_norm, 1.0e-12);
}

double ComputeConeSupportError(const std::vector<DenseContactPoint>& dense_points,
                               const DenseSubpatch& subpatch,
                               const std::vector<ReducedSupportAggregate>& supports,
                               const CompressedContactConfig& cfg) {
    if (subpatch.members.empty() || supports.empty()) {
        return 0.0;
    }

    const chrono::ChVector3d center_W = SubpatchReferenceCenter(dense_points, subpatch);
    const int dir_count = std::max(8, cfg.cone_direction_count);
    double max_error = 0.0;
    for (int dir_index = 0; dir_index < dir_count; ++dir_index) {
        const double theta =
            (2.0 * std::acos(-1.0) * static_cast<double>(dir_index)) / static_cast<double>(dir_count);
        const chrono::ChVector3d dir_W = std::cos(theta) * subpatch.t1_W + std::sin(theta) * subpatch.t2_W;

        double dense_support = -std::numeric_limits<double>::infinity();
        for (const auto point_index : subpatch.members) {
            dense_support = std::max(dense_support, chrono::Vdot(dense_points[point_index].x_W - center_W, dir_W));
        }

        double reduced_support = -std::numeric_limits<double>::infinity();
        for (const auto& support : supports) {
            reduced_support = std::max(reduced_support, chrono::Vdot(support.x_W - center_W, dir_W));
        }

        max_error = std::max(max_error, std::abs(reduced_support - dense_support));
    }

    return max_error / std::max(subpatch.diameter, 1.0e-9);
}

chrono::ChVector3d ComputeReducedForce(const std::vector<ReducedSupportAggregate>& supports) {
    chrono::ChVector3d force_W(0.0, 0.0, 0.0);
    for (const auto& support : supports) {
        force_W += support.allocated_force_W;
    }
    return force_W;
}

chrono::ChVector3d ComputeReducedMoment(const std::vector<ReducedSupportAggregate>& supports,
                                        const chrono::ChVector3d& origin_W) {
    chrono::ChVector3d moment_W(0.0, 0.0, 0.0);
    for (const auto& support : supports) {
        moment_W += chrono::Vcross(support.x_W - origin_W, support.allocated_force_W);
    }
    return moment_W;
}

chrono::ChVector3d ComputeReferenceCoP(const ReferenceWrench& reference, const DenseSubpatch& subpatch) {
    const double normal_load = chrono::Vdot(reference.force_W, subpatch.avg_normal_W);
    if (!(normal_load > 1.0e-12)) {
        return reference.origin_W;
    }
    return reference.origin_W + chrono::Vcross(subpatch.avg_normal_W, reference.moment_W) * (1.0 / normal_load);
}

chrono::ChVector3d ComputeReducedCoP(const std::vector<ReducedSupportAggregate>& supports,
                                     const DenseSubpatch& subpatch) {
    chrono::ChVector3d cop_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto& support : supports) {
        const double normal_load = std::max(0.0, chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W));
        cop_W += normal_load * support.x_W;
        weight_sum += normal_load;
    }
    if (!(weight_sum > 1.0e-12)) {
        return subpatch.centroid_W;
    }
    return cop_W * (1.0 / weight_sum);
}

double ComputeReferenceWrenchError(const std::vector<ReducedSupportAggregate>& supports,
                                   const DenseSubpatch& subpatch,
                                   const ReferenceWrench& reference) {
    if (supports.empty()) {
        return 0.0;
    }

    const chrono::ChVector3d reduced_force_W = ComputeReducedForce(supports);
    const chrono::ChVector3d reduced_moment_W = ComputeReducedMoment(supports, reference.origin_W);
    const double force_scale = std::max(reference.force_W.Length(), 1.0e-12);
    const double radius_scale = std::max(0.5 * subpatch.diameter, 1.0e-6);
    const double moment_scale =
        std::max({reference.moment_W.Length(), reference.force_W.Length() * radius_scale, 1.0e-12});
    const double force_error = (reduced_force_W - reference.force_W).Length() / force_scale;
    const double moment_error = (reduced_moment_W - reference.moment_W).Length() / moment_scale;
    return std::sqrt(force_error * force_error + moment_error * moment_error);
}

double ComputeReferenceCoPError(const std::vector<ReducedSupportAggregate>& supports,
                                const DenseSubpatch& subpatch,
                                const ReferenceWrench& reference) {
    if (supports.empty()) {
        return 0.0;
    }
    return (ComputeReducedCoP(supports, subpatch) - ComputeReferenceCoP(reference, subpatch)).Length();
}

chrono::ChVector3d ChooseStencilAxis(const ReducedSupportAggregate& support, const DenseSubpatch& subpatch) {
    const chrono::ChVector3d tangential_allocated_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const chrono::ChVector3d from_allocated_force =
        ProjectToTangentUnit(chrono::Vcross(subpatch.avg_normal_W, tangential_allocated_force_W),
                             subpatch.avg_normal_W);
    if (from_allocated_force.Length2() > 0.0) {
        return from_allocated_force;
    }

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

chrono::ChVector3d ChooseSecondaryStencilAxis(const ReducedSupportAggregate& support,
                                              const DenseSubpatch& subpatch,
                                              const chrono::ChVector3d& primary_axis_W) {
    chrono::ChVector3d secondary_axis_W = chrono::Vcross(subpatch.avg_normal_W, primary_axis_W);
    secondary_axis_W = ProjectToTangentUnit(secondary_axis_W, subpatch.avg_normal_W);
    if (secondary_axis_W.Length2() > 0.0) {
        return secondary_axis_W;
    }

    const chrono::ChVector3d tangential_velocity_W =
        support.v_rel_W - chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    secondary_axis_W = ProjectToTangentUnit(tangential_velocity_W, subpatch.avg_normal_W);
    if (secondary_axis_W.Length2() > 0.0) {
        if (std::abs(chrono::Vdot(secondary_axis_W, primary_axis_W)) < 0.95) {
            return secondary_axis_W;
        }
    }

    secondary_axis_W = ProjectToTangentUnit(subpatch.t2_W, subpatch.avg_normal_W);
    if (secondary_axis_W.Length2() > 0.0 && std::abs(chrono::Vdot(secondary_axis_W, primary_axis_W)) < 0.98) {
        return secondary_axis_W;
    }

    return chrono::ChVector3d(0.0, 0.0, 0.0);
}

int ChooseEmissionCount(const ReducedSupportAggregate& support,
                        const DenseSubpatch& subpatch,
                        double mean_allocated_load) {
    const chrono::ChVector3d tangential_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const chrono::ChVector3d tangential_velocity_W =
        support.v_rel_W - chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    const double total_force_norm = support.allocated_force_W.Length();
    const double tangential_ratio = tangential_force_W.Length() / std::max(total_force_norm, 1.0e-12);
    const double coverage_ratio = support.coverage_radius / std::max(subpatch.diameter, 1.0e-9);
    const double load_ratio = std::max(0.0, support.allocated_load) / std::max(mean_allocated_load, 1.0e-12);
    const double approach_speed = std::max(0.0, -chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W));
    const double tangential_speed = tangential_velocity_W.Length();
    const double slip_dominance =
        tangential_speed / std::max(tangential_speed + 1.5 * approach_speed, 1.0e-9);

    if (tangential_ratio < 0.12 || coverage_ratio < 0.08 || load_ratio < 0.45) {
        return 1;
    }

    const double expand_score =
        1.40 * tangential_ratio + 0.85 * coverage_ratio + 0.30 * std::min(load_ratio, 2.5) +
        0.35 * slip_dominance;
    if (expand_score > 1.45 && tangential_ratio > 0.38 && coverage_ratio > 0.22 && support.dense_members >= 6 &&
        load_ratio > 0.9 && slip_dominance > 0.72) {
        return 5;
    }
    if (expand_score > 1.18 && tangential_ratio > 0.24 && coverage_ratio > 0.14 && support.dense_members >= 3 &&
        slip_dominance > 0.55) {
        return 3;
    }
    return (expand_score > 0.72 && slip_dominance > 0.22) ? 2 : 1;
}

double ChooseStencilHalfExtent(const ReducedSupportAggregate& support,
                               const DenseSubpatch& subpatch,
                               int emission_count) {
    if (emission_count <= 1) {
        return 0.0;
    }

    const chrono::ChVector3d tangential_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const chrono::ChVector3d tangential_velocity_W =
        support.v_rel_W - chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    const double total_force_norm = support.allocated_force_W.Length();
    const double tangential_ratio = tangential_force_W.Length() / std::max(total_force_norm, 1.0e-12);
    const double approach_speed = std::max(0.0, -chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W));
    const double tangential_speed = tangential_velocity_W.Length();
    const double slip_dominance =
        tangential_speed / std::max(tangential_speed + 1.5 * approach_speed, 1.0e-9);
    const double tangent_scale = std::clamp(0.5 + tangential_ratio, 0.5, 1.0);
    const double motion_scale = std::clamp(0.55 + 0.45 * slip_dominance, 0.55, 1.0);
    double extent = tangent_scale * motion_scale *
                    std::min(0.28 * support.coverage_radius, 0.14 * subpatch.diameter);
    if (emission_count >= 3) {
        extent *= 0.82;
    }
    return extent;
}

double ChooseSecondaryStencilHalfExtent(const ReducedSupportAggregate& support,
                                        const DenseSubpatch& subpatch,
                                        int emission_count) {
    if (emission_count < 5) {
        return 0.0;
    }

    const chrono::ChVector3d tangential_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const double total_force_norm = support.allocated_force_W.Length();
    const double tangential_ratio = tangential_force_W.Length() / std::max(total_force_norm, 1.0e-12);
    const double base_extent = std::min(0.18 * support.coverage_radius, 0.09 * subpatch.diameter);
    return std::clamp((0.45 + 0.35 * tangential_ratio) * base_extent, 0.0, base_extent);
}

void NormalizeStencilWeights(std::array<double, 5>& weights, int emission_count) {
    const int active_slots =
        (emission_count >= 5) ? 5 : ((emission_count >= 3) ? 3 : ((emission_count == 2) ? 2 : 1));
    double sum = 0.0;
    for (int slot = 0; slot < active_slots; ++slot) {
        weights[slot] = std::max(0.0, weights[slot]);
        sum += weights[slot];
    }
    if (!(sum > 1.0e-12)) {
        weights.fill(0.0);
        if (active_slots == 1) {
            weights[0] = 1.0;
        } else if (active_slots == 2) {
            weights[1] = 0.5;
            weights[2] = 0.5;
        } else if (active_slots == 3) {
            weights[0] = 0.34;
            weights[1] = 0.33;
            weights[2] = 0.33;
        } else {
            weights[0] = 0.20;
            weights[1] = 0.25;
            weights[2] = 0.25;
            weights[3] = 0.15;
            weights[4] = 0.15;
        }
        return;
    }
    for (int slot = 0; slot < active_slots; ++slot) {
        weights[slot] /= sum;
    }
    for (int slot = active_slots; slot < 5; ++slot) {
        weights[slot] = 0.0;
    }
}

std::array<double, 5> ComputeNormalStencilWeights(const std::vector<DenseContactPoint>& dense_points,
                                                  const ReducedSupportAggregate& support,
                                                  const DenseSubpatch& subpatch,
                                                  const chrono::ChVector3d& primary_axis_W,
                                                  const chrono::ChVector3d& secondary_axis_W,
                                                  int emission_count) {
    std::array<double, 5> weights{1.0, 0.0, 0.0, 0.0, 0.0};
    if (emission_count <= 1) {
        return weights;
    }

    if (emission_count == 2) {
        weights = {0.0, 0.5, 0.5, 0.0, 0.0};
        return weights;
    }

    const double coverage_ratio = support.coverage_radius / std::max(subpatch.diameter, 1.0e-9);
    const double primary_spread =
        ComputeSupportAxisSpread(dense_points, support, support.x_W, primary_axis_W);
    const double secondary_spread =
        (emission_count >= 5)
            ? ComputeSupportAxisSpread(dense_points, support, support.x_W, secondary_axis_W)
            : 0.0;

    double side_mass = std::clamp(0.24 + 0.40 * coverage_ratio, 0.24, (emission_count >= 5) ? 0.78 : 0.68);
    if (support.dense_members <= 2) {
        side_mass = std::min(side_mass, 0.5);
    }
    const double center_mass = std::max(0.08, 1.0 - side_mass);
    weights[0] = center_mass;

    if (emission_count >= 5) {
        const double spread_total = std::max(primary_spread + secondary_spread, 1.0e-9);
        const double primary_score = 0.35 + 0.65 * (primary_spread / spread_total);
        const double secondary_score = 0.35 + 0.65 * (secondary_spread / spread_total);
        const double score_total = primary_score + secondary_score;
        const double primary_mass = side_mass * primary_score / score_total;
        const double secondary_mass = side_mass * secondary_score / score_total;
        weights[1] = 0.5 * primary_mass;
        weights[2] = 0.5 * primary_mass;
        weights[3] = 0.5 * secondary_mass;
        weights[4] = 0.5 * secondary_mass;
    } else {
        weights[1] = 0.5 * side_mass;
        weights[2] = 0.5 * side_mass;
    }

    NormalizeStencilWeights(weights, emission_count);
    return weights;
}

std::array<double, 5> ComputeTangentialStencilWeights(const std::vector<DenseContactPoint>& dense_points,
                                                      const ReducedSupportAggregate& support,
                                                      const DenseSubpatch& subpatch,
                                                      const chrono::ChVector3d& primary_axis_W,
                                                      const chrono::ChVector3d& secondary_axis_W,
                                                      int emission_count) {
    std::array<double, 5> weights{1.0, 0.0, 0.0, 0.0, 0.0};
    if (emission_count <= 1) {
        return weights;
    }

    const chrono::ChVector3d tangential_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const double tangential_force_norm = tangential_force_W.Length();
    const double total_force_norm = support.allocated_force_W.Length();
    const double tangential_ratio = tangential_force_norm / std::max(total_force_norm, 1.0e-12);
    const double coverage_ratio = support.coverage_radius / std::max(subpatch.diameter, 1.0e-9);

    chrono::ChVector3d tangential_dir_W = ProjectToTangentUnit(tangential_force_W, subpatch.avg_normal_W);
    if (tangential_dir_W.Length2() <= 0.0) {
        const chrono::ChVector3d tangential_velocity_W =
            support.v_rel_W - chrono::Vdot(support.v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
        tangential_dir_W = ProjectToTangentUnit(tangential_velocity_W, subpatch.avg_normal_W);
    }
    if (tangential_dir_W.Length2() <= 0.0) {
        tangential_dir_W = primary_axis_W;
    }

    if (emission_count == 2) {
        weights = {0.0, 0.5, 0.5, 0.0, 0.0};
        return weights;
    }

    const double primary_spread =
        ComputeSupportAxisSpread(dense_points, support, support.x_W, primary_axis_W);
    const double secondary_spread =
        (emission_count >= 5)
            ? ComputeSupportAxisSpread(dense_points, support, support.x_W, secondary_axis_W)
            : 0.0;
    const double primary_align = std::abs(chrono::Vdot(tangential_dir_W, primary_axis_W));
    const double secondary_align = std::abs(chrono::Vdot(tangential_dir_W, secondary_axis_W));

    double side_mass = std::clamp(0.55 + 0.30 * tangential_ratio + 0.15 * coverage_ratio, 0.55,
                                  (emission_count >= 5) ? 0.92 : 0.85);
    if (support.dense_members <= 2) {
        side_mass = std::min(side_mass, 0.7);
    }
    const double center_mass = std::max(0.04, 1.0 - side_mass);
    weights[0] = center_mass;

    if (emission_count >= 5) {
        const double spread_total = std::max(primary_spread + secondary_spread, 1.0e-9);
        const double primary_score =
            0.20 + 0.35 * (primary_spread / spread_total) + 0.75 * primary_align;
        const double secondary_score =
            0.20 + 0.35 * (secondary_spread / spread_total) + 0.75 * secondary_align;
        const double score_total = primary_score + secondary_score;
        const double primary_mass = side_mass * primary_score / score_total;
        const double secondary_mass = side_mass * secondary_score / score_total;
        weights[1] = 0.5 * primary_mass;
        weights[2] = 0.5 * primary_mass;
        weights[3] = 0.5 * secondary_mass;
        weights[4] = 0.5 * secondary_mass;
    } else {
        weights[1] = 0.5 * side_mass;
        weights[2] = 0.5 * side_mass;
    }

    NormalizeStencilWeights(weights, emission_count);
    return weights;
}

void SeedReactionCaches(const std::vector<DenseContactPoint>& dense_points,
                        const DenseSubpatch& subpatch,
                        const ReducedSupportAggregate& support,
                        double step_size,
                        ReducedContactPoint& reduced) {
    reduced.reaction_cache_primary.fill(0.0f);
    reduced.reaction_cache_secondary.fill(0.0f);
    reduced.reaction_cache_tertiary.fill(0.0f);
    reduced.reaction_cache_quaternary.fill(0.0f);
    reduced.reaction_cache_quinary.fill(0.0f);

    const chrono::ChVector3d primary_axis_W = ProjectToTangentUnit(reduced.stencil_axis_W, reduced.n_W);
    chrono::ChVector3d secondary_axis_W =
        ProjectToTangentUnit(reduced.stencil_axis_secondary_W, reduced.n_W);
    if (secondary_axis_W.Length2() <= 0.0) {
        secondary_axis_W = ProjectToTangentUnit(chrono::Vcross(reduced.n_W, primary_axis_W), reduced.n_W);
    }

    const auto normal_weights = ComputeNormalStencilWeights(dense_points, support, subpatch, primary_axis_W,
                                                            secondary_axis_W, reduced.emission_count);
    const auto tangential_weights = ComputeTangentialStencilWeights(dense_points, support, subpatch, primary_axis_W,
                                                                    secondary_axis_W, reduced.emission_count);
    reduced.stencil_gap_offsets.fill(0.0);
    const int active_slots =
        (reduced.emission_count >= 5) ? 5 : ((reduced.emission_count >= 3) ? 3 : ((reduced.emission_count == 2) ? 2 : 1));
    if (active_slots > 1) {
        const double uniform_weight = 1.0 / static_cast<double>(active_slots);
        const double gap_bias_scale = std::min(0.45 * std::abs(reduced.phi_eff), 3.0e-4);
        for (int slot = 0; slot < active_slots; ++slot) {
            reduced.stencil_gap_offsets[slot] = gap_bias_scale * (uniform_weight - normal_weights[slot]);
        }
    }

    const chrono::ChVector3d normal_force_W = chrono::Vdot(reduced.allocated_force_W, reduced.n_W) * reduced.n_W;
    const chrono::ChVector3d tangential_force_W = reduced.allocated_force_W - normal_force_W;
    const double seed_step = std::max(step_size, 0.0);

    auto seed_slot = [&](std::array<float, 6>& reaction_cache, int slot_index) {
        const chrono::ChVector3d seeded_force_W =
            normal_weights[slot_index] * normal_force_W + tangential_weights[slot_index] * tangential_force_W;
        EncodeReactionCacheWorldImpulse(reaction_cache, seed_step * seeded_force_W, reduced.n_W);
    };

    if (reduced.emission_count <= 1) {
        seed_slot(reduced.reaction_cache_primary, 0);
        return;
    }
    if (reduced.emission_count == 2) {
        seed_slot(reduced.reaction_cache_secondary, 1);
        seed_slot(reduced.reaction_cache_tertiary, 2);
        return;
    }

    seed_slot(reduced.reaction_cache_primary, 0);
    seed_slot(reduced.reaction_cache_secondary, 1);
    seed_slot(reduced.reaction_cache_tertiary, 2);
    if (reduced.emission_count >= 5) {
        seed_slot(reduced.reaction_cache_quaternary, 3);
        seed_slot(reduced.reaction_cache_quinary, 4);
    }
}

void BuildSupportBasis(const chrono::ChVector3d& n_W,
                       const chrono::ChVector3d& preferred_t1_W,
                       chrono::ChVector3d& t1_W,
                       chrono::ChVector3d& t2_W) {
    t1_W = preferred_t1_W - chrono::Vdot(preferred_t1_W, n_W) * n_W;
    const double t1_len = t1_W.Length();
    if (!(t1_len > 1.0e-12)) {
        const chrono::ChVector3d seed =
            (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0) : chrono::ChVector3d(1.0, 0.0, 0.0);
        t1_W = chrono::Vcross(seed, n_W);
    }
    const double t1_safe_len = t1_W.Length();
    if (t1_safe_len > 1.0e-12) {
        t1_W *= (1.0 / t1_safe_len);
    } else {
        t1_W = chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    t2_W = chrono::Vcross(n_W, t1_W);
    const double t2_len = t2_W.Length();
    if (t2_len > 1.0e-12) {
        t2_W *= (1.0 / t2_len);
    } else {
        t2_W = chrono::ChVector3d(0.0, 1.0, 0.0);
    }
    t1_W = chrono::Vcross(t2_W, n_W);
    const double t1_final_len = t1_W.Length();
    if (t1_final_len > 1.0e-12) {
        t1_W *= (1.0 / t1_final_len);
    }
}

double SubpatchMatchGate(const CompressedContactConfig& cfg,
                         const DenseSubpatch& current,
                         const TemporalSubpatchState& previous) {
    const double base =
        std::max({cfg.warm_start_match_radius, cfg.max_subpatch_diameter, cfg.max_patch_diameter,
                  0.5 * current.diameter, 0.5 * previous.diameter, 1.0e-6});
    return 1.5 * base;
}

double ComputeSubpatchMatchConfidence(const CompressedContactConfig& cfg,
                                      const DenseSubpatch& current,
                                      const TemporalSubpatchState& previous) {
    const double gate = SubpatchMatchGate(cfg, current, previous);
    const double distance = (current.centroid_W - previous.centroid_W).Length();
    if (distance > gate) {
        return 0.0;
    }

    const double normal_cos = chrono::Vdot(current.avg_normal_W, previous.avg_normal_W);
    if (normal_cos < cfg.normal_cos_min) {
        return 0.0;
    }

    const double distance_conf = Clamp01(1.0 - distance / std::max(gate, 1.0e-9));
    const double normal_conf = Clamp01((normal_cos - cfg.normal_cos_min) / std::max(1.0 - cfg.normal_cos_min, 1.0e-6));
    const double diameter_ratio =
        std::min(current.diameter, previous.diameter) / std::max({current.diameter, previous.diameter, 1.0e-9});
    return Clamp01(distance_conf * (0.25 + 0.75 * normal_conf) * std::sqrt(diameter_ratio));
}

chrono::ChVector3d ComputeAverageRelativeVelocity(const std::vector<DenseContactPoint>& dense_points,
                                                  const DenseSubpatch& subpatch) {
    chrono::ChVector3d avg_v_rel_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = std::max(ProxyLoad(point), std::max(1.0e-12, point.area_weight));
        avg_v_rel_W += weight * point.v_rel_W;
        weight_sum += weight;
    }
    if (!(weight_sum > 1.0e-12)) {
        return avg_v_rel_W;
    }
    return avg_v_rel_W * (1.0 / weight_sum);
}

double ComputeSubpatchMotionBlend(const CompressedContactConfig& cfg,
                                  const std::vector<DenseContactPoint>& dense_points,
                                  const DenseSubpatch& subpatch) {
    const chrono::ChVector3d avg_v_rel_W = ComputeAverageRelativeVelocity(dense_points, subpatch);
    const double approach_speed = std::max(0.0, -chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W));
    const double separation_speed = std::max(0.0, chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W));
    const chrono::ChVector3d tangential_v_W =
        avg_v_rel_W - chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    const double tangential_speed = tangential_v_W.Length();

    const double approach_weight =
        1.0 / (1.0 + approach_speed / std::max(cfg.temporal_approach_velocity_scale, 1.0e-6));
    const double normal_weight =
        1.0 / (1.0 + separation_speed / std::max(cfg.temporal_separation_velocity_scale, 1.0e-6));
    const double tangential_weight =
        1.0 / (1.0 + tangential_speed / std::max(cfg.temporal_slip_velocity_scale, 1.0e-6));
    return Clamp01(std::sqrt(std::max(0.0, approach_weight * normal_weight * tangential_weight)));
}

double ComputeReferenceLoadBlend(const ReferenceWrench& current, const TemporalSubpatchState& previous) {
    const double current_load = std::max(0.0, current.total_load);
    const double previous_load = std::max(0.0, previous.reference_total_load);
    const double load_ratio =
        std::min(current_load, previous_load) / std::max({current_load, previous_load, 1.0e-9});
    return std::sqrt(load_ratio);
}

double ComputeImpulseLoadBlend(const ReferenceWrench& current, const TemporalSubpatchState& previous) {
    if (!previous.has_impulse_wrench) {
        return 0.0;
    }
    const double current_load = std::max(0.0, current.total_load);
    const double previous_load = std::max(0.0, previous.impulse_total_load);
    const double load_ratio =
        std::min(current_load, previous_load) / std::max({current_load, previous_load, 1.0e-9});
    return std::sqrt(load_ratio);
}

ReferenceWrench MakeImpulseTransportReference(const ReferenceWrench& current,
                                              const TemporalSubpatchState& previous) {
    ReferenceWrench transported = current;
    transported.force_W = previous.impulse_force_W;
    transported.moment_W =
        previous.impulse_moment_W + chrono::Vcross(previous.impulse_origin_W - current.origin_W, previous.impulse_force_W);
    transported.total_load = previous.impulse_total_load;
    return transported;
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

bool BuildTransportedImpulseWarmStart(const std::vector<ReducedSupportAggregate>& supports,
                                      const DenseSubpatch& subpatch,
                                      const TemporalSubpatchState* previous,
                                      const ReferenceWrench& current_reference,
                                      const CompressedContactConfig& cfg,
                                      double mu_default,
                                      double transport_alpha,
                                      std::vector<double>& out_loads,
                                      std::vector<chrono::ChVector3d>& out_forces_W) {
    out_loads.clear();
    out_forces_W.clear();
    if (!previous || !previous->has_impulse_wrench || supports.empty() || !(transport_alpha > 0.0)) {
        return false;
    }

    std::vector<SupportWrenchPoint> support_points;
    support_points.reserve(supports.size());
    for (const auto& support : supports) {
        SupportWrenchPoint point;
        point.x_W = support.x_W;
        point.n_W = support.n_W;
        BuildSupportBasis(point.n_W, subpatch.t1_W, point.t1_W, point.t2_W);
        point.mu = mu_default;
        point.initial_load = 0.0;
        point.initial_force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
        support_points.push_back(point);
    }

    const ReferenceWrench transport_reference = MakeImpulseTransportReference(current_reference, *previous);
    WrenchAllocationResult transport_result;
    LocalWrenchAllocator::Allocate(support_points, transport_reference, std::max(1.0e-10, cfg.temporal_load_regularization),
                                   transport_result);
    if (!transport_result.feasible || transport_result.loads.size() != supports.size() ||
        transport_result.forces_W.size() != supports.size()) {
        return false;
    }

    out_loads = std::move(transport_result.loads);
    out_forces_W = std::move(transport_result.forces_W);
    return true;
}

double ComputeSupportMatchConfidence(const ReducedSupportAggregate& support,
                                     const ReducedContactPoint& previous,
                                     double match_radius) {
    const double gate = std::max(2.0 * match_radius, 1.0e-6);
    const double distance = (support.x_W - previous.x_W).Length();
    if (distance > gate) {
        return 0.0;
    }

    const double normal_cos = chrono::Vdot(support.n_W, previous.n_W);
    if (normal_cos < 0.85) {
        return 0.0;
    }

    const double distance_conf = Clamp01(1.0 - distance / gate);
    const double normal_conf = Clamp01((normal_cos - 0.85) / 0.15);
    return Clamp01(distance_conf * (0.25 + 0.75 * normal_conf));
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

SupportTransportSeed LookupPreviousTransportSeed(const ReducedSupportAggregate& support,
                                                 const std::vector<ReducedContactPoint>& previous_contacts,
                                                 double previous_step_size,
                                                 double match_radius,
                                                 std::ptrdiff_t matched_previous = -1) {
    SupportTransportSeed seed;
    seed.load = std::max(0.0, support.allocated_load);
    seed.force_W = support.allocated_force_W;

    const ReducedContactPoint* candidate = nullptr;
    if (matched_previous >= 0 && static_cast<std::size_t>(matched_previous) < previous_contacts.size()) {
        candidate = &previous_contacts[static_cast<std::size_t>(matched_previous)];
    } else {
        double best_distance = std::numeric_limits<double>::infinity();
        for (const auto& previous : previous_contacts) {
            const double distance = (previous.x_W - support.x_W).Length();
            if (distance < best_distance) {
                best_distance = distance;
                candidate = &previous;
            }
        }
    }

    if (!candidate) {
        return seed;
    }

    const double confidence = ComputeSupportMatchConfidence(support, *candidate, match_radius);
    if (!(confidence > 0.0)) {
        return seed;
    }

    const chrono::ChVector3d cached_force_W = DecodeTotalReactionEquivalentForce(*candidate, previous_step_size);
    if (cached_force_W.Length2() > 1.0e-24) {
        seed.force_W = cached_force_W;
        seed.load = std::max(0.0, chrono::Vdot(cached_force_W, support.n_W));
    } else {
        seed.load = std::max(0.0, candidate->allocated_load);
        seed.force_W = candidate->allocated_force_W;
    }
    seed.confidence = confidence;
    return seed;
}

TransportBlendWeights ComputeSupportTransportBlend(const CompressedContactConfig& cfg,
                                                   const ReducedSupportAggregate& support,
                                                   const SupportTransportSeed& previous_seed) {
    TransportBlendWeights blend;
    if (!(previous_seed.confidence > 0.0) || !(cfg.temporal_force_transport_blend > 0.0)) {
        return blend;
    }

    const double base = Clamp01(cfg.temporal_force_transport_blend * previous_seed.confidence);
    const double approach_speed = std::max(0.0, -chrono::Vdot(support.v_rel_W, support.n_W));
    const double separation_speed = std::max(0.0, chrono::Vdot(support.v_rel_W, support.n_W));
    const chrono::ChVector3d tangential_v_W =
        support.v_rel_W - chrono::Vdot(support.v_rel_W, support.n_W) * support.n_W;
    const double tangential_speed = tangential_v_W.Length();
    const double approach_weight =
        1.0 / (1.0 + approach_speed / std::max(cfg.temporal_approach_velocity_scale, 1.0e-6));
    const double normal_weight =
        1.0 / (1.0 + separation_speed / std::max(cfg.temporal_separation_velocity_scale, 1.0e-6));
    double tangential_weight =
        1.0 / (1.0 + tangential_speed / std::max(cfg.temporal_slip_velocity_scale, 1.0e-6));

    const double current_load = std::max(0.0, support.allocated_load);
    const double previous_load = std::max(0.0, previous_seed.load);
    const double load_ratio =
        std::min(current_load, previous_load) / std::max({current_load, previous_load, 1.0e-9});
    const double load_weight = std::sqrt(load_ratio);

    const chrono::ChVector3d previous_tangential_force_W =
        previous_seed.force_W - chrono::Vdot(previous_seed.force_W, support.n_W) * support.n_W;
    if (tangential_speed > 1.0e-12 && previous_tangential_force_W.Length2() > 1.0e-24) {
        const chrono::ChVector3d expected_friction_dir_W = -tangential_v_W * (1.0 / tangential_speed);
        const chrono::ChVector3d previous_force_dir_W =
            previous_tangential_force_W * (1.0 / previous_tangential_force_W.Length());
        tangential_weight *= Clamp01(0.5 * (1.0 + chrono::Vdot(expected_friction_dir_W, previous_force_dir_W)));
    }

    blend.scalar = Clamp01(base * approach_weight * normal_weight * load_weight);
    blend.normal = Clamp01(base * approach_weight * normal_weight * load_weight);
    blend.tangential = Clamp01(base * approach_weight * tangential_weight * load_weight);
    return blend;
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
        std::vector<std::size_t> member_indices;

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
            member_indices.push_back(member_index);

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
        support.allocated_force_W = sum_load * support.n_W;
        support.coverage_radius = max_coverage_radius;
        support.dense_members = count;
        support.member_indices = std::move(member_indices);
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

double ComputeSubpatchGapMismatch(const std::vector<DenseContactPoint>& dense_points,
                                  const DenseSubpatch& subpatch,
                                  const std::vector<ReducedSupportAggregate>& supports) {
    if (subpatch.members.empty() || supports.empty()) {
        return 0.0;
    }

    double dense_worst_gap = 0.0;
    for (const auto point_index : subpatch.members) {
        dense_worst_gap = std::min(dense_worst_gap, dense_points[point_index].phi_eff);
    }

    double reduced_worst_gap = 0.0;
    for (const auto& support : supports) {
        reduced_worst_gap = std::min(reduced_worst_gap, support.phi_eff);
    }

    return std::abs(reduced_worst_gap - dense_worst_gap);
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
    previous_step_size_ = 0.0;
}

void CompressedContactPipeline::SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples) {
    slave_surface_samples_ = std::move(samples);
    dense_sample_bvh_.Build(slave_surface_samples_, cfg_.bvh_leaf_size);
    previous_contacts_.clear();
    previous_subpatches_.clear();
    next_persistent_id_ = 1;
    previous_step_size_ = 0.0;
}

void CompressedContactPipeline::SyncTemporalWarmStart(const std::vector<ReducedContactPoint>& emitted_contacts) const {
    previous_contacts_ = emitted_contacts;
    for (auto& subpatch_state : previous_subpatches_) {
        std::vector<ReducedContactPoint> synced_contacts;
        for (const auto& emitted : emitted_contacts) {
            if (emitted.persistent_id == subpatch_state.persistent_id) {
                synced_contacts.push_back(emitted);
            }
        }
        subpatch_state.contacts = std::move(synced_contacts);

        subpatch_state.has_impulse_wrench = false;
        subpatch_state.impulse_origin_W = subpatch_state.reference_origin_W;
        subpatch_state.impulse_force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
        subpatch_state.impulse_moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
        subpatch_state.impulse_total_load = 0.0;
        for (const auto& contact : subpatch_state.contacts) {
            const chrono::ChVector3d reaction_force_W =
                DecodeTotalReactionEquivalentForce(contact, previous_step_size_);
            if (!(reaction_force_W.Length2() > 1.0e-24)) {
                continue;
            }
            subpatch_state.has_impulse_wrench = true;
            subpatch_state.impulse_force_W += reaction_force_W;
            subpatch_state.impulse_moment_W +=
                chrono::Vcross(contact.x_W - subpatch_state.impulse_origin_W, reaction_force_W);
            subpatch_state.impulse_total_load +=
                std::max(0.0, chrono::Vdot(reaction_force_W, subpatch_state.avg_normal_W));
        }
    }
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
    double max_subpatch_second_moment_error = 0.0;
    double max_subpatch_cone_error = 0.0;
    double max_subpatch_gap_error = 0.0;
    double max_subpatch_force_residual = 0.0;
    double max_subpatch_moment_residual = 0.0;
    double max_subpatch_reference_wrench_error = 0.0;
    double max_subpatch_reference_cop_error = 0.0;
    double max_dense_micro_force_residual = 0.0;
    double max_dense_micro_moment_residual = 0.0;

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

        int target_points = std::min(cfg_.max_reduced_points_per_patch,
                                     (subpatch.members.size() <= 1) ? 1 : ((subpatch.members.size() <= 2) ? 2 : 3));
        DenseMicroReferenceResult dense_micro_reference;
        LocalWrenchAllocator::BuildDenseMicroReference(dense_points, subpatch.members, subpatch.centroid_W,
                                                       mu_default, step_size, dense_micro_reference);
        const ReferenceWrench& dense_reference = dense_micro_reference.reference;
        max_dense_micro_force_residual =
            std::max(max_dense_micro_force_residual, dense_micro_reference.force_residual);
        max_dense_micro_moment_residual =
            std::max(max_dense_micro_moment_residual, dense_micro_reference.moment_residual);
        double temporal_reference_alpha = 0.0;
        if (matched_previous_state) {
            temporal_reference_alpha = cfg_.temporal_reference_blend *
                                       ComputeSubpatchMatchConfidence(cfg_, subpatch, *matched_previous_state) *
                                       ComputeSubpatchMotionBlend(cfg_, dense_points, subpatch) *
                                       ComputeReferenceLoadBlend(dense_reference, *matched_previous_state);
        }
        const ReferenceWrench reference =
            BlendTemporalReference(dense_reference, matched_previous_state, temporal_reference_alpha);
        std::vector<ReducedSupportAggregate> supports;
        std::vector<std::ptrdiff_t> matched_supports;
        std::vector<std::size_t> support_ids;
        std::vector<std::size_t> support_order;
        double subpatch_second_moment_error = 0.0;
        double subpatch_cone_error = 0.0;
        double subpatch_gap_error = 0.0;
        double subpatch_gap_mismatch = 0.0;
        double subpatch_reference_wrench_error = 0.0;
        double subpatch_reference_cop_error = 0.0;
        while (true) {
            supports = BuildReducedSupportsForSubpatch(dense_points, subpatch, selection_previous_contacts,
                                                       cfg_.warm_start_match_radius, target_points);
            matched_supports =
                MatchSupportsToPrevious(supports, previous_local_contacts, cfg_.warm_start_match_radius);
            support_ids = AssignSupportIds(matched_supports, previous_local_contacts);
            support_order.resize(supports.size());
            std::iota(support_order.begin(), support_order.end(), 0);
            std::sort(support_order.begin(), support_order.end(),
                      [&support_ids](std::size_t a, std::size_t b) { return support_ids[a] < support_ids[b]; });

            if (!supports.empty()) {
                const double transport_alpha =
                    (matched_previous_state && matched_previous_state->has_impulse_wrench)
                        ? (cfg_.temporal_force_transport_blend *
                           ComputeSubpatchMatchConfidence(cfg_, subpatch, *matched_previous_state) *
                           ComputeSubpatchMotionBlend(cfg_, dense_points, subpatch) *
                           ComputeImpulseLoadBlend(reference, *matched_previous_state))
                        : 0.0;
                std::vector<double> transported_loads;
                std::vector<chrono::ChVector3d> transported_forces_W;
                const bool has_transport = BuildTransportedImpulseWarmStart(
                    supports, subpatch, matched_previous_state, reference, cfg_, mu_default, transport_alpha,
                    transported_loads, transported_forces_W);

                std::vector<SupportWrenchPoint> support_points;
                support_points.reserve(supports.size());
                for (std::size_t support_index = 0; support_index < supports.size(); ++support_index) {
                    const auto& support = supports[support_index];
                    SupportWrenchPoint point;
                    point.x_W = support.x_W;
                    point.n_W = support.n_W;
                    BuildSupportBasis(point.n_W, subpatch.t1_W, point.t1_W, point.t2_W);
                    point.mu = mu_default;
                    if (has_transport && support_index < transported_loads.size() &&
                        support_index < transported_forces_W.size()) {
                        point.initial_load = (1.0 - transport_alpha) * std::max(0.0, support.allocated_load) +
                                             transport_alpha * std::max(0.0, transported_loads[support_index]);
                        point.initial_force_W =
                            (1.0 - transport_alpha) * support.allocated_force_W +
                            transport_alpha * transported_forces_W[support_index];
                    } else {
                        const auto previous_seed = LookupPreviousTransportSeed(
                            support, previous_local_contacts, previous_step_size_, cfg_.warm_start_match_radius,
                            matched_supports[support_index]);
                        const auto transport_blend = ComputeSupportTransportBlend(cfg_, support, previous_seed);
                        const chrono::ChVector3d current_normal_force_W =
                            chrono::Vdot(support.allocated_force_W, support.n_W) * support.n_W;
                        const chrono::ChVector3d previous_normal_force_W =
                            chrono::Vdot(previous_seed.force_W, support.n_W) * support.n_W;
                        const chrono::ChVector3d current_tangential_force_W =
                            support.allocated_force_W - current_normal_force_W;
                        const chrono::ChVector3d previous_tangential_force_W =
                            previous_seed.force_W - previous_normal_force_W;

                        point.initial_load = (1.0 - transport_blend.scalar) * std::max(0.0, support.allocated_load) +
                                             transport_blend.scalar * std::max(0.0, previous_seed.load);
                        point.initial_force_W =
                            (1.0 - transport_blend.normal) * current_normal_force_W +
                            transport_blend.normal * previous_normal_force_W +
                            (1.0 - transport_blend.tangential) * current_tangential_force_W +
                            transport_blend.tangential * previous_tangential_force_W;
                    }
                    support_points.push_back(point);
                }

                WrenchAllocationResult allocation;
                LocalWrenchAllocator::Allocate(support_points, reference, cfg_.temporal_load_regularization,
                                               allocation);
                max_subpatch_force_residual = std::max(max_subpatch_force_residual, allocation.force_residual);
                max_subpatch_moment_residual = std::max(max_subpatch_moment_residual, allocation.moment_residual);

                if (allocation.feasible && allocation.loads.size() == supports.size()) {
                    for (std::size_t i = 0; i < supports.size(); ++i) {
                        supports[i].allocated_load = allocation.loads[i];
                        if (allocation.forces_W.size() == supports.size()) {
                            supports[i].allocated_force_W = allocation.forces_W[i];
                        } else {
                            supports[i].allocated_force_W = allocation.loads[i] * supports[i].n_W;
                        }
                        if (supports[i].phi_eff < -1.0e-12) {
                            supports[i].support_weight =
                                allocation.loads[i] / std::max(1.0e-12, -supports[i].phi_eff);
                        } else {
                            supports[i].support_weight = supports[i].area_weight;
                        }
                    }
                }
            }

            subpatch_second_moment_error = ComputeSecondMomentError(dense_points, subpatch, supports);
            subpatch_cone_error = ComputeConeSupportError(dense_points, subpatch, supports, cfg_);
            subpatch_gap_error = ComputeSubpatchGapError(subpatch, supports, master_state, sdf, cfg_);
            subpatch_gap_mismatch = ComputeSubpatchGapMismatch(dense_points, subpatch, supports);
            subpatch_reference_wrench_error = ComputeReferenceWrenchError(supports, subpatch, reference);
            subpatch_reference_cop_error = ComputeReferenceCoPError(supports, subpatch, reference);
            const double combined_gap_metric = std::max(subpatch_gap_error, subpatch_gap_mismatch);
            const bool sigma_ok =
                !(cfg_.max_second_moment_error > 0.0) || subpatch_second_moment_error <= cfg_.max_second_moment_error;
            const bool cone_ok = !(cfg_.max_cone_error > 0.0) || subpatch_cone_error <= cfg_.max_cone_error;
            const bool gap_ok = !(cfg_.max_gap_error > 0.0) || combined_gap_metric <= cfg_.max_gap_error;
            const bool wrench_ok =
                !(cfg_.max_wrench_error > 0.0) || subpatch_reference_wrench_error <= cfg_.max_wrench_error;
            const bool cop_ok =
                !(cfg_.max_cop_error > 0.0) || subpatch_reference_cop_error <= cfg_.max_cop_error;
            if ((sigma_ok && cone_ok && gap_ok && wrench_ok && cop_ok) ||
                target_points >= cfg_.max_reduced_points_per_patch) {
                break;
            }
            ++target_points;
        }

        double mean_allocated_load = 0.0;
        for (const auto& support : supports) {
            mean_allocated_load += std::max(0.0, support.allocated_load);
        }
        if (!supports.empty()) {
            mean_allocated_load /= static_cast<double>(supports.size());
        }

        max_subpatch_second_moment_error =
            std::max(max_subpatch_second_moment_error, subpatch_second_moment_error);
        max_subpatch_cone_error = std::max(max_subpatch_cone_error, subpatch_cone_error);
        max_subpatch_gap_error = std::max(max_subpatch_gap_error, std::max(subpatch_gap_error, subpatch_gap_mismatch));
        max_subpatch_reference_wrench_error =
            std::max(max_subpatch_reference_wrench_error, subpatch_reference_wrench_error);
        max_subpatch_reference_cop_error =
            std::max(max_subpatch_reference_cop_error, subpatch_reference_cop_error);

        const std::size_t contacts_begin = out_contacts.size();
        for (const auto support_index : support_order) {
            const auto& support = supports[support_index];
            ReducedContactPoint reduced;
            reduced.persistent_id = persistent_id;
            reduced.patch_id = support.patch_id;
            reduced.subpatch_id = support.subpatch_id;
            reduced.support_id = support_ids[support_index];
            reduced.dense_members = support.dense_members;
            reduced.emission_count = ChooseEmissionCount(support, subpatch, mean_allocated_load);
            reduced.x_W = support.x_W;
            reduced.x_master_M = support.x_master_M;
            reduced.x_master_surface_W = support.x_master_surface_W;
            reduced.n_W = support.n_W;
            reduced.v_rel_W = support.v_rel_W;
            reduced.stencil_axis_W = ChooseStencilAxis(support, subpatch);
            reduced.stencil_axis_secondary_W =
                ChooseSecondaryStencilAxis(support, subpatch, reduced.stencil_axis_W);
            reduced.phi = support.phi;
            reduced.phi_eff = support.phi_eff;
            reduced.area_weight = support.area_weight;
            reduced.support_weight = support.support_weight;
            reduced.allocated_load = support.allocated_load;
            reduced.allocated_force_W = support.allocated_force_W;
            reduced.stencil_half_extent = ChooseStencilHalfExtent(support, subpatch, reduced.emission_count);
            reduced.stencil_half_extent_secondary =
                ChooseSecondaryStencilHalfExtent(support, subpatch, reduced.emission_count);
            reduced.mu = mu_default;
            SeedReactionCaches(dense_points, subpatch, support, step_size, reduced);
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
        temporal_state.impulse_origin_W = reference.origin_W;
        temporal_state.impulse_force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
        temporal_state.impulse_moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
        temporal_state.impulse_total_load = 0.0;
        temporal_state.has_impulse_wrench = false;
        for (std::size_t contact_index = contacts_begin; contact_index < out_contacts.size(); ++contact_index) {
            temporal_state.contacts.push_back(out_contacts[contact_index]);
        }
        current_subpatches.push_back(std::move(temporal_state));
    }

    previous_contacts_ = out_contacts;
    previous_subpatches_ = std::move(current_subpatches);
    previous_step_size_ = step_size;

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
    out_stats->max_subpatch_second_moment_error = max_subpatch_second_moment_error;
    out_stats->max_subpatch_cone_error = max_subpatch_cone_error;
    out_stats->max_subpatch_gap_error = max_subpatch_gap_error;
    out_stats->max_subpatch_force_residual = max_subpatch_force_residual;
    out_stats->max_subpatch_moment_residual = max_subpatch_moment_residual;
    out_stats->max_subpatch_reference_wrench_error = max_subpatch_reference_wrench_error;
    out_stats->max_subpatch_reference_cop_error = max_subpatch_reference_cop_error;
    out_stats->max_dense_micro_force_residual = max_dense_micro_force_residual;
    out_stats->max_dense_micro_moment_residual = max_dense_micro_moment_residual;

    if (dense_points.empty() || out_contacts.empty()) {
        out_stats->epsilon_F = 0.0;
        out_stats->epsilon_M = 0.0;
        out_stats->epsilon_CoP = 0.0;
        out_stats->epsilon_gap = max_subpatch_gap_error;
        out_stats->dense_worst_gap = 0.0;
        out_stats->reduced_worst_gap = 0.0;
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
    double dense_worst_gap = 0.0;
    for (const auto& dense : dense_points) {
        dense_worst_gap = std::min(dense_worst_gap, dense.phi_eff);
    }
    for (const auto& reduced : out_contacts) {
        reduced_worst_gap = std::min(reduced_worst_gap, reduced.phi_eff);
    }
    const double worst_gap_mismatch = std::abs(reduced_worst_gap - dense_worst_gap);
    out_stats->epsilon_gap = std::max(worst_gap_mismatch, max_subpatch_gap_error);
    out_stats->dense_worst_gap = dense_worst_gap;
    out_stats->reduced_worst_gap = reduced_worst_gap;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
