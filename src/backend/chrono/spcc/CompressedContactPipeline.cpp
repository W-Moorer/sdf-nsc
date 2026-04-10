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
    double slip_dominance = 0.0;
    double tangential_heterogeneity = 0.0;
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

struct PlanarStats2D {
    double mean_u = 0.0;
    double mean_v = 0.0;
    double cov_uu = 0.0;
    double cov_uv = 0.0;
    double cov_vv = 0.0;
    double weight_sum = 0.0;
};

struct DirectionalSupportTargets {
    chrono::ChVector3d center_W;
    std::vector<chrono::ChVector3d> directions_W;
    std::vector<double> dense_support;
    std::vector<double> direction_weights;
};

struct ClusterDirectionalEnvelope {
    std::vector<double> positive_support;
};

struct TangentialDistributionStats {
    chrono::ChVector3d dominant_tangent_W;
    double slip_dominance = 0.0;
    double angular_spread = 0.0;
    double leverage_spread = 0.0;
    double heterogeneity = 0.0;
};

double Clamp01(double value) {
    return std::clamp(value, 0.0, 1.0);
}

double ClampScale(double value, double min_value, double max_value) {
    return std::clamp(value, min_value, max_value);
}

chrono::ChVector3d SafeNormalized(const chrono::ChVector3d& v_W, const chrono::ChVector3d& fallback_W) {
    const double len = v_W.Length();
    if (!(len > 1.0e-12)) {
        return fallback_W;
    }
    return v_W * (1.0 / len);
}

chrono::ChVector3d ProjectToTangentUnit(const chrono::ChVector3d& v_W, const chrono::ChVector3d& n_W);

DenseMicroSolverOptions MakeDenseMicroSolverOptions(const CompressedContactConfig& cfg) {
    DenseMicroSolverOptions options;
    options.friction_ray_count = cfg.dense_micro_friction_rays;
    options.normal_response_weight = cfg.dense_micro_normal_response_weight;
    options.tangential_response_weight = cfg.dense_micro_tangential_response_weight;
    options.gap_drive_weight = cfg.dense_micro_gap_drive_weight;
    options.approach_drive_weight = cfg.dense_micro_approach_drive_weight;
    options.slip_drive_weight = cfg.dense_micro_slip_drive_weight;
    options.wrench_coupling_weight = cfg.dense_micro_wrench_coupling_weight;
    options.regularization = cfg.dense_micro_regularization;
    return options;
}

ReducedSolveOptions MakeReducedSolveOptions(const CompressedContactConfig& cfg, double temporal_regularization) {
    ReducedSolveOptions options;
    options.friction_ray_count = cfg.reduced_friction_rays;
    options.temporal_regularization = temporal_regularization;
    return options;
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

chrono::ChVector3d DenseTangentialDirection(const DenseContactPoint& point,
                                            const chrono::ChVector3d& normal_W,
                                            const chrono::ChVector3d& fallback_tangent_W) {
    chrono::ChVector3d tangent_W =
        point.v_rel_W - chrono::Vdot(point.v_rel_W, normal_W) * normal_W;
    const double tangent_len = tangent_W.Length();
    if (tangent_len > 1.0e-12) {
        return tangent_W * (1.0 / tangent_len);
    }
    tangent_W = ProjectToTangentUnit(fallback_tangent_W, normal_W);
    if (tangent_W.Length2() > 0.0) {
        return tangent_W;
    }
    return chrono::ChVector3d(0.0, 0.0, 0.0);
}

TangentialDistributionStats ComputeTangentialDistributionStats(
    const std::vector<DenseContactPoint>& dense_points,
    const std::vector<std::size_t>& member_indices,
    const chrono::ChVector3d& normal_W,
    const chrono::ChVector3d& center_W,
    const chrono::ChVector3d& fallback_tangent_W,
    double diameter) {
    TangentialDistributionStats stats;
    chrono::ChVector3d tangent_sum_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    double tangential_speed_sum = 0.0;
    double approach_speed_sum = 0.0;
    double weighted_mean_proj = 0.0;
    std::vector<std::pair<double, double>> weighted_proj;
    weighted_proj.reserve(member_indices.size());

    for (const auto point_index : member_indices) {
        const auto& point = dense_points[point_index];
        const chrono::ChVector3d tangential_velocity_W =
            point.v_rel_W - chrono::Vdot(point.v_rel_W, normal_W) * normal_W;
        const double tangential_speed = tangential_velocity_W.Length();
        const double approach_speed = std::max(0.0, -chrono::Vdot(point.v_rel_W, normal_W));
        const double weight = std::max(ProxyLoad(point), std::max(1.0e-12, point.area_weight));
        const chrono::ChVector3d tangential_dir_W =
            DenseTangentialDirection(point, normal_W, fallback_tangent_W);
        tangent_sum_W += weight * tangential_dir_W;
        tangential_speed_sum += weight * tangential_speed;
        approach_speed_sum += weight * approach_speed;
        weight_sum += weight;
    }

    if (!(weight_sum > 1.0e-12)) {
        return stats;
    }

    stats.dominant_tangent_W = SafeNormalized(tangent_sum_W, ProjectToTangentUnit(fallback_tangent_W, normal_W));
    stats.slip_dominance =
        tangential_speed_sum / std::max(tangential_speed_sum + 1.5 * approach_speed_sum, 1.0e-12);

    double angular_sum = 0.0;
    for (const auto point_index : member_indices) {
        const auto& point = dense_points[point_index];
        const double weight = std::max(ProxyLoad(point), std::max(1.0e-12, point.area_weight));
        const chrono::ChVector3d tangential_dir_W =
            DenseTangentialDirection(point, normal_W, stats.dominant_tangent_W);
        if (tangential_dir_W.Length2() > 0.0 && stats.dominant_tangent_W.Length2() > 0.0) {
            const double align = std::clamp(chrono::Vdot(tangential_dir_W, stats.dominant_tangent_W), -1.0, 1.0);
            angular_sum += weight * (1.0 - align);
        }
        const double projected = chrono::Vdot(point.x_W - center_W, stats.dominant_tangent_W);
        weighted_mean_proj += weight * projected;
        weighted_proj.emplace_back(weight, projected);
    }
    weighted_mean_proj /= weight_sum;

    double variance_proj = 0.0;
    for (const auto& entry : weighted_proj) {
        const double delta = entry.second - weighted_mean_proj;
        variance_proj += entry.first * delta * delta;
    }
    variance_proj /= weight_sum;

    stats.angular_spread = angular_sum / std::max(weight_sum, 1.0e-12);
    stats.leverage_spread =
        std::sqrt(std::max(variance_proj, 0.0)) / std::max(0.5 * diameter, 1.0e-9);
    stats.heterogeneity =
        Clamp01(0.65 * stats.slip_dominance * std::min(stats.leverage_spread, 1.5) +
                0.35 * std::min(stats.angular_spread, 1.0));
    return stats;
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

void ProjectToSubpatchUV(const DenseSubpatch& subpatch,
                         const chrono::ChVector3d& x_W,
                         double& u,
                         double& v) {
    const chrono::ChVector3d rel = x_W - subpatch.centroid_W;
    u = chrono::Vdot(rel, subpatch.t1_W);
    v = chrono::Vdot(rel, subpatch.t2_W);
}

chrono::ChVector3d LiftSubpatchUV(const DenseSubpatch& subpatch, double u, double v) {
    return subpatch.centroid_W + u * subpatch.t1_W + v * subpatch.t2_W;
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
        DecodeReactionCacheWorldForce(contact.reaction_cache_primary, contact.slot_n_W[0]);
    const chrono::ChVector3d negative_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_secondary, contact.slot_n_W[1]);
    const chrono::ChVector3d positive_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_tertiary, contact.slot_n_W[2]);
    const chrono::ChVector3d secondary_negative_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_quaternary, contact.slot_n_W[3]);
    const chrono::ChVector3d secondary_positive_force_W =
        DecodeReactionCacheWorldForce(contact.reaction_cache_quinary, contact.slot_n_W[4]);

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

std::vector<std::size_t> AssignDensePointsToCenters(const std::vector<DenseContactPoint>& dense_points,
                                                    const DenseSubpatch& subpatch,
                                                    const std::vector<chrono::ChVector3d>& centers_W) {
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
    return assignments;
}

PlanarStats2D ComputeDensePlanarStats(const std::vector<DenseContactPoint>& dense_points,
                                      const DenseSubpatch& subpatch) {
    PlanarStats2D stats;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = DenseMetricWeight(point);
        double u = 0.0;
        double v = 0.0;
        ProjectToSubpatchUV(subpatch, point.x_W, u, v);
        stats.mean_u += weight * u;
        stats.mean_v += weight * v;
        stats.weight_sum += weight;
    }
    if (!(stats.weight_sum > 1.0e-12)) {
        return stats;
    }

    stats.mean_u /= stats.weight_sum;
    stats.mean_v /= stats.weight_sum;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = DenseMetricWeight(point);
        double u = 0.0;
        double v = 0.0;
        ProjectToSubpatchUV(subpatch, point.x_W, u, v);
        const double du = u - stats.mean_u;
        const double dv = v - stats.mean_v;
        stats.cov_uu += weight * du * du;
        stats.cov_uv += weight * du * dv;
        stats.cov_vv += weight * dv * dv;
    }
    stats.cov_uu /= stats.weight_sum;
    stats.cov_uv /= stats.weight_sum;
    stats.cov_vv /= stats.weight_sum;
    return stats;
}

PlanarStats2D ComputeSupportPlanarStats(const DenseSubpatch& subpatch,
                                        const std::vector<chrono::ChVector3d>& centers_W,
                                        const std::vector<std::size_t>& assignments,
                                        const std::vector<DenseContactPoint>& dense_points) {
    PlanarStats2D stats;
    if (centers_W.empty()) {
        return stats;
    }

    std::vector<double> cluster_weights(centers_W.size(), 0.0);
    for (std::size_t local_index = 0; local_index < assignments.size(); ++local_index) {
        const std::size_t cluster = assignments[local_index];
        const auto member_index = subpatch.members[local_index];
        cluster_weights[cluster] += DenseMetricWeight(dense_points[member_index]);
    }

    for (std::size_t cluster = 0; cluster < centers_W.size(); ++cluster) {
        const double weight = std::max(cluster_weights[cluster], 1.0e-12);
        double u = 0.0;
        double v = 0.0;
        ProjectToSubpatchUV(subpatch, centers_W[cluster], u, v);
        stats.mean_u += weight * u;
        stats.mean_v += weight * v;
        stats.weight_sum += weight;
    }
    if (!(stats.weight_sum > 1.0e-12)) {
        return stats;
    }

    stats.mean_u /= stats.weight_sum;
    stats.mean_v /= stats.weight_sum;
    for (std::size_t cluster = 0; cluster < centers_W.size(); ++cluster) {
        const double weight = std::max(cluster_weights[cluster], 1.0e-12);
        double u = 0.0;
        double v = 0.0;
        ProjectToSubpatchUV(subpatch, centers_W[cluster], u, v);
        const double du = u - stats.mean_u;
        const double dv = v - stats.mean_v;
        stats.cov_uu += weight * du * du;
        stats.cov_uv += weight * du * dv;
        stats.cov_vv += weight * dv * dv;
    }
    stats.cov_uu /= stats.weight_sum;
    stats.cov_uv /= stats.weight_sum;
    stats.cov_vv /= stats.weight_sum;
    return stats;
}

double ComputePlanarMomentObjective(const PlanarStats2D& dense_stats, const PlanarStats2D& support_stats) {
    if (!(dense_stats.weight_sum > 1.0e-12) || !(support_stats.weight_sum > 1.0e-12)) {
        return std::numeric_limits<double>::infinity();
    }

    const double mean_du = support_stats.mean_u - dense_stats.mean_u;
    const double mean_dv = support_stats.mean_v - dense_stats.mean_v;
    const double dense_cov_norm = std::sqrt(dense_stats.cov_uu * dense_stats.cov_uu +
                                            2.0 * dense_stats.cov_uv * dense_stats.cov_uv +
                                            dense_stats.cov_vv * dense_stats.cov_vv);
    const double cov_scale = std::max(dense_cov_norm, 1.0e-9);
    const double cov_duu = (support_stats.cov_uu - dense_stats.cov_uu) / cov_scale;
    const double cov_duv = (support_stats.cov_uv - dense_stats.cov_uv) / cov_scale;
    const double cov_dvv = (support_stats.cov_vv - dense_stats.cov_vv) / cov_scale;
    return 4.0 * (mean_du * mean_du + mean_dv * mean_dv) +
           (cov_duu * cov_duu + 2.0 * cov_duv * cov_duv + cov_dvv * cov_dvv);
}

DirectionalSupportTargets BuildDirectionalSupportTargets(const std::vector<DenseContactPoint>& dense_points,
                                                         const DenseSubpatch& subpatch,
                                                         int direction_count) {
    DirectionalSupportTargets targets;
    targets.center_W = SubpatchReferenceCenter(dense_points, subpatch);
    const int dir_count = std::max(8, direction_count);
    targets.directions_W.reserve(dir_count);
    targets.dense_support.reserve(dir_count);

    for (int dir_index = 0; dir_index < dir_count; ++dir_index) {
        const double theta = (2.0 * std::acos(-1.0) * static_cast<double>(dir_index)) / static_cast<double>(dir_count);
        const chrono::ChVector3d dir_W = std::cos(theta) * subpatch.t1_W + std::sin(theta) * subpatch.t2_W;
        double support_value = -std::numeric_limits<double>::infinity();
        for (const auto point_index : subpatch.members) {
            support_value = std::max(support_value, chrono::Vdot(dense_points[point_index].x_W - targets.center_W, dir_W));
        }
        targets.directions_W.push_back(dir_W);
        targets.dense_support.push_back(support_value);
        targets.direction_weights.push_back(1.0);
    }

    return targets;
}

chrono::ChVector3d ComputePrincipalPlanarAxis(const PlanarStats2D& stats, const DenseSubpatch& subpatch) {
    const double theta = 0.5 * std::atan2(2.0 * stats.cov_uv, stats.cov_uu - stats.cov_vv);
    chrono::ChVector3d axis_W = std::cos(theta) * subpatch.t1_W + std::sin(theta) * subpatch.t2_W;
    axis_W = ProjectToTangentUnit(axis_W, subpatch.avg_normal_W);
    if (axis_W.Length2() > 0.0) {
        return axis_W;
    }
    return subpatch.t1_W;
}

DirectionalSupportTargets BuildReinjectionDirectionalTargets(const std::vector<DenseContactPoint>& dense_points,
                                                             const DenseSubpatch& subpatch,
                                                             const PlanarStats2D& dense_stats,
                                                             const chrono::ChVector3d& avg_v_rel_W) {
    DirectionalSupportTargets targets;
    const chrono::ChVector3d tangential_v_W =
        avg_v_rel_W - chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    const double tangential_speed = tangential_v_W.Length();
    if (!(tangential_speed > 1.0e-6)) {
        return targets;
    }

    const double approach_speed = std::max(0.0, -chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W));
    const double slip_dominance =
        tangential_speed / std::max(tangential_speed + 1.25 * approach_speed, 1.0e-9);

    chrono::ChVector3d primary_axis_W = tangential_v_W * (1.0 / tangential_speed);
    primary_axis_W = ProjectToTangentUnit(primary_axis_W, subpatch.avg_normal_W);
    if (primary_axis_W.Length2() <= 0.0) {
        primary_axis_W = ComputePrincipalPlanarAxis(dense_stats, subpatch);
    }
    chrono::ChVector3d secondary_axis_W =
        ProjectToTangentUnit(chrono::Vcross(subpatch.avg_normal_W, primary_axis_W), subpatch.avg_normal_W);
    if (secondary_axis_W.Length2() <= 0.0) {
        secondary_axis_W = ProjectToTangentUnit(subpatch.t2_W, subpatch.avg_normal_W);
    }

    targets.center_W = SubpatchReferenceCenter(dense_points, subpatch);
    const double primary_weight = 1.0 + 0.85 * slip_dominance;
    const double secondary_weight = 0.30 + 0.35 * (1.0 - slip_dominance);
    const std::array<chrono::ChVector3d, 4> directions_W = {
        primary_axis_W, -primary_axis_W, secondary_axis_W, -secondary_axis_W};
    const std::array<double, 4> direction_weights = {
        primary_weight, primary_weight, secondary_weight, secondary_weight};

    for (std::size_t dir_index = 0; dir_index < directions_W.size(); ++dir_index) {
        const auto& dir_W = directions_W[dir_index];
        if (dir_W.Length2() <= 0.0) {
            continue;
        }
        double support_value = -std::numeric_limits<double>::infinity();
        for (const auto point_index : subpatch.members) {
            support_value = std::max(support_value, chrono::Vdot(dense_points[point_index].x_W - targets.center_W, dir_W));
        }
        targets.directions_W.push_back(dir_W);
        targets.dense_support.push_back(support_value);
        targets.direction_weights.push_back(direction_weights[dir_index]);
    }

    return targets;
}

double ComputeCenterConeObjective(const DirectionalSupportTargets& targets,
                                  const std::vector<chrono::ChVector3d>& centers_W,
                                  double diameter) {
    if (centers_W.empty() || targets.directions_W.empty()) {
        return 0.0;
    }

    double sum_sq = 0.0;
    double max_abs = 0.0;
    for (std::size_t dir_index = 0; dir_index < targets.directions_W.size(); ++dir_index) {
        double reduced_support = -std::numeric_limits<double>::infinity();
        for (const auto& center_W : centers_W) {
            reduced_support =
                std::max(reduced_support, chrono::Vdot(center_W - targets.center_W, targets.directions_W[dir_index]));
        }
        const double error = (reduced_support - targets.dense_support[dir_index]) / std::max(diameter, 1.0e-9);
        sum_sq += error * error;
        max_abs = std::max(max_abs, std::abs(error));
    }

    const double rms = std::sqrt(sum_sq / static_cast<double>(targets.directions_W.size()));
    return 0.65 * rms + 0.35 * max_abs;
}

double ComputeCenterSeparationPenalty(const DenseSubpatch& subpatch,
                                      const std::vector<chrono::ChVector3d>& centers_W) {
    if (centers_W.size() <= 1) {
        return 0.0;
    }

    const double min_target = 0.12 * std::max(subpatch.diameter, 1.0e-6);
    double penalty = 0.0;
    std::size_t pair_count = 0;
    for (std::size_t i = 0; i < centers_W.size(); ++i) {
        for (std::size_t j = i + 1; j < centers_W.size(); ++j) {
            const double distance = (centers_W[i] - centers_W[j]).Length();
            if (distance < min_target) {
                const double deficit = (min_target - distance) / min_target;
                penalty += deficit * deficit;
            }
            ++pair_count;
        }
    }
    if (pair_count == 0) {
        return 0.0;
    }
    return penalty / static_cast<double>(pair_count);
}

std::vector<ClusterDirectionalEnvelope> BuildClusterDirectionalEnvelopes(const std::vector<DenseContactPoint>& dense_points,
                                                                         const DenseSubpatch& subpatch,
                                                                         const std::vector<chrono::ChVector3d>& centers_W,
                                                                         const std::vector<std::size_t>& assignments,
                                                                         const DirectionalSupportTargets& targets) {
    std::vector<ClusterDirectionalEnvelope> envelopes(centers_W.size());
    for (auto& envelope : envelopes) {
        envelope.positive_support.assign(targets.directions_W.size(), 0.0);
    }

    for (std::size_t local_index = 0; local_index < assignments.size(); ++local_index) {
        const std::size_t cluster = assignments[local_index];
        const auto point_index = subpatch.members[local_index];
        const auto& point = dense_points[point_index];
        for (std::size_t dir_index = 0; dir_index < targets.directions_W.size(); ++dir_index) {
            const double projected = chrono::Vdot(point.x_W - centers_W[cluster], targets.directions_W[dir_index]);
            envelopes[cluster].positive_support[dir_index] =
                std::max(envelopes[cluster].positive_support[dir_index], projected);
        }
    }

    return envelopes;
}

double ComputeReinjectionDirectionalReach(const DenseSubpatch& subpatch,
                                          const chrono::ChVector3d& reference_center_W,
                                          const chrono::ChVector3d& center_W,
                                          const chrono::ChVector3d& direction_W,
                                          double local_support_radius,
                                          std::size_t center_count) {
    const double global_cap = std::max(0.18 * subpatch.diameter, 1.0e-6);
    const double local_cap = std::min(global_cap, 0.82 * local_support_radius);
    const double center_scale = (center_count <= 1) ? 0.95 : ((center_count == 2) ? 0.82 : 0.72);
    const double tangential_extent = center_scale * local_cap;
    return chrono::Vdot(center_W - reference_center_W, direction_W) + tangential_extent;
}

double ComputeReinjectionObjective(const DenseSubpatch& subpatch,
                                   const DirectionalSupportTargets& targets,
                                   const std::vector<chrono::ChVector3d>& centers_W,
                                   const std::vector<ClusterDirectionalEnvelope>& envelopes) {
    if (centers_W.empty() || targets.directions_W.empty() || envelopes.size() != centers_W.size()) {
        return 0.0;
    }

    double weighted_sum = 0.0;
    double weight_sum = 0.0;
    double max_under = 0.0;
    for (std::size_t dir_index = 0; dir_index < targets.directions_W.size(); ++dir_index) {
        double reduced_reach = -std::numeric_limits<double>::infinity();
        for (std::size_t center_index = 0; center_index < centers_W.size(); ++center_index) {
            reduced_reach = std::max(reduced_reach,
                                     ComputeReinjectionDirectionalReach(subpatch, targets.center_W,
                                                                        centers_W[center_index],
                                                                        targets.directions_W[dir_index],
                                                                        envelopes[center_index].positive_support[dir_index],
                                                                        centers_W.size()));
        }

        const double normalized_error =
            (targets.dense_support[dir_index] - reduced_reach) / std::max(subpatch.diameter, 1.0e-9);
        const double under = std::max(0.0, normalized_error);
        const double over = std::max(0.0, -normalized_error);
        const double direction_weight =
            (dir_index < targets.direction_weights.size()) ? targets.direction_weights[dir_index] : 1.0;
        weighted_sum += direction_weight * (2.5 * under * under + 0.35 * over * over);
        weight_sum += direction_weight;
        max_under = std::max(max_under, under);
    }

    const double rms = std::sqrt(weighted_sum / std::max(weight_sum, 1.0e-12));
    return 0.7 * rms + 0.3 * max_under;
}

double ComputeSupportGeometryObjective(const PlanarStats2D& dense_stats,
                                       const PlanarStats2D& support_stats,
                                       const DirectionalSupportTargets& targets,
                                       const DenseSubpatch& subpatch,
                                       const std::vector<chrono::ChVector3d>& centers_W) {
    const double cone_objective = ComputeCenterConeObjective(targets, centers_W, subpatch.diameter);
    const double moment_objective = ComputePlanarMomentObjective(dense_stats, support_stats);
    const double separation_penalty = ComputeCenterSeparationPenalty(subpatch, centers_W);
    return 2.75 * cone_objective + 0.55 * moment_objective + 0.20 * separation_penalty;
}

void ApplyConeDrivenCenterUpdate(const DirectionalSupportTargets& targets,
                                 const DenseSubpatch& subpatch,
                                 double blend,
                                 std::vector<chrono::ChVector3d>& centers_W) {
    if (!(blend > 0.0) || centers_W.size() <= 1 || targets.directions_W.empty()) {
        return;
    }

    std::vector<chrono::ChVector3d> deltas_W(centers_W.size(), chrono::ChVector3d(0.0, 0.0, 0.0));
    std::vector<double> weights(centers_W.size(), 0.0);
    for (std::size_t dir_index = 0; dir_index < targets.directions_W.size(); ++dir_index) {
        std::size_t best_center = 0;
        double best_projection = -std::numeric_limits<double>::infinity();
        for (std::size_t center_index = 0; center_index < centers_W.size(); ++center_index) {
            const double projection =
                chrono::Vdot(centers_W[center_index] - targets.center_W, targets.directions_W[dir_index]);
            if (projection > best_projection) {
                best_projection = projection;
                best_center = center_index;
            }
        }

        const double support_error = targets.dense_support[dir_index] - best_projection;
        deltas_W[best_center] += support_error * targets.directions_W[dir_index];
        weights[best_center] += 1.0;
    }

    for (std::size_t center_index = 0; center_index < centers_W.size(); ++center_index) {
        if (!(weights[center_index] > 0.0)) {
            continue;
        }
        centers_W[center_index] += (blend / weights[center_index]) * deltas_W[center_index];
    }
}

void ApplyReinjectionDrivenCenterUpdate(const DirectionalSupportTargets& targets,
                                        const DenseSubpatch& subpatch,
                                        double blend,
                                        const std::vector<ClusterDirectionalEnvelope>& envelopes,
                                        std::vector<chrono::ChVector3d>& centers_W) {
    if (!(blend > 0.0) || centers_W.size() <= 1 || targets.directions_W.empty() ||
        envelopes.size() != centers_W.size()) {
        return;
    }

    std::vector<chrono::ChVector3d> deltas_W(centers_W.size(), chrono::ChVector3d(0.0, 0.0, 0.0));
    std::vector<double> weights(centers_W.size(), 0.0);
    for (std::size_t dir_index = 0; dir_index < targets.directions_W.size(); ++dir_index) {
        std::size_t best_center = 0;
        double best_reach = -std::numeric_limits<double>::infinity();
        for (std::size_t center_index = 0; center_index < centers_W.size(); ++center_index) {
            const double reach =
                ComputeReinjectionDirectionalReach(subpatch, targets.center_W, centers_W[center_index],
                                                   targets.directions_W[dir_index],
                                                   envelopes[center_index].positive_support[dir_index], centers_W.size());
            if (reach > best_reach) {
                best_reach = reach;
                best_center = center_index;
            }
        }

        const double reach_error = targets.dense_support[dir_index] - best_reach;
        const double weight = (dir_index < targets.direction_weights.size()) ? targets.direction_weights[dir_index] : 1.0;
        deltas_W[best_center] += weight * reach_error * targets.directions_W[dir_index];
        weights[best_center] += weight;
    }

    for (std::size_t center_index = 0; center_index < centers_W.size(); ++center_index) {
        if (!(weights[center_index] > 0.0)) {
            continue;
        }
        centers_W[center_index] += (blend / weights[center_index]) * deltas_W[center_index];
    }
}

void ComputeSubpatchBoundsUV(const std::vector<DenseContactPoint>& dense_points,
                             const DenseSubpatch& subpatch,
                             double& min_u,
                             double& max_u,
                             double& min_v,
                             double& max_v) {
    min_u = std::numeric_limits<double>::infinity();
    max_u = -std::numeric_limits<double>::infinity();
    min_v = std::numeric_limits<double>::infinity();
    max_v = -std::numeric_limits<double>::infinity();
    for (const auto point_index : subpatch.members) {
        double u = 0.0;
        double v = 0.0;
        ProjectToSubpatchUV(subpatch, dense_points[point_index].x_W, u, v);
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
    }
    if (!std::isfinite(min_u) || !std::isfinite(max_u) || !std::isfinite(min_v) || !std::isfinite(max_v)) {
        min_u = max_u = min_v = max_v = 0.0;
    }
}

void OptimizeSupportCenters(const std::vector<DenseContactPoint>& dense_points,
                            const DenseSubpatch& subpatch,
                            const CompressedContactConfig& cfg,
                            bool has_temporal_history,
                            std::vector<chrono::ChVector3d>& centers_W) {
    if (centers_W.empty()) {
        return;
    }

    chrono::ChVector3d avg_v_rel_W(0.0, 0.0, 0.0);
    double v_weight_sum = 0.0;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double weight = std::max(ProxyLoad(point), std::max(1.0e-12, point.area_weight));
        avg_v_rel_W += weight * point.v_rel_W;
        v_weight_sum += weight;
    }
    if (v_weight_sum > 1.0e-12) {
        avg_v_rel_W *= (1.0 / v_weight_sum);
    }
    const double approach_speed = std::max(0.0, -chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W));
    const chrono::ChVector3d tangential_v_W =
        avg_v_rel_W - chrono::Vdot(avg_v_rel_W, subpatch.avg_normal_W) * subpatch.avg_normal_W;
    const double tangential_speed = tangential_v_W.Length();
    const double motion_gate =
        1.0 / (1.0 + 2.5 * approach_speed / std::max(cfg.temporal_approach_velocity_scale, 1.0e-6) +
               tangential_speed / std::max(cfg.temporal_slip_velocity_scale, 1.0e-6));
    if (motion_gate < 0.55) {
        return;
    }

    auto assignments = AssignDensePointsToCenters(dense_points, subpatch, centers_W);
    const auto dense_stats = ComputeDensePlanarStats(dense_points, subpatch);
    const auto cone_targets =
        cfg.enable_cone_objective ? BuildDirectionalSupportTargets(dense_points, subpatch, cfg.cone_direction_count)
                                  : DirectionalSupportTargets{};
    const auto reinjection_targets =
        cfg.enable_reinjection_acceptance
            ? BuildReinjectionDirectionalTargets(dense_points, subpatch, dense_stats, avg_v_rel_W)
            : DirectionalSupportTargets{};
    if (!(dense_stats.weight_sum > 1.0e-12)) {
        return;
    }
    const auto original_centers_W = centers_W;
    const auto original_stats = ComputeSupportPlanarStats(subpatch, centers_W, assignments, dense_points);
    const auto original_envelopes = cfg.enable_reinjection_acceptance
                                        ? BuildClusterDirectionalEnvelopes(dense_points, subpatch, centers_W,
                                                                           assignments, reinjection_targets)
                                        : std::vector<ClusterDirectionalEnvelope>{};
    const double original_objective =
        ComputeSupportGeometryObjective(dense_stats, original_stats, cone_targets, subpatch, centers_W);
    const double original_reinjection_objective = cfg.enable_reinjection_acceptance
                                                      ? ComputeReinjectionObjective(subpatch, reinjection_targets,
                                                                                    centers_W, original_envelopes)
                                                      : 0.0;
    const double reinjection_weight = cfg.enable_reinjection_acceptance ? (has_temporal_history ? 2.2 : 2.8) : 0.0;
    const double original_combined_objective = original_objective + reinjection_weight * original_reinjection_objective;

    double min_u = 0.0;
    double max_u = 0.0;
    double min_v = 0.0;
    double max_v = 0.0;
    ComputeSubpatchBoundsUV(dense_points, subpatch, min_u, max_u, min_v, max_v);

    const int max_iters = has_temporal_history ? 1 : 2;
    for (int iter = 0; iter < max_iters; ++iter) {
        const auto support_stats = ComputeSupportPlanarStats(subpatch, centers_W, assignments, dense_points);
        if (!(support_stats.weight_sum > 1.0e-12)) {
            break;
        }
        const double dense_std_u = std::sqrt(std::max(dense_stats.cov_uu, 1.0e-12));
        const double dense_std_v = std::sqrt(std::max(dense_stats.cov_vv, 1.0e-12));
        const double support_std_u = std::sqrt(std::max(support_stats.cov_uu, 1.0e-12));
        const double support_std_v = std::sqrt(std::max(support_stats.cov_vv, 1.0e-12));
        const double scale_u = ClampScale(dense_std_u / support_std_u, 0.65, 1.65);
        const double scale_v = ClampScale(dense_std_v / support_std_v, 0.65, 1.65);
        const double cone_blend =
            cfg.enable_cone_objective
                ? motion_gate * (has_temporal_history ? 0.10 : ((centers_W.size() <= 2) ? 0.28 : 0.22))
                : 0.0;
        const auto current_envelopes =
            cfg.enable_reinjection_acceptance
                ? BuildClusterDirectionalEnvelopes(dense_points, subpatch, centers_W, assignments, reinjection_targets)
                : std::vector<ClusterDirectionalEnvelope>{};
        const double reinjection_blend =
            cfg.enable_reinjection_acceptance
                ? motion_gate * (has_temporal_history ? 0.14 : ((centers_W.size() <= 2) ? 0.26 : 0.20))
                : 0.0;
        const double moment_blend = motion_gate * (has_temporal_history ? 0.06 : ((centers_W.size() <= 2) ? 0.18 : 0.14));

        ApplyConeDrivenCenterUpdate(cone_targets, subpatch, cone_blend, centers_W);
        ApplyReinjectionDrivenCenterUpdate(reinjection_targets, subpatch, reinjection_blend, current_envelopes, centers_W);

        for (auto& center_W : centers_W) {
            double u = 0.0;
            double v = 0.0;
            ProjectToSubpatchUV(subpatch, center_W, u, v);
            const double target_u = dense_stats.mean_u + scale_u * (u - support_stats.mean_u);
            const double target_v = dense_stats.mean_v + scale_v * (v - support_stats.mean_v);
            const double blended_u = (1.0 - moment_blend) * u + moment_blend * target_u;
            const double blended_v = (1.0 - moment_blend) * v + moment_blend * target_v;
            center_W = LiftSubpatchUV(subpatch, std::clamp(blended_u, min_u, max_u),
                                      std::clamp(blended_v, min_v, max_v));
        }

        assignments = AssignDensePointsToCenters(dense_points, subpatch, centers_W);
    }

    const auto optimized_stats = ComputeSupportPlanarStats(subpatch, centers_W, assignments, dense_points);
    const auto optimized_envelopes = cfg.enable_reinjection_acceptance
                                         ? BuildClusterDirectionalEnvelopes(dense_points, subpatch, centers_W,
                                                                            assignments, reinjection_targets)
                                         : std::vector<ClusterDirectionalEnvelope>{};
    const double optimized_objective =
        ComputeSupportGeometryObjective(dense_stats, optimized_stats, cone_targets, subpatch, centers_W);
    const double optimized_reinjection_objective = cfg.enable_reinjection_acceptance
                                                       ? ComputeReinjectionObjective(subpatch, reinjection_targets,
                                                                                     centers_W, optimized_envelopes)
                                                       : 0.0;
    const double optimized_combined_objective =
        optimized_objective + reinjection_weight * optimized_reinjection_objective;
    const double primary_improvement_tol = std::max(1.0e-4, 0.01 * original_objective);
    const double primary_nonworse_tol = std::max(2.5e-4, 0.03 * original_objective);
    const double combined_improvement_tol = std::max(2.0e-4, 0.01 * original_combined_objective);
    const double reinjection_nonworse_tol =
        std::max(2.0e-4, (has_temporal_history ? 0.03 : 0.05) * original_reinjection_objective);
    const double reinjection_improve_tol =
        std::max(2.0e-4, 0.08 * original_reinjection_objective);
    const bool primary_improved = optimized_objective + primary_improvement_tol < original_objective;
    const bool combined_improved = optimized_combined_objective + combined_improvement_tol < original_combined_objective;
    const bool reinjection_nonworse =
        optimized_reinjection_objective <= original_reinjection_objective + reinjection_nonworse_tol;
    const bool primary_nonworse = optimized_objective <= original_objective + primary_nonworse_tol;
    const bool reinjection_improved =
        optimized_reinjection_objective + reinjection_improve_tol < original_reinjection_objective;
    if (cfg.enable_reinjection_acceptance) {
        if (!(combined_improved || (primary_improved && reinjection_nonworse) ||
              (primary_nonworse && reinjection_improved))) {
            centers_W = original_centers_W;
        }
    } else if (!primary_improved) {
        centers_W = original_centers_W;
    }
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

ReferenceWrench BuildProxyReference(const std::vector<DenseContactPoint>& dense_points,
                                    const DenseSubpatch& subpatch,
                                    const chrono::ChVector3d& origin_W) {
    ReferenceWrench reference;
    reference.origin_W = origin_W;
    reference.force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.total_load = 0.0;
    for (const auto point_index : subpatch.members) {
        const auto& point = dense_points[point_index];
        const double load = ProxyLoad(point);
        const chrono::ChVector3d force_W = load * point.n_W;
        reference.force_W += force_W;
        reference.moment_W += chrono::Vcross(point.x_W - origin_W, force_W);
        reference.total_load += std::max(0.0, chrono::Vdot(force_W, subpatch.avg_normal_W));
    }
    return reference;
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
                        double mean_allocated_load,
                        const ReducedContactPoint* previous_contact,
                        const CompressedContactConfig& cfg) {
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
        if (previous_contact && previous_contact->emission_count > 1 && support.slip_dominance > 0.7 &&
            support.tangential_heterogeneity > 0.5 * cfg.tangential_emission_threshold && load_ratio > 0.35) {
            return std::max(1, previous_contact->emission_count - 1);
        }
        return 1;
    }

    const double expand_score =
        1.40 * tangential_ratio + 0.85 * coverage_ratio + 0.30 * std::min(load_ratio, 2.5) +
        0.35 * slip_dominance + 0.65 * support.tangential_heterogeneity;
    int emission_count = 1;
    if (expand_score > 1.45 && tangential_ratio > 0.38 && coverage_ratio > 0.22 && support.dense_members >= 6 &&
        load_ratio > 0.9 && slip_dominance > 0.72) {
        emission_count = 5;
    } else if (expand_score > 1.18 && tangential_ratio > 0.24 && coverage_ratio > 0.14 && support.dense_members >= 3 &&
               slip_dominance > 0.55) {
        emission_count = 3;
    } else {
        emission_count = (expand_score > 0.72 && slip_dominance > 0.22) ? 2 : 1;
    }

    if (support.tangential_heterogeneity > cfg.tangential_emission_threshold && slip_dominance > 0.6 &&
        support.dense_members >= 4 && emission_count < 5) {
        emission_count = (support.tangential_heterogeneity > 1.5 * cfg.tangential_emission_threshold) ? 5 : 3;
    }

    if (previous_contact && previous_contact->emission_count > emission_count && slip_dominance > 0.45 &&
        tangential_ratio > 0.18 && load_ratio > 0.55) {
        emission_count = std::max(emission_count, previous_contact->emission_count);
    }

    return emission_count;
}

double ChooseStencilHalfExtent(const ReducedSupportAggregate& support,
                               const DenseSubpatch& subpatch,
                               int emission_count,
                               const ReducedContactPoint* previous_contact,
                               const CompressedContactConfig& cfg) {
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
    extent *= std::clamp(1.0 + 0.30 * support.tangential_heterogeneity, 1.0, 1.25);
    if (previous_contact && previous_contact->stencil_half_extent > 0.0 && support.slip_dominance > 0.45) {
        const double blend = std::clamp(cfg.temporal_stencil_blend * support.slip_dominance, 0.0, 0.45);
        extent = (1.0 - blend) * extent + blend * previous_contact->stencil_half_extent;
    }
    return extent;
}

double ChooseSecondaryStencilHalfExtent(const ReducedSupportAggregate& support,
                                        const DenseSubpatch& subpatch,
                                        int emission_count,
                                        const ReducedContactPoint* previous_contact,
                                        const CompressedContactConfig& cfg) {
    if (emission_count < 5) {
        return 0.0;
    }

    const chrono::ChVector3d tangential_force_W =
        support.allocated_force_W - chrono::Vdot(support.allocated_force_W, subpatch.avg_normal_W) *
                                        subpatch.avg_normal_W;
    const double total_force_norm = support.allocated_force_W.Length();
    const double tangential_ratio = tangential_force_W.Length() / std::max(total_force_norm, 1.0e-12);
    const double base_extent = std::min(0.18 * support.coverage_radius, 0.09 * subpatch.diameter);
    double extent =
        std::clamp((0.45 + 0.35 * tangential_ratio + 0.25 * support.tangential_heterogeneity) * base_extent, 0.0,
                   1.15 * base_extent);
    if (previous_contact && previous_contact->stencil_half_extent_secondary > 0.0 && support.slip_dominance > 0.55) {
        const double blend = std::clamp(cfg.temporal_stencil_blend * support.slip_dominance, 0.0, 0.45);
        extent = (1.0 - blend) * extent + blend * previous_contact->stencil_half_extent_secondary;
    }
    return extent;
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

constexpr double kReinjectionSlipVelocityScale = 1.0e-3;

int ActiveStencilSlotCount(int emission_count) {
    if (emission_count >= 5) {
        return 5;
    }
    if (emission_count >= 3) {
        return 3;
    }
    return (emission_count == 2) ? 2 : 1;
}

std::vector<int> ActiveStencilSlots(int emission_count) {
    if (emission_count >= 5) {
        return {0, 1, 2, 3, 4};
    }
    if (emission_count >= 3) {
        return {0, 1, 2};
    }
    if (emission_count == 2) {
        return {1, 2};
    }
    return {0};
}

std::array<chrono::ChVector3d, 5> BuildStencilSlotOffsets(const ReducedContactPoint& reduced,
                                                          const chrono::ChVector3d& primary_axis_W,
                                                          const chrono::ChVector3d& secondary_axis_W) {
    std::array<chrono::ChVector3d, 5> offsets{
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0),
    };
    offsets[1] = -reduced.stencil_half_extent * primary_axis_W;
    offsets[2] = reduced.stencil_half_extent * primary_axis_W;
    offsets[3] = -reduced.stencil_half_extent_secondary * secondary_axis_W;
    offsets[4] = reduced.stencil_half_extent_secondary * secondary_axis_W;
    return offsets;
}

chrono::ChVector3d BuildDenseSeedProxyForce(const DenseContactPoint& point,
                                            const chrono::ChVector3d& fallback_tangent_W,
                                            double mu_default) {
    const double normal_load = ProxyLoad(point);
    chrono::ChVector3d force_W = normal_load * point.n_W;
    if (!(mu_default > 0.0) || !(normal_load > 0.0)) {
        return force_W;
    }

    chrono::ChVector3d tangential_dir_W = ProjectToTangentUnit(point.v_rel_W, point.n_W);
    if (tangential_dir_W.Length2() <= 0.0) {
        tangential_dir_W = ProjectToTangentUnit(fallback_tangent_W, point.n_W);
    }
    if (tangential_dir_W.Length2() <= 0.0) {
        return force_W;
    }

    const chrono::ChVector3d tangential_velocity_W =
        point.v_rel_W - chrono::Vdot(point.v_rel_W, point.n_W) * point.n_W;
    const double tangential_speed = tangential_velocity_W.Length();
    const double slip_ratio = tangential_speed / std::max(tangential_speed + kReinjectionSlipVelocityScale, 1.0e-9);
    force_W -= (mu_default * slip_ratio * normal_load) * tangential_dir_W;
    return force_W;
}

double DotStencilColumn(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

double EvaluateStencilScalarObjective(const std::vector<std::vector<double>>& columns,
                                      const std::vector<double>& target,
                                      const std::vector<double>& initial,
                                      double regularization,
                                      const std::vector<double>& weights) {
    std::vector<double> residual(target.size(), 0.0);
    for (std::size_t row = 0; row < target.size(); ++row) {
        residual[row] = -target[row];
    }
    for (std::size_t column_index = 0; column_index < columns.size(); ++column_index) {
        const double weight = weights[column_index];
        if (!(std::abs(weight) > 0.0)) {
            continue;
        }
        for (std::size_t row = 0; row < target.size(); ++row) {
            residual[row] += columns[column_index][row] * weight;
        }
    }

    double objective = 0.5 * DotStencilColumn(residual, residual);
    if (regularization > 0.0) {
        for (std::size_t i = 0; i < weights.size(); ++i) {
            const double delta = weights[i] - initial[i];
            objective += 0.5 * regularization * delta * delta;
        }
    }
    return objective;
}

double EstimateStencilScalarLipschitz(const std::vector<std::vector<double>>& columns, double regularization) {
    double trace_bound = regularization;
    for (const auto& column : columns) {
        trace_bound += DotStencilColumn(column, column);
    }
    return std::max(trace_bound, 1.0e-6);
}

void ComputeStencilScalarGradient(const std::vector<std::vector<double>>& columns,
                                  const std::vector<double>& target,
                                  const std::vector<double>& initial,
                                  double regularization,
                                  const std::vector<double>& weights,
                                  std::vector<double>& gradient) {
    std::vector<double> residual(target.size(), 0.0);
    for (std::size_t row = 0; row < target.size(); ++row) {
        residual[row] = -target[row];
    }
    for (std::size_t column_index = 0; column_index < columns.size(); ++column_index) {
        const double weight = weights[column_index];
        if (!(std::abs(weight) > 0.0)) {
            continue;
        }
        for (std::size_t row = 0; row < target.size(); ++row) {
            residual[row] += columns[column_index][row] * weight;
        }
    }

    gradient.assign(columns.size(), 0.0);
    for (std::size_t column_index = 0; column_index < columns.size(); ++column_index) {
        double value = DotStencilColumn(columns[column_index], residual);
        if (regularization > 0.0) {
            value += regularization * (weights[column_index] - initial[column_index]);
        }
        gradient[column_index] = value;
    }
}

bool SolveNonnegativeStencilScalars(const std::vector<std::vector<double>>& columns,
                                    const std::vector<double>& target,
                                    const std::vector<double>& initial,
                                    double regularization,
                                    std::vector<double>& out_weights) {
    out_weights = initial;
    if (columns.empty()) {
        return true;
    }

    std::vector<double> x(columns.size(), 0.0);
    for (std::size_t i = 0; i < columns.size(); ++i) {
        x[i] = std::max(0.0, initial[i]);
    }
    std::vector<double> y = x;
    std::vector<double> x_next(columns.size(), 0.0);
    std::vector<double> gradient;
    const double step = 1.0 / EstimateStencilScalarLipschitz(columns, regularization);
    double t = 1.0;

    constexpr int kMaxIterations = 512;
    constexpr double kTolerance = 1.0e-11;
    for (int iter = 0; iter < kMaxIterations; ++iter) {
        ComputeStencilScalarGradient(columns, target, initial, regularization, y, gradient);

        double delta_sq = 0.0;
        for (std::size_t i = 0; i < columns.size(); ++i) {
            x_next[i] = std::max(0.0, y[i] - step * gradient[i]);
            const double delta = x_next[i] - x[i];
            delta_sq += delta * delta;
        }
        if (delta_sq <= kTolerance * kTolerance) {
            x = x_next;
            break;
        }

        const double t_next = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t * t));
        const double momentum = (t - 1.0) / t_next;
        for (std::size_t i = 0; i < columns.size(); ++i) {
            y[i] = x_next[i] + momentum * (x_next[i] - x[i]);
        }
        x.swap(x_next);
        t = t_next;
    }

    const double initial_objective =
        EvaluateStencilScalarObjective(columns, target, initial, regularization, initial);
    const double solved_objective = EvaluateStencilScalarObjective(columns, target, initial, regularization, x);
    if (solved_objective <= initial_objective + 1.0e-12) {
        out_weights = std::move(x);
    } else {
        out_weights = initial;
    }
    return true;
}

void SeedReactionCaches(const std::vector<DenseContactPoint>& dense_points,
                        const DenseSubpatch& subpatch,
                        const ReducedSupportAggregate& support,
                        const CompressedContactConfig& cfg,
                        double step_size,
                        const ReducedContactPoint* previous_contact,
                        ReducedContactPoint& reduced) {
    reduced.reaction_cache_primary.fill(0.0f);
    reduced.reaction_cache_secondary.fill(0.0f);
    reduced.reaction_cache_tertiary.fill(0.0f);
    reduced.reaction_cache_quaternary.fill(0.0f);
    reduced.reaction_cache_quinary.fill(0.0f);
    reduced.stencil_gap_offsets.fill(0.0);
    reduced.slot_weight.fill(0.0);
    for (int slot = 0; slot < 5; ++slot) {
        reduced.slot_x_W[slot] = reduced.x_W;
        reduced.slot_x_master_surface_W[slot] = reduced.x_master_surface_W;
        reduced.slot_n_W[slot] = reduced.n_W;
    }

    chrono::ChVector3d primary_axis_W = ProjectToTangentUnit(reduced.stencil_axis_W, reduced.n_W);
    chrono::ChVector3d secondary_axis_W =
        ProjectToTangentUnit(reduced.stencil_axis_secondary_W, reduced.n_W);
    if (secondary_axis_W.Length2() <= 0.0) {
        secondary_axis_W = ProjectToTangentUnit(chrono::Vcross(reduced.n_W, primary_axis_W), reduced.n_W);
    }
    if (primary_axis_W.Length2() <= 0.0) {
        primary_axis_W = subpatch.t1_W;
    }
    if (secondary_axis_W.Length2() <= 0.0) {
        secondary_axis_W = subpatch.t2_W;
    }

    const auto normal_weights = ComputeNormalStencilWeights(dense_points, support, subpatch, primary_axis_W,
                                                            secondary_axis_W, reduced.emission_count);
    const auto tangential_weights = ComputeTangentialStencilWeights(dense_points, support, subpatch, primary_axis_W,
                                                                    secondary_axis_W, reduced.emission_count);
    const auto slot_ids = ActiveStencilSlots(reduced.emission_count);
    const auto slot_offsets_W = BuildStencilSlotOffsets(reduced, primary_axis_W, secondary_axis_W);
    const int active_slots = ActiveStencilSlotCount(reduced.emission_count);

    const chrono::ChVector3d normal_force_W = chrono::Vdot(reduced.allocated_force_W, reduced.n_W) * reduced.n_W;
    const double normal_force_mag = std::max(0.0, chrono::Vdot(reduced.allocated_force_W, reduced.n_W));
    const chrono::ChVector3d tangential_force_W = reduced.allocated_force_W - normal_force_W;
    chrono::ChVector3d tangential_dir_W = ProjectToTangentUnit(tangential_force_W, reduced.n_W);
    if (tangential_dir_W.Length2() <= 0.0) {
        tangential_dir_W = primary_axis_W;
    }
    const double tangential_force_mag =
        (tangential_dir_W.Length2() > 0.0) ? chrono::Vdot(tangential_force_W, tangential_dir_W) : 0.0;

    chrono::ChVector3d proxy_normal_moment_W(0.0, 0.0, 0.0);
    double proxy_normal_sum = 0.0;
    chrono::ChVector3d proxy_tangential_moment_W(0.0, 0.0, 0.0);
    double proxy_tangential_sum = 0.0;
    for (const auto dense_index : support.member_indices) {
        const auto& point = dense_points[dense_index];
        const chrono::ChVector3d proxy_force_W =
            BuildDenseSeedProxyForce(point, tangential_dir_W, std::max(0.0, reduced.mu));
        const chrono::ChVector3d rel_W = point.x_W - reduced.x_W;

        const double proxy_normal_load = std::max(0.0, chrono::Vdot(proxy_force_W, reduced.n_W));
        proxy_normal_sum += proxy_normal_load;
        proxy_normal_moment_W += chrono::Vcross(rel_W, proxy_normal_load * reduced.n_W);

        const chrono::ChVector3d tangential_proxy_force_W = proxy_force_W - proxy_normal_load * reduced.n_W;
        const double proxy_tangential_load = tangential_proxy_force_W.Length();
        if (proxy_tangential_load > 1.0e-12) {
            proxy_tangential_sum += proxy_tangential_load;
            proxy_tangential_moment_W += chrono::Vcross(rel_W, tangential_proxy_force_W);
        }
    }

    std::vector<double> normal_initial(slot_ids.size(), 0.0);
    std::vector<double> tangential_initial(slot_ids.size(), 0.0);
    for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
        normal_initial[local_index] = normal_force_mag * normal_weights[slot_ids[local_index]];
        tangential_initial[local_index] =
            std::max(0.0, tangential_force_mag) * tangential_weights[slot_ids[local_index]];
    }

    std::vector<double> normal_solution = normal_initial;
    if (normal_force_mag > 1.0e-12) {
        std::vector<std::vector<double>> columns(slot_ids.size(), std::vector<double>{1.0});
        std::vector<double> target{normal_force_mag};
        const chrono::ChVector3d target_normal_moment_W =
            (proxy_normal_sum > 1.0e-12) ? ((normal_force_mag / proxy_normal_sum) * proxy_normal_moment_W)
                                         : chrono::ChVector3d(0.0, 0.0, 0.0);
        if (active_slots > 1 && cfg.reinjection_normal_moment_weight > 0.0) {
            const double moment_weight = cfg.reinjection_normal_moment_weight;
            target.push_back(moment_weight * chrono::Vdot(target_normal_moment_W, primary_axis_W));
            for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
                const auto& offset_W = slot_offsets_W[slot_ids[local_index]];
                columns[local_index].push_back(
                    moment_weight * chrono::Vdot(chrono::Vcross(offset_W, reduced.n_W), primary_axis_W));
            }
            if (reduced.emission_count >= 5) {
                target.push_back(moment_weight * chrono::Vdot(target_normal_moment_W, secondary_axis_W));
                for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
                    const auto& offset_W = slot_offsets_W[slot_ids[local_index]];
                    columns[local_index].push_back(
                        moment_weight * chrono::Vdot(chrono::Vcross(offset_W, reduced.n_W), secondary_axis_W));
                }
            }
        }
        SolveNonnegativeStencilScalars(columns, target, normal_initial, std::max(0.0, cfg.reinjection_seed_regularization),
                                       normal_solution);
        const double solved_sum = std::accumulate(normal_solution.begin(), normal_solution.end(), 0.0);
        if (solved_sum > 1.0e-12) {
            const double rescale = normal_force_mag / solved_sum;
            for (auto& value : normal_solution) {
                value *= rescale;
            }
        }
    }

    std::vector<double> tangential_solution = tangential_initial;
    if (tangential_dir_W.Length2() > 0.0 && tangential_force_mag > 1.0e-12) {
        std::vector<std::vector<double>> columns(slot_ids.size(), std::vector<double>{1.0});
        std::vector<double> target{tangential_force_mag};
        if (active_slots > 1 && cfg.reinjection_tangential_moment_weight > 0.0) {
            const double moment_weight = cfg.reinjection_tangential_moment_weight;
            const double target_torsion =
                (proxy_tangential_sum > 1.0e-12)
                    ? chrono::Vdot(reduced.n_W,
                                   (std::abs(tangential_force_mag) / proxy_tangential_sum) * proxy_tangential_moment_W)
                    : 0.0;
            target.push_back(moment_weight * target_torsion);
            for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
                const auto& offset_W = slot_offsets_W[slot_ids[local_index]];
                columns[local_index].push_back(
                    moment_weight * chrono::Vdot(reduced.n_W, chrono::Vcross(offset_W, tangential_dir_W)));
            }
        }
        SolveNonnegativeStencilScalars(columns, target, tangential_initial,
                                       std::max(0.0, cfg.reinjection_seed_regularization), tangential_solution);
        const double solved_sum = std::accumulate(tangential_solution.begin(), tangential_solution.end(), 0.0);
        if (solved_sum > 1.0e-12) {
            const double rescale = tangential_force_mag / solved_sum;
            for (auto& value : tangential_solution) {
                value *= rescale;
            }
        }
    }

    std::array<chrono::ChVector3d, 5> slot_forces_W{
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0),
    };
    for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
        const int slot = slot_ids[local_index];
        slot_forces_W[slot] = normal_solution[local_index] * reduced.n_W;
        if (tangential_dir_W.Length2() > 0.0) {
            slot_forces_W[slot] += tangential_solution[local_index] * tangential_dir_W;
        }
    }

    std::array<chrono::ChVector3d, 5> slot_default_x_W{
        reduced.x_W, reduced.x_W + slot_offsets_W[1], reduced.x_W + slot_offsets_W[2], reduced.x_W + slot_offsets_W[3],
        reduced.x_W + slot_offsets_W[4],
    };
    std::array<chrono::ChVector3d, 5> slot_default_master_W{
        reduced.x_master_surface_W,
        reduced.x_master_surface_W + slot_offsets_W[1],
        reduced.x_master_surface_W + slot_offsets_W[2],
        reduced.x_master_surface_W + slot_offsets_W[3],
        reduced.x_master_surface_W + slot_offsets_W[4],
    };

    std::array<double, 5> cluster_weight_sum{0.0, 0.0, 0.0, 0.0, 0.0};
    std::array<chrono::ChVector3d, 5> cluster_x_sum{
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
    };
    std::array<chrono::ChVector3d, 5> cluster_master_sum{
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
    };
    std::array<chrono::ChVector3d, 5> cluster_normal_sum{
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 0.0, 0.0),
    };

    for (const auto dense_index : support.member_indices) {
        const auto& point = dense_points[dense_index];
        const double weight = DenseMetricWeight(point);
        int best_slot = slot_ids.front();
        double best_distance = std::numeric_limits<double>::infinity();
        for (const int slot : slot_ids) {
            const double distance = TangentialRadiusSquared(point.x_W, slot_default_x_W[slot], reduced.n_W);
            if (distance < best_distance) {
                best_distance = distance;
                best_slot = slot;
            }
        }
        cluster_weight_sum[best_slot] += weight;
        cluster_x_sum[best_slot] += weight * point.x_W;
        cluster_master_sum[best_slot] += weight * point.x_master_surface_W;
        cluster_normal_sum[best_slot] += weight * point.n_W;
    }

    for (int slot = 0; slot < 5; ++slot) {
        reduced.slot_x_W[slot] = slot_default_x_W[slot];
        reduced.slot_x_master_surface_W[slot] = slot_default_master_W[slot];
        reduced.slot_n_W[slot] = reduced.n_W;
    }

    double cluster_total = 0.0;
    for (const int slot : slot_ids) {
        cluster_total += cluster_weight_sum[slot];
    }
    for (const int slot : slot_ids) {
        const double weight = cluster_weight_sum[slot];
        if (weight > 1.0e-12) {
            const chrono::ChVector3d clustered_x_W = cluster_x_sum[slot] * (1.0 / weight);
            const chrono::ChVector3d clustered_master_W = cluster_master_sum[slot] * (1.0 / weight);
            const chrono::ChVector3d clustered_normal_W =
                SafeNormalized(cluster_normal_sum[slot], reduced.n_W);
            const double blend = std::clamp(0.45 + 0.35 * (weight / std::max(cluster_total, 1.0e-12)), 0.45, 0.85);
            reduced.slot_x_W[slot] = (1.0 - blend) * slot_default_x_W[slot] + blend * clustered_x_W;
            reduced.slot_x_master_surface_W[slot] =
                (1.0 - blend) * slot_default_master_W[slot] + blend * clustered_master_W;
            reduced.slot_n_W[slot] = SafeNormalized((1.0 - blend) * reduced.n_W + blend * clustered_normal_W, reduced.n_W);
            reduced.slot_weight[slot] = weight;
        }
    }

    if (previous_contact && previous_contact->emission_count >= reduced.emission_count && support.slip_dominance > 0.45) {
        const double approach_speed = std::max(0.0, -chrono::Vdot(support.v_rel_W, reduced.n_W));
        const double motion_gate =
            support.slip_dominance /
            (1.0 + approach_speed / std::max(cfg.temporal_approach_velocity_scale, 1.0e-6));
        const double blend = std::clamp(cfg.temporal_stencil_blend * motion_gate, 0.0, 0.45);
        for (const int slot : slot_ids) {
            reduced.slot_x_W[slot] =
                (1.0 - blend) * reduced.slot_x_W[slot] + blend * previous_contact->slot_x_W[slot];
            reduced.slot_x_master_surface_W[slot] =
                (1.0 - blend) * reduced.slot_x_master_surface_W[slot] +
                blend * previous_contact->slot_x_master_surface_W[slot];
            reduced.slot_n_W[slot] = SafeNormalized((1.0 - blend) * reduced.slot_n_W[slot] +
                                                        blend * previous_contact->slot_n_W[slot],
                                                    reduced.slot_n_W[slot]);
            reduced.slot_weight[slot] =
                (1.0 - blend) * reduced.slot_weight[slot] + blend * previous_contact->slot_weight[slot];
        }
    }

    const auto build_slot_basis = [](const chrono::ChVector3d& n_W,
                                     const chrono::ChVector3d& preferred_t1_W,
                                     chrono::ChVector3d& t1_W,
                                     chrono::ChVector3d& t2_W) {
        t1_W = preferred_t1_W - chrono::Vdot(preferred_t1_W, n_W) * n_W;
        double t1_len = t1_W.Length();
        if (!(t1_len > 1.0e-12)) {
            const chrono::ChVector3d seed =
                (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0) : chrono::ChVector3d(1.0, 0.0, 0.0);
            t1_W = chrono::Vcross(seed, n_W);
            t1_len = t1_W.Length();
        }
        if (t1_len > 1.0e-12) {
            t1_W *= (1.0 / t1_len);
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
        t1_W = SafeNormalized(chrono::Vcross(t2_W, n_W), t1_W);
    };

    std::vector<SupportWrenchPoint> slot_supports;
    slot_supports.reserve(slot_ids.size());
    for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
        const int slot = slot_ids[local_index];
        SupportWrenchPoint point;
        point.x_W = reduced.slot_x_W[slot];
        point.n_W = SafeNormalized(reduced.slot_n_W[slot], reduced.n_W);
        const chrono::ChVector3d preferred_axis_W =
            (tangential_dir_W.Length2() > 0.0) ? tangential_dir_W : primary_axis_W;
        build_slot_basis(point.n_W, preferred_axis_W, point.t1_W, point.t2_W);
        point.mu = std::max(0.0, reduced.mu);
        point.initial_force_W = slot_forces_W[slot];
        point.initial_load = std::max(0.0, chrono::Vdot(slot_forces_W[slot], point.n_W));
        slot_supports.push_back(point);
    }

    ReferenceWrench reinjection_reference;
    reinjection_reference.origin_W = reduced.x_W;
    reinjection_reference.force_W = reduced.allocated_force_W;
    reinjection_reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reinjection_reference.total_load = normal_force_mag;
    if (proxy_normal_sum > 1.0e-12 && normal_force_mag > 1.0e-12) {
        reinjection_reference.moment_W += (normal_force_mag / proxy_normal_sum) * proxy_normal_moment_W;
    }
    if (proxy_tangential_sum > 1.0e-12 && tangential_force_W.Length() > 1.0e-12) {
        reinjection_reference.moment_W +=
            (tangential_force_W.Length() / proxy_tangential_sum) * proxy_tangential_moment_W;
    }

    WrenchAllocationResult slot_allocation;
    ReducedSolveOptions slot_options;
    slot_options.friction_ray_count = std::max(8, cfg.reduced_friction_rays);
    slot_options.temporal_regularization = std::max(1.0e-10, cfg.reinjection_seed_regularization);
    LocalWrenchAllocator::Allocate(slot_supports, reinjection_reference, slot_options, slot_allocation);
    if (slot_allocation.feasible && slot_allocation.forces_W.size() == slot_supports.size()) {
        for (std::size_t local_index = 0; local_index < slot_ids.size(); ++local_index) {
            slot_forces_W[slot_ids[local_index]] = slot_allocation.forces_W[local_index];
        }
    }

    if (active_slots > 1) {
        const double uniform_weight = 1.0 / static_cast<double>(active_slots);
        const double gap_bias_scale = std::min(0.45 * std::abs(reduced.phi_eff), 3.0e-4);
        double normal_sum = 0.0;
        for (const int slot : slot_ids) {
            normal_sum += std::max(0.0, chrono::Vdot(slot_forces_W[slot], reduced.slot_n_W[slot]));
        }
        normal_sum = std::max(normal_sum, 1.0e-12);
        for (const int slot : slot_ids) {
            const double weight = std::max(0.0, chrono::Vdot(slot_forces_W[slot], reduced.slot_n_W[slot])) / normal_sum;
            reduced.stencil_gap_offsets[slot] = gap_bias_scale * (uniform_weight - weight);
            reduced.slot_weight[slot] = weight;
        }
    }
    const double seed_step = std::max(step_size, 0.0);
    auto seed_slot = [&](std::array<float, 6>& reaction_cache, int slot_index) {
        EncodeReactionCacheWorldImpulse(reaction_cache, seed_step * slot_forces_W[slot_index], reduced.slot_n_W[slot_index]);
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
    if (!cfg.enable_impulse_transport || !previous || !previous->has_impulse_wrench || supports.empty() ||
        !(transport_alpha > 0.0)) {
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
    LocalWrenchAllocator::Allocate(
        support_points, transport_reference,
        MakeReducedSolveOptions(cfg, std::max(1.0e-10, cfg.temporal_load_regularization)), transport_result);
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
    if (!cfg.enable_impulse_transport || !(previous_seed.confidence > 0.0) ||
        !(cfg.temporal_force_transport_blend > 0.0)) {
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
                                                                     const CompressedContactConfig& cfg,
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

    OptimizeSupportCenters(dense_points, subpatch, cfg, !previous_contacts.empty(), centers_W);
    std::vector<std::size_t> assignments = AssignDensePointsToCenters(dense_points, subpatch, centers_W);

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

        const auto tangential_stats =
            ComputeTangentialDistributionStats(dense_points, member_indices, subpatch.avg_normal_W, avg_x_W,
                                               subpatch.t1_W, subpatch.diameter);

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
        support.slip_dominance = tangential_stats.slip_dominance;
        support.tangential_heterogeneity = tangential_stats.heterogeneity;
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
    if (!cfg.enable_sentinel_monitor || supports.empty()) {
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

        const auto subpatch_tangential_stats =
            ComputeTangentialDistributionStats(dense_points, subpatch.members, subpatch.avg_normal_W,
                                               subpatch.centroid_W, subpatch.t1_W, subpatch.diameter);
        const int max_dynamic_points =
            std::max(cfg_.max_reduced_points_per_patch, cfg_.max_dynamic_reduced_points_per_patch);
        int target_point_cap = cfg_.max_reduced_points_per_patch;
        if (subpatch_tangential_stats.heterogeneity > cfg_.tangential_heterogeneity_threshold &&
            subpatch_tangential_stats.slip_dominance > 0.45) {
            target_point_cap = max_dynamic_points;
        }

        int target_points = std::min(target_point_cap,
                                     (subpatch.members.size() <= 1)
                                         ? 1
                                         : ((subpatch.members.size() <= 2) ? 2
                                                                           : ((subpatch_tangential_stats.heterogeneity >
                                                                               cfg_.tangential_heterogeneity_threshold)
                                                                                  ? 4
                                                                                  : 3)));
        DenseMicroReferenceResult dense_micro_reference;
        if (cfg_.enable_dense_micro_solver) {
            const auto dense_micro_options = MakeDenseMicroSolverOptions(cfg_);
            LocalWrenchAllocator::BuildDenseMicroReference(dense_points, subpatch.members, subpatch.centroid_W,
                                                           mu_default, step_size, dense_micro_options,
                                                           dense_micro_reference);
        } else {
            dense_micro_reference.reference = BuildProxyReference(dense_points, subpatch, subpatch.centroid_W);
            dense_micro_reference.force_residual = 0.0;
            dense_micro_reference.moment_residual = 0.0;
            dense_micro_reference.feasible = true;
        }
        const ReferenceWrench& dense_reference = dense_micro_reference.reference;
        max_dense_micro_force_residual =
            std::max(max_dense_micro_force_residual, dense_micro_reference.force_residual);
        max_dense_micro_moment_residual =
            std::max(max_dense_micro_moment_residual, dense_micro_reference.moment_residual);
        double temporal_reference_alpha = 0.0;
        if (cfg_.enable_impulse_transport && matched_previous_state) {
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
            supports = BuildReducedSupportsForSubpatch(dense_points, subpatch, cfg_, selection_previous_contacts,
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
                    (cfg_.enable_impulse_transport && matched_previous_state && matched_previous_state->has_impulse_wrench)
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
                LocalWrenchAllocator::Allocate(
                    support_points, reference,
                    MakeReducedSolveOptions(cfg_, cfg_.temporal_load_regularization), allocation);
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
            const bool require_extra_tangential_capacity =
                subpatch_tangential_stats.heterogeneity > cfg_.tangential_heterogeneity_threshold &&
                subpatch_tangential_stats.slip_dominance > 0.45 && target_points < target_point_cap;
            if (((sigma_ok && cone_ok && gap_ok && wrench_ok && cop_ok) && !require_extra_tangential_capacity) ||
                target_points >= target_point_cap) {
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
            const ReducedContactPoint* previous_support_contact = nullptr;
            for (const auto& previous : previous_local_contacts) {
                if (previous.support_id == support_ids[support_index]) {
                    previous_support_contact = &previous;
                    break;
                }
            }
            ReducedContactPoint reduced;
            reduced.persistent_id = persistent_id;
            reduced.patch_id = support.patch_id;
            reduced.subpatch_id = support.subpatch_id;
            reduced.support_id = support_ids[support_index];
            reduced.dense_members = support.dense_members;
            reduced.emission_count =
                ChooseEmissionCount(support, subpatch, mean_allocated_load, previous_support_contact, cfg_);
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
            reduced.stencil_half_extent =
                ChooseStencilHalfExtent(support, subpatch, reduced.emission_count, previous_support_contact, cfg_);
            reduced.stencil_half_extent_secondary =
                ChooseSecondaryStencilHalfExtent(support, subpatch, reduced.emission_count, previous_support_contact,
                                                cfg_);
            reduced.mu = mu_default;
            SeedReactionCaches(dense_points, subpatch, support, cfg_, step_size, previous_support_contact, reduced);
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
