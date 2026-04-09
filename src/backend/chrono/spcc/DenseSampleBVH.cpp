#include "platform/backend/spcc/DenseSampleBVH.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

namespace platform {
namespace backend {
namespace spcc {

namespace {

chrono::ChVector3d SafeNormalized(const chrono::ChVector3d& v, const chrono::ChVector3d& fallback) {
    const double len = v.Length();
    if (!(len > 1.0e-12) || !std::isfinite(len)) {
        return fallback;
    }
    return v * (1.0 / len);
}

double Dot(const chrono::ChVector3d& a, const chrono::ChVector3d& b) {
    return chrono::Vdot(a, b);
}

chrono::ChVector3d MulT(const chrono::ChMatrix33<>& A, const chrono::ChVector3d& v) {
    return A.transpose() * v;
}

double Clamp(double value, double lo, double hi) {
    return std::max(lo, std::min(value, hi));
}

chrono::ChMatrix33<> OuterProduct(const chrono::ChVector3d& a, const chrono::ChVector3d& b) {
    chrono::ChMatrix33<> out(0.0);
    out(0, 0) = a.x() * b.x();
    out(0, 1) = a.x() * b.y();
    out(0, 2) = a.x() * b.z();
    out(1, 0) = a.y() * b.x();
    out(1, 1) = a.y() * b.y();
    out(1, 2) = a.y() * b.z();
    out(2, 0) = a.z() * b.x();
    out(2, 1) = a.z() * b.y();
    out(2, 2) = a.z() * b.z();
    return out;
}

chrono::ChVector3d Mul(const chrono::ChMatrix33<>& A, const chrono::ChVector3d& v) {
    return A * v;
}

chrono::ChVector3d PowerIterateSymmetric(const chrono::ChMatrix33<>& A,
                                         const chrono::ChVector3d& seed,
                                         const chrono::ChVector3d& fallback) {
    chrono::ChVector3d v = SafeNormalized(seed, fallback);
    for (int iter = 0; iter < 16; ++iter) {
        const chrono::ChVector3d next = Mul(A, v);
        v = SafeNormalized(next, v);
    }
    return SafeNormalized(v, fallback);
}

chrono::ChVector3d Orthogonalized(const chrono::ChVector3d& v,
                                  const chrono::ChVector3d& axis,
                                  const chrono::ChVector3d& fallback) {
    const chrono::ChVector3d projected = v - Dot(v, axis) * axis;
    return SafeNormalized(projected, fallback);
}

std::array<chrono::ChVector3d, 3> ComputeAxes(const std::vector<DenseSurfaceSample>& samples,
                                              const std::vector<std::size_t>& ordered_indices,
                                              uint32_t first,
                                              uint32_t count) {
    chrono::ChVector3d mean(0.0, 0.0, 0.0);
    for (uint32_t offset = 0; offset < count; ++offset) {
        mean += samples[ordered_indices[first + offset]].xi_slave_S;
    }
    mean *= (1.0 / static_cast<double>(count));

    chrono::ChMatrix33<> covariance(0.0);
    for (uint32_t offset = 0; offset < count; ++offset) {
        const chrono::ChVector3d rel = samples[ordered_indices[first + offset]].xi_slave_S - mean;
        covariance += OuterProduct(rel, rel);
    }

    const chrono::ChVector3d axis0 =
        PowerIterateSymmetric(covariance, chrono::ChVector3d(1.0, 0.0, 0.0), chrono::ChVector3d(1.0, 0.0, 0.0));
    const double lambda0 = Dot(axis0, Mul(covariance, axis0));
    chrono::ChMatrix33<> deflated = covariance - lambda0 * OuterProduct(axis0, axis0);
    chrono::ChVector3d axis1 =
        PowerIterateSymmetric(deflated, chrono::ChVector3d(0.0, 1.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    axis1 = Orthogonalized(axis1, axis0, chrono::ChVector3d(0.0, 1.0, 0.0));
    chrono::ChVector3d axis2 = SafeNormalized(chrono::Vcross(axis0, axis1), chrono::ChVector3d(0.0, 0.0, 1.0));
    axis1 = SafeNormalized(chrono::Vcross(axis2, axis0), chrono::ChVector3d(0.0, 1.0, 0.0));

    return {axis0, axis1, axis2};
}

GeometryOBB ComputeOBB(const std::vector<DenseSurfaceSample>& samples,
                       const std::vector<std::size_t>& ordered_indices,
                       uint32_t first,
                       uint32_t count) {
    GeometryOBB obb;
    if (count == 0) {
        return obb;
    }

    const auto axes = ComputeAxes(samples, ordered_indices, first, count);
    chrono::ChMatrix33<> R_SO(1.0);
    for (int row = 0; row < 3; ++row) {
        R_SO(row, 0) = axes[0][row];
        R_SO(row, 1) = axes[1][row];
        R_SO(row, 2) = axes[2][row];
    }

    chrono::ChVector3d min_proj(std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity());
    chrono::ChVector3d max_proj(-std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity());
    for (uint32_t offset = 0; offset < count; ++offset) {
        const chrono::ChVector3d p_S = samples[ordered_indices[first + offset]].xi_slave_S;
        const chrono::ChVector3d p_O = MulT(R_SO, p_S);
        min_proj.x() = std::min(min_proj.x(), p_O.x());
        min_proj.y() = std::min(min_proj.y(), p_O.y());
        min_proj.z() = std::min(min_proj.z(), p_O.z());
        max_proj.x() = std::max(max_proj.x(), p_O.x());
        max_proj.y() = std::max(max_proj.y(), p_O.y());
        max_proj.z() = std::max(max_proj.z(), p_O.z());
    }

    const chrono::ChVector3d center_O = 0.5 * (min_proj + max_proj);
    obb.center_S = R_SO * center_O;
    obb.R_SO = R_SO;
    obb.half_extent = 0.5 * (max_proj - min_proj);
    obb.outer_radius = obb.half_extent.Length();
    return obb;
}

bool IntersectsSphere(const GeometryOBB& obb_S, const chrono::ChVector3d& center_S, double radius) {
    const chrono::ChVector3d rel_O = MulT(obb_S.R_SO, center_S - obb_S.center_S);
    const chrono::ChVector3d closest_O(Clamp(rel_O.x(), -obb_S.half_extent.x(), obb_S.half_extent.x()),
                                       Clamp(rel_O.y(), -obb_S.half_extent.y(), obb_S.half_extent.y()),
                                       Clamp(rel_O.z(), -obb_S.half_extent.z(), obb_S.half_extent.z()));
    return (closest_O - rel_O).Length2() <= radius * radius;
}

double ComputePredictiveSpeedBound(const RigidBodyStateW& master_state,
                                   const RigidBodyStateW& slave_state,
                                   const chrono::ChVector3d& node_center_W,
                                   double node_radius) {
    const double translational = (slave_state.v_com_W - master_state.v_com_W).Length();
    const double master_rot =
        master_state.w_W.Length() * ((node_center_W - master_state.x_com_W).Length() + node_radius);
    const double slave_rot =
        slave_state.w_W.Length() * ((node_center_W - slave_state.x_com_W).Length() + node_radius);
    return translational + master_rot + slave_rot;
}

}  // namespace

void DenseSampleBVH::Build(const std::vector<DenseSurfaceSample>& samples, int leaf_size) {
    leaf_size_ = std::max(4, leaf_size);
    root_outer_radius_ = 0.0;
    ordered_indices_.resize(samples.size());
    std::iota(ordered_indices_.begin(), ordered_indices_.end(), std::size_t(0));
    nodes_.clear();

    if (samples.empty()) {
        return;
    }

    const int root = BuildRecursive(samples, 0, static_cast<uint32_t>(samples.size()));
    if (root >= 0 && !nodes_.empty()) {
        root_outer_radius_ = nodes_[static_cast<std::size_t>(root)].obb_S.outer_radius;
    }
}

bool DenseSampleBVH::Empty() const {
    return nodes_.empty();
}

std::size_t DenseSampleBVH::SampleCount() const {
    return ordered_indices_.size();
}

double DenseSampleBVH::RootOuterRadius() const {
    return root_outer_radius_;
}

void DenseSampleBVH::CollectCandidateSampleIndices(const RigidBodyStateW& master_state,
                                                   const RigidBodyStateW& slave_state,
                                                   const FirstOrderSDF& sdf,
                                                   const CompressedContactConfig& cfg,
                                                   double step_size,
                                                   std::vector<std::size_t>& out_indices,
                                                   DenseSampleBVHQueryStats* out_stats) const {
    out_indices.clear();
    DenseSampleBVHQueryStats stats;

    if (nodes_.empty()) {
        if (out_stats) {
            *out_stats = stats;
        }
        return;
    }

    chrono::ChVector3d master_center_M;
    double master_radius = 0.0;
    const bool has_master_sphere = sdf.GetBoundingSphereM(master_center_M, master_radius);

    chrono::ChVector3d query_center_S(0.0, 0.0, 0.0);
    double query_radius = 0.0;
    if (has_master_sphere) {
        const chrono::ChVector3d master_center_W = master_state.x_ref_W + master_state.R_WRef * master_center_M;
        query_center_S = slave_state.R_WRef.transpose() * (master_center_W - slave_state.x_ref_W);

        const double activation_margin = std::max(cfg.delta_on, cfg.delta_off) + std::max(0.0, cfg.bvh_query_margin);
        query_radius = master_radius + activation_margin;
        if (cfg.predictive_gap) {
            query_radius +=
                step_size * ComputePredictiveSpeedBound(master_state, slave_state, master_center_W, master_radius);
        }
    }

    const double activation_margin = std::max(cfg.delta_on, cfg.delta_off) + std::max(0.0, cfg.bvh_query_margin);
    std::vector<int> stack;
    stack.push_back(0);
    while (!stack.empty()) {
        const int node_index = stack.back();
        stack.pop_back();

        const auto& node = nodes_[static_cast<std::size_t>(node_index)];
        ++stats.nodes_visited;

        if (has_master_sphere && !IntersectsSphere(node.obb_S, query_center_S, query_radius)) {
            ++stats.nodes_pruned_obb;
            continue;
        }

        if (cfg.bvh_enable_sdf_node_bound) {
            const chrono::ChVector3d node_center_W = slave_state.x_ref_W + slave_state.R_WRef * node.obb_S.center_S;
            const chrono::ChVector3d node_center_M =
                master_state.R_WRef.transpose() * (node_center_W - master_state.x_ref_W);

            double phi_center = 0.0;
            if (sdf.QueryPhiM(node_center_M, phi_center)) {
                double lower_bound = phi_center - node.obb_S.outer_radius;
                if (cfg.predictive_gap) {
                    lower_bound -= step_size * cfg.bvh_velocity_bound_scale *
                                   ComputePredictiveSpeedBound(master_state, slave_state, node_center_W,
                                                               node.obb_S.outer_radius);
                }
                if (lower_bound > activation_margin) {
                    ++stats.nodes_pruned_sdf;
                    continue;
                }
            }
        }

        if (node.is_leaf) {
            ++stats.leaf_nodes_visited;
            stats.leaf_samples_tested += node.count;
            for (uint32_t offset = 0; offset < node.count; ++offset) {
                out_indices.push_back(ordered_indices_[node.first + offset]);
            }
            continue;
        }

        if (node.right >= 0) {
            stack.push_back(node.right);
        }
        if (node.left >= 0) {
            stack.push_back(node.left);
        }
    }

    stats.candidate_samples = out_indices.size();
    if (out_stats) {
        *out_stats = stats;
    }
}

int DenseSampleBVH::BuildRecursive(const std::vector<DenseSurfaceSample>& samples, uint32_t first, uint32_t count) {
    if (count == 0) {
        return -1;
    }

    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{});
    nodes_[static_cast<std::size_t>(node_index)].first = first;
    nodes_[static_cast<std::size_t>(node_index)].count = count;
    nodes_[static_cast<std::size_t>(node_index)].obb_S = ComputeOBB(samples, ordered_indices_, first, count);
    nodes_[static_cast<std::size_t>(node_index)].is_leaf = (count <= static_cast<uint32_t>(leaf_size_));
    if (nodes_[static_cast<std::size_t>(node_index)].is_leaf) {
        return node_index;
    }

    int split_axis = 0;
    const auto& node_obb = nodes_[static_cast<std::size_t>(node_index)].obb_S;
    if (node_obb.half_extent.y() > node_obb.half_extent.x()) {
        split_axis = 1;
    }
    const double current_extent = (split_axis == 0) ? node_obb.half_extent.x() : node_obb.half_extent.y();
    if (node_obb.half_extent.z() > current_extent) {
        split_axis = 2;
    }

    const chrono::ChVector3d axis_S(node_obb.R_SO(0, split_axis),
                                    node_obb.R_SO(1, split_axis),
                                    node_obb.R_SO(2, split_axis));
    const uint32_t mid = first + count / 2;
    std::nth_element(
        ordered_indices_.begin() + first,
        ordered_indices_.begin() + mid,
        ordered_indices_.begin() + first + count,
        [&](std::size_t a, std::size_t b) {
            return Dot(samples[a].xi_slave_S, axis_S) < Dot(samples[b].xi_slave_S, axis_S);
        });

    if (mid == first || mid == first + count) {
        nodes_[static_cast<std::size_t>(node_index)].is_leaf = true;
        return node_index;
    }

    const int left = BuildRecursive(samples, first, mid - first);
    const int right = BuildRecursive(samples, mid, first + count - mid);
    nodes_[static_cast<std::size_t>(node_index)].left = left;
    nodes_[static_cast<std::size_t>(node_index)].right = right;
    nodes_[static_cast<std::size_t>(node_index)].is_leaf = (left < 0 || right < 0);
    return node_index;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
