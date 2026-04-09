#include "platform/backend/spcc/SampleBVH.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include <chrono/core/ChMatrix33.h>

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
    if (!std::isfinite(n) || n <= 1.0e-12) {
        return false;
    }
    v *= (1.0 / n);
    return true;
}

chrono::ChVector3d MultiplySymmetric3x3(const double A[3][3], const chrono::ChVector3d& x) {
    return chrono::ChVector3d(
        A[0][0] * x.x() + A[0][1] * x.y() + A[0][2] * x.z(),
        A[1][0] * x.x() + A[1][1] * x.y() + A[1][2] * x.z(),
        A[2][0] * x.x() + A[2][1] * x.y() + A[2][2] * x.z());
}

chrono::ChVector3d EstimatePrincipalAxis(const double A[3][3]) {
    chrono::ChVector3d axis(1.0, 0.0, 0.0);
    for (int iter = 0; iter < 8; ++iter) {
        chrono::ChVector3d next = MultiplySymmetric3x3(A, axis);
        if (!NormalizeOrReject(next)) {
            break;
        }
        axis = next;
    }
    if (!NormalizeOrReject(axis)) {
        return chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    return axis;
}

std::array<chrono::ChVector3d, 3> BuildOrthonormalBasis(const chrono::ChVector3d& primary) {
    chrono::ChVector3d e0 = primary;
    if (!NormalizeOrReject(e0)) {
        e0 = chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    chrono::ChVector3d seed = (std::abs(e0.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0)
                                                        : chrono::ChVector3d(0.0, 1.0, 0.0);
    chrono::ChVector3d e1 = chrono::Vcross(seed, e0);
    if (!NormalizeOrReject(e1)) {
        seed = chrono::ChVector3d(1.0, 0.0, 0.0);
        e1 = chrono::Vcross(seed, e0);
        if (!NormalizeOrReject(e1)) {
            e1 = chrono::ChVector3d(0.0, 1.0, 0.0);
        }
    }
    chrono::ChVector3d e2 = chrono::Vcross(e0, e1);
    if (!NormalizeOrReject(e2)) {
        e2 = chrono::ChVector3d(0.0, 0.0, 1.0);
    }
    return {e0, e1, e2};
}

}  // namespace

void SampleBVH::Build(const std::vector<chrono::ChVector3d>& samples_S, int leaf_size) {
    samples_S_ = &samples_S;
    perm_.resize(samples_S.size());
    std::iota(perm_.begin(), perm_.end(), std::uint32_t{0});
    nodes_.clear();
    leaf_size_ = std::max(1, leaf_size);
    if (!samples_S.empty()) {
        BuildRecursive(0, static_cast<std::uint32_t>(samples_S.size()));
    }
}

int SampleBVH::BuildRecursive(std::uint32_t begin, std::uint32_t end) {
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back(SampleBVHNode{});
    auto& node = nodes_[static_cast<std::size_t>(node_index)];
    node.begin = begin;
    node.end = end;
    node.obb_S = FitOBB(begin, end);
    node.probe_S = node.obb_S.center_S;
    node.probe_radius = node.obb_S.bounding_radius;

    if (samples_S_ != nullptr && begin < end) {
        double best_dist_sq = std::numeric_limits<double>::infinity();
        std::uint32_t best_sample = perm_[begin];
        for (std::uint32_t i = begin; i < end; ++i) {
            const std::uint32_t sample_index = perm_[i];
            const auto& p = (*samples_S_)[sample_index];
            const double dist_sq = (p - node.obb_S.center_S).Length2();
            if (dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best_sample = sample_index;
            }
        }

        node.probe_S = (*samples_S_)[best_sample];
        double max_dist_sq = 0.0;
        for (std::uint32_t i = begin; i < end; ++i) {
            const auto& p = (*samples_S_)[perm_[i]];
            max_dist_sq = std::max(max_dist_sq, (p - node.probe_S).Length2());
        }
        node.probe_radius = std::sqrt(max_dist_sq);
    }

    if (end <= begin || static_cast<int>(end - begin) <= leaf_size_) {
        node.leaf = true;
        return node_index;
    }

    const int axis_index =
        (node.obb_S.half_extents.x() >= node.obb_S.half_extents.y() &&
         node.obb_S.half_extents.x() >= node.obb_S.half_extents.z())
            ? 0
            : ((node.obb_S.half_extents.y() >= node.obb_S.half_extents.z()) ? 1 : 2);
    const chrono::ChVector3d split_axis = node.obb_S.axes_S[static_cast<std::size_t>(axis_index)];
    const chrono::ChVector3d center_S = node.obb_S.center_S;
    const std::uint32_t mid = begin + (end - begin) / 2;

    auto proj = [&](std::uint32_t idx) {
        return chrono::Vdot((*samples_S_)[idx] - center_S, split_axis);
    };

    std::nth_element(perm_.begin() + begin, perm_.begin() + mid, perm_.begin() + end,
                     [&](std::uint32_t a, std::uint32_t b) { return proj(a) < proj(b); });

    if (mid == begin || mid == end) {
        node.leaf = true;
        return node_index;
    }

    const int left = BuildRecursive(begin, mid);
    const int right = BuildRecursive(mid, end);
    nodes_[static_cast<std::size_t>(node_index)].left = left;
    nodes_[static_cast<std::size_t>(node_index)].right = right;
    return node_index;
}

OBBLocal SampleBVH::FitOBB(std::uint32_t begin, std::uint32_t end) const {
    OBBLocal obb;
    if (!samples_S_ || begin >= end) {
        return obb;
    }

    chrono::ChVector3d mean(0.0, 0.0, 0.0);
    double count = 0.0;
    for (std::uint32_t i = begin; i < end; ++i) {
        const auto& p = (*samples_S_)[perm_[i]];
        if (!IsFiniteVec(p)) {
            continue;
        }
        mean += p;
        count += 1.0;
    }
    if (!(count > 0.0)) {
        return obb;
    }
    mean *= (1.0 / count);

    double cov[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for (std::uint32_t i = begin; i < end; ++i) {
        const auto& p = (*samples_S_)[perm_[i]];
        if (!IsFiniteVec(p)) {
            continue;
        }
        const chrono::ChVector3d d = p - mean;
        cov[0][0] += d.x() * d.x();
        cov[0][1] += d.x() * d.y();
        cov[0][2] += d.x() * d.z();
        cov[1][1] += d.y() * d.y();
        cov[1][2] += d.y() * d.z();
        cov[2][2] += d.z() * d.z();
    }
    cov[1][0] = cov[0][1];
    cov[2][0] = cov[0][2];
    cov[2][1] = cov[1][2];

    const chrono::ChVector3d primary = EstimatePrincipalAxis(cov);
    obb.axes_S = BuildOrthonormalBasis(primary);

    double min_u[3] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity()};
    double max_u[3] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()};

    for (std::uint32_t i = begin; i < end; ++i) {
        const auto& p = (*samples_S_)[perm_[i]];
        if (!IsFiniteVec(p)) {
            continue;
        }
        const chrono::ChVector3d d = p - mean;
        for (int axis = 0; axis < 3; ++axis) {
            const double u = chrono::Vdot(d, obb.axes_S[static_cast<std::size_t>(axis)]);
            min_u[axis] = std::min(min_u[axis], u);
            max_u[axis] = std::max(max_u[axis], u);
        }
    }

    obb.center_S = mean;
    for (int axis = 0; axis < 3; ++axis) {
        const double mid = 0.5 * (min_u[axis] + max_u[axis]);
        const double half = 0.5 * (max_u[axis] - min_u[axis]);
        obb.center_S += obb.axes_S[static_cast<std::size_t>(axis)] * mid;
        if (axis == 0) {
            obb.half_extents.x() = half;
        } else if (axis == 1) {
            obb.half_extents.y() = half;
        } else {
            obb.half_extents.z() = half;
        }
    }
    obb.bounding_radius = std::sqrt(obb.half_extents.x() * obb.half_extents.x() +
                                    obb.half_extents.y() * obb.half_extents.y() +
                                    obb.half_extents.z() * obb.half_extents.z());
    return obb;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
