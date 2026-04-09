#include "platform/backend/spcc/SubpatchRefiner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

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

void BuildOrthonormalBasis(const chrono::ChVector3d& n_W,
                           chrono::ChVector3d& t1_W,
                           chrono::ChVector3d& t2_W) {
    const chrono::ChVector3d seed = (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0)
                                                               : chrono::ChVector3d(1.0, 0.0, 0.0);
    t1_W = SafeNormalized(chrono::Vcross(seed, n_W), chrono::ChVector3d(1.0, 0.0, 0.0));
    t2_W = SafeNormalized(chrono::Vcross(n_W, t1_W), chrono::ChVector3d(0.0, 1.0, 0.0));
}

double TangentialDistance(const chrono::ChVector3d& a_W,
                          const chrono::ChVector3d& b_W,
                          const chrono::ChVector3d& n_W) {
    const chrono::ChVector3d rel = a_W - b_W;
    const chrono::ChVector3d tangential = rel - chrono::Vdot(rel, n_W) * n_W;
    return tangential.Length();
}

DensePatch MakePatchGeometry(const std::vector<DenseContactPoint>& dense_points,
                             std::size_t patch_id,
                             std::vector<std::size_t> members) {
    DensePatch patch;
    patch.patch_id = patch_id;
    patch.members = std::move(members);
    if (patch.members.empty()) {
        return patch;
    }

    chrono::ChVector3d centroid_W(0.0, 0.0, 0.0);
    chrono::ChVector3d avg_normal_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto member_index : patch.members) {
        const auto& point = dense_points[member_index];
        const double weight = std::max(1.0e-12, point.area_weight);
        centroid_W += weight * point.x_W;
        avg_normal_W += weight * point.n_W;
        weight_sum += weight;
    }

    if (weight_sum > 1.0e-12) {
        centroid_W *= (1.0 / weight_sum);
    } else {
        centroid_W = dense_points[patch.members.front()].x_W;
    }

    patch.centroid_W = centroid_W;
    patch.avg_normal_W = SafeNormalized(avg_normal_W, dense_points[patch.members.front()].n_W);
    BuildOrthonormalBasis(patch.avg_normal_W, patch.t1_W, patch.t2_W);
    return patch;
}

struct GridKey {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const GridKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct GridKeyHash {
    std::size_t operator()(const GridKey& key) const {
        std::size_t seed = 1469598103934665603ull;
        seed ^= static_cast<std::size_t>(key.x + 0x9e3779b9);
        seed *= 1099511628211ull;
        seed ^= static_cast<std::size_t>(key.y + 0x9e3779b9);
        seed *= 1099511628211ull;
        seed ^= static_cast<std::size_t>(key.z + 0x9e3779b9);
        seed *= 1099511628211ull;
        return seed;
    }
};

double PatchAdjacencyRadius(const CompressedContactConfig& cfg) {
    if (cfg.patch_radius > 0.0) {
        return cfg.patch_radius;
    }
    if (cfg.max_patch_diameter > 0.0) {
        return cfg.max_patch_diameter;
    }
    return 0.0;
}

GridKey MakeGridKey(const chrono::ChVector3d& x_W, double cell_size) {
    return GridKey{static_cast<int>(std::floor(x_W.x() / cell_size)),
                   static_cast<int>(std::floor(x_W.y() / cell_size)),
                   static_cast<int>(std::floor(x_W.z() / cell_size))};
}

bool ArePatchAdjacent(const DenseContactPoint& a,
                      const DenseContactPoint& b,
                      const CompressedContactConfig& cfg,
                      double adjacency_radius) {
    if (chrono::Vdot(a.n_W, b.n_W) < cfg.normal_cos_min) {
        return false;
    }

    if (adjacency_radius <= 0.0) {
        return true;
    }

    chrono::ChVector3d avg_normal = SafeNormalized(a.n_W + b.n_W, a.n_W);
    const double tangential = TangentialDistance(a.x_W, b.x_W, avg_normal);
    if (tangential > adjacency_radius) {
        return false;
    }

    // Guard against linking points that are close tangentially but far apart along the surface normal.
    const double normal_sep = std::abs(chrono::Vdot(a.x_W - b.x_W, avg_normal));
    return normal_sep <= std::max(cfg.delta_off, cfg.delta_on) + adjacency_radius;
}

DenseSubpatch MakeSubpatchGeometry(const std::vector<DenseContactPoint>& dense_points,
                                   std::size_t patch_id,
                                   std::size_t subpatch_id,
                                   std::size_t depth,
                                   std::vector<std::size_t> members,
                                   const CompressedContactConfig& cfg) {
    DenseSubpatch subpatch;
    subpatch.patch_id = patch_id;
    subpatch.subpatch_id = subpatch_id;
    subpatch.depth = depth;
    subpatch.members = std::move(members);

    if (subpatch.members.empty()) {
        return subpatch;
    }

    chrono::ChVector3d centroid_W(0.0, 0.0, 0.0);
    chrono::ChVector3d avg_normal_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    for (const auto member_index : subpatch.members) {
        const auto& point = dense_points[member_index];
        const double weight = std::max(1.0e-12, point.area_weight);
        centroid_W += weight * point.x_W;
        avg_normal_W += weight * point.n_W;
        weight_sum += weight;
    }

    if (weight_sum > 1.0e-12) {
        centroid_W *= (1.0 / weight_sum);
    } else {
        centroid_W = dense_points[subpatch.members.front()].x_W;
    }

    avg_normal_W = SafeNormalized(avg_normal_W, dense_points[subpatch.members.front()].n_W);
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildOrthonormalBasis(avg_normal_W, t1_W, t2_W);

    double max_tangential = 0.0;
    double max_plane_error = 0.0;
    for (const auto member_index : subpatch.members) {
        const auto& point = dense_points[member_index];
        max_tangential = std::max(max_tangential, TangentialDistance(point.x_W, centroid_W, avg_normal_W));
        max_plane_error = std::max(max_plane_error, std::abs(chrono::Vdot(point.x_W - centroid_W, avg_normal_W)));
    }

    subpatch.centroid_W = centroid_W;
    subpatch.avg_normal_W = avg_normal_W;
    subpatch.t1_W = t1_W;
    subpatch.t2_W = t2_W;
    subpatch.diameter = 2.0 * max_tangential;
    subpatch.plane_error = max_plane_error;

    const double spacing = cfg.sentinel_spacing;
    if (spacing > 0.0) {
        double u_min = std::numeric_limits<double>::infinity();
        double u_max = -std::numeric_limits<double>::infinity();
        double v_min = std::numeric_limits<double>::infinity();
        double v_max = -std::numeric_limits<double>::infinity();
        for (const auto member_index : subpatch.members) {
            const chrono::ChVector3d rel = dense_points[member_index].x_W - centroid_W;
            const double u = chrono::Vdot(rel, t1_W);
            const double v = chrono::Vdot(rel, t2_W);
            u_min = std::min(u_min, u);
            u_max = std::max(u_max, u);
            v_min = std::min(v_min, v);
            v_max = std::max(v_max, v);
        }

        const double accept_radius = std::max(cfg.sentinel_margin + 0.75 * spacing, 1.0e-8);
        for (double u = u_min; u <= u_max + 1.0e-12; u += spacing) {
            for (double v = v_min; v <= v_max + 1.0e-12; v += spacing) {
                const chrono::ChVector3d candidate_W = centroid_W + u * t1_W + v * t2_W;
                double nearest = std::numeric_limits<double>::infinity();
                for (const auto member_index : subpatch.members) {
                    const auto& point = dense_points[member_index];
                    nearest = std::min(nearest, TangentialDistance(point.x_W, candidate_W, avg_normal_W));
                }
                if (nearest <= accept_radius) {
                    subpatch.sentinel_W.push_back(candidate_W);
                }
            }
        }
    }

    if (subpatch.sentinel_W.empty()) {
        subpatch.sentinel_W.push_back(centroid_W);
    }

    return subpatch;
}

bool ShouldSplit(const DenseSubpatch& subpatch,
                 const CompressedContactConfig& cfg) {
    if (cfg.max_subpatch_depth > 0 && static_cast<int>(subpatch.depth) >= cfg.max_subpatch_depth) {
        return false;
    }
    if (cfg.min_dense_points_per_subpatch > 0 &&
        static_cast<int>(subpatch.members.size()) < 2 * cfg.min_dense_points_per_subpatch) {
        return false;
    }

    bool diameter_bad = false;
    if (cfg.max_subpatch_diameter > 0.0) {
        diameter_bad = subpatch.diameter > cfg.max_subpatch_diameter;
    }

    bool plane_bad = false;
    if (cfg.max_plane_error > 0.0) {
        plane_bad = subpatch.plane_error > cfg.max_plane_error;
    }

    return diameter_bad || plane_bad;
}

void SplitMembers(const std::vector<DenseContactPoint>& dense_points,
                  const DenseSubpatch& subpatch,
                  std::vector<std::size_t>& left_members,
                  std::vector<std::size_t>& right_members) {
    std::vector<std::pair<double, std::size_t>> projected;
    projected.reserve(subpatch.members.size());

    double u_min = std::numeric_limits<double>::infinity();
    double u_max = -std::numeric_limits<double>::infinity();
    double v_min = std::numeric_limits<double>::infinity();
    double v_max = -std::numeric_limits<double>::infinity();
    for (const auto member_index : subpatch.members) {
        const chrono::ChVector3d rel = dense_points[member_index].x_W - subpatch.centroid_W;
        const double u = chrono::Vdot(rel, subpatch.t1_W);
        const double v = chrono::Vdot(rel, subpatch.t2_W);
        u_min = std::min(u_min, u);
        u_max = std::max(u_max, u);
        v_min = std::min(v_min, v);
        v_max = std::max(v_max, v);
    }

    const bool split_u = (u_max - u_min) >= (v_max - v_min);
    for (const auto member_index : subpatch.members) {
        const chrono::ChVector3d rel = dense_points[member_index].x_W - subpatch.centroid_W;
        const double key = split_u ? chrono::Vdot(rel, subpatch.t1_W) : chrono::Vdot(rel, subpatch.t2_W);
        projected.emplace_back(key, member_index);
    }

    const std::size_t mid = projected.size() / 2;
    std::nth_element(projected.begin(),
                     projected.begin() + static_cast<std::ptrdiff_t>(mid),
                     projected.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });

    left_members.clear();
    right_members.clear();
    left_members.reserve(mid);
    right_members.reserve(projected.size() - mid);
    for (std::size_t i = 0; i < projected.size(); ++i) {
        if (i < mid) {
            left_members.push_back(projected[i].second);
        } else {
            right_members.push_back(projected[i].second);
        }
    }
}

void BuildSubpatchesRecursive(const std::vector<DenseContactPoint>& dense_points,
                              const CompressedContactConfig& cfg,
                              std::size_t patch_id,
                              std::size_t depth,
                              std::vector<std::size_t> members,
                              std::size_t& next_subpatch_id,
                              std::vector<DenseSubpatch>& out_subpatches) {
    DenseSubpatch subpatch =
        MakeSubpatchGeometry(dense_points, patch_id, next_subpatch_id, depth, std::move(members), cfg);

    if (!ShouldSplit(subpatch, cfg)) {
        out_subpatches.push_back(std::move(subpatch));
        ++next_subpatch_id;
        return;
    }

    std::vector<std::size_t> left_members;
    std::vector<std::size_t> right_members;
    SplitMembers(dense_points, subpatch, left_members, right_members);
    if (left_members.empty() || right_members.empty()) {
        out_subpatches.push_back(std::move(subpatch));
        ++next_subpatch_id;
        return;
    }

    BuildSubpatchesRecursive(dense_points, cfg, patch_id, depth + 1, std::move(left_members),
                             next_subpatch_id, out_subpatches);
    BuildSubpatchesRecursive(dense_points, cfg, patch_id, depth + 1, std::move(right_members),
                             next_subpatch_id, out_subpatches);
}

}  // namespace

void SubpatchRefiner::BuildPatches(const std::vector<DenseContactPoint>& dense_points,
                                   const CompressedContactConfig& cfg,
                                   std::vector<DensePatch>& out_patches) {
    out_patches.clear();
    if (dense_points.empty()) {
        return;
    }

    const double adjacency_radius = PatchAdjacencyRadius(cfg);
    if (!(adjacency_radius > 0.0)) {
        std::vector<std::size_t> members(dense_points.size());
        for (std::size_t i = 0; i < dense_points.size(); ++i) {
            members[i] = i;
        }
        out_patches.push_back(MakePatchGeometry(dense_points, 0, std::move(members)));
        return;
    }

    std::unordered_map<GridKey, std::vector<std::size_t>, GridKeyHash> grid;
    grid.reserve(dense_points.size());
    for (std::size_t i = 0; i < dense_points.size(); ++i) {
        grid[MakeGridKey(dense_points[i].x_W, adjacency_radius)].push_back(i);
    }

    std::vector<char> visited(dense_points.size(), 0);

    for (std::size_t seed_index = 0; seed_index < dense_points.size(); ++seed_index) {
        if (visited[seed_index]) {
            continue;
        }

        std::vector<std::size_t> members;
        std::queue<std::size_t> frontier;
        frontier.push(seed_index);
        visited[seed_index] = 1;

        while (!frontier.empty()) {
            const std::size_t point_index = frontier.front();
            frontier.pop();
            members.push_back(point_index);

            const auto base_key = MakeGridKey(dense_points[point_index].x_W, adjacency_radius);
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const GridKey neighbor_key{base_key.x + dx, base_key.y + dy, base_key.z + dz};
                        const auto it = grid.find(neighbor_key);
                        if (it == grid.end()) {
                            continue;
                        }

                        for (const auto neighbor_index : it->second) {
                            if (visited[neighbor_index] || neighbor_index == point_index) {
                                continue;
                            }
                            if (!ArePatchAdjacent(dense_points[point_index], dense_points[neighbor_index], cfg,
                                                  adjacency_radius)) {
                                continue;
                            }
                            visited[neighbor_index] = 1;
                            frontier.push(neighbor_index);
                        }
                    }
                }
            }
        }

        out_patches.push_back(
            MakePatchGeometry(dense_points, out_patches.size(), std::move(members)));
    }
}

void SubpatchRefiner::BuildSubpatches(const std::vector<DenseContactPoint>& dense_points,
                                      const std::vector<DensePatch>& patches,
                                      const CompressedContactConfig& cfg,
                                      std::vector<DenseSubpatch>& out_subpatches) {
    out_subpatches.clear();
    std::size_t next_subpatch_id = 0;
    for (const auto& patch : patches) {
        BuildSubpatchesRecursive(dense_points, cfg, patch.patch_id, 0, patch.members, next_subpatch_id,
                                 out_subpatches);
    }
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
