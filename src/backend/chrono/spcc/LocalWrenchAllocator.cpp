#include "platform/backend/spcc/LocalWrenchAllocator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace platform {
namespace backend {
namespace spcc {

namespace {

constexpr int kFrictionPyramidEdges = 4;
constexpr double kSlipVelocityScale = 1.0e-3;

struct RayColumn {
    std::size_t support_index = 0;
    chrono::ChVector3d force_W;
    std::vector<double> column;
};

double ProxyLoad(const DenseContactPoint& point) {
    return std::max(0.0, -point.phi_eff) * std::max(0.0, point.area_weight);
}

double DotVector(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

bool SolveLinearSystem(std::vector<double> A,
                       std::vector<double> b,
                       std::vector<double>& x,
                       std::size_t n) {
    x.assign(n, 0.0);
    if (n == 0) {
        return true;
    }

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        double pivot_abs = std::abs(A[k * n + k]);
        for (std::size_t i = k + 1; i < n; ++i) {
            const double value_abs = std::abs(A[i * n + k]);
            if (value_abs > pivot_abs) {
                pivot_abs = value_abs;
                pivot = i;
            }
        }

        if (!(pivot_abs > 1.0e-12)) {
            return false;
        }

        if (pivot != k) {
            for (std::size_t j = k; j < n; ++j) {
                std::swap(A[k * n + j], A[pivot * n + j]);
            }
            std::swap(b[k], b[pivot]);
        }

        const double diag = A[k * n + k];
        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A[i * n + k] / diag;
            if (!(std::abs(factor) > 0.0)) {
                continue;
            }
            for (std::size_t j = k; j < n; ++j) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    for (std::size_t row = n; row-- > 0;) {
        double rhs = b[row];
        for (std::size_t col = row + 1; col < n; ++col) {
            rhs -= A[row * n + col] * x[col];
        }
        const double diag = A[row * n + row];
        if (!(std::abs(diag) > 1.0e-12)) {
            return false;
        }
        x[row] = rhs / diag;
    }

    return true;
}

chrono::ChVector3d SafeNormalized(const chrono::ChVector3d& v, const chrono::ChVector3d& fallback) {
    const double len = v.Length();
    if (!(len > 1.0e-12)) {
        return fallback;
    }
    return v * (1.0 / len);
}

void BuildSupportBasis(const SupportWrenchPoint& support,
                       chrono::ChVector3d& t1_W,
                       chrono::ChVector3d& t2_W) {
    t1_W = support.t1_W - chrono::Vdot(support.t1_W, support.n_W) * support.n_W;
    t1_W = SafeNormalized(t1_W, chrono::ChVector3d(1.0, 0.0, 0.0));
    if (std::abs(chrono::Vdot(t1_W, support.n_W)) > 1.0e-3) {
        const chrono::ChVector3d seed = (std::abs(support.n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0)
                                                                           : chrono::ChVector3d(1.0, 0.0, 0.0);
        t1_W = SafeNormalized(chrono::Vcross(seed, support.n_W), chrono::ChVector3d(1.0, 0.0, 0.0));
    }
    t2_W = SafeNormalized(chrono::Vcross(support.n_W, t1_W), chrono::ChVector3d(0.0, 1.0, 0.0));
    t1_W = SafeNormalized(chrono::Vcross(t2_W, support.n_W), t1_W);
}

std::vector<chrono::ChVector3d> BuildRayDirections(const SupportWrenchPoint& support) {
    std::vector<chrono::ChVector3d> rays;
    if (!(support.mu > 0.0)) {
        rays.push_back(support.n_W);
        return rays;
    }

    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildSupportBasis(support, t1_W, t2_W);
    rays.reserve(kFrictionPyramidEdges);
    for (int edge = 0; edge < kFrictionPyramidEdges; ++edge) {
        const double theta =
            (2.0 * std::acos(-1.0) * static_cast<double>(edge)) / static_cast<double>(kFrictionPyramidEdges);
        const chrono::ChVector3d tangent_dir_W = std::cos(theta) * t1_W + std::sin(theta) * t2_W;
        rays.push_back(support.n_W + support.mu * tangent_dir_W);
    }
    return rays;
}

std::vector<double> BuildScaledColumn(const chrono::ChVector3d& x_W,
                                      const chrono::ChVector3d& force_W,
                                      const chrono::ChVector3d& origin_W,
                                      double moment_scale) {
    const chrono::ChVector3d moment = chrono::Vcross(x_W - origin_W, force_W);
    return {
        force_W.x(),
        force_W.y(),
        force_W.z(),
        moment.x() / moment_scale,
        moment.y() / moment_scale,
        moment.z() / moment_scale,
    };
}

std::vector<double> BuildScaledTarget(const ReferenceWrench& reference, double moment_scale) {
    return {
        reference.force_W.x(),
        reference.force_W.y(),
        reference.force_W.z(),
        reference.moment_W.x() / moment_scale,
        reference.moment_W.y() / moment_scale,
        reference.moment_W.z() / moment_scale,
    };
}

double ComputeMomentScale(const std::vector<SupportWrenchPoint>& supports,
                          const ReferenceWrench& reference) {
    double radius = 0.0;
    for (const auto& support : supports) {
        radius = std::max(radius, (support.x_W - reference.origin_W).Length());
    }
    radius = std::max(radius, 1.0e-6);
    const double force_norm = std::max(reference.force_W.Length(), std::abs(reference.total_load));
    return std::max(1.0e-6, radius * std::max(force_norm, 1.0));
}

double EvaluateObjective(const std::vector<RayColumn>& rays,
                         const std::vector<double>& target,
                         const std::vector<double>& initial,
                         const std::vector<double>& ray_weights,
                         double regularization) {
    std::vector<double> residual(target.size(), -target[0]);
    for (std::size_t row = 0; row < target.size(); ++row) {
        residual[row] = -target[row];
    }
    for (std::size_t col = 0; col < rays.size(); ++col) {
        for (std::size_t row = 0; row < target.size(); ++row) {
            residual[row] += rays[col].column[row] * ray_weights[col];
        }
    }

    double objective = DotVector(residual, residual);
    if (regularization > 0.0) {
        for (std::size_t i = 0; i < ray_weights.size(); ++i) {
            const double delta = ray_weights[i] - initial[i];
            objective += regularization * delta * delta;
        }
    }
    return objective;
}

std::size_t CountBits(std::size_t mask) {
    std::size_t count = 0;
    while (mask != 0) {
        count += (mask & 1u);
        mask >>= 1u;
    }
    return count;
}

chrono::ChVector3d BuildTangentialProxyForce(const DenseContactPoint& point, double load, double mu_default) {
    if (!(mu_default > 0.0) || !(load > 0.0)) {
        return chrono::ChVector3d(0.0, 0.0, 0.0);
    }

    const chrono::ChVector3d tangential_v_W = point.v_rel_W - chrono::Vdot(point.v_rel_W, point.n_W) * point.n_W;
    const double tangential_speed = tangential_v_W.Length();
    if (!(tangential_speed > 1.0e-12)) {
        return chrono::ChVector3d(0.0, 0.0, 0.0);
    }

    const double slip_ratio = tangential_speed / (tangential_speed + kSlipVelocityScale);
    const chrono::ChVector3d tangential_dir_W = tangential_v_W * (1.0 / tangential_speed);
    return -(mu_default * slip_ratio * load) * tangential_dir_W;
}

}  // namespace

ReferenceWrench LocalWrenchAllocator::BuildDenseReference(const std::vector<DenseContactPoint>& dense_points,
                                                          const std::vector<std::size_t>& member_indices,
                                                          const chrono::ChVector3d& origin_W,
                                                          double mu_default) {
    ReferenceWrench reference;
    reference.origin_W = origin_W;
    reference.force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.total_load = 0.0;

    for (const auto member_index : member_indices) {
        const auto& point = dense_points[member_index];
        const double load = ProxyLoad(point);
        const chrono::ChVector3d normal_force_W = load * point.n_W;
        const chrono::ChVector3d tangential_force_W = BuildTangentialProxyForce(point, load, mu_default);
        const chrono::ChVector3d contact_force_W = normal_force_W + tangential_force_W;
        reference.force_W += contact_force_W;
        reference.moment_W += chrono::Vcross(point.x_W - origin_W, contact_force_W);
        reference.total_load += load;
    }

    return reference;
}

void LocalWrenchAllocator::Allocate(const std::vector<SupportWrenchPoint>& supports,
                                    const ReferenceWrench& reference,
                                    double temporal_regularization,
                                    WrenchAllocationResult& out_result) {
    out_result = WrenchAllocationResult{};
    out_result.loads.assign(supports.size(), 0.0);
    out_result.forces_W.assign(supports.size(), chrono::ChVector3d(0.0, 0.0, 0.0));
    if (supports.empty()) {
        return;
    }

    const double moment_scale = ComputeMomentScale(supports, reference);
    const auto target = BuildScaledTarget(reference, moment_scale);

    std::vector<RayColumn> rays;
    std::vector<double> initial_ray_weights;
    for (std::size_t support_index = 0; support_index < supports.size(); ++support_index) {
        const auto& support = supports[support_index];
        const auto ray_directions = BuildRayDirections(support);
        const double initial_load = std::max(0.0, support.initial_load);
        const double initial_ray_weight = initial_load / static_cast<double>(ray_directions.size());
        for (const auto& ray_force_W : ray_directions) {
            RayColumn ray;
            ray.support_index = support_index;
            ray.force_W = ray_force_W;
            ray.column = BuildScaledColumn(support.x_W, ray_force_W, reference.origin_W, moment_scale);
            rays.push_back(ray);
            initial_ray_weights.push_back(initial_ray_weight);
        }
    }
    out_result.rays_per_support = supports.empty() ? 0 : static_cast<int>(rays.size() / supports.size());

    const double regularization = std::max(0.0, temporal_regularization);
    const std::size_t n = rays.size();
    const std::size_t mask_limit = static_cast<std::size_t>(1) << n;
    double best_objective = std::numeric_limits<double>::infinity();
    std::vector<double> best_ray_weights(n, 0.0);
    bool found = false;

    for (std::size_t mask = 0; mask < mask_limit; ++mask) {
        const std::size_t active_count = CountBits(mask);
        if (active_count == 0) {
            const double objective =
                EvaluateObjective(rays, target, initial_ray_weights, best_ray_weights, regularization);
            if (objective < best_objective) {
                best_objective = objective;
                best_ray_weights.assign(n, 0.0);
                found = true;
            }
            continue;
        }

        std::vector<std::size_t> active_indices;
        active_indices.reserve(active_count);
        for (std::size_t i = 0; i < n; ++i) {
            if ((mask & (static_cast<std::size_t>(1) << i)) != 0) {
                active_indices.push_back(i);
            }
        }

        std::vector<double> normal_matrix(active_count * active_count, 0.0);
        std::vector<double> rhs(active_count, 0.0);
        for (std::size_t row = 0; row < active_count; ++row) {
            for (std::size_t col = 0; col < active_count; ++col) {
                normal_matrix[row * active_count + col] =
                    DotVector(rays[active_indices[row]].column, rays[active_indices[col]].column) +
                    (row == col ? regularization : 0.0);
            }
            rhs[row] = DotVector(rays[active_indices[row]].column, target) +
                       regularization * initial_ray_weights[active_indices[row]];
        }

        std::vector<double> active_loads;
        if (!SolveLinearSystem(normal_matrix, rhs, active_loads, active_count)) {
            continue;
        }

        bool nonnegative = true;
        std::vector<double> candidate_ray_weights(n, 0.0);
        for (std::size_t i = 0; i < active_count; ++i) {
            if (active_loads[i] < -1.0e-10) {
                nonnegative = false;
                break;
            }
            candidate_ray_weights[active_indices[i]] = std::max(0.0, active_loads[i]);
        }
        if (!nonnegative) {
            continue;
        }

        const double objective =
            EvaluateObjective(rays, target, initial_ray_weights, candidate_ray_weights, regularization);
        if (objective < best_objective) {
            best_objective = objective;
            best_ray_weights = std::move(candidate_ray_weights);
            found = true;
        }
    }

    out_result.ray_weights = best_ray_weights;
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const auto& ray = rays[ray_index];
        const double weight = best_ray_weights[ray_index];
        out_result.loads[ray.support_index] += weight;
        out_result.forces_W[ray.support_index] += weight * ray.force_W;
    }
    out_result.objective = found ? best_objective : 0.0;
    out_result.feasible = found;

    chrono::ChVector3d force_W(0.0, 0.0, 0.0);
    chrono::ChVector3d moment_W(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < supports.size(); ++i) {
        force_W += out_result.forces_W[i];
        moment_W += chrono::Vcross(supports[i].x_W - reference.origin_W, out_result.forces_W[i]);
    }

    const double force_scale = std::max(1.0e-12, reference.force_W.Length());
    const double moment_ref_norm = reference.moment_W.Length();
    const double moment_scale_norm = std::max(1.0e-12, moment_ref_norm);
    out_result.force_residual = (force_W - reference.force_W).Length() / force_scale;
    out_result.moment_residual = (moment_W - reference.moment_W).Length() / moment_scale_norm;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
