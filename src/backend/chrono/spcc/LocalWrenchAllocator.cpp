#include "platform/backend/spcc/LocalWrenchAllocator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace platform {
namespace backend {
namespace spcc {

namespace {

constexpr double kSlipVelocityScale = 1.0e-3;

struct RayColumn {
    std::size_t support_index = 0;
    chrono::ChVector3d force_W;
    std::vector<double> column;
};

struct DenseMicroDriveTarget {
    chrono::ChVector3d desired_force_W;
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    chrono::ChVector3d n_W;
    double normal_metric_weight = 0.0;
    double tangential_metric_weight = 0.0;
};

bool SolveProjectedRayQP(const std::vector<RayColumn>& rays,
                         const std::vector<double>& target,
                         const std::vector<double>& initial_ray_weights,
                         double regularization,
                         std::vector<double>& out_ray_weights,
                         double& out_objective);

bool RefineActiveRayWeights(const std::vector<RayColumn>& rays,
                           const std::vector<double>& target,
                           const std::vector<double>& initial_ray_weights,
                           double regularization,
                           std::vector<double>& io_ray_weights,
                           double& io_objective);

double DotVector(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
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

std::vector<chrono::ChVector3d> BuildRayDirections(const SupportWrenchPoint& support, int friction_ray_count) {
    std::vector<chrono::ChVector3d> rays;
    if (!(support.mu > 0.0)) {
        rays.push_back(support.n_W);
        return rays;
    }

    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildSupportBasis(support, t1_W, t2_W);
    const int edge_count = std::max(4, friction_ray_count);
    rays.reserve(edge_count);
    for (int edge = 0; edge < edge_count; ++edge) {
        const double theta =
            (2.0 * std::acos(-1.0) * static_cast<double>(edge)) / static_cast<double>(edge_count);
        const chrono::ChVector3d tangent_dir_W = std::cos(theta) * t1_W + std::sin(theta) * t2_W;
        rays.push_back(support.n_W + support.mu * tangent_dir_W);
    }
    return rays;
}

std::vector<double> BuildInitialRayWeights(const std::vector<chrono::ChVector3d>& ray_directions,
                                           const SupportWrenchPoint& support) {
    std::vector<double> ray_weights(ray_directions.size(), 0.0);
    if (ray_directions.empty()) {
        return ray_weights;
    }

    const chrono::ChVector3d initial_force_W = support.initial_force_W;
    if (initial_force_W.Length2() <= 1.0e-24) {
        const double initial_load = std::max(0.0, support.initial_load);
        const double initial_ray_weight = initial_load / static_cast<double>(ray_directions.size());
        std::fill(ray_weights.begin(), ray_weights.end(), initial_ray_weight);
        return ray_weights;
    }

    std::vector<RayColumn> force_rays;
    force_rays.reserve(ray_directions.size());
    for (const auto& ray_direction_W : ray_directions) {
        RayColumn ray;
        ray.force_W = ray_direction_W;
        ray.column = {ray_direction_W.x(), ray_direction_W.y(), ray_direction_W.z()};
        force_rays.push_back(std::move(ray));
    }

    const std::vector<double> target = {initial_force_W.x(), initial_force_W.y(), initial_force_W.z()};
    std::vector<double> seed(ray_directions.size(), 0.0);
    const double initial_load = std::max(0.0, support.initial_load);
    const double initial_ray_weight = initial_load / static_cast<double>(ray_directions.size());
    std::fill(seed.begin(), seed.end(), initial_ray_weight);

    double objective = 0.0;
    if (!SolveProjectedRayQP(force_rays, target, seed, 0.0, ray_weights, objective)) {
        return seed;
    }
    RefineActiveRayWeights(force_rays, target, seed, 0.0, ray_weights, objective);
    return ray_weights;
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

bool SolveDenseLinearSystem(std::vector<std::vector<double>> matrix,
                            std::vector<double> rhs,
                            std::vector<double>& solution) {
    const std::size_t n = matrix.size();
    if (rhs.size() != n) {
        return false;
    }

    for (std::size_t pivot = 0; pivot < n; ++pivot) {
        std::size_t best_row = pivot;
        double best_abs = std::abs(matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row) {
            const double value_abs = std::abs(matrix[row][pivot]);
            if (value_abs > best_abs) {
                best_abs = value_abs;
                best_row = row;
            }
        }
        if (!(best_abs > 1.0e-14)) {
            return false;
        }

        if (best_row != pivot) {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        const double diag = matrix[pivot][pivot];
        for (std::size_t col = pivot; col < n; ++col) {
            matrix[pivot][col] /= diag;
        }
        rhs[pivot] /= diag;

        for (std::size_t row = pivot + 1; row < n; ++row) {
            const double factor = matrix[row][pivot];
            if (!(std::abs(factor) > 0.0)) {
                continue;
            }
            for (std::size_t col = pivot; col < n; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    solution.assign(n, 0.0);
    for (std::size_t row = n; row-- > 0;) {
        double value = rhs[row];
        for (std::size_t col = row + 1; col < n; ++col) {
            value -= matrix[row][col] * solution[col];
        }
        solution[row] = value;
    }
    return true;
}

double RayColumnNormSquared(const RayColumn& ray) {
    return DotVector(ray.column, ray.column);
}

double EstimateLipschitz(const std::vector<RayColumn>& rays, double regularization) {
    double trace_bound = regularization;
    for (const auto& ray : rays) {
        trace_bound += RayColumnNormSquared(ray);
    }
    return std::max(trace_bound, 1.0e-6);
}

void ComputeResidual(const std::vector<RayColumn>& rays,
                     const std::vector<double>& target,
                     const std::vector<double>& ray_weights,
                     std::vector<double>& residual) {
    residual.assign(target.size(), 0.0);
    for (std::size_t row = 0; row < target.size(); ++row) {
        residual[row] = -target[row];
    }
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const double weight = ray_weights[ray_index];
        if (!(std::abs(weight) > 0.0)) {
            continue;
        }
        for (std::size_t row = 0; row < residual.size(); ++row) {
            residual[row] += rays[ray_index].column[row] * weight;
        }
    }
}

void ComputeGradient(const std::vector<RayColumn>& rays,
                     const std::vector<double>& residual,
                     const std::vector<double>& ray_weights,
                     const std::vector<double>& initial_ray_weights,
                     double regularization,
                     std::vector<double>& gradient) {
    gradient.assign(ray_weights.size(), 0.0);
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        gradient[ray_index] = DotVector(rays[ray_index].column, residual);
        if (regularization > 0.0) {
            gradient[ray_index] += regularization * (ray_weights[ray_index] - initial_ray_weights[ray_index]);
        }
    }
}

bool SolveProjectedRayQP(const std::vector<RayColumn>& rays,
                         const std::vector<double>& target,
                         const std::vector<double>& initial_ray_weights,
                         double regularization,
                         std::vector<double>& out_ray_weights,
                         double& out_objective) {
    const std::size_t n = rays.size();
    out_ray_weights.assign(n, 0.0);
    if (n == 0) {
        out_objective = 0.0;
        return true;
    }

    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::max(0.0, initial_ray_weights[i]);
    }
    std::vector<double> y = x;
    std::vector<double> x_next(n, 0.0);
    std::vector<double> residual;
    std::vector<double> gradient;
    const double step = 1.0 / EstimateLipschitz(rays, regularization);
    double t = 1.0;

    constexpr int kMaxIterations = 4000;
    constexpr double kTolerance = 1.0e-12;
    for (int iter = 0; iter < kMaxIterations; ++iter) {
        ComputeResidual(rays, target, y, residual);
        ComputeGradient(rays, residual, y, initial_ray_weights, regularization, gradient);

        double delta_sq = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
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
        for (std::size_t i = 0; i < n; ++i) {
            y[i] = x_next[i] + momentum * (x_next[i] - x[i]);
        }
        x.swap(x_next);
        t = t_next;
    }

    out_ray_weights = x;
    out_objective = EvaluateObjective(rays, target, initial_ray_weights, out_ray_weights, regularization);
    return true;
}

bool RefineActiveRayWeights(const std::vector<RayColumn>& rays,
                           const std::vector<double>& target,
                           const std::vector<double>& initial_ray_weights,
                           double regularization,
                           std::vector<double>& io_ray_weights,
                           double& io_objective) {
    constexpr double kActiveThreshold = 1.0e-10;
    constexpr double kNegativeTolerance = 1.0e-12;

    std::vector<std::size_t> active_indices;
    active_indices.reserve(io_ray_weights.size());
    for (std::size_t i = 0; i < io_ray_weights.size(); ++i) {
        if (io_ray_weights[i] > kActiveThreshold) {
            active_indices.push_back(i);
        }
    }
    if (active_indices.empty()) {
        return false;
    }

    while (!active_indices.empty()) {
        const std::size_t active_count = active_indices.size();
        std::vector<std::vector<double>> normal_matrix(active_count, std::vector<double>(active_count, 0.0));
        std::vector<double> rhs(active_count, 0.0);

        for (std::size_t i = 0; i < active_count; ++i) {
            const auto active_i = active_indices[i];
            rhs[i] = DotVector(rays[active_i].column, target);
            if (regularization > 0.0) {
                rhs[i] += regularization * initial_ray_weights[active_i];
            }
            for (std::size_t j = i; j < active_count; ++j) {
                const auto active_j = active_indices[j];
                double value = DotVector(rays[active_i].column, rays[active_j].column);
                if (regularization > 0.0 && i == j) {
                    value += regularization;
                }
                normal_matrix[i][j] = value;
                normal_matrix[j][i] = value;
            }
        }

        std::vector<double> active_solution;
        if (!SolveDenseLinearSystem(std::move(normal_matrix), std::move(rhs), active_solution)) {
            return false;
        }

        std::vector<std::size_t> next_active;
        next_active.reserve(active_count);
        bool all_nonnegative = true;
        std::vector<double> candidate(io_ray_weights.size(), 0.0);
        for (std::size_t i = 0; i < active_count; ++i) {
            const auto active_i = active_indices[i];
            const double value = active_solution[i];
            if (value > kNegativeTolerance) {
                candidate[active_i] = value;
                next_active.push_back(active_i);
            } else {
                all_nonnegative = false;
            }
        }

        if (all_nonnegative) {
            const double refined_objective =
                EvaluateObjective(rays, target, initial_ray_weights, candidate, regularization);
            if (refined_objective <= io_objective + 1.0e-12) {
                io_ray_weights = std::move(candidate);
                io_objective = refined_objective;
                return true;
            }
            return false;
        }

        active_indices.swap(next_active);
    }

    return false;
}

double EvaluateDenseMicroObjective(const std::vector<RayColumn>& rays,
                                   const std::vector<DenseMicroDriveTarget>& drive_targets,
                                   const std::vector<double>& drive_target,
                                   const std::vector<double>& initial_ray_weights,
                                   const DenseMicroSolverOptions& options,
                                   const std::vector<double>& ray_weights) {
    std::vector<chrono::ChVector3d> support_forces_W(drive_targets.size(), chrono::ChVector3d(0.0, 0.0, 0.0));
    std::vector<double> wrench_residual(drive_target.size(), 0.0);
    for (std::size_t row = 0; row < drive_target.size(); ++row) {
        wrench_residual[row] = -drive_target[row];
    }
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const double weight = ray_weights[ray_index];
        if (!(std::abs(weight) > 0.0)) {
            continue;
        }
        support_forces_W[rays[ray_index].support_index] += weight * rays[ray_index].force_W;
        for (std::size_t row = 0; row < wrench_residual.size(); ++row) {
            wrench_residual[row] += weight * rays[ray_index].column[row];
        }
    }

    double objective = 0.0;
    for (std::size_t support_index = 0; support_index < drive_targets.size(); ++support_index) {
        const auto& target = drive_targets[support_index];
        const chrono::ChVector3d residual_W = support_forces_W[support_index] - target.desired_force_W;
        const double residual_t1 = chrono::Vdot(residual_W, target.t1_W);
        const double residual_t2 = chrono::Vdot(residual_W, target.t2_W);
        const double residual_n = chrono::Vdot(residual_W, target.n_W);
        objective += 0.5 * target.tangential_metric_weight *
                     (residual_t1 * residual_t1 + residual_t2 * residual_t2);
        objective += 0.5 * target.normal_metric_weight * residual_n * residual_n;
    }
    if (options.wrench_coupling_weight > 0.0) {
        objective += 0.5 * options.wrench_coupling_weight * DotVector(wrench_residual, wrench_residual);
    }
    if (options.regularization > 0.0) {
        for (std::size_t i = 0; i < ray_weights.size(); ++i) {
            const double delta = ray_weights[i] - initial_ray_weights[i];
            objective += 0.5 * options.regularization * delta * delta;
        }
    }
    return objective;
}

double EstimateDenseMicroLipschitz(const std::vector<RayColumn>& rays,
                                   const std::vector<std::size_t>& support_offsets,
                                   const std::vector<DenseMicroDriveTarget>& drive_targets,
                                   const DenseMicroSolverOptions& options) {
    double max_support_trace = 0.0;
    for (std::size_t support_index = 0; support_index + 1 < support_offsets.size(); ++support_index) {
        double support_trace = 0.0;
        const double metric_scale = drive_targets[support_index].normal_metric_weight +
                                    2.0 * drive_targets[support_index].tangential_metric_weight;
        for (std::size_t ray_index = support_offsets[support_index];
             ray_index < support_offsets[support_index + 1];
             ++ray_index) {
            support_trace += metric_scale * rays[ray_index].force_W.Length2();
        }
        max_support_trace = std::max(max_support_trace, support_trace);
    }

    double wrench_trace = 0.0;
    if (options.wrench_coupling_weight > 0.0) {
        for (const auto& ray : rays) {
            wrench_trace += DotVector(ray.column, ray.column);
        }
        wrench_trace *= options.wrench_coupling_weight;
    }
    return std::max(1.0e-6, max_support_trace + wrench_trace + options.regularization);
}

void ComputeDenseMicroGradient(const std::vector<RayColumn>& rays,
                               const std::vector<DenseMicroDriveTarget>& drive_targets,
                               const std::vector<double>& drive_target,
                               const std::vector<double>& initial_ray_weights,
                               const DenseMicroSolverOptions& options,
                               const std::vector<double>& ray_weights,
                               std::vector<double>& out_gradient) {
    std::vector<chrono::ChVector3d> support_residuals_W(drive_targets.size(), chrono::ChVector3d(0.0, 0.0, 0.0));
    std::vector<double> wrench_residual(drive_target.size(), 0.0);
    for (std::size_t row = 0; row < drive_target.size(); ++row) {
        wrench_residual[row] = -drive_target[row];
    }
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const double weight = ray_weights[ray_index];
        if (!(std::abs(weight) > 0.0)) {
            continue;
        }
        support_residuals_W[rays[ray_index].support_index] += weight * rays[ray_index].force_W;
        for (std::size_t row = 0; row < wrench_residual.size(); ++row) {
            wrench_residual[row] += weight * rays[ray_index].column[row];
        }
    }
    for (std::size_t support_index = 0; support_index < drive_targets.size(); ++support_index) {
        support_residuals_W[support_index] -= drive_targets[support_index].desired_force_W;
    }

    out_gradient.assign(ray_weights.size(), 0.0);
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const auto& target = drive_targets[rays[ray_index].support_index];
        const auto& residual_W = support_residuals_W[rays[ray_index].support_index];
        double value =
            target.tangential_metric_weight *
                (chrono::Vdot(rays[ray_index].force_W, target.t1_W) * chrono::Vdot(residual_W, target.t1_W) +
                 chrono::Vdot(rays[ray_index].force_W, target.t2_W) * chrono::Vdot(residual_W, target.t2_W)) +
            target.normal_metric_weight *
                chrono::Vdot(rays[ray_index].force_W, target.n_W) * chrono::Vdot(residual_W, target.n_W);
        if (options.wrench_coupling_weight > 0.0) {
            value += options.wrench_coupling_weight * DotVector(rays[ray_index].column, wrench_residual);
        }
        if (options.regularization > 0.0) {
            value += options.regularization * (ray_weights[ray_index] - initial_ray_weights[ray_index]);
        }
        out_gradient[ray_index] = value;
    }
}

bool SolveDenseMicroRayQP(const std::vector<RayColumn>& rays,
                          const std::vector<std::size_t>& support_offsets,
                          const std::vector<DenseMicroDriveTarget>& drive_targets,
                          const std::vector<double>& drive_target,
                          const std::vector<double>& initial_ray_weights,
                          const DenseMicroSolverOptions& options,
                          std::vector<double>& out_ray_weights,
                          double& out_objective) {
    const std::size_t n = rays.size();
    out_ray_weights.assign(n, 0.0);
    if (n == 0) {
        out_objective = 0.0;
        return true;
    }

    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        x[i] = std::max(0.0, initial_ray_weights[i]);
    }
    std::vector<double> y = x;
    std::vector<double> x_next(n, 0.0);
    std::vector<double> gradient;
    const double step = 1.0 / EstimateDenseMicroLipschitz(rays, support_offsets, drive_targets, options);
    double t = 1.0;

    constexpr int kMaxIterations = 2000;
    constexpr double kTolerance = 1.0e-10;
    for (int iter = 0; iter < kMaxIterations; ++iter) {
        ComputeDenseMicroGradient(rays, drive_targets, drive_target, initial_ray_weights, options, y, gradient);

        double delta_sq = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
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
        for (std::size_t i = 0; i < n; ++i) {
            y[i] = x_next[i] + momentum * (x_next[i] - x[i]);
        }
        x.swap(x_next);
        t = t_next;
    }

    out_ray_weights = x;
    out_objective = EvaluateDenseMicroObjective(rays, drive_targets, drive_target, initial_ray_weights, options,
                                                out_ray_weights);
    return true;
}

DenseMicroDriveTarget BuildDenseMicroDriveTarget(const DenseContactPoint& point,
                                                 double mu_default,
                                                 double step_size,
                                                 const DenseMicroSolverOptions& options) {
    DenseMicroDriveTarget target;
    target.n_W = point.n_W;

    const chrono::ChVector3d tangential_v_W = point.v_rel_W - chrono::Vdot(point.v_rel_W, point.n_W) * point.n_W;
    const double tangential_speed = tangential_v_W.Length();
    const chrono::ChVector3d tangential_dir_W =
        (tangential_speed > 1.0e-12) ? (tangential_v_W * (1.0 / tangential_speed))
                                     : chrono::ChVector3d(0.0, 0.0, 0.0);

    target.t1_W = (tangential_dir_W.Length2() > 0.0) ? tangential_dir_W : point.v_rel_W;
    if (!(target.t1_W.Length2() > 0.0)) {
        target.t1_W = chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    target.t1_W = target.t1_W - chrono::Vdot(target.t1_W, point.n_W) * point.n_W;
    const double t1_len = target.t1_W.Length();
    if (t1_len > 1.0e-12) {
        target.t1_W *= (1.0 / t1_len);
    } else {
        const chrono::ChVector3d seed =
            (std::abs(point.n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0)
                                            : chrono::ChVector3d(1.0, 0.0, 0.0);
        target.t1_W = SafeNormalized(chrono::Vcross(seed, point.n_W), chrono::ChVector3d(1.0, 0.0, 0.0));
    }
    target.t2_W = SafeNormalized(chrono::Vcross(point.n_W, target.t1_W), chrono::ChVector3d(0.0, 1.0, 0.0));

    const double area = std::max(point.area_weight, 1.0e-12);
    const double approach_speed = std::max(0.0, -chrono::Vdot(point.v_rel_W, point.n_W));
    const double slip_ratio = tangential_speed / (tangential_speed + kSlipVelocityScale);
    const double normal_drive =
        area * (options.gap_drive_weight * std::max(0.0, -point.phi_eff) +
                options.approach_drive_weight * step_size * approach_speed);
    const chrono::ChVector3d tangential_drive_W =
        -(options.slip_drive_weight * mu_default * slip_ratio * normal_drive) * target.t1_W;
    target.desired_force_W = normal_drive * point.n_W + tangential_drive_W;
    target.normal_metric_weight = std::max(1.0e-8, options.normal_response_weight * area);
    target.tangential_metric_weight = std::max(1.0e-8, options.tangential_response_weight * area);
    return target;
}

}  // namespace

void LocalWrenchAllocator::BuildDenseMicroReference(const std::vector<DenseContactPoint>& dense_points,
                                                    const std::vector<std::size_t>& member_indices,
                                                    const chrono::ChVector3d& origin_W,
                                                    double mu_default,
                                                    double step_size,
                                                    const DenseMicroSolverOptions& options,
                                                    DenseMicroReferenceResult& out_result) {
    out_result = DenseMicroReferenceResult{};
    out_result.reference.origin_W = origin_W;

    std::vector<SupportWrenchPoint> dense_supports;
    std::vector<DenseMicroDriveTarget> drive_targets;
    dense_supports.reserve(member_indices.size());
    drive_targets.reserve(member_indices.size());

    ReferenceWrench drive_reference;
    drive_reference.origin_W = origin_W;
    drive_reference.force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    drive_reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    drive_reference.total_load = 0.0;

    for (const auto member_index : member_indices) {
        const auto& point = dense_points[member_index];
        const auto drive_target = BuildDenseMicroDriveTarget(point, mu_default, step_size, options);
        drive_reference.force_W += drive_target.desired_force_W;
        drive_reference.moment_W += chrono::Vcross(point.x_W - origin_W, drive_target.desired_force_W);
        drive_reference.total_load += std::max(0.0, chrono::Vdot(drive_target.desired_force_W, point.n_W));

        SupportWrenchPoint support;
        support.x_W = point.x_W;
        support.n_W = point.n_W;
        support.t1_W = drive_target.t1_W;
        if (!(support.t1_W.Length2() > 0.0)) {
            support.t1_W = point.v_rel_W;
        }
        support.t2_W = chrono::Vcross(point.n_W, support.t1_W);
        support.mu = mu_default;
        support.initial_load = std::max(0.0, chrono::Vdot(drive_target.desired_force_W, point.n_W));
        support.initial_force_W = drive_target.desired_force_W;
        dense_supports.push_back(support);
        drive_targets.push_back(drive_target);
    }

    if (dense_supports.empty()) {
        return;
    }

    const double moment_scale = ComputeMomentScale(dense_supports, drive_reference);
    const auto drive_target = BuildScaledTarget(drive_reference, moment_scale);
    std::vector<RayColumn> rays;
    std::vector<std::size_t> support_offsets;
    std::vector<double> initial_ray_weights;
    support_offsets.reserve(dense_supports.size() + 1);
    support_offsets.push_back(0);
    for (std::size_t support_index = 0; support_index < dense_supports.size(); ++support_index) {
        const auto ray_directions = BuildRayDirections(dense_supports[support_index], options.friction_ray_count);
        const auto support_initial_ray_weights = BuildInitialRayWeights(ray_directions, dense_supports[support_index]);
        for (const auto& ray_force_W : ray_directions) {
            RayColumn ray;
            ray.support_index = support_index;
            ray.force_W = ray_force_W;
            ray.column = BuildScaledColumn(dense_supports[support_index].x_W, ray_force_W, origin_W, moment_scale);
            rays.push_back(std::move(ray));
        }
        initial_ray_weights.insert(initial_ray_weights.end(), support_initial_ray_weights.begin(),
                                   support_initial_ray_weights.end());
        if (support_initial_ray_weights.size() != ray_directions.size()) {
            const double initial_load = std::max(0.0, dense_supports[support_index].initial_load);
            const double initial_ray_weight = initial_load / static_cast<double>(ray_directions.size());
            for (std::size_t ray_index = support_initial_ray_weights.size(); ray_index < ray_directions.size();
                 ++ray_index) {
                initial_ray_weights.push_back(initial_ray_weight);
            }
        }
        support_offsets.push_back(rays.size());
    }

    std::vector<double> dense_ray_weights;
    double dense_objective = 0.0;
    const bool solved = SolveDenseMicroRayQP(rays, support_offsets, drive_targets, drive_target, initial_ray_weights,
                                             options,
                                             dense_ray_weights, dense_objective);
    out_result.feasible = solved && dense_ray_weights.size() == rays.size();
    if (!out_result.feasible) {
        return;
    }

    ReferenceWrench reference;
    reference.origin_W = origin_W;
    reference.force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.total_load = 0.0;
    for (std::size_t ray_index = 0; ray_index < rays.size(); ++ray_index) {
        const double weight = dense_ray_weights[ray_index];
        if (!(weight > 0.0)) {
            continue;
        }
        const auto& ray = rays[ray_index];
        const chrono::ChVector3d force_W = weight * ray.force_W;
        reference.force_W += force_W;
        reference.moment_W += chrono::Vcross(dense_supports[ray.support_index].x_W - origin_W, force_W);
        reference.total_load += std::max(0.0, chrono::Vdot(force_W, dense_supports[ray.support_index].n_W));
    }
    out_result.reference = reference;

    const double force_scale = std::max(drive_reference.force_W.Length(), 1.0e-12);
    const double radius_scale = std::max(moment_scale, 1.0e-6);
    const double moment_scale_norm =
        std::max({drive_reference.moment_W.Length(), drive_reference.force_W.Length() * radius_scale, 1.0e-12});
    out_result.force_residual = (reference.force_W - drive_reference.force_W).Length() / force_scale;
    out_result.moment_residual = (reference.moment_W - drive_reference.moment_W).Length() / moment_scale_norm;
}

void LocalWrenchAllocator::Allocate(const std::vector<SupportWrenchPoint>& supports,
                                    const ReferenceWrench& reference,
                                    const ReducedSolveOptions& options,
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
        const auto ray_directions = BuildRayDirections(support, options.friction_ray_count);
        const auto support_initial_ray_weights = BuildInitialRayWeights(ray_directions, support);
        for (const auto& ray_force_W : ray_directions) {
            RayColumn ray;
            ray.support_index = support_index;
            ray.force_W = ray_force_W;
            ray.column = BuildScaledColumn(support.x_W, ray_force_W, reference.origin_W, moment_scale);
            rays.push_back(ray);
        }
        initial_ray_weights.insert(initial_ray_weights.end(), support_initial_ray_weights.begin(),
                                   support_initial_ray_weights.end());
        if (support_initial_ray_weights.size() != ray_directions.size()) {
            const double initial_load = std::max(0.0, support.initial_load);
            const double initial_ray_weight = initial_load / static_cast<double>(ray_directions.size());
            for (std::size_t ray_index = support_initial_ray_weights.size(); ray_index < ray_directions.size();
                 ++ray_index) {
                initial_ray_weights.push_back(initial_ray_weight);
            }
        }
    }
    out_result.rays_per_support = supports.empty() ? 0 : static_cast<int>(rays.size() / supports.size());

    const double regularization = std::max(0.0, options.temporal_regularization);
    double best_objective = 0.0;
    std::vector<double> best_ray_weights;
    const bool found =
        SolveProjectedRayQP(rays, target, initial_ray_weights, regularization, best_ray_weights, best_objective);
    if (found) {
        RefineActiveRayWeights(rays, target, initial_ray_weights, regularization, best_ray_weights, best_objective);
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
