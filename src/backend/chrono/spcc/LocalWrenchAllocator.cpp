#include "platform/backend/spcc/LocalWrenchAllocator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace platform {
namespace backend {
namespace spcc {

namespace {

double ProxyLoad(const DenseContactPoint& point) {
    return std::max(0.0, -point.phi_eff) * std::max(0.0, point.area_weight);
}

double Dot6(const std::vector<double>& a, const std::vector<double>& b) {
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

std::vector<double> BuildScaledColumn(const SupportWrenchPoint& support,
                                      const chrono::ChVector3d& origin_W,
                                      double moment_scale) {
    const chrono::ChVector3d moment = chrono::Vcross(support.x_W - origin_W, support.n_W);
    return {
        support.n_W.x(),
        support.n_W.y(),
        support.n_W.z(),
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

double EvaluateObjective(const std::vector<std::vector<double>>& columns,
                         const std::vector<double>& target,
                         const std::vector<double>& initial,
                         const std::vector<double>& loads,
                         double regularization) {
    std::vector<double> residual(target.size(), -target[0]);
    for (std::size_t row = 0; row < target.size(); ++row) {
        residual[row] = -target[row];
    }
    for (std::size_t col = 0; col < columns.size(); ++col) {
        for (std::size_t row = 0; row < target.size(); ++row) {
            residual[row] += columns[col][row] * loads[col];
        }
    }

    double objective = Dot6(residual, residual);
    if (regularization > 0.0) {
        for (std::size_t i = 0; i < loads.size(); ++i) {
            const double delta = loads[i] - initial[i];
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

}  // namespace

ReferenceWrench LocalWrenchAllocator::BuildDenseReference(const std::vector<DenseContactPoint>& dense_points,
                                                          const std::vector<std::size_t>& member_indices,
                                                          const chrono::ChVector3d& origin_W) {
    ReferenceWrench reference;
    reference.origin_W = origin_W;
    reference.force_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.moment_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    reference.total_load = 0.0;

    for (const auto member_index : member_indices) {
        const auto& point = dense_points[member_index];
        const double load = ProxyLoad(point);
        reference.force_W += load * point.n_W;
        reference.moment_W += chrono::Vcross(point.x_W - origin_W, load * point.n_W);
        reference.total_load += load;
    }

    return reference;
}

void LocalWrenchAllocator::Allocate(const std::vector<SupportWrenchPoint>& supports,
                                    const ReferenceWrench& reference,
                                    WrenchAllocationResult& out_result) {
    out_result = WrenchAllocationResult{};
    out_result.loads.assign(supports.size(), 0.0);
    if (supports.empty()) {
        return;
    }

    const double moment_scale = ComputeMomentScale(supports, reference);
    const auto target = BuildScaledTarget(reference, moment_scale);

    std::vector<std::vector<double>> columns;
    columns.reserve(supports.size());
    std::vector<double> initial_loads;
    initial_loads.reserve(supports.size());
    for (const auto& support : supports) {
        columns.push_back(BuildScaledColumn(support, reference.origin_W, moment_scale));
        initial_loads.push_back(std::max(0.0, support.initial_load));
    }

    const double regularization = 1.0e-10;
    const std::size_t n = supports.size();
    const std::size_t mask_limit = static_cast<std::size_t>(1) << n;
    double best_objective = std::numeric_limits<double>::infinity();
    std::vector<double> best_loads(n, 0.0);
    bool found = false;

    for (std::size_t mask = 0; mask < mask_limit; ++mask) {
        const std::size_t active_count = CountBits(mask);
        if (active_count == 0) {
            const double objective = EvaluateObjective(columns, target, initial_loads, best_loads, regularization);
            if (objective < best_objective) {
                best_objective = objective;
                best_loads.assign(n, 0.0);
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
                    Dot6(columns[active_indices[row]], columns[active_indices[col]]) +
                    (row == col ? regularization : 0.0);
            }
            rhs[row] = Dot6(columns[active_indices[row]], target) +
                       regularization * initial_loads[active_indices[row]];
        }

        std::vector<double> active_loads;
        if (!SolveLinearSystem(normal_matrix, rhs, active_loads, active_count)) {
            continue;
        }

        bool nonnegative = true;
        std::vector<double> candidate_loads(n, 0.0);
        for (std::size_t i = 0; i < active_count; ++i) {
            if (active_loads[i] < -1.0e-10) {
                nonnegative = false;
                break;
            }
            candidate_loads[active_indices[i]] = std::max(0.0, active_loads[i]);
        }
        if (!nonnegative) {
            continue;
        }

        const double objective =
            EvaluateObjective(columns, target, initial_loads, candidate_loads, regularization);
        if (objective < best_objective) {
            best_objective = objective;
            best_loads = std::move(candidate_loads);
            found = true;
        }
    }

    out_result.loads = std::move(best_loads);
    out_result.objective = found ? best_objective : 0.0;
    out_result.feasible = found;

    chrono::ChVector3d force_W(0.0, 0.0, 0.0);
    chrono::ChVector3d moment_W(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < supports.size(); ++i) {
        force_W += out_result.loads[i] * supports[i].n_W;
        moment_W += chrono::Vcross(supports[i].x_W - reference.origin_W, out_result.loads[i] * supports[i].n_W);
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
