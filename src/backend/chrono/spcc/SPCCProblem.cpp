#include "platform/backend/spcc/SPCCProblem.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace platform {
namespace backend {
namespace spcc {

namespace {

bool IsFiniteScalar(double v) {
    return std::isfinite(v);
}

bool IsFiniteVec(const chrono::ChVector3d& v) {
    return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

chrono::ChMatrix33<> BuildProjectorFromNormal(const chrono::ChVector3d& n_W) {
    chrono::ChMatrix33<> P(1);
    P(0, 0) -= n_W.x() * n_W.x();
    P(0, 1) -= n_W.x() * n_W.y();
    P(0, 2) -= n_W.x() * n_W.z();
    P(1, 0) -= n_W.y() * n_W.x();
    P(1, 1) -= n_W.y() * n_W.y();
    P(1, 2) -= n_W.y() * n_W.z();
    P(2, 0) -= n_W.z() * n_W.x();
    P(2, 1) -= n_W.z() * n_W.y();
    P(2, 2) -= n_W.z() * n_W.z();
    return P;
}

}  // namespace

chrono::ChVector3d SPCCProblem::ProjectToTangent(const chrono::ChMatrix33<>& P_W,
                                                 const chrono::ChVector3d& v_W) {
    return P_W * v_W;
}

bool SPCCProblem::Build(const SPCCProblemInput& input) {
    data_ = SPCCProblemData{};

    if (input.h <= 0.0 || !IsFiniteScalar(input.h) || !IsFiniteScalar(input.rho_n) ||
        !IsFiniteScalar(input.rho_t)) {
        return false;
    }

    data_.h = input.h;
    data_.rho_n = input.rho_n;
    data_.rho_t = input.rho_t;
    data_.master_pred = input.master_pred;
    data_.slave_pred = input.slave_pred;

    std::size_t use_count = input.active_samples.size();
    if (input.max_active_used > 0) {
        use_count = std::min(use_count, input.max_active_used);
    }

    data_.contacts.clear();
    data_.contacts.reserve(use_count);

    for (std::size_t i = 0; i < use_count; ++i) {
        const auto& in_c = input.active_samples[i];

        if (!IsFiniteScalar(in_c.phi) || !IsFiniteScalar(in_c.mu) || in_c.mu < 0.0 ||
            !IsFiniteScalar(in_c.u_n_pred) || !IsFiniteVec(in_c.n_W) || !IsFiniteVec(in_c.u_pred_W) ||
            !IsFiniteVec(in_c.u_tau_pred_W) || !IsFiniteVec(in_c.x_W) || !IsFiniteVec(in_c.x_master_M) ||
            !IsFiniteVec(in_c.rA_W) || !IsFiniteVec(in_c.rB_W)) {
            continue;
        }

        const double n_len = in_c.n_W.Length();
        if (!IsFiniteScalar(n_len) || n_len <= 1e-12) {
            continue;
        }

        ActiveContactSample c = in_c;
        c.n_W = in_c.n_W * (1.0 / n_len);
        c.P_W = BuildProjectorFromNormal(c.n_W);

        data_.contacts.push_back(c);
    }

    data_.m = data_.contacts.size();

    data_.b_n.resize(data_.m, 0.0);
    data_.b_tau_W.resize(data_.m, chrono::ChVector3d(0, 0, 0));

    for (std::size_t i = 0; i < data_.m; ++i) {
        const auto& c = data_.contacts[i];
        data_.b_n[i] = c.phi + data_.h * c.u_n_pred;
        data_.b_tau_W[i] = c.u_tau_pred_W;
    }

    data_.valid = true;
    return true;
}

void SPCCProblem::EvaluateFnFt(const SPCCLambda& lambda,
                               std::vector<double>& Fn,
                               std::vector<chrono::ChVector3d>& Ft_W) const {
    if (!data_.valid) {
        Fn.clear();
        Ft_W.clear();
        return;
    }

    const std::size_t m = data_.m;
    if (!lambda.IsSizeConsistent(m)) {
        throw std::invalid_argument("SPCCProblem::EvaluateFnFt lambda size mismatch");
    }

    Fn.assign(m, 0.0);
    Ft_W.assign(m, chrono::ChVector3d(0, 0, 0));

    chrono::ChVector3d P_A_W(0, 0, 0);
    chrono::ChVector3d L_A_W(0, 0, 0);
    chrono::ChVector3d P_B_W(0, 0, 0);
    chrono::ChVector3d L_B_W(0, 0, 0);

    for (std::size_t i = 0; i < m; ++i) {
        const auto& c = data_.contacts[i];
        const double lambda_n_i = lambda.lambda_n[i];
        const chrono::ChVector3d lambda_tau_tan_W = ProjectToTangent(c.P_W, lambda.lambda_tau_W[i]);

        const chrono::ChVector3d w_i_W = c.n_W * lambda_n_i + lambda_tau_tan_W;

        P_B_W += w_i_W;
        P_A_W -= w_i_W;

        L_B_W += chrono::Vcross(c.rB_W, w_i_W);
        L_A_W -= chrono::Vcross(c.rA_W, w_i_W);
    }

    const chrono::ChVector3d Delta_vA_W = data_.master_pred.inv_mass * P_A_W;
    const chrono::ChVector3d Delta_vB_W = data_.slave_pred.inv_mass * P_B_W;

    const chrono::ChVector3d Delta_wA_W = data_.master_pred.I_inv_W * L_A_W;
    const chrono::ChVector3d Delta_wB_W = data_.slave_pred.I_inv_W * L_B_W;

    for (std::size_t i = 0; i < m; ++i) {
        const auto& c = data_.contacts[i];

        const chrono::ChVector3d delta_vA_at_i_W = Delta_vA_W + chrono::Vcross(Delta_wA_W, c.rA_W);
        const chrono::ChVector3d delta_vB_at_i_W = Delta_vB_W + chrono::Vcross(Delta_wB_W, c.rB_W);
        const chrono::ChVector3d delta_u_i_W = delta_vB_at_i_W - delta_vA_at_i_W;

        const chrono::ChVector3d lambda_tau_tan_W = ProjectToTangent(c.P_W, lambda.lambda_tau_W[i]);

        Fn[i] = data_.b_n[i] + data_.h * chrono::Vdot(c.n_W, delta_u_i_W) + data_.rho_n * lambda.lambda_n[i];

        Ft_W[i] = data_.b_tau_W[i] + c.P_W * delta_u_i_W + data_.rho_t * lambda_tau_tan_W;
    }
}

void SPCCProblem::ComputeResiduals(const SPCCLambda& lambda,
                                   double& comp_res_inf,
                                   double& cone_res_inf,
                                   double& tangency_res_inf) const {
    comp_res_inf = 0.0;
    cone_res_inf = 0.0;
    tangency_res_inf = 0.0;

    if (!data_.valid || data_.m == 0) {
        return;
    }

    std::vector<double> Fn;
    std::vector<chrono::ChVector3d> Ft_W;
    EvaluateFnFt(lambda, Fn, Ft_W);

    for (std::size_t i = 0; i < data_.m; ++i) {
        const auto& c = data_.contacts[i];
        const double lambda_n_i = lambda.lambda_n[i];
        const chrono::ChVector3d lambda_tau_raw_W = lambda.lambda_tau_W[i];

        const double r_nonneg_lambda = std::max(0.0, -lambda_n_i);
        const double r_nonneg_fn = std::max(0.0, -Fn[i]);
        const double r_comp_prod = std::abs(lambda_n_i * Fn[i]);
        const double comp_i = std::max({r_nonneg_lambda, r_nonneg_fn, r_comp_prod});

        const double cone_i = std::max(0.0, lambda_tau_raw_W.Length() - c.mu * std::max(0.0, lambda_n_i));
        const double tangency_i = std::abs(chrono::Vdot(c.n_W, lambda_tau_raw_W));

        comp_res_inf = std::max(comp_res_inf, comp_i);
        cone_res_inf = std::max(cone_res_inf, cone_i);
        tangency_res_inf = std::max(tangency_res_inf, tangency_i);
    }
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
