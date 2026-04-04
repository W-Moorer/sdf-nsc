#include "platform/backend/spcc/SPCCProblem.h"

#include <cmath>

#include <gtest/gtest.h>

namespace {

platform::backend::spcc::SPCCProblemInput MakeBaseInput() {
    using platform::backend::spcc::ActiveContactSample;
    using platform::backend::spcc::SPCCProblemInput;

    SPCCProblemInput input;
    input.h = 1e-3;
    input.rho_n = 1e-3;
    input.rho_t = 1e-3;

    input.master_pred.mass = 2.0;
    input.master_pred.inv_mass = 0.5;
    input.master_pred.I_inv_W = chrono::ChMatrix33<>(1);

    input.slave_pred.mass = 1.0;
    input.slave_pred.inv_mass = 1.0;
    input.slave_pred.I_inv_W = chrono::ChMatrix33<>(1);

    ActiveContactSample c;
    c.sample_id = 0;
    c.n_W = chrono::ChVector3d(1.0, 0.0, 0.0);

    c.P_W = chrono::ChMatrix33<>(1);
    c.P_W(0, 0) = 0.0;
    c.P_W(0, 1) = 0.0;
    c.P_W(0, 2) = 0.0;
    c.P_W(1, 0) = 0.0;
    c.P_W(2, 0) = 0.0;

    c.rA_W = chrono::ChVector3d(0.0, 0.02, 0.0);
    c.rB_W = chrono::ChVector3d(0.0, -0.02, 0.0);

    c.mu = 0.30;
    c.phi = -1.1e-4;
    c.u_n_pred = 0.0;
    c.u_pred_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    c.u_tau_pred_W = chrono::ChVector3d(0.0, 0.02, 0.0);

    input.active_samples.push_back(c);
    return input;
}

TEST(SPCCProblemTest, EvaluateFnFtSingleSample) {
    using platform::backend::spcc::SPCCLambda;
    using platform::backend::spcc::SPCCProblem;
    const auto input = MakeBaseInput();
    const auto& c = input.active_samples[0];

    SPCCProblem problem;
    ASSERT_TRUE(problem.Build(input));

    SPCCLambda lambda;
    lambda.Resize(1);
    lambda.lambda_n[0] = 0.5;
    lambda.lambda_tau_W[0] = chrono::ChVector3d(0.0, 0.1, 0.0);

    std::vector<double> Fn;
    std::vector<chrono::ChVector3d> Ft_W;
    problem.EvaluateFnFt(lambda, Fn, Ft_W);

    ASSERT_EQ(Fn.size(), 1u);
    ASSERT_EQ(Ft_W.size(), 1u);

    EXPECT_NEAR(Ft_W[0].x(), 0.0, 1e-10);
    EXPECT_TRUE(std::isfinite(Fn[0]));
    EXPECT_GT(Fn[0], 0.0);

    const chrono::ChVector3d lambda_tau_tan = c.P_W * lambda.lambda_tau_W[0];
    const double tangency = chrono::Vdot(c.n_W, lambda_tau_tan);
    EXPECT_NEAR(tangency, 0.0, 1e-12);
}

TEST(SPCCProblemTest, ComputeResidualsDetectsNormalComponentInLambdaTau) {
    using platform::backend::spcc::SPCCLambda;
    using platform::backend::spcc::SPCCProblem;

    auto input = MakeBaseInput();
    input.active_samples[0].mu = 1.0;

    SPCCProblem problem;
    ASSERT_TRUE(problem.Build(input));

    SPCCLambda lambda;
    lambda.Resize(1);
    lambda.lambda_n[0] = 0.5;
    // Inject a clear normal component in lambda_tau (x-component).
    lambda.lambda_tau_W[0] = chrono::ChVector3d(0.2, 0.1, 0.0);

    double comp_res_inf = 0.0;
    double cone_res_inf = 0.0;
    double tangency_res_inf = 0.0;
    problem.ComputeResiduals(lambda, comp_res_inf, cone_res_inf, tangency_res_inf);

    EXPECT_GT(tangency_res_inf, 1e-3);
    EXPECT_GE(cone_res_inf, 0.0);
}

TEST(SPCCProblemTest, EvaluateFnFtPureNormalImpulseKeepsTangentialBehavior) {
    using platform::backend::spcc::SPCCLambda;
    using platform::backend::spcc::SPCCProblem;

    const auto input = MakeBaseInput();
    SPCCProblem problem;
    ASSERT_TRUE(problem.Build(input));

    SPCCLambda lambda;
    lambda.Resize(1);
    lambda.lambda_n[0] = 0.5;
    lambda.lambda_tau_W[0] = chrono::ChVector3d(0.0, 0.0, 0.0);

    std::vector<double> Fn;
    std::vector<chrono::ChVector3d> Ft_W;
    problem.EvaluateFnFt(lambda, Fn, Ft_W);

    ASSERT_EQ(Fn.size(), 1u);
    ASSERT_EQ(Ft_W.size(), 1u);
    EXPECT_TRUE(std::isfinite(Fn[0]));

    const chrono::ChVector3d n = problem.Data().contacts[0].n_W;
    const double ft_normal = chrono::Vdot(n, Ft_W[0]);
    EXPECT_NEAR(ft_normal, 0.0, 1e-10);
    EXPECT_NEAR(Ft_W[0].x(), 0.0, 1e-10);
}

}  // namespace
