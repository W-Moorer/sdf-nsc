#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChMatrix33.h>
#include <chrono/core/ChVector3.h>

namespace platform {
namespace backend {
namespace spcc {

struct RigidBodyStateW {
    int body_id = -1;

    chrono::ChVector3d x_com_W;
    chrono::ChMatrix33<> R_WL = chrono::ChMatrix33<>(1);

    chrono::ChVector3d x_ref_W;
    chrono::ChMatrix33<> R_WRef = chrono::ChMatrix33<>(1);

    chrono::ChVector3d v_com_W;
    chrono::ChVector3d w_W;

    double mass = 0.0;
    double inv_mass = 0.0;

    chrono::ChMatrix33<> I_inv_L = chrono::ChMatrix33<>(1);
    chrono::ChMatrix33<> I_inv_W = chrono::ChMatrix33<>(1);
};

struct ActiveContactSample {
    std::size_t sample_id = 0;

    chrono::ChVector3d xi_slave_S;
    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_deepest_master_M;

    double phi = 0.0;

    chrono::ChVector3d grad_M;
    chrono::ChMatrix33<> hessian_M = chrono::ChMatrix33<>(0);
    chrono::ChMatrix33<> hessian_W = chrono::ChMatrix33<>(0);

    chrono::ChVector3d n_W;
    chrono::ChMatrix33<> P_W = chrono::ChMatrix33<>(1);

    chrono::ChVector3d rA_W;
    chrono::ChVector3d rB_W;

    double mu = 0.0;

    chrono::ChVector3d u_pred_W;
    double u_n_pred = 0.0;
    chrono::ChVector3d u_tau_pred_W;

    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d v_rel_W;
    double curvature_term = 0.0;
    double phi_eff = 0.0;
    double curvature_gate = 1.0;
    bool curvature_tangential_only = false;
    double curvature_term_abs_max = 0.0;
    double curvature_term_ratio_max = 0.0;
    double curvature_gap_floor = 0.0;
    double queried_normal_alignment = 1.0;
    double hessian_frobenius = 0.0;

    int active_age = 0;
    int cluster_size = 1;
    std::size_t manifold_id = 0;
    bool manifold_matched = false;
    double patch_weight_sum = 1.0;
};

struct SPCCProblemInput {
    double h = 0.0;
    double rho_n = 0.0;
    double rho_t = 0.0;

    RigidBodyStateW master_pred;
    RigidBodyStateW slave_pred;

    std::vector<ActiveContactSample> active_samples;
    std::size_t max_active_used = 0;
};

struct SPCCProblemData {
    double h = 0.0;
    double rho_n = 0.0;
    double rho_t = 0.0;

    RigidBodyStateW master_pred;
    RigidBodyStateW slave_pred;

    std::vector<ActiveContactSample> contacts;
    std::vector<double> b_n;
    std::vector<chrono::ChVector3d> b_tau_W;

    std::size_t m = 0;
    bool valid = false;
};

struct SPCCLambda {
    std::vector<double> lambda_n;
    std::vector<chrono::ChVector3d> lambda_tau_W;

    void Resize(std::size_t m) {
        lambda_n.assign(m, 0.0);
        lambda_tau_W.assign(m, chrono::ChVector3d(0, 0, 0));
    }

    bool IsSizeConsistent(std::size_t m) const {
        return lambda_n.size() == m && lambda_tau_W.size() == m;
    }
};

class SPCCProblem {
  public:
    bool Build(const SPCCProblemInput& input);

    void EvaluateFnFt(const SPCCLambda& lambda,
                      std::vector<double>& Fn,
                      std::vector<chrono::ChVector3d>& Ft_W) const;

    void ComputeResiduals(const SPCCLambda& lambda,
                          double& comp_res_inf,
                          double& cone_res_inf,
                          double& tangency_res_inf) const;

    const SPCCProblemData& Data() const { return data_; }

  private:
    static chrono::ChVector3d ProjectToTangent(const chrono::ChMatrix33<>& P_W,
                                                const chrono::ChVector3d& v_W);

    SPCCProblemData data_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
