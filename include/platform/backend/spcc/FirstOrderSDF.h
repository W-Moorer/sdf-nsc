#pragma once

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

class FirstOrderSDF {
  public:
    virtual ~FirstOrderSDF() = default;

    virtual bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                               double& phi,
                               chrono::ChVector3d& grad_M) const = 0;

    virtual bool QueryPhiM(const chrono::ChVector3d& x_M, double& phi) const {
        chrono::ChVector3d grad_M(0.0, 0.0, 0.0);
        return QueryPhiGradM(x_M, phi, grad_M);
    }

    virtual bool GetBoundingSphereM(chrono::ChVector3d& center_M, double& radius) const {
        center_M = chrono::ChVector3d(0.0, 0.0, 0.0);
        radius = 0.0;
        return false;
    }
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
