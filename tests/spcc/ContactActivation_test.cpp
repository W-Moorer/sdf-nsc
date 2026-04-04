#include "platform/backend/spcc/ContactActivation.h"

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

namespace {

class PlaneSDFField final : public platform::backend::spcc::SDFField {
  public:
    explicit PlaneSDFField(double bias) : bias_(bias) {}

    void SetBias(double bias) { bias_ = bias; }

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        phi = x_M.x() + bias_;
        grad_M = chrono::ChVector3d(1.0, 0.0, 0.0);
        return true;
    }

  private:
    double bias_ = 0.0;
};

platform::backend::spcc::RigidBodyStateW MakeMasterState() {
    platform::backend::spcc::RigidBodyStateW s;
    s.body_id = 1;
    s.x_com_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    s.R_WL = chrono::ChMatrix33<>(1);
    s.v_com_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    s.w_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    return s;
}

platform::backend::spcc::RigidBodyStateW MakeSlaveState() {
    platform::backend::spcc::RigidBodyStateW s;
    s.body_id = 2;
    s.x_com_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    s.R_WL = chrono::ChMatrix33<>(1);
    s.v_com_W = chrono::ChVector3d(2.0, 1.0, 0.0);
    s.w_W = chrono::ChVector3d(0.0, 0.0, 0.0);
    return s;
}

TEST(ContactActivationTest, BuildActiveSetConstructsKinematicsAndProjection) {
    using platform::backend::spcc::ContactActivation;

    ContactActivation activation;
    activation.Configure(0.0, 0.002, 0, 0);

    PlaneSDFField sdf(0.0);

    const auto master = MakeMasterState();
    const auto slave = MakeSlaveState();

    std::vector<chrono::ChVector3d> local_samples_S = {
        chrono::ChVector3d(0.0, 0.0, 0.0),
    };
    std::vector<chrono::ChVector3d> world_samples_W = {
        chrono::ChVector3d(-0.001, 0.0, 0.0),
    };

    std::vector<platform::backend::spcc::ActiveContactSample> out_active;
    activation.BuildActiveSet(master, slave, sdf, local_samples_S, world_samples_W, 0.3, out_active);

    ASSERT_EQ(out_active.size(), 1u);
    const auto& c = out_active[0];

    EXPECT_NEAR(c.n_W.Length(), 1.0, 1e-12);

    const chrono::ChVector3d p_times_n = c.P_W * c.n_W;
    EXPECT_NEAR(p_times_n.x(), 0.0, 1e-10);
    EXPECT_NEAR(p_times_n.y(), 0.0, 1e-10);
    EXPECT_NEAR(p_times_n.z(), 0.0, 1e-10);

    EXPECT_GT(c.u_n_pred, 0.0);

    const double tangency = chrono::Vdot(c.n_W, c.u_tau_pred_W);
    EXPECT_NEAR(tangency, 0.0, 1e-10);
}

TEST(ContactActivationTest, HysteresisKeepsActiveBetweenOnAndOffThresholds) {
    using platform::backend::spcc::ContactActivation;

    ContactActivation activation;
    activation.Configure(0.0, 0.002, 0, 0);

    PlaneSDFField sdf(0.0);

    const auto master = MakeMasterState();
    const auto slave = MakeSlaveState();

    std::vector<chrono::ChVector3d> local_samples_S = {
        chrono::ChVector3d(0.0, 0.0, 0.0),
    };

    std::vector<chrono::ChVector3d> world_step1 = {
        chrono::ChVector3d(-0.001, 0.0, 0.0),  // phi=-0.001 <= delta_on
    };

    std::vector<platform::backend::spcc::ActiveContactSample> out_active;
    activation.BuildActiveSet(master, slave, sdf, local_samples_S, world_step1, 0.3, out_active);
    ASSERT_EQ(out_active.size(), 1u);
    EXPECT_EQ(out_active[0].active_age, 1);

    std::vector<chrono::ChVector3d> world_step2 = {
        chrono::ChVector3d(0.001, 0.0, 0.0),  // phi=0.001 > delta_on && <= delta_off
    };

    activation.BuildActiveSet(master, slave, sdf, local_samples_S, world_step2, 0.3, out_active);
    ASSERT_EQ(out_active.size(), 1u);
    EXPECT_GE(out_active[0].active_age, 2);
}

}  // namespace
