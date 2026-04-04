#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/SPCCProblem.h"

namespace platform {
namespace backend {
namespace spcc {

// Minimal SDF query interface for ContactActivation.
// Real SDF implementations can inherit from this interface.
class SDFField {
  public:
    virtual ~SDFField() = default;

    // Input : x_M (master local point)
    // Output: phi and grad_M
    virtual bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                               double& phi,
                               chrono::ChVector3d& grad_M) const = 0;

    // Output: phi, grad_M, and hessian_M
    virtual bool QueryPhiGradHessianM(const chrono::ChVector3d& x_M,
                                      double& phi,
                                      chrono::ChVector3d& grad_M,
                                      chrono::ChMatrix33<>& hessian_M) const {
        // default fallback
        hessian_M.setZero();
        return QueryPhiGradM(x_M, phi, grad_M);
    }
};

class ContactActivation {
  public:
    struct Stats {
        std::size_t queried = 0;
        std::size_t accepted_before_cap = 0;
        std::size_t accepted_after_cap = 0;
        std::size_t rejected_invalid = 0;
    };

    void Configure(double delta_on, double delta_off, int hold_steps, std::size_t max_active_keep);

    void Reset(std::size_t sample_count);

    // first-pass assumption: local_samples_S.size() == world_samples_W.size()
    void BuildActiveSet(const RigidBodyStateW& master_pred,
                        const RigidBodyStateW& slave_pred,
                        const SDFField& sdf,
                        const std::vector<chrono::ChVector3d>& local_samples_S,
                        const std::vector<chrono::ChVector3d>& world_samples_W,
                        double mu_default,
                        std::vector<ActiveContactSample>& out_active);

    const Stats& GetStats() const { return stats_; }

  private:
    double delta_on_ = 0.0;
    double delta_off_ = 0.0;
    int hold_steps_ = 0;
    std::size_t max_active_keep_ = 0;

    std::vector<uint8_t> prev_active_;
    std::vector<int> hold_counter_;
    std::vector<int> active_age_;
    Stats stats_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
