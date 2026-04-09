#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/DenseSurfaceSampler.h"
#include "platform/backend/spcc/FirstOrderSDF.h"

namespace platform {
namespace backend {
namespace spcc {

struct ReducedContactPoint {
    std::size_t patch_id = 0;
    std::size_t support_id = 0;
    std::size_t dense_members = 0;

    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d n_W;
    chrono::ChVector3d v_rel_W;

    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
    double support_weight = 0.0;
    double mu = 0.0;
};

struct CompressionStats {
    std::size_t dense_count = 0;
    std::size_t reduced_count = 0;
    std::size_t patch_count = 0;
    double epsilon_F = 0.0;
    double epsilon_M = 0.0;
    double epsilon_CoP = 0.0;
    double epsilon_gap = 0.0;
};

class CompressedContactPipeline {
  public:
    void Configure(const CompressedContactConfig& cfg);
    void SetSlaveSurfaceSamples(std::vector<DenseSurfaceSample> samples);

    void BuildReducedContacts(const RigidBodyStateW& master_state,
                              const RigidBodyStateW& slave_state,
                              const FirstOrderSDF& sdf,
                              double mu_default,
                              double step_size,
                              std::vector<ReducedContactPoint>& out_contacts,
                              CompressionStats* out_stats = nullptr) const;

  private:
    CompressedContactConfig cfg_;
    std::vector<DenseSurfaceSample> slave_surface_samples_;
    mutable std::vector<ReducedContactPoint> previous_contacts_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
