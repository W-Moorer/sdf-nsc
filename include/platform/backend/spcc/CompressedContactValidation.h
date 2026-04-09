#pragma once

#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/CompressedContactPipeline.h"

namespace platform {
namespace backend {
namespace spcc {

struct DenseValidationContact {
    std::size_t sample_id = 0;

    chrono::ChVector3d x_W;
    chrono::ChVector3d x_master_M;
    chrono::ChVector3d x_master_surface_W;
    chrono::ChVector3d n_W;
    chrono::ChVector3d v_rel_W;

    double phi = 0.0;
    double phi_eff = 0.0;
    double area_weight = 0.0;
};

struct CompressionWrenchSummary {
    chrono::ChVector3d force_W;
    chrono::ChVector3d moment_at_origin_W;
    chrono::ChVector3d cop_W;
    double total_proxy_load = 0.0;
    double worst_gap = 0.0;
};

struct CompressionValidationReport {
    CompressionStats stats;
    CompressionWrenchSummary dense_wrench;
    CompressionWrenchSummary reduced_wrench;
    std::vector<DenseValidationContact> dense_contacts;
    std::vector<ReducedContactPoint> reduced_contacts;
};

class CompressedContactValidation {
  public:
    static void Validate(const CompressedContactConfig& cfg,
                         const std::vector<DenseSurfaceSample>& slave_surface_samples,
                         const RigidBodyStateW& master_state,
                         const RigidBodyStateW& slave_state,
                         const FirstOrderSDF& sdf,
                         double mu_default,
                         double step_size,
                         CompressionValidationReport& out_report);
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
