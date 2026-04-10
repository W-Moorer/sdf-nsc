#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/DenseSampleBVH.h"
#include "platform/backend/spcc/DenseSurfaceSampler.h"
#include "platform/backend/spcc/FirstOrderSDF.h"

namespace platform {
namespace backend {
namespace spcc {

struct DenseContactPoint {
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

struct DenseContactCloudStats {
    std::size_t total_samples = 0;
    std::size_t candidate_samples = 0;
    std::size_t phi_prefilter_samples = 0;
    std::size_t exact_samples = 0;
    std::size_t active_samples = 0;
    DenseSampleBVHQueryStats bvh;
};

class DenseContactCloudBuilder {
  public:
    static void Build(const CompressedContactConfig& cfg,
                      const std::vector<DenseSurfaceSample>& slave_surface_samples,
                      const DenseSampleBVH& dense_sample_bvh,
                      const RigidBodyStateW& master_state,
                      const RigidBodyStateW& slave_state,
                      const FirstOrderSDF& sdf,
                      double step_size,
                      std::vector<DenseContactPoint>& out_dense_points,
                      DenseContactCloudStats* out_stats = nullptr);
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
