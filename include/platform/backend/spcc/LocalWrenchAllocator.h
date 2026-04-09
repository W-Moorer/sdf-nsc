#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/DenseContactCloud.h"

namespace platform {
namespace backend {
namespace spcc {

struct SupportWrenchPoint {
    chrono::ChVector3d x_W;
    chrono::ChVector3d n_W;
    double initial_load = 0.0;
};

struct ReferenceWrench {
    chrono::ChVector3d origin_W;
    chrono::ChVector3d force_W;
    chrono::ChVector3d moment_W;
    double total_load = 0.0;
};

struct WrenchAllocationResult {
    std::vector<double> loads;
    double objective = 0.0;
    double force_residual = 0.0;
    double moment_residual = 0.0;
    bool feasible = false;
};

class LocalWrenchAllocator {
  public:
    static ReferenceWrench BuildDenseReference(const std::vector<DenseContactPoint>& dense_points,
                                              const std::vector<std::size_t>& member_indices,
                                              const chrono::ChVector3d& origin_W);

    static void Allocate(const std::vector<SupportWrenchPoint>& supports,
                         const ReferenceWrench& reference,
                         double temporal_regularization,
                         WrenchAllocationResult& out_result);
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
