#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>

namespace platform {
namespace backend {
namespace spcc {

struct DenseSurfaceSample {
    std::size_t source_face_index = 0;
    chrono::ChVector3d xi_slave_S;
    chrono::ChVector3d normal_slave_S;
    double area_weight = 0.0;
};

class DenseSurfaceSampler {
  public:
    static std::vector<DenseSurfaceSample> BuildFromMesh(const chrono::ChTriangleMeshConnected& mesh,
                                                         std::size_t max_samples,
                                                         double surface_res);
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
