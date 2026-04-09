#pragma once

#include <cstddef>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/DenseContactCloud.h"

namespace platform {
namespace backend {
namespace spcc {

struct DensePatch {
    std::size_t patch_id = 0;
    std::vector<std::size_t> members;
    chrono::ChVector3d centroid_W;
    chrono::ChVector3d avg_normal_W;
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
};

struct DenseSubpatch {
    std::size_t patch_id = 0;
    std::size_t subpatch_id = 0;
    std::size_t depth = 0;
    std::vector<std::size_t> members;
    chrono::ChVector3d centroid_W;
    chrono::ChVector3d avg_normal_W;
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    double diameter = 0.0;
    double plane_error = 0.0;
    std::vector<chrono::ChVector3d> sentinel_W;
};

class SubpatchRefiner {
  public:
    static void BuildPatches(const std::vector<DenseContactPoint>& dense_points,
                             const CompressedContactConfig& cfg,
                             std::vector<DensePatch>& out_patches);

    static void BuildSubpatches(const std::vector<DenseContactPoint>& dense_points,
                                const std::vector<DensePatch>& patches,
                                const CompressedContactConfig& cfg,
                                std::vector<DenseSubpatch>& out_subpatches);
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
