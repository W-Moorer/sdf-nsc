#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/DenseSurfaceSampler.h"
#include "platform/backend/spcc/FirstOrderSDF.h"
#include "platform/backend/spcc/GeometryOBB.h"

namespace platform {
namespace backend {
namespace spcc {

struct DenseSampleBVHQueryStats {
    std::size_t nodes_visited = 0;
    std::size_t nodes_pruned_obb = 0;
    std::size_t nodes_pruned_sdf = 0;
    std::size_t leaf_nodes_visited = 0;
    std::size_t leaf_samples_tested = 0;
    std::size_t candidate_samples = 0;
};

class DenseSampleBVH {
  public:
    void Build(const std::vector<DenseSurfaceSample>& samples, int leaf_size);

    bool Empty() const;
    std::size_t SampleCount() const;
    double RootOuterRadius() const;

    void CollectCandidateSampleIndices(const RigidBodyStateW& master_state,
                                       const RigidBodyStateW& slave_state,
                                       const FirstOrderSDF& sdf,
                                       const CompressedContactConfig& cfg,
                                       double step_size,
                                       std::vector<std::size_t>& out_indices,
                                       DenseSampleBVHQueryStats* out_stats = nullptr) const;

  private:
    struct Node {
        GeometryOBB obb_S;
        uint32_t first = 0;
        uint32_t count = 0;
        int left = -1;
        int right = -1;
        bool is_leaf = false;
    };

    int BuildRecursive(const std::vector<DenseSurfaceSample>& samples, uint32_t first, uint32_t count);

    int leaf_size_ = 24;
    double root_outer_radius_ = 0.0;
    std::vector<std::size_t> ordered_indices_;
    std::vector<Node> nodes_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
