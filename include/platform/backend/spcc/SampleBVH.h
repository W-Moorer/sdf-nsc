#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include <chrono/core/ChVector3.h>

namespace platform {
namespace backend {
namespace spcc {

struct OBBLocal {
    chrono::ChVector3d center_S;
    std::array<chrono::ChVector3d, 3> axes_S{
        chrono::ChVector3d(1.0, 0.0, 0.0),
        chrono::ChVector3d(0.0, 1.0, 0.0),
        chrono::ChVector3d(0.0, 0.0, 1.0)};
    chrono::ChVector3d half_extents{0.0, 0.0, 0.0};
    double bounding_radius = 0.0;
};

struct SampleBVHNode {
    OBBLocal obb_S;
    chrono::ChVector3d probe_S{0.0, 0.0, 0.0};
    double probe_radius = 0.0;
    std::uint32_t begin = 0;
    std::uint32_t end = 0;
    int left = -1;
    int right = -1;
    bool leaf = false;
};

class SampleBVH {
  public:
    void Build(const std::vector<chrono::ChVector3d>& samples_S, int leaf_size = 32);

    const std::vector<SampleBVHNode>& Nodes() const { return nodes_; }
    const std::vector<std::uint32_t>& Permutation() const { return perm_; }
    bool Empty() const { return nodes_.empty(); }

  private:
    int BuildRecursive(std::uint32_t begin, std::uint32_t end);
    OBBLocal FitOBB(std::uint32_t begin, std::uint32_t end) const;

    const std::vector<chrono::ChVector3d>* samples_S_ = nullptr;
    std::vector<std::uint32_t> perm_;
    std::vector<SampleBVHNode> nodes_;
    int leaf_size_ = 32;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
