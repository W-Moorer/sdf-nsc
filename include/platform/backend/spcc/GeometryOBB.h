#pragma once

#include <chrono/core/ChMatrix33.h>
#include <chrono/core/ChVector3.h>

namespace platform {
namespace backend {
namespace spcc {

struct GeometryOBB {
    chrono::ChVector3d center_S;
    chrono::ChMatrix33<> R_SO = chrono::ChMatrix33<>(1);
    chrono::ChVector3d half_extent;
    double outer_radius = 0.0;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
