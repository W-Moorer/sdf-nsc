#pragma once

#include "platform/validation/ReferenceCurve.h"

namespace platform {
namespace validation {

struct MetricsResult {
    double pos_y_rmse = 0.0;
    double pos_y_max_err = 0.0;
    double vel_y_rmse = 0.0;
    double vel_y_max_err = 0.0;
    bool valid = false;
};

class ErrorMetrics {
public:
    // Compare simulation results (sim_curve) with reference (ref_curve)
    // Points are matched by evaluating ref_curve at sim_curve's time points.
    static MetricsResult Compute(const ReferenceCurve& sim_curve, const ReferenceCurve& ref_curve);
};

} // namespace validation
} // namespace platform
