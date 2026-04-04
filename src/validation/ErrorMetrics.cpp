#include "platform/validation/ErrorMetrics.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace platform {
namespace validation {

MetricsResult ErrorMetrics::Compute(const ReferenceCurve& sim_curve, const ReferenceCurve& ref_curve) {
    MetricsResult res;
    
    const auto& sim_data = sim_curve.GetData();
    if (sim_data.empty()) {
        std::cerr << "[ErrorMetrics] Simulation curve is empty." << std::endl;
        return res;
    }

    int count = 0;
    double pos_sq_err = 0.0;
    double vel_sq_err = 0.0;

    for (const auto& sim_pt : sim_data) {
        DataPoint ref_pt;
        if (ref_curve.GetInterpolated(sim_pt.time, ref_pt)) {
            double p_err = std::abs(sim_pt.pos_y - ref_pt.pos_y);
            double v_err = std::abs(sim_pt.vel_y - ref_pt.vel_y);

            res.pos_y_max_err = std::max(res.pos_y_max_err, p_err);
            res.vel_y_max_err = std::max(res.vel_y_max_err, v_err);

            pos_sq_err += p_err * p_err;
            vel_sq_err += v_err * v_err;
            count++;
        }
    }

    if (count > 0) {
        res.pos_y_rmse = std::sqrt(pos_sq_err / count);
        res.vel_y_rmse = std::sqrt(vel_sq_err / count);
        res.valid = true;
    } else {
        std::cerr << "[ErrorMetrics] No overlapping time points found between sim and reference." << std::endl;
    }

    return res;
}

} // namespace validation
} // namespace platform
