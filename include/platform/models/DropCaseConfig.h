#pragma once

#include <string>

namespace platform {
namespace models {

struct DropCaseConfig {
    // Ball parameters
    double ball_radius = 0.5;
    double ball_density = 1000.0;
    double ball_initial_height = 5.0;

    // Ground parameters
    double ground_size_x = 10.0;
    double ground_size_y = 1.0;
    double ground_size_z = 10.0;

    // Physical material
    double friction = 0.5;
    double restitution = 0.0; // Keep it 0.0 to stay aligned with the baseline template
    double gravity = -9.81;

    // Simulation settings
    double step_size = 0.005;
    double total_time = 1.5;

    // Output
    std::string output_csv_path = "data/outputs/baseline_drop_nsc.csv";
};

} // namespace models
} // namespace platform
