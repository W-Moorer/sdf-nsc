#pragma once

#include "platform/common/ContactAlgorithm.h"

#include <string>
#include <vector>

namespace platform {
namespace models {

struct CamCaseConfig {
    // Meshes
    std::string cam_mesh_path = "assets/cam/models/cam_body1.obj";
    std::string follower_mesh_path = "assets/cam/models/cam_body2.obj";

    // Mesh offsets
    double cam_init_pos[3] = {-0.320462025491936, -0.325583626629385, -0.116};
    double follower_init_pos[3] = {-0.1662, 0.2338, -0.0662 + 0.05};

    // Material
    double density = 7800.0;
    double friction = 0.2;
    double restitution = 0.0;

    // Dynamics
    double motor_speed = -3.0; // rad/s (Matches Expression 1 in RecurDyn exactly)
    double gravity_y = -9.81;

    // Simulation settings
    double step_size = 0.001; // finer step size for mesh contact
    double total_time = 3.0;
    platform::common::ContactAlgorithm contact_algorithm = platform::common::ContactAlgorithm::SdfSecondOrder;

    // Output
    std::string output_csv_path = "data/outputs/baseline_cam_nsc.csv";
    std::string snapshot_output_path;
    std::vector<double> snapshot_times;
    std::string vtk_output_dir;
    int vtk_stride = 1;
};

} // namespace models
} // namespace platform
