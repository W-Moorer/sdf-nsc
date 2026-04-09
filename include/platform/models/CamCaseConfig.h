#pragma once

#include "platform/common/ContactAlgorithm.h"
#include "platform/backend/spcc/ContactTuning.h"

#include <string>
#include <vector>

namespace platform {
namespace models {

struct FollowerPreloadConfig {
    bool enabled = false;
    double anchor_pos[3] = {0.0, 0.0, 0.0};
    double rest_length = 0.0;
    double stiffness = 0.0;
    double damping = 0.0;
};

struct CamCaseConfig {
    // Meshes
    std::string cam_mesh_path = "assets/cam/models/cam_body1.obj";
    std::string follower_mesh_path = "assets/cam/models/cam_body2.obj";

    // Mesh offsets
    double cam_init_pos[3] = {-0.320462025491936, -0.325583626629385, -0.116};
    double follower_init_pos[3] = {-0.1662, 0.2338, -0.0662 + 0.05};
    double motor_joint_pos[3] = {0.0, 0.0, 0.05};

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
    int dynamics_substeps = 4;
    std::string env_prefix = "SPCC_CAM";
    platform::common::ContactAlgorithm contact_algorithm = platform::common::ContactAlgorithm::SdfFirstOrder;

    // Case-specific SPCC tuning.
    platform::backend::spcc::SdfBuildTuning sdf_build = platform::backend::spcc::MakeCamSdfBuildDefaults();
    platform::backend::spcc::SurfaceSampleTuning sample_tuning =
        platform::backend::spcc::MakeCamSurfaceSampleDefaults();
    platform::backend::spcc::CompressedContactConfig contact_regime =
        platform::backend::spcc::MakeCamCompressedDefaults();
    FollowerPreloadConfig follower_preload;

    // Output
    std::string output_csv_path = "data/outputs/baseline_cam_nsc.csv";
    std::string snapshot_output_path;
    std::vector<double> snapshot_times;
    std::string vtk_output_dir;
    int vtk_stride = 1;
};

} // namespace models
} // namespace platform
