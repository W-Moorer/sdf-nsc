#include "platform/models/BaselineCamCase.h"
#include "platform/common/ContactAlgorithm.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<double> ParseSnapshotTimes(const std::string& csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char ch) { return std::isspace(ch); }),
                    token.end());
        if (!token.empty()) {
            out.push_back(std::stod(token));
        }
    }
    return out;
}

}  // namespace

int main(int argc, char* argv[]) {
    platform::models::CamCaseConfig config;
    config.cam_mesh_path = "assets/onset_stress/models/onset_cam.obj";
    config.follower_mesh_path = "assets/onset_stress/models/roller_follower.obj";

    config.cam_init_pos[0] = 0.0;
    config.cam_init_pos[1] = 0.0;
    config.cam_init_pos[2] = 0.0;

    config.follower_init_pos[0] = 0.0;
    config.follower_init_pos[1] = 0.04136029035991934;
    config.follower_init_pos[2] = 0.0;

    config.motor_joint_pos[0] = 0.0;
    config.motor_joint_pos[1] = 0.0;
    config.motor_joint_pos[2] = 0.0;

    config.density = 7800.0;
    config.friction = 0.0;
    config.restitution = 0.0;
    config.gravity_y = 0.0;
    config.motor_speed = -2.0;
    config.step_size = 0.001;
    config.total_time = 0.45;
    config.dynamics_substeps = 2;
    config.env_prefix = "SPCC_ONSET";
    config.contact_algorithm = platform::common::ContactAlgorithm::SdfFirstOrder;
    config.output_csv_path = "data/outputs/onset_stress_sdf1.csv";

    config.sdf_build.voxel_size = 2.0e-4;
    config.sdf_build.half_band_width_voxels = 8.0;
    config.sample_tuning.surface_res = 5.0e-4;
    config.sample_tuning.max_samples = 20000;

    config.contact_regime.delta_on = 4.0e-4;
    config.contact_regime.delta_off = 1.0e-3;
    config.contact_regime.max_active_dense = 48;
    config.contact_regime.patch_radius = 3.0e-3;
    config.contact_regime.normal_cos_min = 0.93;
    config.contact_regime.max_patch_diameter = 6.0e-3;
    config.contact_regime.max_reduced_points_per_patch = 4;
    config.contact_regime.warm_start_match_radius = 2.0e-3;
    config.contact_regime.max_wrench_error = 0.04;
    config.contact_regime.max_cop_error = 5.0e-4;
    config.contact_regime.max_gap_error = 5.0e-4;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--dt" && i + 1 < argc) {
            config.step_size = std::stod(argv[++i]);
        } else if (arg == "--T" && i + 1 < argc) {
            config.total_time = std::stod(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            config.output_csv_path = argv[++i];
        } else if (arg == "--speed" && i + 1 < argc) {
            config.motor_speed = std::stod(argv[++i]);
        } else if (arg == "--snapshot-out" && i + 1 < argc) {
            config.snapshot_output_path = argv[++i];
        } else if (arg == "--snapshot-times" && i + 1 < argc) {
            config.snapshot_times = ParseSnapshotTimes(argv[++i]);
        } else if (arg == "--vtk-dir" && i + 1 < argc) {
            config.vtk_output_dir = argv[++i];
        } else if (arg == "--vtk-stride" && i + 1 < argc) {
            config.vtk_stride = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--contact-algorithm" && i + 1 < argc) {
            const std::string value = argv[++i];
            if (!platform::common::ParseContactAlgorithm(value, config.contact_algorithm)) {
                std::cerr << "Unknown contact algorithm: " << value << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0]
                      << " [--dt <step_size>] [--T <total_time>] [--output <csv_path>] [--speed <rad_s>]"
                      << " [--snapshot-out <json_path>] [--snapshot-times <t1,t2,...>]"
                      << " [--vtk-dir <dir>] [--vtk-stride <N>]"
                      << " [--contact-algorithm <mesh|sdf_1st>]" << std::endl;
            return 1;
        }
    }

    platform::models::BaselineCamCase test_case(config);
    auto t_start = std::chrono::high_resolution_clock::now();
    test_case.Run();
    auto t_end = std::chrono::high_resolution_clock::now();
    double exec_time = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "[PERF] Total Run Time: " << exec_time << " s" << std::endl;

    return 0;
}
