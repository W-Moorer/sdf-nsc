#include "platform/models/BaselineCamCase.h"
#include "platform/common/ContactAlgorithm.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <chrono>
#include <sstream>
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
    config.contact_regime.curvature.small_step_dt_threshold = 0.005;
    config.contact_regime.curvature.small_step_max_hessian_frobenius = 70.0;

    // Basic CLI parsing
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
                      << " [--contact-algorithm <mesh|sdf_1st|sdf_2nd>]" << std::endl;
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
