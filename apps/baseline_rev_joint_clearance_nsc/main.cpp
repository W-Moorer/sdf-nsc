#include "platform/common/ContactAlgorithm.h"
#include "platform/models/BaselineRevoluteClearanceCase.h"

#include <chrono>
#include <iostream>
#include <string>

namespace {

std::string DefaultOutputFor(platform::common::ContactAlgorithm algorithm) {
    switch (algorithm) {
        case platform::common::ContactAlgorithm::NativeMesh:
            return "data/outputs/rev_joint_clearance_mesh.csv";
        case platform::common::ContactAlgorithm::SdfFirstOrder:
            return "data/outputs/rev_joint_clearance_sdf1.csv";
        case platform::common::ContactAlgorithm::SdfSecondOrder:
            return "data/outputs/rev_joint_clearance_sdf2.csv";
        default:
            return "data/outputs/rev_joint_clearance.csv";
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    platform::models::RevoluteClearanceCaseConfig config;
    bool output_overridden = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--dt" && i + 1 < argc) {
            config.step_size = std::stod(argv[++i]);
        } else if (arg == "--T" && i + 1 < argc) {
            config.total_time = std::stod(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            config.output_csv_path = argv[++i];
            output_overridden = true;
        } else if (arg == "--contact-algorithm" && i + 1 < argc) {
            const std::string value = argv[++i];
            if (!platform::common::ParseContactAlgorithm(value, config.contact_algorithm)) {
                std::cerr << "Unknown contact algorithm: " << value << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0]
                      << " [--dt <step_size>] [--T <total_time>] [--output <csv_path>]"
                      << " [--contact-algorithm <mesh|sdf_1st|sdf_2nd>]" << std::endl;
            return 1;
        }
    }

    if (!output_overridden) {
        config.output_csv_path = DefaultOutputFor(config.contact_algorithm);
    }

    platform::models::BaselineRevoluteClearanceCase test_case(config);
    const auto t_start = std::chrono::high_resolution_clock::now();
    test_case.Run();
    const auto t_end = std::chrono::high_resolution_clock::now();
    const double exec_time = std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "[PERF] Total Run Time: " << exec_time << " s" << std::endl;
    return 0;
}
