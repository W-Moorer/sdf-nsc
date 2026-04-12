#include "platform/models/BaselineSimpleGearCase.h"
#include "platform/common/ContactAlgorithm.h"
#include <algorithm>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    platform::models::SimpleGearCaseConfig config;

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
        } else if (arg == "--sdf-type" && i + 1 < argc) {
            config.sdf_type = std::stoi(argv[++i]);
        } else if (arg == "--contact-algorithm" && i + 1 < argc) {
            platform::common::ContactAlgorithm algorithm = platform::common::ContactAlgorithm::SdfSecondOrder;
            const std::string value = argv[++i];
            if (!platform::common::ParseContactAlgorithm(value, algorithm)) {
                std::cerr << "Unknown contact algorithm: " << value << std::endl;
                return 1;
            }
            if (algorithm == platform::common::ContactAlgorithm::NativeMesh) {
                config.sdf_type = 0;
            } else if (algorithm == platform::common::ContactAlgorithm::SdfFirstOrder) {
                config.sdf_type = 1;
            } else {
                config.sdf_type = 2;
            }
        } else if (arg == "--vtk-dir" && i + 1 < argc) {
            config.vtk_output_dir = argv[++i];
        } else if (arg == "--vtk-stride" && i + 1 < argc) {
            config.vtk_stride = std::max(1, std::stoi(argv[++i]));
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0]
                      << " [--dt <step_size>] [--T <total_time>] [--output <csv_path>] [--speed <rad_s>]"
                      << " [--sdf-type <0|1|2>] [--contact-algorithm <mesh|sdf_1st|sdf_2nd>]"
                      << " [--vtk-dir <dir>] [--vtk-stride <N>]"
                      << std::endl;
            return 1;
        }
    }

    platform::models::BaselineSimpleGearCase test_case(config);
    test_case.Run();
    return 0;
}
