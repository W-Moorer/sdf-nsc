#include "platform/common/ContactAlgorithm.h"
#include "platform/models/BaselineSimpleGearCase.h"
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
        } else if (arg == "--contact-algorithm" && i + 1 < argc) {
            const std::string value = argv[++i];
            if (!platform::common::ParseContactAlgorithm(value, config.contact_algorithm)) {
                std::cerr << "Unknown contact algorithm: " << value << std::endl;
                return 1;
            }
        } else if (arg == "--sdf-type" && i + 1 < argc) {
            const int legacy_type = std::stoi(argv[++i]);
            config.contact_algorithm = (legacy_type == 0) ? platform::common::ContactAlgorithm::NativeMesh
                                                          : platform::common::ContactAlgorithm::SdfFirstOrder;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0]
                      << " [--dt <step_size>] [--T <total_time>] [--output <csv_path>] [--speed <rad_s>]"
                      << " [--contact-algorithm <mesh|sdf_1st>] [--sdf-type <0|1|2>]"
                      << std::endl;
            return 1;
        }
    }

    platform::models::BaselineSimpleGearCase test_case(config);
    test_case.Run();
    return 0;
}
