#include "platform/models/BaselineDropCase.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    platform::models::DropCaseConfig config;

    // Basic CLI parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--dt" && i + 1 < argc) {
            config.step_size = std::stod(argv[++i]);
        } else if (arg == "--T" && i + 1 < argc) {
            config.total_time = std::stod(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            config.output_csv_path = argv[++i];
        } else if (arg == "--ball-height" && i + 1 < argc) {
            config.ball_initial_height = std::stod(argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            std::cerr << "Usage: " << argv[0] << " [--dt <step_size>] [--T <total_time>] [--output <csv_path>] [--ball-height <height>]" << std::endl;
            return 1;
        }
    }

    platform::models::BaselineDropCase test_case(config);
    test_case.Run();

    return 0;
}
