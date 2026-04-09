#pragma once

#include "platform/common/ContactAlgorithm.h"
#include "platform/backend/spcc/ContactTuning.h"

#include <string>

namespace platform {
namespace models {

struct HeadOnSphereCaseConfig {
    // Assets
    std::string sphere_a_mesh_path = "assets/headon_spheres/models/ball_sphere.obj";
    std::string sphere_b_mesh_path = "assets/headon_spheres/models/ball_sphere.obj";

    // Geometry / mass
    double sphere_radius = 0.05;
    double sphere_a_density = 1000.0;
    double sphere_b_density = 1000.0;

    // Initial state
    double sphere_a_init_pos[3] = {-0.15, 0.0, 0.0};
    double sphere_b_init_pos[3] = {0.15, 0.0, 0.0};
    double sphere_a_init_vel[3] = {1.0, 0.0, 0.0};
    double sphere_b_init_vel[3] = {0.0, 0.0, 0.0};

    // Material
    double friction = 0.0;
    double restitution = 1.0;
    double gravity_y = 0.0;

    // Simulation
    double step_size = 5.0e-4;
    double total_time = 0.5;
    int dynamics_substeps = 4;
    std::string env_prefix = "SPCC_HEADON";
    platform::common::ContactAlgorithm contact_algorithm = platform::common::ContactAlgorithm::SdfSecondOrder;

    // SPCC tuning
    platform::backend::spcc::SdfBuildTuning sdf_build =
        platform::backend::spcc::MakeHeadOnSphereSdfBuildDefaults();
    platform::backend::spcc::SurfaceSampleTuning sample_tuning =
        platform::backend::spcc::MakeHeadOnSphereSurfaceSampleDefaults();
    platform::backend::spcc::ContactRegimeConfig contact_regime =
        platform::backend::spcc::MakeHeadOnSphereDefaults();

    // Output
    std::string output_csv_path = "data/outputs/headon_spheres_sdf2.csv";
};

}  // namespace models
}  // namespace platform
