#pragma once

#include "platform/common/ContactAlgorithm.h"
#include "platform/backend/spcc/ContactTuning.h"

#include <string>

namespace platform {
namespace models {

struct SimpleGearCaseConfig {
    // Meshes extracted from RecurDyn RMD (in mm); loaded with mesh_scale = 1e-3.
    std::string gear1_mesh_path = "assets/simple_gear/model/GEAR21.obj";
    std::string gear2_mesh_path = "assets/simple_gear/model/GEAR22.obj";
    double mesh_scale = 1e-3;

    // Body reference positions from BaseGSurfacePatchRefMarker (converted to meters).
    double gear1_ref_pos[3] = {0.01746, -0.001366543853054, -0.000478799197013};
    double gear2_ref_pos[3] = {0.01746, -0.001366699711717, -0.003209194558867};

    // Joint centers (Ground.Marker5 and Ground.Marker7 with locked backlash), meters.
    double joint1_pos[3] = {0.0209136192988016, 0.000503589451286037, 0.00136499444260163};
    double joint2_pos[3] = {0.0209136192988086, 0.000503589459392722, -0.00139500556929607};

    // Material and dynamics from RecurDyn model/locked parameters.
    double density = 7800.0;
    double friction = 0.16;
    double restitution = 0.0;
    double gravity_y = 0.0;
    double motor_speed = -1.0;  // Sign aligned with Gear22 Y:Vel_RX reference convention

    // RecurDyn masses/inertia converted to SI (kg, kg*m^2).
    double gear1_mass = 0.000133136344904928;
    double gear1_inertia_xx[3] = {1.01373364551386e-10, 2.28361045813721e-10, 2.28351850267242e-10};
    double gear2_mass = 0.000133136344903739;
    double gear2_inertia_xx[3] = {1.01373395205820e-10, 2.28361098940409e-10, 2.28351827787442e-10};

    // Simulation settings.
    double step_size = 0.001;
    double total_time = 1.0;
    platform::common::ContactAlgorithm contact_algorithm = platform::common::ContactAlgorithm::SdfFirstOrder;

    // Case-specific SPCC tuning.
    platform::backend::spcc::SdfBuildTuning sdf_build = platform::backend::spcc::MakeGearSdfBuildDefaults();
    platform::backend::spcc::SurfaceSampleTuning sample_tuning =
        platform::backend::spcc::MakeGearSurfaceSampleDefaults();
    platform::backend::spcc::CompressedContactConfig contact_regime =
        platform::backend::spcc::MakeGearCompressedDefaults();

    // Output.
    std::string output_csv_path = "data/outputs/baseline_simple_gear_nsc.csv";
};

} // namespace models
} // namespace platform
