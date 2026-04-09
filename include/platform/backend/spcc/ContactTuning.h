#pragma once

namespace platform {
namespace backend {
namespace spcc {

struct SdfBuildTuning {
    double voxel_size = 5.0e-4;
    double half_band_width_voxels = 3.0;
};

struct SurfaceSampleTuning {
    double surface_res = 5.0e-4;
    int max_samples = 0;
};

struct CompressedContactConfig {
    double delta_on = 0.0;
    double delta_off = 0.0;
    int max_active_dense = 0;
    double patch_radius = 0.0;
    double normal_cos_min = 0.9;
    double max_patch_diameter = 0.0;
    int max_reduced_points_per_patch = 4;
    double warm_start_match_radius = 0.0;
    double max_wrench_error = 0.05;
    double max_cop_error = 0.001;
    double max_gap_error = 0.001;
    bool predictive_gap = true;
    bool analytic_sphere_toi_contact = false;
};

inline SdfBuildTuning MakeCamSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 4.0e-4;
    tuning.half_band_width_voxels = 12.0;
    return tuning;
}

inline SurfaceSampleTuning MakeCamSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 1.0e-3;
    tuning.max_samples = 60000;
    return tuning;
}

inline CompressedContactConfig MakeCamCompressedDefaults() {
    CompressedContactConfig cfg;
    cfg.delta_on = 4.0e-3;
    cfg.delta_off = 5.0e-3;
    cfg.max_active_dense = 160;
    cfg.patch_radius = 8.0e-3;
    cfg.normal_cos_min = 0.93;
    cfg.max_patch_diameter = 1.6e-2;
    cfg.max_reduced_points_per_patch = 4;
    cfg.warm_start_match_radius = 4.0e-3;
    cfg.max_wrench_error = 0.05;
    cfg.max_cop_error = 1.0e-3;
    cfg.max_gap_error = 1.0e-3;
    cfg.predictive_gap = true;
    return cfg;
}

inline SdfBuildTuning MakeEccentricRollerSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 2.0e-4;
    tuning.half_band_width_voxels = 8.0;
    return tuning;
}

inline SurfaceSampleTuning MakeEccentricRollerSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 5.0e-4;
    tuning.max_samples = 20000;
    return tuning;
}

inline CompressedContactConfig MakeEccentricRollerDefaults() {
    CompressedContactConfig cfg = MakeCamCompressedDefaults();
    cfg.max_active_dense = 64;
    cfg.patch_radius = 5.0e-3;
    cfg.max_patch_diameter = 1.0e-2;
    return cfg;
}

inline SdfBuildTuning MakeHeadOnSphereSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 1.0e-3;
    tuning.half_band_width_voxels = 8.0;
    return tuning;
}

inline SurfaceSampleTuning MakeHeadOnSphereSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 3.0e-3;
    tuning.max_samples = 12000;
    return tuning;
}

inline CompressedContactConfig MakeHeadOnSphereCompressedDefaults() {
    CompressedContactConfig cfg;
    cfg.delta_on = 2.0e-3;
    cfg.delta_off = 3.0e-3;
    cfg.max_active_dense = 16;
    cfg.patch_radius = 1.0e-3;
    cfg.normal_cos_min = 0.98;
    cfg.max_patch_diameter = 2.0e-3;
    cfg.max_reduced_points_per_patch = 1;
    cfg.warm_start_match_radius = 1.0e-3;
    cfg.max_wrench_error = 0.01;
    cfg.max_cop_error = 5.0e-4;
    cfg.max_gap_error = 5.0e-4;
    cfg.predictive_gap = true;
    cfg.analytic_sphere_toi_contact = true;
    return cfg;
}

inline SdfBuildTuning MakeGearSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 2.0e-5;
    tuning.half_band_width_voxels = 20.0;
    return tuning;
}

inline SurfaceSampleTuning MakeGearSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 5.0e-5;
    tuning.max_samples = 25000;
    return tuning;
}

inline CompressedContactConfig MakeGearCompressedDefaults() {
    CompressedContactConfig cfg;
    cfg.delta_on = 1.0e-4;
    cfg.delta_off = 5.0e-4;
    cfg.max_active_dense = 192;
    cfg.patch_radius = 8.0e-4;
    cfg.normal_cos_min = 0.90;
    cfg.max_patch_diameter = 1.2e-3;
    cfg.max_reduced_points_per_patch = 4;
    cfg.warm_start_match_radius = 2.0e-4;
    cfg.max_wrench_error = 0.04;
    cfg.max_cop_error = 2.0e-4;
    cfg.max_gap_error = 2.0e-4;
    cfg.predictive_gap = true;
    return cfg;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
