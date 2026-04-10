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
    int bvh_leaf_size = 24;
    double bvh_query_margin = 0.0;
    double bvh_velocity_bound_scale = 1.0;
    bool bvh_enable_sdf_node_bound = true;
    double patch_radius = 0.0;
    double normal_cos_min = 0.9;
    double max_patch_diameter = 0.0;
    double max_subpatch_diameter = 0.0;
    double max_plane_error = 0.0;
    double max_second_moment_error = 0.0;
    double max_cone_error = 0.0;
    int cone_direction_count = 16;
    double sentinel_spacing = 0.0;
    double sentinel_margin = 0.0;
    int max_subpatch_depth = 0;
    int min_dense_points_per_subpatch = 0;
    int max_reduced_points_per_patch = 4;
    double warm_start_match_radius = 0.0;
    double temporal_load_regularization = 1.0e-10;
    double temporal_reference_blend = 0.0;
    double temporal_force_transport_blend = 0.5;
    double temporal_slip_velocity_scale = 0.05;
    double temporal_approach_velocity_scale = 0.05;
    double temporal_separation_velocity_scale = 0.02;
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
    cfg.bvh_leaf_size = 24;
    cfg.bvh_query_margin = 2.0e-3;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 8.0e-3;
    cfg.normal_cos_min = 0.93;
    cfg.max_patch_diameter = 1.6e-2;
    cfg.max_subpatch_diameter = 8.0e-3;
    cfg.max_plane_error = 5.0e-4;
    cfg.max_second_moment_error = 8.0e-2;
    cfg.max_cone_error = 8.0e-2;
    cfg.cone_direction_count = 24;
    cfg.sentinel_spacing = 1.5e-3;
    cfg.sentinel_margin = 7.5e-4;
    cfg.max_subpatch_depth = 3;
    cfg.min_dense_points_per_subpatch = 12;
    cfg.max_reduced_points_per_patch = 6;
    cfg.warm_start_match_radius = 4.0e-3;
    cfg.temporal_load_regularization = 1.0e-8;
    cfg.temporal_reference_blend = 0.05;
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
    cfg.bvh_leaf_size = 16;
    cfg.bvh_query_margin = 1.0e-3;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 1.0e-3;
    cfg.normal_cos_min = 0.98;
    cfg.max_patch_diameter = 2.0e-3;
    cfg.max_subpatch_diameter = 1.0e-3;
    cfg.max_plane_error = 2.0e-4;
    cfg.max_second_moment_error = 2.0e-2;
    cfg.max_cone_error = 2.0e-2;
    cfg.cone_direction_count = 12;
    cfg.sentinel_spacing = 7.5e-4;
    cfg.sentinel_margin = 2.5e-4;
    cfg.max_subpatch_depth = 2;
    cfg.min_dense_points_per_subpatch = 8;
    cfg.max_reduced_points_per_patch = 1;
    cfg.warm_start_match_radius = 1.0e-3;
    cfg.temporal_load_regularization = 1.0e-10;
    cfg.temporal_reference_blend = 0.0;
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
    cfg.bvh_leaf_size = 24;
    cfg.bvh_query_margin = 2.5e-4;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 8.0e-4;
    cfg.normal_cos_min = 0.90;
    cfg.max_patch_diameter = 1.2e-3;
    cfg.max_subpatch_diameter = 6.0e-4;
    cfg.max_plane_error = 2.0e-4;
    cfg.max_second_moment_error = 5.0e-2;
    cfg.max_cone_error = 5.0e-2;
    cfg.cone_direction_count = 24;
    cfg.sentinel_spacing = 1.5e-4;
    cfg.sentinel_margin = 5.0e-5;
    cfg.max_subpatch_depth = 4;
    cfg.min_dense_points_per_subpatch = 10;
    cfg.max_reduced_points_per_patch = 6;
    cfg.warm_start_match_radius = 2.0e-4;
    cfg.temporal_load_regularization = 1.0e-8;
    cfg.temporal_reference_blend = 0.05;
    cfg.max_wrench_error = 0.04;
    cfg.max_cop_error = 2.0e-4;
    cfg.max_gap_error = 2.0e-4;
    cfg.predictive_gap = true;
    return cfg;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
