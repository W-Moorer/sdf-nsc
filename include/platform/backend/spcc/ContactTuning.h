#pragma once

namespace platform {
namespace backend {
namespace spcc {

enum class ContactRegimeType {
    CompactDiscrete,
    SlidingPatch,
};

enum class PatchGeometryMode {
    DeepestPoint,
    RepresentativeQuery,
};

enum class ContactManifoldKind {
    ImpactPair,
    CompactPoint,
    SlidingPatch,
    ArcSliding,
};

struct SdfBuildTuning {
    double voxel_size = 5.0e-4;
    double half_band_width_voxels = 3.0;
    bool direct_phi_hessian = false;
};

struct SurfaceSampleTuning {
    double surface_res = 5.0e-4;
    int max_samples = 0;
};

struct ContactActivationTuning {
    double delta_on = 0.0;
    double delta_off = 0.0;
    int hold_steps = 0;
    int max_active_keep = 0;

    int full_scan_period = 1;
    double local_scan_radius = 0.0;

    double cluster_angle_deg = 30.0;
    double separating_cutoff = 1.0e-4;
    double cluster_radius = 0.0;
    bool avg_point = false;

    double persistent_match_radius = 0.0;
    double persistent_normal_cos_min = -1.0;
    double persistent_blend_alpha = 1.0;
    int persistent_path_samples = 1;
    double coverage_spacing_radius = 0.0;
    int onset_refine_steps = 0;
    int onset_refine_path_samples = 1;
    double onset_refine_backtrack_scale = 1.0;
    int local_fit_onset_steps = 0;
    int local_fit_min_cluster_size = 3;
    int local_fit_path_samples = 1;
    double local_fit_max_shift_ratio = 0.5;
    double local_fit_blend = 1.0;
    double local_fit_reject_positive_phi = -1.0;
    double onset_gate_current_phi_max = -1.0;
    int single_point_local_fit_path_samples = 0;
    double single_point_local_fit_backtrack_scale = 1.0;
    bool use_sample_bvh = false;
    int sample_bvh_leaf_size = 32;
    double sample_bvh_margin_scale = 1.0;
    bool sample_bvh_use_persistent_seeds_only = false;
    int sample_bvh_warmup_steps = 0;
};

struct CurvatureGateTuning {
    bool enabled = false;
    bool tangential_only = false;
    double normal_alignment_cos_min = -1.0;
    double max_hessian_frobenius = 0.0;
    double small_step_dt_threshold = 0.0;
    double small_step_max_hessian_frobenius = 0.0;
    double max_curvature_term_abs = 0.0;
    double max_curvature_term_ratio = 0.0;
    double gap_floor = 0.0;
    int ramp_steps = 0;
};

struct ManifoldQuadratureTuning {
    ContactManifoldKind manifold_kind = ContactManifoldKind::CompactPoint;
    int target_contacts = 1;
    double span_scale = 1.0;
    double min_half_span = 0.0;
};

struct ContactRegimeConfig {
    ContactRegimeType regime = ContactRegimeType::CompactDiscrete;
    PatchGeometryMode patch_geometry_mode = PatchGeometryMode::DeepestPoint;
    ContactActivationTuning activation;
    CurvatureGateTuning curvature;
    ManifoldQuadratureTuning quadrature;
};

inline SdfBuildTuning MakeCamSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 4.0e-4;
    tuning.half_band_width_voxels = 12.0;
    tuning.direct_phi_hessian = false;
    return tuning;
}

inline SurfaceSampleTuning MakeCamSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 1.0e-3;
    tuning.max_samples = 60000;
    return tuning;
}

inline ContactRegimeConfig MakeCamSlidingPatchDefaults() {
    ContactRegimeConfig cfg;
    cfg.regime = ContactRegimeType::SlidingPatch;
    cfg.patch_geometry_mode = PatchGeometryMode::RepresentativeQuery;
    cfg.quadrature.manifold_kind = ContactManifoldKind::SlidingPatch;
    cfg.quadrature.target_contacts = 1;
    cfg.activation.delta_on = 4.0e-3;
    cfg.activation.delta_off = 5.0e-3;
    cfg.activation.hold_steps = 2;
    cfg.activation.max_active_keep = 8;
    cfg.activation.full_scan_period = 1;
    cfg.activation.local_scan_radius = 1.0e-2;
    cfg.activation.cluster_angle_deg = 25.0;
    cfg.activation.separating_cutoff = 1.0e-3;
    cfg.activation.cluster_radius = 2.5e-3;
    cfg.activation.avg_point = true;
    cfg.activation.persistent_match_radius = 6.0e-3;
    cfg.activation.persistent_normal_cos_min = 0.92;
    cfg.activation.persistent_blend_alpha = 0.60;
    cfg.activation.persistent_path_samples = 1;
    cfg.activation.coverage_spacing_radius = 4.0e-3;
    cfg.activation.onset_refine_steps = 0;
    cfg.activation.onset_refine_path_samples = 5;
    cfg.activation.onset_refine_backtrack_scale = 1.0;
    cfg.activation.local_fit_onset_steps = 1;
    cfg.activation.local_fit_min_cluster_size = 2;
    cfg.activation.local_fit_path_samples = 1;
    cfg.activation.local_fit_max_shift_ratio = 1.0;
    cfg.activation.local_fit_blend = 1.0;
    cfg.activation.local_fit_reject_positive_phi = 2.5e-3;
    cfg.activation.onset_gate_current_phi_max = -1.0;
    cfg.activation.single_point_local_fit_path_samples = 0;
    cfg.activation.single_point_local_fit_backtrack_scale = 1.0;
    cfg.activation.use_sample_bvh = true;
    cfg.curvature.enabled = true;
    cfg.curvature.tangential_only = true;
    cfg.curvature.normal_alignment_cos_min = 0.99;
    cfg.curvature.max_hessian_frobenius = 120.0;
    cfg.curvature.max_curvature_term_abs = 6.0e-4;
    cfg.curvature.max_curvature_term_ratio = 0.20;
    cfg.curvature.gap_floor = 4.0e-3;
    cfg.curvature.ramp_steps = 0;
    return cfg;
}

inline SdfBuildTuning MakeEccentricRollerSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 2.0e-4;
    tuning.half_band_width_voxels = 8.0;
    tuning.direct_phi_hessian = false;
    return tuning;
}

inline SurfaceSampleTuning MakeEccentricRollerSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 5.0e-4;
    tuning.max_samples = 20000;
    return tuning;
}

inline ContactRegimeConfig MakeEccentricRollerDefaults() {
    ContactRegimeConfig cfg = MakeCamSlidingPatchDefaults();
    cfg.activation.max_active_keep = 1;
    cfg.activation.local_scan_radius = 6.0e-3;
    cfg.activation.cluster_radius = 1.5e-3;
    cfg.activation.avg_point = false;
    cfg.activation.use_sample_bvh = true;
    cfg.activation.sample_bvh_use_persistent_seeds_only = true;
    cfg.activation.local_fit_onset_steps = 0;
    cfg.activation.onset_refine_steps = 0;
    return cfg;
}

inline SdfBuildTuning MakeHeadOnSphereSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 1.0e-3;
    tuning.half_band_width_voxels = 8.0;
    tuning.direct_phi_hessian = false;
    return tuning;
}

inline SurfaceSampleTuning MakeHeadOnSphereSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 3.0e-3;
    tuning.max_samples = 12000;
    return tuning;
}

inline ContactRegimeConfig MakeHeadOnSphereDefaults() {
    ContactRegimeConfig cfg;
    cfg.regime = ContactRegimeType::CompactDiscrete;
    cfg.patch_geometry_mode = PatchGeometryMode::DeepestPoint;
    cfg.quadrature.manifold_kind = ContactManifoldKind::ImpactPair;
    cfg.quadrature.target_contacts = 1;
    cfg.activation.delta_on = 2.0e-3;
    cfg.activation.delta_off = 3.0e-3;
    cfg.activation.hold_steps = 0;
    cfg.activation.max_active_keep = 1;
    cfg.activation.full_scan_period = 1;
    cfg.activation.local_scan_radius = 0.0;
    cfg.activation.cluster_angle_deg = 20.0;
    cfg.activation.separating_cutoff = 1.0e-5;
    cfg.activation.cluster_radius = 0.0;
    cfg.activation.avg_point = false;
    cfg.activation.use_sample_bvh = false;
    cfg.curvature.max_hessian_frobenius = 90.0;
    cfg.curvature.max_curvature_term_abs = 3.0e-4;
    cfg.curvature.max_curvature_term_ratio = 0.10;
    cfg.curvature.gap_floor = 2.0e-3;
    return cfg;
}

inline SdfBuildTuning MakeGearSdfBuildDefaults() {
    SdfBuildTuning tuning;
    tuning.voxel_size = 2.0e-5;
    tuning.half_band_width_voxels = 20.0;
    tuning.direct_phi_hessian = false;
    return tuning;
}

inline SurfaceSampleTuning MakeGearSurfaceSampleDefaults() {
    SurfaceSampleTuning tuning;
    tuning.surface_res = 5.0e-5;
    tuning.max_samples = 25000;
    return tuning;
}

inline ContactRegimeConfig MakeGearCompactDefaults() {
    ContactRegimeConfig cfg;
    cfg.regime = ContactRegimeType::CompactDiscrete;
    cfg.patch_geometry_mode = PatchGeometryMode::DeepestPoint;
    cfg.quadrature.manifold_kind = ContactManifoldKind::CompactPoint;
    cfg.quadrature.target_contacts = 1;
    cfg.activation.delta_on = 1.0e-4;
    cfg.activation.delta_off = 5.0e-4;
    cfg.activation.hold_steps = 2;
    cfg.activation.max_active_keep = 100;
    cfg.activation.full_scan_period = 10;
    cfg.activation.local_scan_radius = 8.0e-4;
    cfg.activation.cluster_angle_deg = 30.0;
    cfg.activation.separating_cutoff = 1.0e-4;
    cfg.activation.cluster_radius = 3.0e-4;
    cfg.activation.avg_point = false;
    cfg.activation.persistent_match_radius = 0.0;
    cfg.activation.persistent_normal_cos_min = -1.0;
    cfg.activation.persistent_blend_alpha = 1.0;
    cfg.activation.persistent_path_samples = 1;
    cfg.activation.coverage_spacing_radius = 0.0;
    cfg.activation.onset_refine_steps = 0;
    cfg.activation.onset_refine_path_samples = 1;
    cfg.activation.onset_refine_backtrack_scale = 1.0;
    cfg.activation.local_fit_onset_steps = 0;
    cfg.activation.local_fit_min_cluster_size = 3;
    cfg.activation.local_fit_path_samples = 1;
    cfg.activation.local_fit_max_shift_ratio = 0.5;
    cfg.activation.local_fit_blend = 1.0;
    cfg.activation.local_fit_reject_positive_phi = -1.0;
    cfg.activation.onset_gate_current_phi_max = -1.0;
    cfg.activation.single_point_local_fit_path_samples = 0;
    cfg.activation.single_point_local_fit_backtrack_scale = 1.0;
    return cfg;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
