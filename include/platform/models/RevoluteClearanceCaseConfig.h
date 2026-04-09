#pragma once

#include "platform/common/ContactAlgorithm.h"
#include "platform/backend/spcc/ContactTuning.h"

#include <string>

namespace platform {
namespace models {

struct RevoluteClearanceCaseConfig {
    std::string body1_mesh_path = "assets/rev_joint_clearance/models/body1_subtract1_centered.obj";
    std::string body3_mesh_path = "assets/rev_joint_clearance/models/body3_cylinder1_centered.obj";

    double body1_init_pos[3] = {0.0, 0.0, 0.0};
    double body3_init_pos[3] = {0.0, 0.0, 0.0};
    double body2_cm_offset[3] = {1.64388630568497e-06, 8.35588256883334e-05, 2.65983598394373};

    double body3_mass = 12984.280977103;
    double body3_inertia_xx[3] = {1314.65844893167, 7971.80750823384, 7971.80750823384};
    double body3_inertia_xy[3] = {0.0, 3.02128692434659e-13, 0.0};

    double body2_mass = 271.0;
    double body2_inertia_xx[3] = {202.0, 359.0, 203.0};
    double body2_inertia_xy[3] = {0.244945783962547, -0.0295154838928402, -1.39749450062735};

    double friction = 0.0;
    double restitution = 0.0;
    double gravity_y = -9.80665;
    double contact_compliance = 1.0e-9;
    double contact_compliance_t = 1.0e-9;
    double contact_damping_f = 1.0e-5;
    double collision_envelope = 2.0e-3;

    double step_size = 0.001;
    double total_time = 3.0;
    int dynamics_substeps = 1;
    std::string env_prefix = "SPCC_REVCLR";
    platform::common::ContactAlgorithm contact_algorithm = platform::common::ContactAlgorithm::SdfSecondOrder;

    platform::backend::spcc::SdfBuildTuning sdf_build = [] {
        auto t = platform::backend::spcc::MakeCamSdfBuildDefaults();
        t.voxel_size = 2.5e-3;
        t.half_band_width_voxels = 10.0;
        return t;
    }();
    platform::backend::spcc::SurfaceSampleTuning sample_tuning = [] {
        platform::backend::spcc::SurfaceSampleTuning t;
        t.surface_res = 5.0e-2;
        t.max_samples = 4000;
        return t;
    }();
    platform::backend::spcc::ContactRegimeConfig contact_regime = [] {
        auto cfg = platform::backend::spcc::MakeCamSlidingPatchDefaults();
        cfg.quadrature.manifold_kind = platform::backend::spcc::ContactManifoldKind::ArcSliding;
        cfg.quadrature.target_contacts = 3;
        cfg.quadrature.span_scale = 0.35;
        cfg.quadrature.min_half_span = 6.0e-3;
        cfg.activation.delta_on = 1.5e-2;
        cfg.activation.delta_off = 2.5e-2;
        cfg.activation.max_active_keep = 8;
        cfg.activation.full_scan_period = 1;
        cfg.activation.local_scan_radius = 5.0e-2;
        cfg.activation.cluster_radius = 1.0e-2;
        cfg.activation.cluster_angle_deg = 20.0;
        cfg.activation.separating_cutoff = 1.0e-3;
        cfg.activation.persistent_match_radius = 2.0e-2;
        cfg.activation.persistent_normal_cos_min = 0.90;
        cfg.activation.coverage_spacing_radius = 1.0e-2;
        cfg.activation.use_sample_bvh = true;
        cfg.activation.sample_bvh_use_persistent_seeds_only = true;
        cfg.curvature.enabled = true;
        cfg.curvature.tangential_only = true;
        cfg.curvature.normal_alignment_cos_min = 0.98;
        cfg.curvature.max_hessian_frobenius = 80.0;
        cfg.curvature.max_curvature_term_abs = 8.0e-4;
        cfg.curvature.max_curvature_term_ratio = 0.06;
        cfg.curvature.gap_floor = 1.5e-3;
        return cfg;
    }();
    bool use_lcp_manifold_quadrature = false;
    int manifold_quadrature_contacts = 3;
    double manifold_quadrature_span_scale = 0.75;
    double manifold_quadrature_min_half_span = 1.0e-2;

    std::string output_csv_path = "data/outputs/rev_joint_clearance_sdf2.csv";
};

}  // namespace models
}  // namespace platform
