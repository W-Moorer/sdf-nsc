#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <chrono/collision/ChCollisionShapeSphere.h>
#include <chrono/core/ChQuaternion.h>
#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChContactMaterialNSC.h>
#include <chrono/physics/ChSystemNSC.h>

#include "platform/backend/spcc/CompressedContactValidation.h"

namespace fs = std::filesystem;

namespace {

using platform::backend::spcc::CompressedContactConfig;
using platform::backend::spcc::CompressedContactPipeline;
using platform::backend::spcc::CompressedContactValidation;
using platform::backend::spcc::CompressionStats;
using platform::backend::spcc::DenseSurfaceSample;
using platform::backend::spcc::DenseValidationContact;
using platform::backend::spcc::FirstOrderSDF;
using platform::backend::spcc::ReducedContactPoint;
using platform::backend::spcc::RigidBodyStateW;

enum class ContactMode {
    DenseReference,
    ReducedCompressed,
};

class PlaneSDF final : public FirstOrderSDF {
  public:
    PlaneSDF(const chrono::ChVector3d& origin_W, const chrono::ChVector3d& normal_W)
        : origin_W_(origin_W), normal_W_(normal_W.GetNormalized()) {}

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        grad_M = normal_W_;
        phi = chrono::Vdot(normal_W_, x_M - origin_W_);
        return true;
    }

  private:
    chrono::ChVector3d origin_W_;
    chrono::ChVector3d normal_W_;
};

struct DynamicsScenario {
    std::string name;
    std::string description;
    CompressedContactConfig cfg;
    chrono::ChVector3d gravity_W{0.0, -9.81, 0.0};
    double step_size = 5.0e-4;
    int steps = 240;
    double friction = 0.25;
    double restitution = 0.0;
    double mass = 1.0;
    chrono::ChMatrix33<> inertia_body = chrono::ChMatrix33<>(1);
    chrono::ChVector3d initial_pos_W;
    chrono::ChQuaternion<> initial_rot = chrono::QUNIT;
    chrono::ChVector3d initial_vel_W;
    chrono::ChVector3d initial_w_W;
    std::vector<DenseSurfaceSample> samples;
    std::unique_ptr<FirstOrderSDF> sdf;
};

struct BodyKinematics {
    chrono::ChVector3d pos_W;
    chrono::ChVector3d vel_W;
    chrono::ChVector3d w_W;
    chrono::ChMatrix33<> R_WL = chrono::ChMatrix33<>(1);
    chrono::ChMatrix33<> I_W = chrono::ChMatrix33<>(1);
    double mass = 0.0;
};

struct StepMetrics {
    double time = 0.0;
    double epsF = 0.0;
    double epsM = 0.0;
    double epsCoP = 0.0;
    double epsGap = 0.0;
    double pos_error = 0.0;
    double vel_error = 0.0;
    double ang_vel_error = 0.0;
    double linear_impulse_error = 0.0;
    double angular_impulse_error = 0.0;
    double energy_dense = 0.0;
    double energy_reduced = 0.0;
    double energy_drift_diff = 0.0;
    double temporal_hausdorff = 0.0;
    double temporal_mean_drift = 0.0;
    double support_churn = 0.0;
    std::size_t patch_count = 0;
    std::size_t subpatch_count = 0;
    std::size_t dense_contacts = 0;
    std::size_t reduced_contacts = 0;
    double max_subpatch_plane_error = 0.0;
    double max_subpatch_second_moment_error = 0.0;
    double max_subpatch_cone_error = 0.0;
    double max_subpatch_gap_error = 0.0;
    double max_subpatch_force_residual = 0.0;
    double max_subpatch_moment_residual = 0.0;
    double max_subpatch_reference_wrench_error = 0.0;
    double max_subpatch_reference_cop_error = 0.0;
    double max_dense_micro_force_residual = 0.0;
    double max_dense_micro_moment_residual = 0.0;
};

struct ScenarioSummary {
    double max_epsF = 0.0;
    double max_epsM = 0.0;
    double max_epsCoP = 0.0;
    double max_epsGap = 0.0;
    double max_pos_error = 0.0;
    double max_vel_error = 0.0;
    double max_ang_vel_error = 0.0;
    double max_linear_impulse_error = 0.0;
    double max_angular_impulse_error = 0.0;
    double max_energy_drift_diff = 0.0;
    double max_temporal_hausdorff = 0.0;
    double max_temporal_mean_drift = 0.0;
    double max_support_churn = 0.0;
    double max_contact_pos_error = 0.0;
    double max_contact_vel_error = 0.0;
    double max_contact_ang_vel_error = 0.0;
    double max_contact_linear_impulse_error = 0.0;
    double max_contact_angular_impulse_error = 0.0;
    double max_contact_energy_drift_diff = 0.0;
    double max_contact_temporal_hausdorff = 0.0;
    double max_contact_temporal_mean_drift = 0.0;
    double max_contact_support_churn = 0.0;
    int contact_steps = 0;
};

struct VariantDefinition {
    std::string name;
    std::string description;
};

struct TemporalCoherenceMetrics {
    double hausdorff = 0.0;
    double mean_drift = 0.0;
    double support_churn = 0.0;
};

chrono::ChVector3d Normalized(const chrono::ChVector3d& v) {
    const double len = v.Length();
    if (!(len > 1.0e-12)) {
        return chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    return v * (1.0 / len);
}

void EmitReducedContactStencil(chrono::ChSystem* sys,
                               const std::shared_ptr<chrono::ChBody>& master,
                               const std::shared_ptr<chrono::ChBody>& slave,
                               const std::shared_ptr<chrono::ChContactMaterial>& material,
                               ReducedContactPoint& contact,
                               std::size_t& emitted_contacts) {
    auto emit_one = [&](const chrono::ChVector3d& vpA_W,
                        const chrono::ChVector3d& vpB_W,
                        const chrono::ChVector3d& normal_W,
                        float* reaction_cache,
                        double gap_offset) {
        chrono::ChCollisionInfo cinfo;
        cinfo.modelA = master->GetCollisionModel().get();
        cinfo.modelB = slave->GetCollisionModel().get();
        cinfo.shapeA = nullptr;
        cinfo.shapeB = nullptr;
        cinfo.vN = normal_W;
        cinfo.vpA = vpA_W;
        cinfo.vpB = vpB_W;
        cinfo.distance = contact.phi_eff + gap_offset;
        cinfo.reaction_cache = reaction_cache;
        sys->GetContactContainer()->AddContact(cinfo, material, material);
        ++emitted_contacts;
    };

    auto emit_slot = [&](int slot_index,
                         const chrono::ChVector3d& fallback_vpA_W,
                         const chrono::ChVector3d& fallback_vpB_W,
                         float* reaction_cache) {
        const chrono::ChVector3d normal_W =
            (contact.slot_n_W[slot_index].Length2() > 0.0) ? contact.slot_n_W[slot_index] : contact.n_W;
        const chrono::ChVector3d vpA_W =
            (contact.slot_x_master_surface_W[slot_index].Length2() > 0.0) ? contact.slot_x_master_surface_W[slot_index]
                                                                           : fallback_vpA_W;
        const chrono::ChVector3d vpB_W =
            (contact.slot_x_W[slot_index].Length2() > 0.0) ? contact.slot_x_W[slot_index] : fallback_vpB_W;
        emit_one(vpA_W, vpB_W, normal_W, reaction_cache, contact.stencil_gap_offsets[slot_index]);
    };

    if (contact.emission_count <= 1 || !(contact.stencil_half_extent > 1.0e-8)) {
        contact.emission_count = 1;
        emit_slot(0, contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        return;
    }

    const chrono::ChVector3d axis_W = Normalized(contact.stencil_axis_W);
    if (axis_W.Length2() <= 0.0) {
        contact.emission_count = 1;
        emit_slot(0, contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        return;
    }

    const chrono::ChVector3d offset_W = contact.stencil_half_extent * axis_W;
    if (contact.emission_count == 2) {
        emit_slot(1, contact.x_master_surface_W - offset_W, contact.x_W - offset_W, contact.reaction_cache_secondary.data());
        emit_slot(2, contact.x_master_surface_W + offset_W, contact.x_W + offset_W, contact.reaction_cache_tertiary.data());
        return;
    }

    if (contact.emission_count >= 5) {
        const chrono::ChVector3d secondary_axis_W = Normalized(contact.stencil_axis_secondary_W);
        if (!(secondary_axis_W.Length2() > 0.0) || !(contact.stencil_half_extent_secondary > 1.0e-8)) {
            contact.emission_count = 3;
            emit_slot(0, contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
            emit_slot(1, contact.x_master_surface_W - offset_W, contact.x_W - offset_W,
                      contact.reaction_cache_secondary.data());
            emit_slot(2, contact.x_master_surface_W + offset_W, contact.x_W + offset_W,
                      contact.reaction_cache_tertiary.data());
            return;
        }

        const chrono::ChVector3d secondary_offset_W = contact.stencil_half_extent_secondary * secondary_axis_W;
        contact.emission_count = 5;
        emit_slot(0, contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        emit_slot(1, contact.x_master_surface_W - offset_W, contact.x_W - offset_W, contact.reaction_cache_secondary.data());
        emit_slot(2, contact.x_master_surface_W + offset_W, contact.x_W + offset_W, contact.reaction_cache_tertiary.data());
        emit_slot(3, contact.x_master_surface_W - secondary_offset_W, contact.x_W - secondary_offset_W,
                  contact.reaction_cache_quaternary.data());
        emit_slot(4, contact.x_master_surface_W + secondary_offset_W, contact.x_W + secondary_offset_W,
                  contact.reaction_cache_quinary.data());
        return;
    }

    contact.emission_count = 3;
    emit_slot(0, contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
    emit_slot(1, contact.x_master_surface_W - offset_W, contact.x_W - offset_W, contact.reaction_cache_secondary.data());
    emit_slot(2, contact.x_master_surface_W + offset_W, contact.x_W + offset_W, contact.reaction_cache_tertiary.data());
}

void BuildBasis(const chrono::ChVector3d& n_W,
                chrono::ChVector3d& t1_W,
                chrono::ChVector3d& t2_W) {
    const chrono::ChVector3d seed =
        (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0) : chrono::ChVector3d(1.0, 0.0, 0.0);
    t1_W = Normalized(chrono::Vcross(seed, n_W));
    t2_W = Normalized(chrono::Vcross(n_W, t1_W));
}

std::vector<DenseSurfaceSample> BuildRectSamples(const chrono::ChVector3d& center_S,
                                                 const chrono::ChVector3d& axis_u_S,
                                                 const chrono::ChVector3d& axis_v_S,
                                                 double size_u,
                                                 double size_v,
                                                 int nu,
                                                 int nv) {
    std::vector<DenseSurfaceSample> samples;
    if (nu <= 0 || nv <= 0) {
        return samples;
    }

    const chrono::ChVector3d u_S = Normalized(axis_u_S);
    const chrono::ChVector3d v_S = Normalized(axis_v_S);
    const chrono::ChVector3d n_S = Normalized(chrono::Vcross(u_S, v_S));
    const double area_weight = (size_u * size_v) / static_cast<double>(nu * nv);

    samples.reserve(static_cast<std::size_t>(nu * nv));
    for (int iu = 0; iu < nu; ++iu) {
        for (int iv = 0; iv < nv; ++iv) {
            const double u = ((static_cast<double>(iu) + 0.5) / static_cast<double>(nu) - 0.5) * size_u;
            const double v = ((static_cast<double>(iv) + 0.5) / static_cast<double>(nv) - 0.5) * size_v;

            DenseSurfaceSample sample;
            sample.xi_slave_S = center_S + u * u_S + v * v_S;
            sample.normal_slave_S = n_S;
            sample.area_weight = area_weight;
            samples.push_back(sample);
        }
    }
    return samples;
}

CompressedContactConfig MakeDynamicConfig() {
    CompressedContactConfig cfg;
    cfg.delta_on = 4.0e-3;
    cfg.delta_off = 5.0e-3;
    cfg.max_active_dense = 0;
    cfg.bvh_leaf_size = 24;
    cfg.bvh_query_margin = 2.0e-3;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 4.0e-2;
    cfg.normal_cos_min = 0.93;
    cfg.max_patch_diameter = 8.0e-2;
    cfg.max_subpatch_diameter = 3.0e-2;
    cfg.max_plane_error = 1.0e-3;
    cfg.max_second_moment_error = 0.26;
    cfg.max_cone_error = 0.25;
    cfg.cone_direction_count = 24;
    cfg.sentinel_spacing = 8.0e-3;
    cfg.sentinel_margin = 2.0e-3;
    cfg.max_subpatch_depth = 2;
    cfg.min_dense_points_per_subpatch = 12;
    cfg.max_reduced_points_per_patch = 6;
    cfg.max_dynamic_reduced_points_per_patch = 8;
    cfg.warm_start_match_radius = 6.0e-3;
    cfg.temporal_load_regularization = 1.0e-6;
    cfg.temporal_reference_blend = 0.03;
    cfg.temporal_force_transport_blend = 0.6;
    cfg.temporal_stencil_blend = 0.28;
    cfg.temporal_slip_velocity_scale = 0.12;
    cfg.temporal_approach_velocity_scale = 0.06;
    cfg.temporal_separation_velocity_scale = 0.04;
    cfg.tangential_heterogeneity_threshold = 0.14;
    cfg.tangential_emission_threshold = 0.10;
    cfg.max_wrench_error = 0.08;
    cfg.max_cop_error = 2.5e-3;
    cfg.max_gap_error = 0.05;
    cfg.predictive_gap = true;
    return cfg;
}

std::vector<VariantDefinition> BuildVariants() {
    return {
        {"full", "full compressed contact pipeline"},
        {"fixed4", "fixed four-point reduction without subpatch refinement"},
        {"single_patch", "single-patch adaptive reduction without subpatch refinement"},
        {"no_dense_micro", "ablation without dense micro-solver"},
        {"no_eC", "ablation without cone/support-function objective"},
        {"no_sentinel", "ablation without sentinel gap monitoring"},
        {"no_impulse_transport", "ablation without temporal impulse transport"},
        {"no_reinj_accept", "ablation without reinjection-aware Q acceptance"},
    };
}

CompressedContactConfig ApplyVariant(const CompressedContactConfig& base_cfg, const std::string& variant_name) {
    CompressedContactConfig cfg = base_cfg;
    if (variant_name == "full") {
        return cfg;
    }
    if (variant_name == "fixed4") {
        cfg.max_subpatch_depth = 0;
        cfg.max_subpatch_diameter = std::max(cfg.max_patch_diameter, 1.0);
        cfg.max_reduced_points_per_patch = 4;
        return cfg;
    }
    if (variant_name == "single_patch") {
        cfg.max_subpatch_depth = 0;
        cfg.max_subpatch_diameter = std::max(cfg.max_patch_diameter, 1.0);
        return cfg;
    }
    if (variant_name == "no_dense_micro") {
        cfg.enable_dense_micro_solver = false;
        return cfg;
    }
    if (variant_name == "no_eC") {
        cfg.enable_cone_objective = false;
        cfg.max_cone_error = 0.0;
        cfg.cone_direction_count = 0;
        return cfg;
    }
    if (variant_name == "no_sentinel") {
        cfg.enable_sentinel_monitor = false;
        cfg.sentinel_spacing = 0.0;
        cfg.sentinel_margin = 0.0;
        cfg.max_gap_error = 0.0;
        return cfg;
    }
    if (variant_name == "no_impulse_transport") {
        cfg.enable_impulse_transport = false;
        cfg.temporal_reference_blend = 0.0;
        cfg.temporal_force_transport_blend = 0.0;
        return cfg;
    }
    if (variant_name == "no_reinj_accept") {
        cfg.enable_reinjection_acceptance = false;
        return cfg;
    }
    return cfg;
}

chrono::ChMatrix33<> MakeBoxInertia(double mass, const chrono::ChVector3d& half_extents) {
    const double sx = 2.0 * half_extents.x();
    const double sy = 2.0 * half_extents.y();
    const double sz = 2.0 * half_extents.z();
    chrono::ChMatrix33<> I(0);
    I(0, 0) = mass * (sy * sy + sz * sz) / 12.0;
    I(1, 1) = mass * (sx * sx + sz * sz) / 12.0;
    I(2, 2) = mass * (sx * sx + sy * sy) / 12.0;
    return I;
}

DynamicsScenario MakeTiltedPlateImpactScenario() {
    DynamicsScenario s;
    s.name = "tilted_plate_impact";
    s.description = "tilted plate impacts plane under gravity";
    s.cfg = MakeDynamicConfig();
    const chrono::ChVector3d half_extents(3.0e-2, 8.0e-3, 2.0e-2);
    s.mass = 0.8;
    s.inertia_body = MakeBoxInertia(s.mass, half_extents);
    s.initial_pos_W = chrono::ChVector3d(0.0, 1.3e-2, 0.0);
    s.initial_rot = chrono::QuatFromAngleZ(0.23) * chrono::QuatFromAngleX(0.12);
    s.initial_vel_W = chrono::ChVector3d(0.08, -0.12, 0.03);
    s.initial_w_W = chrono::ChVector3d(5.0, -1.0, 2.5);
    s.samples = BuildRectSamples(chrono::ChVector3d(0.0, -half_extents.y(), 0.0),
                                 chrono::ChVector3d(1.0, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 2.0 * half_extents.x(), 2.0 * half_extents.z(), 26, 18);
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    return s;
}

DynamicsScenario MakeTiltedPlateSlideScenario() {
    DynamicsScenario s;
    s.name = "tilted_plate_slide";
    s.description = "tilted plate slides and spins on plane";
    s.cfg = MakeDynamicConfig();
    const chrono::ChVector3d half_extents(4.0e-2, 7.0e-3, 1.6e-2);
    s.mass = 0.9;
    s.inertia_body = MakeBoxInertia(s.mass, half_extents);
    s.initial_pos_W = chrono::ChVector3d(-1.0e-2, 1.05e-2, 5.0e-3);
    s.initial_rot = chrono::QuatFromAngleY(0.16) * chrono::QuatFromAngleZ(-0.18);
    s.initial_vel_W = chrono::ChVector3d(0.25, -0.02, -0.08);
    s.initial_w_W = chrono::ChVector3d(-2.0, 1.0, 6.0);
    s.samples = BuildRectSamples(chrono::ChVector3d(0.0, -half_extents.y(), 0.0),
                                 chrono::ChVector3d(1.0, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 2.0 * half_extents.x(), 2.0 * half_extents.z(), 28, 14);
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    return s;
}

RigidBodyStateW MakeRigidBodyStateW(const std::shared_ptr<chrono::ChBody>& body, int body_id) {
    RigidBodyStateW state;
    state.body_id = body_id;
    if (!body) {
        return state;
    }

    state.x_com_W = body->GetPos();
    state.R_WL = body->GetRotMat();
    state.x_ref_W = state.x_com_W;
    state.R_WRef = state.R_WL;
    state.v_com_W = body->GetPosDt();
    state.w_W = body->GetAngVelParent();
    state.mass = body->GetMass();
    state.inv_mass = (state.mass > 0.0) ? (1.0 / state.mass) : 0.0;
    state.I_inv_L = body->GetInvInertia();
    state.I_inv_W = state.R_WL * state.I_inv_L * state.R_WL.transpose();
    return state;
}

BodyKinematics CaptureKinematics(const std::shared_ptr<chrono::ChBody>& body) {
    BodyKinematics out;
    out.pos_W = body->GetPos();
    out.vel_W = body->GetPosDt();
    out.w_W = body->GetAngVelParent();
    out.R_WL = body->GetRotMat();
    out.mass = body->GetMass();
    out.I_W = out.R_WL * body->GetInertia() * out.R_WL.transpose();
    return out;
}

double ComputeEnergy(const BodyKinematics& state, const chrono::ChVector3d& gravity_W) {
    const double translational = 0.5 * state.mass * state.vel_W.Length2();
    const chrono::ChVector3d Iw = state.I_W * state.w_W;
    const double rotational = 0.5 * chrono::Vdot(state.w_W, Iw);
    const double potential = -state.mass * chrono::Vdot(gravity_W, state.pos_W);
    return translational + rotational + potential;
}

chrono::ChVector3d ComputeLinearImpulse(const BodyKinematics& before,
                                        const BodyKinematics& after,
                                        const chrono::ChVector3d& gravity_W,
                                        double dt) {
    return before.mass * (after.vel_W - before.vel_W) - dt * before.mass * gravity_W;
}

chrono::ChVector3d ComputeAngularImpulse(const BodyKinematics& before,
                                         const BodyKinematics& after) {
    const chrono::ChVector3d L_before = before.I_W * before.w_W;
    const chrono::ChVector3d L_after = after.I_W * after.w_W;
    return L_after - L_before;
}

double NearestContactDistance(const ReducedContactPoint& query,
                              const std::vector<ReducedContactPoint>& supports) {
    if (supports.empty()) {
        return 0.0;
    }

    double best = std::numeric_limits<double>::infinity();
    for (const auto& support : supports) {
        best = std::min(best, (query.x_W - support.x_W).Length());
    }
    return std::isfinite(best) ? best : 0.0;
}

double DirectedHausdorffDistance(const std::vector<ReducedContactPoint>& from,
                                 const std::vector<ReducedContactPoint>& to) {
    if (from.empty() || to.empty()) {
        return 0.0;
    }

    double directed = 0.0;
    for (const auto& point : from) {
        directed = std::max(directed, NearestContactDistance(point, to));
    }
    return directed;
}

TemporalCoherenceMetrics ComputeTemporalCoherence(const std::vector<ReducedContactPoint>& previous_supports,
                                                  const std::vector<ReducedContactPoint>& current_supports,
                                                  double match_radius) {
    TemporalCoherenceMetrics metrics;
    if (previous_supports.empty() && current_supports.empty()) {
        return metrics;
    }

    if (previous_supports.empty() || current_supports.empty()) {
        metrics.support_churn = 1.0;
        return metrics;
    }

    metrics.hausdorff = std::max(DirectedHausdorffDistance(previous_supports, current_supports),
                                 DirectedHausdorffDistance(current_supports, previous_supports));

    bool has_persistent_ids = false;
    for (const auto& point : previous_supports) {
        if (point.persistent_id != 0) {
            has_persistent_ids = true;
            break;
        }
    }
    if (!has_persistent_ids) {
        for (const auto& point : current_supports) {
            if (point.persistent_id != 0) {
                has_persistent_ids = true;
                break;
            }
        }
    }

    if (has_persistent_ids) {
        double sum_drift = 0.0;
        std::size_t matched = 0;
        for (const auto& point : current_supports) {
            for (const auto& previous : previous_supports) {
                if (point.persistent_id == previous.persistent_id && point.support_id == previous.support_id) {
                    sum_drift += (point.x_W - previous.x_W).Length();
                    ++matched;
                    break;
                }
            }
        }
        if (matched > 0) {
            metrics.mean_drift = sum_drift / static_cast<double>(matched);
        }
        metrics.support_churn =
            1.0 - static_cast<double>(matched) /
                      static_cast<double>(std::max(previous_supports.size(), current_supports.size()));
        return metrics;
    }

    double sum_drift = 0.0;
    std::size_t matched = 0;
    for (const auto& point : current_supports) {
        const double drift = NearestContactDistance(point, previous_supports);
        sum_drift += drift;
        if (match_radius > 0.0 && drift <= match_radius) {
            ++matched;
        }
    }
    metrics.mean_drift = sum_drift / static_cast<double>(current_supports.size());
    metrics.support_churn =
        1.0 - static_cast<double>(matched) /
                  static_cast<double>(std::max(previous_supports.size(), current_supports.size()));
    return metrics;
}

class ValidationCollisionCallback : public chrono::ChSystem::CustomCollisionCallback {
  public:
    ValidationCollisionCallback(std::shared_ptr<chrono::ChBody> master,
                                std::shared_ptr<chrono::ChBody> slave,
                                const FirstOrderSDF* sdf,
                                const std::vector<DenseSurfaceSample>* samples,
                                const CompressedContactConfig& cfg,
                                std::shared_ptr<chrono::ChContactMaterial> material,
                                ContactMode mode)
        : master_(std::move(master)),
          slave_(std::move(slave)),
          sdf_(sdf),
          samples_(samples),
          cfg_(cfg),
          material_(std::move(material)),
          mode_(mode) {
        pipeline_.Configure(cfg_);
        if (samples_) {
            pipeline_.SetSlaveSurfaceSamples(*samples_);
        }
    }

    std::size_t last_contact_count = 0;
    CompressionStats last_stats;
    std::vector<ReducedContactPoint> last_reduced_contacts;
    std::vector<ReducedContactPoint> previous_reduced_contacts;

    void OnCustomCollision(chrono::ChSystem* sys) override {
        last_contact_count = 0;
        last_stats = CompressionStats{};
        previous_reduced_contacts = last_reduced_contacts;
        last_reduced_contacts.clear();
        if (!master_ || !slave_ || !sdf_ || !samples_) {
            return;
        }

        const auto master_state = MakeRigidBodyStateW(master_, 1);
        const auto slave_state = MakeRigidBodyStateW(slave_, 2);
        const double mu = static_cast<double>(
            std::static_pointer_cast<chrono::ChContactMaterialNSC>(material_)->GetSlidingFriction());

        if (mode_ == ContactMode::DenseReference) {
            std::vector<DenseValidationContact> dense_contacts;
            CompressedContactValidation::BuildDenseContacts(cfg_, *samples_, master_state, slave_state, *sdf_,
                                                            sys->GetStep(), dense_contacts);
            for (const auto& contact : dense_contacts) {
                chrono::ChCollisionInfo cinfo;
                cinfo.modelA = master_->GetCollisionModel().get();
                cinfo.modelB = slave_->GetCollisionModel().get();
                cinfo.shapeA = nullptr;
                cinfo.shapeB = nullptr;
                cinfo.vN = contact.n_W;
                cinfo.vpA = contact.x_master_surface_W;
                cinfo.vpB = contact.x_W;
                cinfo.distance = contact.phi_eff;
                sys->GetContactContainer()->AddContact(cinfo, material_, material_);
                ++last_contact_count;
            }
            last_stats.dense_count = dense_contacts.size();
            last_stats.reduced_count = dense_contacts.size();
            last_stats.patch_count = dense_contacts.empty() ? 0 : 1;
            return;
        }

        std::vector<ReducedContactPoint> reduced_contacts;
        pipeline_.SyncTemporalWarmStart(previous_reduced_contacts);
        pipeline_.BuildReducedContacts(master_state, slave_state, *sdf_, mu, sys->GetStep(), reduced_contacts,
                                       &last_stats);
        last_reduced_contacts = std::move(reduced_contacts);
        for (auto& contact : last_reduced_contacts) {
            EmitReducedContactStencil(sys, master_, slave_, material_, contact, last_contact_count);
        }
    }

  private:
    std::shared_ptr<chrono::ChBody> master_;
    std::shared_ptr<chrono::ChBody> slave_;
    const FirstOrderSDF* sdf_ = nullptr;
    const std::vector<DenseSurfaceSample>* samples_ = nullptr;
    CompressedContactConfig cfg_;
    std::shared_ptr<chrono::ChContactMaterial> material_;
    ContactMode mode_ = ContactMode::ReducedCompressed;
    CompressedContactPipeline pipeline_;
};

struct SimulationInstance {
    std::unique_ptr<chrono::ChSystemNSC> system;
    std::shared_ptr<chrono::ChBody> master;
    std::shared_ptr<chrono::ChBody> slave;
    std::shared_ptr<ValidationCollisionCallback> callback;
};

SimulationInstance MakeSimulationInstance(const DynamicsScenario& scenario, ContactMode mode) {
    SimulationInstance sim;
    sim.system = std::make_unique<chrono::ChSystemNSC>();
    sim.system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    sim.system->SetGravitationalAcceleration(scenario.gravity_W);

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction(static_cast<float>(scenario.friction));
    material->SetRestitution(static_cast<float>(scenario.restitution));

    sim.master = std::make_shared<chrono::ChBody>();
    sim.master->SetFixed(true);
    sim.master->EnableCollision(true);
    sim.master->AddCollisionShape(chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6));
    sim.system->AddBody(sim.master);

    sim.slave = std::make_shared<chrono::ChBody>();
    sim.slave->SetFixed(false);
    sim.slave->SetMass(scenario.mass);
    sim.slave->SetInertia(scenario.inertia_body);
    sim.slave->SetPos(scenario.initial_pos_W);
    sim.slave->SetRot(scenario.initial_rot);
    sim.slave->SetPosDt(scenario.initial_vel_W);
    sim.slave->SetAngVelParent(scenario.initial_w_W);
    sim.slave->SetSleepingAllowed(false);
    sim.slave->EnableCollision(true);
    sim.slave->AddCollisionShape(chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6));
    sim.system->AddBody(sim.slave);

    sim.callback = std::make_shared<ValidationCollisionCallback>(sim.master, sim.slave, scenario.sdf.get(),
                                                                 &scenario.samples, scenario.cfg, material, mode);
    sim.system->RegisterCustomCollisionCallback(sim.callback);
    return sim;
}

void WriteCsv(const std::string& path,
              const std::string& scenario_name,
              const std::string& variant_name,
              const std::vector<StepMetrics>& steps) {
    if (path.empty()) {
        return;
    }

    const fs::path out_path(path);
    if (out_path.has_parent_path()) {
        fs::create_directories(out_path.parent_path());
    }

    const bool append = fs::exists(out_path);
    std::ofstream out(path, append ? std::ios::app : std::ios::trunc);
    if (!append) {
        out << "scenario,variant,time,epsF,epsM,epsCoP,epsGap,pos_error,vel_error,ang_vel_error,linear_impulse_error,"
               "angular_impulse_error,energy_dense,energy_reduced,energy_drift_diff,temporal_hausdorff,"
               "temporal_mean_drift,support_churn,patch_count,subpatch_count,dense_contacts,reduced_contacts,"
               "max_subpatch_plane_error,max_subpatch_second_moment_error,max_subpatch_cone_error,"
               "max_subpatch_gap_error,max_subpatch_force_residual,max_subpatch_moment_residual,"
               "max_subpatch_reference_wrench_error,max_subpatch_reference_cop_error,"
               "max_dense_micro_force_residual,max_dense_micro_moment_residual\n";
    }
    for (const auto& step : steps) {
        out << scenario_name << ','
            << variant_name << ','
            << step.time << ','
            << step.epsF << ','
            << step.epsM << ','
            << step.epsCoP << ','
            << step.epsGap << ','
            << step.pos_error << ','
            << step.vel_error << ','
            << step.ang_vel_error << ','
            << step.linear_impulse_error << ','
            << step.angular_impulse_error << ','
            << step.energy_dense << ','
            << step.energy_reduced << ','
            << step.energy_drift_diff << ','
            << step.temporal_hausdorff << ','
            << step.temporal_mean_drift << ','
            << step.support_churn << ','
            << step.patch_count << ','
            << step.subpatch_count << ','
            << step.dense_contacts << ','
            << step.reduced_contacts << ','
            << step.max_subpatch_plane_error << ','
            << step.max_subpatch_second_moment_error << ','
            << step.max_subpatch_cone_error << ','
            << step.max_subpatch_gap_error << ','
            << step.max_subpatch_force_residual << ','
            << step.max_subpatch_moment_residual << ','
            << step.max_subpatch_reference_wrench_error << ','
            << step.max_subpatch_reference_cop_error << ','
            << step.max_dense_micro_force_residual << ','
            << step.max_dense_micro_moment_residual << '\n';
    }
}

ScenarioSummary RunScenario(const DynamicsScenario& scenario,
                            const std::string& variant_name,
                            const std::string& csv_path) {
    auto dense_sim = MakeSimulationInstance(scenario, ContactMode::DenseReference);
    auto reduced_sim = MakeSimulationInstance(scenario, ContactMode::ReducedCompressed);

    const double dense_initial_energy = ComputeEnergy(CaptureKinematics(dense_sim.slave), scenario.gravity_W);
    const double reduced_initial_energy = ComputeEnergy(CaptureKinematics(reduced_sim.slave), scenario.gravity_W);

    std::vector<StepMetrics> steps;
    steps.reserve(static_cast<std::size_t>(scenario.steps));
    ScenarioSummary summary;
    std::vector<ReducedContactPoint> previous_reduced_supports;

    for (int step_index = 0; step_index < scenario.steps; ++step_index) {
        const bool has_temporal_reference = !previous_reduced_supports.empty();
        const auto dense_before = CaptureKinematics(dense_sim.slave);
        const auto reduced_before = CaptureKinematics(reduced_sim.slave);

        dense_sim.system->DoStepDynamics(scenario.step_size);
        reduced_sim.system->DoStepDynamics(scenario.step_size);

        const auto dense_after = CaptureKinematics(dense_sim.slave);
        const auto reduced_after = CaptureKinematics(reduced_sim.slave);

        StepMetrics step;
        step.time = (step_index + 1) * scenario.step_size;
        step.epsF = reduced_sim.callback->last_stats.epsilon_F;
        step.epsM = reduced_sim.callback->last_stats.epsilon_M;
        step.epsCoP = reduced_sim.callback->last_stats.epsilon_CoP;
        step.epsGap = reduced_sim.callback->last_stats.epsilon_gap;
        step.pos_error = (dense_after.pos_W - reduced_after.pos_W).Length();
        step.vel_error = (dense_after.vel_W - reduced_after.vel_W).Length();
        step.ang_vel_error = (dense_after.w_W - reduced_after.w_W).Length();
        step.linear_impulse_error =
            (ComputeLinearImpulse(dense_before, dense_after, scenario.gravity_W, scenario.step_size) -
             ComputeLinearImpulse(reduced_before, reduced_after, scenario.gravity_W, scenario.step_size))
                .Length();
        step.angular_impulse_error =
            (ComputeAngularImpulse(dense_before, dense_after) - ComputeAngularImpulse(reduced_before, reduced_after))
                .Length();
        step.energy_dense = ComputeEnergy(dense_after, scenario.gravity_W);
        step.energy_reduced = ComputeEnergy(reduced_after, scenario.gravity_W);
        step.energy_drift_diff =
            std::abs((step.energy_dense - dense_initial_energy) - (step.energy_reduced - reduced_initial_energy));
        const auto temporal_metrics = ComputeTemporalCoherence(previous_reduced_supports,
                                                               reduced_sim.callback->last_reduced_contacts,
                                                               scenario.cfg.warm_start_match_radius);
        step.temporal_hausdorff = temporal_metrics.hausdorff;
        step.temporal_mean_drift = temporal_metrics.mean_drift;
        step.support_churn = temporal_metrics.support_churn;
        step.patch_count = reduced_sim.callback->last_stats.patch_count;
        step.subpatch_count = reduced_sim.callback->last_stats.subpatch_count;
        step.dense_contacts = dense_sim.callback->last_contact_count;
        step.reduced_contacts = reduced_sim.callback->last_contact_count;
        step.max_subpatch_plane_error = reduced_sim.callback->last_stats.max_subpatch_plane_error;
        step.max_subpatch_second_moment_error =
            reduced_sim.callback->last_stats.max_subpatch_second_moment_error;
        step.max_subpatch_cone_error = reduced_sim.callback->last_stats.max_subpatch_cone_error;
        step.max_subpatch_gap_error = reduced_sim.callback->last_stats.max_subpatch_gap_error;
        step.max_subpatch_force_residual = reduced_sim.callback->last_stats.max_subpatch_force_residual;
        step.max_subpatch_moment_residual = reduced_sim.callback->last_stats.max_subpatch_moment_residual;
        step.max_subpatch_reference_wrench_error =
            reduced_sim.callback->last_stats.max_subpatch_reference_wrench_error;
        step.max_subpatch_reference_cop_error = reduced_sim.callback->last_stats.max_subpatch_reference_cop_error;
        step.max_dense_micro_force_residual = reduced_sim.callback->last_stats.max_dense_micro_force_residual;
        step.max_dense_micro_moment_residual = reduced_sim.callback->last_stats.max_dense_micro_moment_residual;
        steps.push_back(step);
        previous_reduced_supports = reduced_sim.callback->last_reduced_contacts;

        summary.max_epsF = std::max(summary.max_epsF, step.epsF);
        summary.max_epsM = std::max(summary.max_epsM, step.epsM);
        summary.max_epsCoP = std::max(summary.max_epsCoP, step.epsCoP);
        summary.max_epsGap = std::max(summary.max_epsGap, step.epsGap);
        summary.max_pos_error = std::max(summary.max_pos_error, step.pos_error);
        summary.max_vel_error = std::max(summary.max_vel_error, step.vel_error);
        summary.max_ang_vel_error = std::max(summary.max_ang_vel_error, step.ang_vel_error);
        summary.max_linear_impulse_error =
            std::max(summary.max_linear_impulse_error, step.linear_impulse_error);
        summary.max_angular_impulse_error =
            std::max(summary.max_angular_impulse_error, step.angular_impulse_error);
        summary.max_energy_drift_diff =
            std::max(summary.max_energy_drift_diff, step.energy_drift_diff);
        summary.max_temporal_hausdorff =
            std::max(summary.max_temporal_hausdorff, step.temporal_hausdorff);
        summary.max_temporal_mean_drift =
            std::max(summary.max_temporal_mean_drift, step.temporal_mean_drift);
        if (has_temporal_reference) {
            summary.max_support_churn = std::max(summary.max_support_churn, step.support_churn);
        }

        if (step.dense_contacts > 0 || step.reduced_contacts > 0) {
            ++summary.contact_steps;
            summary.max_contact_pos_error = std::max(summary.max_contact_pos_error, step.pos_error);
            summary.max_contact_vel_error = std::max(summary.max_contact_vel_error, step.vel_error);
            summary.max_contact_ang_vel_error =
                std::max(summary.max_contact_ang_vel_error, step.ang_vel_error);
            summary.max_contact_linear_impulse_error =
                std::max(summary.max_contact_linear_impulse_error, step.linear_impulse_error);
            summary.max_contact_angular_impulse_error =
                std::max(summary.max_contact_angular_impulse_error, step.angular_impulse_error);
            summary.max_contact_energy_drift_diff =
                std::max(summary.max_contact_energy_drift_diff, step.energy_drift_diff);
            summary.max_contact_temporal_hausdorff =
                std::max(summary.max_contact_temporal_hausdorff, step.temporal_hausdorff);
            summary.max_contact_temporal_mean_drift =
                std::max(summary.max_contact_temporal_mean_drift, step.temporal_mean_drift);
            if (has_temporal_reference) {
                summary.max_contact_support_churn =
                    std::max(summary.max_contact_support_churn, step.support_churn);
            }
        }
    }

    WriteCsv(csv_path, scenario.name, variant_name, steps);
    return summary;
}

std::vector<DynamicsScenario> BuildScenarios() {
    std::vector<DynamicsScenario> scenarios;
    scenarios.push_back(MakeTiltedPlateImpactScenario());
    scenarios.push_back(MakeTiltedPlateSlideScenario());
    return scenarios;
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string csv_path = "data/generated/compressed_contact_dynamics.csv";
    std::string scenario_filter = "all";
    std::string variant_filter = "full";

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--csv" && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (arg == "--scenario" && i + 1 < argc) {
            scenario_filter = argv[++i];
        } else if (arg == "--variant" && i + 1 < argc) {
            variant_filter = argv[++i];
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--scenario all|tilted_plate_impact|tilted_plate_slide]"
                      << " [--variant full|all|fixed4|single_patch|no_dense_micro|no_eC|no_sentinel|"
                         "no_impulse_transport|no_reinj_accept] [--csv output.csv]\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    if (!csv_path.empty()) {
        const fs::path out_path(csv_path);
        if (fs::exists(out_path)) {
            fs::remove(out_path);
        }
    }

    const auto variants = BuildVariants();
    bool ran_any = false;
    for (const auto& variant : variants) {
        if (variant_filter != "all" && variant_filter != variant.name) {
            continue;
        }
        auto scenarios = BuildScenarios();
        for (auto& scenario : scenarios) {
            if (scenario_filter != "all" && scenario_filter != scenario.name) {
                continue;
            }
            scenario.cfg = ApplyVariant(scenario.cfg, variant.name);
            ran_any = true;
            const auto summary = RunScenario(scenario, variant.name, csv_path);
            std::cout << scenario.name << " [" << variant.name << "]\n";
            std::cout << "  variant               : " << variant.description << '\n';
            std::cout << "  desc                  : " << scenario.description << '\n';
            std::cout << "  max_epsF              : " << std::setprecision(8) << summary.max_epsF << '\n';
            std::cout << "  max_epsM              : " << summary.max_epsM << '\n';
            std::cout << "  max_epsCoP            : " << summary.max_epsCoP << '\n';
            std::cout << "  max_epsGap            : " << summary.max_epsGap << '\n';
            std::cout << "  max_pos_error         : " << std::setprecision(8) << summary.max_pos_error << '\n';
            std::cout << "  max_vel_error         : " << summary.max_vel_error << '\n';
            std::cout << "  max_ang_vel_error     : " << summary.max_ang_vel_error << '\n';
            std::cout << "  max_linear_imp_error  : " << summary.max_linear_impulse_error << '\n';
            std::cout << "  max_angular_imp_error : " << summary.max_angular_impulse_error << '\n';
            std::cout << "  max_energy_drift_diff : " << summary.max_energy_drift_diff << '\n';
            std::cout << "  max_temporal_H        : " << summary.max_temporal_hausdorff << '\n';
            std::cout << "  max_temporal_mean     : " << summary.max_temporal_mean_drift << '\n';
            std::cout << "  max_support_churn     : " << summary.max_support_churn << '\n';
            std::cout << "  contact_steps         : " << summary.contact_steps << '\n';
            std::cout << "  contact_max_pos_error : " << summary.max_contact_pos_error << '\n';
            std::cout << "  contact_max_vel_error : " << summary.max_contact_vel_error << '\n';
            std::cout << "  contact_max_ang_error : " << summary.max_contact_ang_vel_error << '\n';
            std::cout << "  contact_max_lin_imp   : " << summary.max_contact_linear_impulse_error << '\n';
            std::cout << "  contact_max_ang_imp   : " << summary.max_contact_angular_impulse_error << '\n';
            std::cout << "  contact_max_Ediff     : " << summary.max_contact_energy_drift_diff << '\n';
            std::cout << "  contact_max_H         : " << summary.max_contact_temporal_hausdorff << '\n';
            std::cout << "  contact_max_mean      : " << summary.max_contact_temporal_mean_drift << '\n';
            std::cout << "  contact_max_churn     : " << summary.max_contact_support_churn << '\n';
        }
    }

    if (!ran_any) {
        std::cerr << "No scenario/variant matched filter: scenario=" << scenario_filter
                  << " variant=" << variant_filter << '\n';
        return 1;
    }

    std::cout << "dynamics_csv : " << csv_path << '\n';
    return 0;
}
