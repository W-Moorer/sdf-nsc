#include "platform/backend/ChronoRigidSystemNSC.h"
#include "platform/models/CamCaseConfig.h"
#include <chrono/physics/ChContactMaterialNSC.h>
#include <chrono/physics/ChLinkMotorRotationSpeed.h>
#include <chrono/physics/ChLinkLock.h>
#include <chrono/physics/ChLinkTSDA.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <chrono/collision/ChCollisionShapeSphere.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#if defined(SPCC_ENABLE_VDB)
#include "platform/backend/spcc/CompressedContactPipeline.h"
#include "platform/backend/spcc/DenseSurfaceSampler.h"
#include "platform/backend/spcc/VDBSDFField.h"
#endif

namespace platform {
namespace backend {

#if defined(SPCC_ENABLE_VDB)
namespace {

double GetEnvDouble(const char* name, double fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || !std::isfinite(parsed)) {
        return fallback;
    }
    return parsed;
}

int GetEnvInt(const char* name, int fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value) {
        return fallback;
    }
    return static_cast<int>(parsed);
}

std::string MakeScopedEnvName(const std::string& prefix, const char* suffix) {
    if (prefix.empty()) {
        return std::string(suffix);
    }
    return prefix + "_" + suffix;
}

double GetScopedEnvDouble(const std::string& prefix, const char* suffix, double fallback) {
    return GetEnvDouble(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

int GetScopedEnvInt(const std::string& prefix, const char* suffix, int fallback) {
    return GetEnvInt(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

spcc::CompressedContactConfig ResolveCompressedContactConfig(const std::string& env_prefix,
                                                            const spcc::CompressedContactConfig& defaults) {
    spcc::CompressedContactConfig cfg = defaults;
    cfg.delta_on = GetScopedEnvDouble(env_prefix, "DELTA_ON", cfg.delta_on);
    cfg.delta_off = GetScopedEnvDouble(env_prefix, "DELTA_OFF", cfg.delta_off);
    cfg.max_active_dense = GetScopedEnvInt(env_prefix, "MAX_ACTIVE_DENSE", cfg.max_active_dense);
    cfg.bvh_leaf_size = GetScopedEnvInt(env_prefix, "BVH_LEAF_SIZE", cfg.bvh_leaf_size);
    cfg.bvh_query_margin = GetScopedEnvDouble(env_prefix, "BVH_QUERY_MARGIN", cfg.bvh_query_margin);
    cfg.bvh_velocity_bound_scale =
        GetScopedEnvDouble(env_prefix, "BVH_VELOCITY_BOUND_SCALE", cfg.bvh_velocity_bound_scale);
    cfg.bvh_enable_sdf_node_bound =
        (GetScopedEnvDouble(env_prefix, "BVH_ENABLE_SDF_NODE_BOUND", cfg.bvh_enable_sdf_node_bound ? 1.0 : 0.0) >
         0.5);
    cfg.patch_radius = GetScopedEnvDouble(env_prefix, "PATCH_RADIUS", cfg.patch_radius);
    cfg.normal_cos_min = GetScopedEnvDouble(env_prefix, "NORMAL_COS_MIN", cfg.normal_cos_min);
    cfg.max_patch_diameter = GetScopedEnvDouble(env_prefix, "MAX_PATCH_DIAMETER", cfg.max_patch_diameter);
    cfg.max_subpatch_diameter =
        GetScopedEnvDouble(env_prefix, "MAX_SUBPATCH_DIAMETER", cfg.max_subpatch_diameter);
    cfg.max_plane_error = GetScopedEnvDouble(env_prefix, "MAX_PLANE_ERROR", cfg.max_plane_error);
    cfg.max_second_moment_error =
        GetScopedEnvDouble(env_prefix, "MAX_SECOND_MOMENT_ERROR", cfg.max_second_moment_error);
    cfg.max_cone_error = GetScopedEnvDouble(env_prefix, "MAX_CONE_ERROR", cfg.max_cone_error);
    cfg.cone_direction_count = GetScopedEnvInt(env_prefix, "CONE_DIRECTION_COUNT", cfg.cone_direction_count);
    cfg.sentinel_spacing = GetScopedEnvDouble(env_prefix, "SENTINEL_SPACING", cfg.sentinel_spacing);
    cfg.sentinel_margin = GetScopedEnvDouble(env_prefix, "SENTINEL_MARGIN", cfg.sentinel_margin);
    cfg.max_subpatch_depth = GetScopedEnvInt(env_prefix, "MAX_SUBPATCH_DEPTH", cfg.max_subpatch_depth);
    cfg.min_dense_points_per_subpatch =
        GetScopedEnvInt(env_prefix, "MIN_DENSE_POINTS_PER_SUBPATCH", cfg.min_dense_points_per_subpatch);
    cfg.max_reduced_points_per_patch =
        GetScopedEnvInt(env_prefix, "MAX_REDUCED_POINTS_PER_PATCH", cfg.max_reduced_points_per_patch);
    cfg.warm_start_match_radius =
        GetScopedEnvDouble(env_prefix, "WARM_START_MATCH_RADIUS", cfg.warm_start_match_radius);
    cfg.temporal_load_regularization =
        GetScopedEnvDouble(env_prefix, "TEMPORAL_LOAD_REGULARIZATION", cfg.temporal_load_regularization);
    cfg.temporal_reference_blend =
        GetScopedEnvDouble(env_prefix, "TEMPORAL_REFERENCE_BLEND", cfg.temporal_reference_blend);
    cfg.max_wrench_error = GetScopedEnvDouble(env_prefix, "MAX_WRENCH_ERROR", cfg.max_wrench_error);
    cfg.max_cop_error = GetScopedEnvDouble(env_prefix, "MAX_COP_ERROR", cfg.max_cop_error);
    cfg.max_gap_error = GetScopedEnvDouble(env_prefix, "MAX_GAP_ERROR", cfg.max_gap_error);
    cfg.predictive_gap = (GetScopedEnvDouble(env_prefix, "PREDICTIVE_GAP", cfg.predictive_gap ? 1.0 : 0.0) > 0.5);
    cfg.analytic_sphere_toi_contact =
        (GetScopedEnvDouble(env_prefix, "USE_ANALYTIC_TOI_CONTACT", cfg.analytic_sphere_toi_contact ? 1.0 : 0.0) >
         0.5);
    return cfg;
}

chrono::ChVector3d NormalizeOrZero(const chrono::ChVector3d& v) {
    const double len = v.Length();
    if (!(len > 1.0e-12)) {
        return chrono::ChVector3d(0.0, 0.0, 0.0);
    }
    return v * (1.0 / len);
}

void TransferReducedReactionCaches(const std::vector<spcc::ReducedContactPoint>& previous_contacts,
                                   std::vector<spcc::ReducedContactPoint>& current_contacts) {
    for (auto& current : current_contacts) {
        for (const auto& previous : previous_contacts) {
            if (current.persistent_id == previous.persistent_id && current.support_id == previous.support_id) {
                current.reaction_cache_primary = previous.reaction_cache_primary;
                current.reaction_cache_secondary = previous.reaction_cache_secondary;
                current.reaction_cache_tertiary = previous.reaction_cache_tertiary;
                current.reaction_cache_quaternary = previous.reaction_cache_quaternary;
                current.reaction_cache_quinary = previous.reaction_cache_quinary;
                break;
            }
        }
    }
}

void EmitReducedContactStencil(chrono::ChSystem* sys,
                               const std::shared_ptr<chrono::ChBody>& master,
                               const std::shared_ptr<chrono::ChBody>& slave,
                               const std::shared_ptr<chrono::ChContactMaterial>& material,
                               spcc::ReducedContactPoint& contact,
                               bool flip_contact_normal,
                               int& emitted_contacts) {
    auto emit_one = [&](const chrono::ChVector3d& vpA_W, const chrono::ChVector3d& vpB_W, float* reaction_cache) {
        chrono::ChCollisionInfo cinfo;
        cinfo.modelA = master->GetCollisionModel().get();
        cinfo.modelB = slave->GetCollisionModel().get();
        cinfo.shapeA = nullptr;
        cinfo.shapeB = nullptr;
        cinfo.vN = flip_contact_normal ? (-contact.n_W) : contact.n_W;
        cinfo.vpB = vpB_W;
        cinfo.vpA = vpA_W;
        cinfo.distance = contact.phi_eff;
        cinfo.reaction_cache = reaction_cache;
        sys->GetContactContainer()->AddContact(cinfo, material, material);
        ++emitted_contacts;
    };

    if (contact.emission_count <= 1 || !(contact.stencil_half_extent > 1.0e-8)) {
        contact.emission_count = 1;
        emit_one(contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        return;
    }

    const chrono::ChVector3d axis_W = NormalizeOrZero(contact.stencil_axis_W);
    if (axis_W.Length2() <= 0.0) {
        contact.emission_count = 1;
        emit_one(contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        return;
    }

    const chrono::ChVector3d offset_W = contact.stencil_half_extent * axis_W;
    if (contact.emission_count == 2) {
        emit_one(contact.x_master_surface_W - offset_W, contact.x_W - offset_W,
                 contact.reaction_cache_secondary.data());
        emit_one(contact.x_master_surface_W + offset_W, contact.x_W + offset_W,
                 contact.reaction_cache_tertiary.data());
        return;
    }

    if (contact.emission_count >= 5) {
        const chrono::ChVector3d secondary_axis_W = NormalizeOrZero(contact.stencil_axis_secondary_W);
        if (!(secondary_axis_W.Length2() > 0.0) || !(contact.stencil_half_extent_secondary > 1.0e-8)) {
            contact.emission_count = 3;
            emit_one(contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
            emit_one(contact.x_master_surface_W - offset_W, contact.x_W - offset_W,
                     contact.reaction_cache_secondary.data());
            emit_one(contact.x_master_surface_W + offset_W, contact.x_W + offset_W,
                     contact.reaction_cache_tertiary.data());
            return;
        }

        const chrono::ChVector3d secondary_offset_W = contact.stencil_half_extent_secondary * secondary_axis_W;
        contact.emission_count = 5;
        emit_one(contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
        emit_one(contact.x_master_surface_W - offset_W, contact.x_W - offset_W, contact.reaction_cache_secondary.data());
        emit_one(contact.x_master_surface_W + offset_W, contact.x_W + offset_W, contact.reaction_cache_tertiary.data());
        emit_one(contact.x_master_surface_W - secondary_offset_W, contact.x_W - secondary_offset_W,
                 contact.reaction_cache_quaternary.data());
        emit_one(contact.x_master_surface_W + secondary_offset_W, contact.x_W + secondary_offset_W,
                 contact.reaction_cache_quinary.data());
        return;
    }

    contact.emission_count = 3;
    emit_one(contact.x_master_surface_W, contact.x_W, contact.reaction_cache_primary.data());
    emit_one(contact.x_master_surface_W - offset_W, contact.x_W - offset_W, contact.reaction_cache_secondary.data());
    emit_one(contact.x_master_surface_W + offset_W, contact.x_W + offset_W, contact.reaction_cache_tertiary.data());
}

spcc::RigidBodyStateW MakeRigidBodyStateW(const std::shared_ptr<chrono::ChBody>& body,
                                          int body_id) {
    spcc::RigidBodyStateW state;
    state.body_id = body_id;
    if (!body) {
        return state;
    }

    state.x_com_W = body->GetPos();
    state.R_WL = body->GetRotMat();
    state.v_com_W = body->GetPosDt();
    state.w_W = body->GetAngVelParent();
    
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        state.x_ref_W = aux_ref->GetFrameRefToAbs().GetPos();
        state.R_WRef  = aux_ref->GetFrameRefToAbs().GetRotMat();
    } else {
        state.x_ref_W = state.x_com_W;
        state.R_WRef  = state.R_WL;
    }

    state.mass = body->GetMass();
    state.inv_mass = (state.mass > 0.0) ? (1.0 / state.mass) : 0.0;

    state.I_inv_L = body->GetInvInertia();
    state.I_inv_W = state.R_WL * state.I_inv_L * state.R_WL.transpose();
    return state;
}

BodyDebugSnapshot MakeBodyDebugSnapshot(const std::shared_ptr<chrono::ChBody>& body) {
    BodyDebugSnapshot out;
    if (!body) {
        return out;
    }

    out.valid = true;
    out.x_com_W = body->GetPos();
    out.q_WL = body->GetRot();
    out.v_com_W = body->GetPosDt();
    out.w_W = body->GetAngVelParent();

    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        const auto frame_ref = aux_ref->GetFrameRefToAbs();
        out.x_ref_W = frame_ref.GetPos();
        out.q_WRef = frame_ref.GetRot();
    } else {
        out.x_ref_W = out.x_com_W;
        out.q_WRef = out.q_WL;
    }

    return out;
}

std::vector<chrono::ChVector3d> BuildWorldSamplesFromLocal(const std::shared_ptr<chrono::ChBody>& slave_body,
                                                            const std::vector<chrono::ChVector3d>& local_samples_S) {
    std::vector<chrono::ChVector3d> world_samples_W;
    if (!slave_body || local_samples_S.empty()) {
        return world_samples_W;
    }

    world_samples_W.reserve(local_samples_S.size());
    // Crucial fix: ChBodyEasyMesh places the actual mesh at the REF frame, not COG frame.
    // If the body has a REF frame offset, we must apply it.
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(slave_body);
    chrono::ChFrame<> frame_W = aux_ref ? aux_ref->GetFrameRefToAbs() : slave_body->GetFrameCOMToAbs();

    for (const auto& xi_S : local_samples_S) {
        world_samples_W.push_back(frame_W.TransformPointLocalToParent(xi_S));
    }
    return world_samples_W;
}

std::shared_ptr<chrono::ChTriangleMeshConnected> BuildUvSphereMesh(double radius, int slices, int stacks) {
    auto mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
    auto& vertices = mesh->GetCoordsVertices();
    auto& faces = mesh->GetIndicesVertexes();

    slices = std::max(8, slices);
    stacks = std::max(4, stacks);
    vertices.reserve(static_cast<std::size_t>((stacks + 1) * slices));
    faces.reserve(static_cast<std::size_t>(2 * slices * (stacks - 1)));

    for (int stack = 0; stack <= stacks; ++stack) {
        const double v = static_cast<double>(stack) / static_cast<double>(stacks);
        const double phi = chrono::CH_PI * v;
        const double sin_phi = std::sin(phi);
        const double cos_phi = std::cos(phi);
        for (int slice = 0; slice < slices; ++slice) {
            const double u = static_cast<double>(slice) / static_cast<double>(slices);
            const double theta = chrono::CH_2PI * u;
            vertices.emplace_back(radius * sin_phi * std::cos(theta), radius * cos_phi,
                                  radius * sin_phi * std::sin(theta));
        }
    }

    auto vid = [slices](int stack, int slice) { return stack * slices + (slice % slices); };
    for (int stack = 0; stack < stacks; ++stack) {
        for (int slice = 0; slice < slices; ++slice) {
            const int a = vid(stack, slice);
            const int b = vid(stack, slice + 1);
            const int c = vid(stack + 1, slice);
            const int d = vid(stack + 1, slice + 1);
            if (stack != 0) {
                faces.emplace_back(a, c, b);
            }
            if (stack != stacks - 1) {
                faces.emplace_back(b, c, d);
            }
        }
    }

    return mesh;
}

void SetBodyFrameOrPos(const std::shared_ptr<chrono::ChBody>& body, const chrono::ChVector3d& pos_W) {
    if (!body) {
        return;
    }
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        aux_ref->SetFrameRefToAbs(chrono::ChFrame<>(pos_W, chrono::QUNIT));
    } else {
        body->SetPos(pos_W);
    }
}

class SDFCollisionCallback : public chrono::ChSystem::CustomCollisionCallback {
public:
    std::shared_ptr<chrono::ChBody> m_master;
    std::shared_ptr<chrono::ChBody> m_slave;
    spcc::VDBSDFField* m_sdf;
    spcc::CompressedContactPipeline* m_pipeline;
    std::shared_ptr<chrono::ChContactMaterial> m_material;
    std::vector<spcc::ReducedContactPoint>* m_reduced_contacts_out;
    std::vector<spcc::ReducedContactPoint> m_emitted_contacts;
    bool m_flip_contact_normal;
    bool m_use_analytic_sphere_toi_contact;
    double m_master_sphere_radius;
    double m_slave_sphere_radius;

    SDFCollisionCallback(std::shared_ptr<chrono::ChBody> master,
                         std::shared_ptr<chrono::ChBody> slave,
                         spcc::VDBSDFField* sdf,
                         spcc::CompressedContactPipeline* pipeline,
                         std::shared_ptr<chrono::ChContactMaterial> mat,
                         std::vector<spcc::ReducedContactPoint>* reduced_contacts_out,
                         bool flip_contact_normal = false,
                         bool use_analytic_sphere_toi_contact = false,
                         double master_sphere_radius = 0.0,
                         double slave_sphere_radius = 0.0)
        : m_master(master), m_slave(slave), m_sdf(sdf), m_pipeline(pipeline), m_material(mat),
          m_reduced_contacts_out(reduced_contacts_out),
          m_flip_contact_normal(flip_contact_normal),
          m_use_analytic_sphere_toi_contact(use_analytic_sphere_toi_contact),
          m_master_sphere_radius(master_sphere_radius),
          m_slave_sphere_radius(slave_sphere_radius) {}

    void AddAnalyticSphereTOIContact(chrono::ChSystem* sys,
                                     const spcc::RigidBodyStateW& master_state,
                                     const spcc::RigidBodyStateW& slave_state) {
        if (m_reduced_contacts_out) {
            m_reduced_contacts_out->clear();
        }

        const double radius_sum = m_master_sphere_radius + m_slave_sphere_radius;
        if (!(radius_sum > 0.0) || !std::isfinite(radius_sum)) {
            return;
        }

        chrono::ChVector3d delta_W = slave_state.x_com_W - master_state.x_com_W;
        double center_dist = delta_W.Length();
        chrono::ChVector3d n_W;
        if (center_dist > 1.0e-12 && std::isfinite(center_dist)) {
            n_W = delta_W * (1.0 / center_dist);
        } else {
            n_W = slave_state.v_com_W - master_state.v_com_W;
            const double n_len = n_W.Length();
            if (!(n_len > 1.0e-12) || !std::isfinite(n_len)) {
                n_W = chrono::ChVector3d(1.0, 0.0, 0.0);
            } else {
                n_W *= (1.0 / n_len);
            }
            center_dist = 0.0;
        }

        const chrono::ChVector3d v_rel = slave_state.v_com_W - master_state.v_com_W;
        const double vn = chrono::Vdot(n_W, v_rel);
        const double gap = center_dist - radius_sum;
        const double dt = sys->GetStep();

        const bool overlapping = gap <= 0.0;
        const bool crosses_this_step = (gap > 0.0) && (vn < -1.0e-12) && ((gap + dt * vn) <= 0.0);
        if ((!overlapping && !crosses_this_step) || !(vn < -1.0e-12)) {
            return;
        }

        double tau = 0.0;
        if (crosses_this_step) {
            tau = std::clamp(gap / (-vn), 0.0, dt);
        }

        const chrono::ChVector3d cA_toi_W = master_state.x_com_W + master_state.v_com_W * tau;
        const chrono::ChVector3d cB_toi_W = slave_state.x_com_W + slave_state.v_com_W * tau;
        const double restitution = static_cast<double>(
            std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetRestitution());
        const double denom = master_state.inv_mass + slave_state.inv_mass;
        if (!(denom > 0.0) || !std::isfinite(denom)) {
            return;
        }

        const double j = -(1.0 + restitution) * vn / denom;
        const chrono::ChVector3d vA_after_W = master_state.v_com_W - j * master_state.inv_mass * n_W;
        const chrono::ChVector3d vB_after_W = slave_state.v_com_W + j * slave_state.inv_mass * n_W;
        const double remain = std::max(0.0, dt - tau);
        const chrono::ChVector3d xA_after_W = cA_toi_W + vA_after_W * remain;
        const chrono::ChVector3d xB_after_W = cB_toi_W + vB_after_W * remain;

        SetBodyFrameOrPos(m_master, xA_after_W);
        SetBodyFrameOrPos(m_slave, xB_after_W);
        m_master->SetPosDt(vA_after_W);
        m_slave->SetPosDt(vB_after_W);
        m_master->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_slave->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_master->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
        m_slave->SetAngVelParent(chrono::ChVector3d(0, 0, 0));

        const chrono::ChVector3d vpA_W = cA_toi_W + m_master_sphere_radius * n_W;
        const chrono::ChVector3d vpB_W = cB_toi_W - m_slave_sphere_radius * n_W;
        const chrono::ChVector3d contact_normal_W = m_flip_contact_normal ? (-n_W) : n_W;

        if (m_reduced_contacts_out) {
            spcc::ReducedContactPoint sample;
            sample.patch_id = 0;
            sample.support_id = 0;
            sample.dense_members = 1;
            sample.x_W = 0.5 * (vpA_W + vpB_W);
            sample.x_master_surface_W = vpA_W;
            sample.x_master_M = master_state.R_WRef.transpose() * (vpA_W - master_state.x_ref_W);
            sample.phi = 0.0;
            sample.phi_eff = 0.0;
            sample.n_W = contact_normal_W;
            sample.v_rel_W = vB_after_W - vA_after_W;
            sample.area_weight = 1.0;
            sample.support_weight = 1.0;
            sample.mu = static_cast<double>(
                std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetSlidingFriction());
            m_reduced_contacts_out->push_back(sample);
        }
    }

    void OnCustomCollision(chrono::ChSystem* sys) override {
        if (!m_master || !m_slave || !m_sdf || !m_pipeline) {
            if (m_reduced_contacts_out) {
                m_reduced_contacts_out->clear();
            }
            return;
        }

        // 1. Build kinematic states
        auto master_state = MakeRigidBodyStateW(m_master, 1);
        auto slave_state = MakeRigidBodyStateW(m_slave, 2);

        if (m_use_analytic_sphere_toi_contact) {
            AddAnalyticSphereTOIContact(sys, master_state, slave_state);
            return;
        }

        std::vector<spcc::ReducedContactPoint> reduced_contacts;
        spcc::CompressionStats stats;
        const double friction = static_cast<double>(
            std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetSlidingFriction());
        m_pipeline->SyncTemporalWarmStart(m_emitted_contacts);
        m_pipeline->BuildReducedContacts(master_state, slave_state, *m_sdf, friction, sys->GetStep(),
                                         reduced_contacts, &stats);
        TransferReducedReactionCaches(m_emitted_contacts, reduced_contacts);
        m_emitted_contacts = std::move(reduced_contacts);

        int penetration_count = 0;
        for (auto& contact : m_emitted_contacts) {
            EmitReducedContactStencil(sys, m_master, m_slave, m_material, contact, m_flip_contact_normal,
                                      penetration_count);
        }

        if (m_reduced_contacts_out) {
            *m_reduced_contacts_out = m_emitted_contacts;
        }
        
        static int step_id = 0;
        if (step_id % 10 == 0) {
             std::cout << "[SDF] Step " << step_id << " reduced contacts: " << penetration_count;
             if (GetEnvDouble("SPCC_DEBUG_CONTACT_STATS", 0.0) > 0.5) {
                 std::cout << " dense=" << stats.dense_count
                           << " total=" << stats.total_samples
                           << " candidates=" << stats.candidate_count
                           << " reduced=" << stats.reduced_count
                           << " patches=" << stats.patch_count
                           << " bvhVisited=" << stats.bvh_nodes_visited
                           << " bvhPrunedObb=" << stats.bvh_nodes_pruned_obb
                           << " bvhPrunedSdf=" << stats.bvh_nodes_pruned_sdf
                           << " epsF=" << stats.epsilon_F
                           << " epsM=" << stats.epsilon_M
                           << " epsCoP=" << stats.epsilon_CoP
                           << " epsGap=" << stats.epsilon_gap;
             }
             std::cout << std::endl;
        }
        step_id++;
    }
};

}  // namespace
#endif

ChronoRigidSystemNSC::ChronoRigidSystemNSC() {
    m_system = std::make_unique<chrono::ChSystemNSC>();
}

ChronoRigidSystemNSC::~ChronoRigidSystemNSC() = default;

void ChronoRigidSystemNSC::Initialize() {
    // Default call with baseline hardcoded params if Interface is strictly used
    Initialize(0.5, 1000.0, 5.0, 10.0, 1.0, 10.0, 0.5, 0.0, -9.81);
}

void ChronoRigidSystemNSC::Initialize(
    double ball_radius, double ball_density, double ball_height,
    double ground_x, double ground_y, double ground_z,
    double friction, double restitution, double gravity_y
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_reduced_contacts.clear();
    m_contact_surface_samples.clear();
#endif

    // 0. Initialize Collision System (Crucial since Chrono 9+)
    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);

    // 1. Setup Gravity
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    // 2. Setup Base Contact Material
    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution); 

    // 3. Fixed Ground Body
    auto ground = std::make_shared<chrono::ChBodyEasyBox>(ground_x, ground_y, ground_z, 1000.0, true, true, material);
    // Align ground top surface to Y=0
    ground->SetPos(chrono::ChVector3d(0, -ground_y / 2.0, 0)); 
    ground->SetFixed(true);
    m_system->AddBody(ground);

    // 4. Dynamic Sphere
    m_sphere = std::make_shared<chrono::ChBodyEasySphere>(ball_radius, ball_density, true, true, material);
    m_sphere->SetPos(chrono::ChVector3d(0, ball_height, 0));
    m_sphere->SetFixed(false);
    m_system->AddBody(m_sphere);
}

void ChronoRigidSystemNSC::StepDynamics(double step_size) {
    double headon_xA0 = 0.0;
    double headon_xB0 = 0.0;
    double headon_vA0 = 0.0;
    double headon_vB0 = 0.0;
    bool use_headon_exact_step = m_headon_manual_elastic_resolution && m_headon_a && m_headon_b;
    if (use_headon_exact_step) {
        headon_xA0 = m_headon_a->GetPos().x();
        headon_xB0 = m_headon_b->GetPos().x();
        headon_vA0 = m_headon_a->GetPosDt().x();
        headon_vB0 = m_headon_b->GetPosDt().x();
    }

    const int substeps = std::max(1, m_dynamics_substeps);
    const double sub_dt = step_size / static_cast<double>(substeps);
    for (int i = 0; i < substeps; ++i) {
        m_system->DoStepDynamics(sub_dt);
    }

    if (use_headon_exact_step) {
        const double mA = std::max(0.0, m_headon_a->GetMass());
        const double mB = std::max(0.0, m_headon_b->GetMass());
        const double radius_sum = 2.0 * m_headon_sphere_radius;
        const double gap0 = (headon_xB0 - headon_xA0) - radius_sum;
        const double closing_speed = headon_vA0 - headon_vB0;

        double xA1 = headon_xA0 + headon_vA0 * step_size;
        double xB1 = headon_xB0 + headon_vB0 * step_size;
        double vA1 = headon_vA0;
        double vB1 = headon_vB0;

        if (mA > 0.0 && mB > 0.0 && closing_speed > 1.0e-12) {
            const double gap1 = gap0 - closing_speed * step_size;
            if (gap0 <= 0.0 || gap1 <= 0.0) {
                double tau = 0.0;
                if (gap0 > 0.0) {
                    tau = std::clamp(gap0 / closing_speed, 0.0, step_size);
                }
                const double xA_toi = headon_xA0 + headon_vA0 * tau;
                const double xB_toi = headon_xB0 + headon_vB0 * tau;
                const double e = std::clamp(m_headon_restitution, 0.0, 1.0);
                vA1 = ((mA - e * mB) * headon_vA0 + (1.0 + e) * mB * headon_vB0) / (mA + mB);
                vB1 = ((mB - e * mA) * headon_vB0 + (1.0 + e) * mA * headon_vA0) / (mA + mB);
                const double remain = std::max(0.0, step_size - tau);
                xA1 = xA_toi + vA1 * remain;
                xB1 = xB_toi + vB1 * remain;
            }
        }

        auto posA = m_headon_a->GetPos();
        auto posB = m_headon_b->GetPos();
        posA.x() = xA1;
        posB.x() = xB1;
        m_headon_a->SetPos(posA);
        m_headon_b->SetPos(posB);

        auto velA = m_headon_a->GetPosDt();
        auto velB = m_headon_b->GetPosDt();
        velA.x() = vA1;
        velB.x() = vB1;
        velA.y() = 0.0;
        velA.z() = 0.0;
        velB.y() = 0.0;
        velB.z() = 0.0;
        m_headon_a->SetPosDt(velA);
        m_headon_b->SetPosDt(velB);
        m_headon_a->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_headon_b->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_headon_a->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
        m_headon_b->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
    }
}

double ChronoRigidSystemNSC::GetTime() const {
    return m_system->GetChTime();
}

unsigned int ChronoRigidSystemNSC::GetNumContacts() const {
    return m_system->GetNumContacts();
}

double ChronoRigidSystemNSC::GetDynamicSpherePosY() const {
    return m_sphere->GetPos().y();
}

double ChronoRigidSystemNSC::GetDynamicSphereVelY() const {
    return m_sphere->GetPosDt().y();
}

double ChronoRigidSystemNSC::GetHeadOnSphereAPosX() const {
    return m_headon_a ? m_headon_a->GetPos().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereAVelX() const {
    return m_headon_a ? m_headon_a->GetPosDt().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereBPosX() const {
    return m_headon_b ? m_headon_b->GetPos().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereBVelX() const {
    return m_headon_b ? m_headon_b->GetPosDt().x() : 0.0;
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody2Pos() const {
    return m_clearance_body2 ? m_clearance_body2->GetPos() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody2Vel() const {
    return m_clearance_body2 ? m_clearance_body2->GetPosDt() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3Pos() const {
    return m_clearance_body3 ? m_clearance_body3->GetPos() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3Vel() const {
    return m_clearance_body3 ? m_clearance_body3->GetPosDt() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3AngVel() const {
    return m_clearance_body3 ? m_clearance_body3->GetAngVelParent() : chrono::ChVector3d(0, 0, 0);
}

void ChronoRigidSystemNSC::InitializeCamCase(
    const std::string& cam_obj, const std::string& follower_obj,
    const double* cam_pos, const double* follower_pos, const double* motor_joint_pos,
    double density, double friction, double restitution,
    double gravity_y, double motor_speed,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::CompressedContactConfig& contact_regime,
    const platform::models::FollowerPreloadConfig& follower_preload
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_reduced_contacts.clear();
    m_contact_surface_samples.clear();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution);

    // Ground anchor for joints
    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
#endif

#if defined(SPCC_ENABLE_VDB)
    std::shared_ptr<chrono::ChTriangleMeshConnected> cam_mesh;
    std::shared_ptr<chrono::ChTriangleMeshConnected> follower_mesh;
    if (use_vdb) {
        cam_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        cam_mesh->LoadWavefrontMesh(cam_obj, true, false);
        follower_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        follower_mesh->LoadWavefrontMesh(follower_obj, true, false);
    }
#endif

    // Load Cam Body
    std::shared_ptr<chrono::ChBodyEasyMesh> cam_body;
    if (use_vdb) {
        cam_body = std::make_shared<chrono::ChBodyEasyMesh>(
            cam_obj, density, true, true, false, material, 0.000001);
        auto dummy_sphere_c = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        cam_body->AddCollisionShape(dummy_sphere_c);
        cam_body->EnableCollision(true);
    } else {
        cam_body = std::make_shared<chrono::ChBodyEasyMesh>(
            cam_obj, density, true, true, true, material, 0.000001);
    }

    // Move frame to align with provided global starting pos
    cam_body->SetFrameRefToAbs(chrono::ChFrame<>(chrono::ChVector3d(cam_pos[0], cam_pos[1], cam_pos[2]), chrono::QUNIT));
    m_system->AddBody(cam_body);
    m_cam_body = cam_body;

    // Motor for Cam to Ground (Revolute around Z axis approx)
    auto motor = std::make_shared<chrono::ChLinkMotorRotationSpeed>();      
    // The revolute motor frame
    motor->Initialize(cam_body, ground,
                      chrono::ChFrame<>(chrono::ChVector3d(motor_joint_pos[0], motor_joint_pos[1], motor_joint_pos[2]),
                                        chrono::QuatFromAngleX(0)));
    auto mfun = std::make_shared<chrono::ChFunctionConst>(motor_speed);     
    motor->SetSpeedFunction(mfun);
    m_system->AddLink(motor);

    // Load Follower Body
    if (use_vdb) {
        m_follower = std::make_shared<chrono::ChBodyEasyMesh>(
            follower_obj, density, true, true, false, material, 0.000001);
        auto dummy_sphere_f = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        m_follower->AddCollisionShape(dummy_sphere_f);
        m_follower->EnableCollision(true);
    } else {
        m_follower = std::make_shared<chrono::ChBodyEasyMesh>(
            follower_obj, density, true, true, true, material, 0.000001);       
    }

    auto follower_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_follower);
    if(follower_aux) {
        follower_aux->SetFrameRefToAbs(chrono::ChFrame<>(chrono::ChVector3d(follower_pos[0], follower_pos[1], follower_pos[2]), chrono::QUNIT));
    }
    m_system->AddBody(m_follower);

    // Prismatic joint: Follower to Ground (along Y axis)
    auto prismatic = std::make_shared<chrono::ChLinkLockPrismatic>();       
    // Initialize passing the center of the joint at follower's current position to allow sliding along Y
    prismatic->Initialize(m_follower, ground, chrono::ChFrame<>(m_follower->GetPos(), chrono::QuatFromAngleX(chrono::CH_PI_2)));
    m_system->AddLink(prismatic);

    if (follower_preload.enabled) {
        auto preload = std::make_shared<chrono::ChLinkTSDA>();
        preload->Initialize(m_follower, ground, false, m_follower->GetPos(),
                            chrono::ChVector3d(follower_preload.anchor_pos[0], follower_preload.anchor_pos[1],
                                               follower_preload.anchor_pos[2]));
        preload->SetRestLength(follower_preload.rest_length);
        preload->SetSpringCoefficient(follower_preload.stiffness);
        preload->SetDampingCoefficient(follower_preload.damping);
        preload->SetActuatorForce(0.0);
        m_system->AddLink(preload);
    }

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));
        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);

        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*cam_mesh, sdf_options);

        const double cam_surface_res = GetScopedEnvDouble(env_prefix, "SURFACE_RES", sample_tuning.surface_res);
        const int cam_max_samples = GetScopedEnvInt(env_prefix, "MAX_SAMPLES", sample_tuning.max_samples);
        m_contact_surface_samples = spcc::DenseSurfaceSampler::BuildFromMesh(
            *follower_mesh, static_cast<std::size_t>(std::max(0, cam_max_samples)), cam_surface_res);
        m_contact_pipeline.Configure(ResolveCompressedContactConfig(env_prefix, contact_regime));
        m_contact_pipeline.SetSlaveSurfaceSamples(m_contact_surface_samples);
        m_contact_mu_default = friction;

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            cam_body, m_follower, m_gear1_sdf.get(), &m_contact_pipeline, material, &m_reduced_contacts);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif
}

double ChronoRigidSystemNSC::GetFollowerPosY() const {
    if (m_follower) return m_follower->GetPos().y();
    return 0.0;
}

double ChronoRigidSystemNSC::GetFollowerVelY() const {
    if (m_follower) return m_follower->GetPosDt().y();
    return 0.0;
}

bool ChronoRigidSystemNSC::GetCamBodySnapshot(BodyDebugSnapshot& out) const {
    out = MakeBodyDebugSnapshot(m_cam_body);
    return out.valid;
}

bool ChronoRigidSystemNSC::GetFollowerBodySnapshot(BodyDebugSnapshot& out) const {
    out = MakeBodyDebugSnapshot(m_follower);
    return out.valid;
}

#if defined(SPCC_ENABLE_VDB)
const std::vector<platform::backend::spcc::ReducedContactPoint>& ChronoRigidSystemNSC::GetReducedContacts() const {
    return m_reduced_contacts;
}
#endif

void ChronoRigidSystemNSC::InitializeSimpleGearCase(
    const std::string& gear1_obj, const std::string& gear2_obj,
    const double* gear1_ref_pos, const double* gear2_ref_pos,
    const double* joint1_pos, const double* joint2_pos,
    double mesh_scale, double density,
    double friction, double restitution,
    double gravity_y, double motor_speed,
    double gear1_mass, const double* gear1_inertia_xx,
    double gear2_mass, const double* gear2_inertia_xx,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::CompressedContactConfig& contact_regime
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_reduced_contacts.clear();
    m_contact_surface_samples.clear();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution);

    auto load_scaled_mesh = [mesh_scale](const std::string& mesh_file) {
        auto mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        mesh->LoadWavefrontMesh(mesh_file, true, false);

        chrono::ChMatrix33<> scale_matrix(1);
        scale_matrix(0, 0) = mesh_scale;
        scale_matrix(1, 1) = mesh_scale;
        scale_matrix(2, 2) = mesh_scale;
        mesh->Transform(chrono::ChVector3d(0, 0, 0), scale_matrix);
        return mesh;
    };

    // Ground anchor for joints
    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    std::cout << "[DEBUG] Loading gear meshes..." << std::endl;
    auto mesh1 = load_scaled_mesh(gear1_obj);
    auto mesh2 = load_scaled_mesh(gear2_obj);
    std::cout << "[DEBUG] Loaded gear meshes" << std::endl;

#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
#endif

    if (use_vdb) {
        m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(mesh1, density, true, true, false, material, 0.000001);
        m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(mesh2, density, true, true, false, material, 0.000001);

        auto dummy_shape1 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        auto dummy_shape2 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        m_gear1->AddCollisionShape(dummy_shape1);
        m_gear1->EnableCollision(true);
        m_gear2->AddCollisionShape(dummy_shape2);
        m_gear2->EnableCollision(true);
    } else {
#if defined(SPCC_ENABLE_VDB)
        m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(mesh1, density, true, true, true, material, 0.000001);
        m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(mesh2, density, true, true, true, material, 0.000001);
#else
    m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh1, density, true, true, true, material, 0.000001);
    m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh2, density, true, true, true, material, 0.000001);
#endif
    }

    auto gear1_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_gear1);
    if (gear1_aux) {
        gear1_aux->SetFrameRefToAbs(
            chrono::ChFrame<>(
                chrono::ChVector3d(gear1_ref_pos[0], gear1_ref_pos[1], gear1_ref_pos[2]),
                chrono::QUNIT));
    }

    auto gear2_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_gear2);
    if (gear2_aux) {
        gear2_aux->SetFrameRefToAbs(
            chrono::ChFrame<>(
                chrono::ChVector3d(gear2_ref_pos[0], gear2_ref_pos[1], gear2_ref_pos[2]),
                chrono::QUNIT));
    }

    // Overwrite mass/inertia with values parsed from RecurDyn model.
    m_gear1->SetMass(gear1_mass);
    m_gear1->SetInertiaXX(
        chrono::ChVector3d(gear1_inertia_xx[0], gear1_inertia_xx[1], gear1_inertia_xx[2]));
    m_gear2->SetMass(gear2_mass);
    m_gear2->SetInertiaXX(
        chrono::ChVector3d(gear2_inertia_xx[0], gear2_inertia_xx[1], gear2_inertia_xx[2]));

    m_gear1->SetSleepingAllowed(false);
    m_gear2->SetSleepingAllowed(false);

    m_system->AddBody(m_gear1);
    m_system->AddBody(m_gear2);

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        m_dynamics_substeps = std::max(1, GetEnvInt("SPCC_GEAR_DYNAMICS_SUBSTEPS", 1));
        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetEnvDouble("SPCC_GEAR_VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetEnvDouble("SPCC_GEAR_HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*mesh1, sdf_options);

        const double gear_surface_res = GetEnvDouble("SPCC_GEAR_SURFACE_RES", sample_tuning.surface_res);
        const int gear_max_samples = GetEnvInt("SPCC_GEAR_MAX_SAMPLES", sample_tuning.max_samples);
        m_contact_surface_samples = spcc::DenseSurfaceSampler::BuildFromMesh(
            *mesh2, static_cast<std::size_t>(std::max(0, gear_max_samples)), gear_surface_res);
        m_contact_pipeline.Configure(ResolveCompressedContactConfig("SPCC_GEAR", contact_regime));
        m_contact_pipeline.SetSlaveSurfaceSamples(m_contact_surface_samples);
        m_contact_mu_default = friction;

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            m_gear1, m_gear2, m_gear1_sdf.get(), &m_contact_pipeline, material, &m_reduced_contacts);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif

    // RecurDyn results indicate RX is dominant, so use joint axis along global X.
    auto x_axis_revolute_rot = chrono::QuatFromAngleY(-chrono::CH_PI_2);

    auto motor = std::make_shared<chrono::ChLinkMotorRotationSpeed>();
    motor->Initialize(
        m_gear1,
        ground,
        chrono::ChFrame<>(
            chrono::ChVector3d(joint1_pos[0], joint1_pos[1], joint1_pos[2]),
            x_axis_revolute_rot));
    auto mfun = std::make_shared<chrono::ChFunctionConst>(motor_speed);
    motor->SetSpeedFunction(mfun);
    m_system->AddLink(motor);

    auto rev2 = std::make_shared<chrono::ChLinkLockRevolute>();
    rev2->Initialize(
        m_gear2,
        ground,
        chrono::ChFrame<>(
            chrono::ChVector3d(joint2_pos[0], joint2_pos[1], joint2_pos[2]),
            x_axis_revolute_rot));
    m_system->AddLink(rev2);
}

void ChronoRigidSystemNSC::InitializeRevoluteClearanceCase(
    const std::string& body1_obj,
    const std::string& body3_obj,
    const double* body1_pos,
    const double* body3_pos,
    const double* body2_cm_offset,
    double body3_mass,
    const double* body3_inertia_xx,
    const double* body3_inertia_xy,
    double body2_mass,
    const double* body2_inertia_xx,
    const double* body2_inertia_xy,
    double friction,
    double restitution,
    double gravity_y,
    double contact_compliance,
    double contact_compliance_t,
    double contact_damping_f,
    double collision_envelope,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::CompressedContactConfig& contact_regime) {
    m_sphere.reset();
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_reduced_contacts.clear();
    m_contact_surface_samples.clear();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction(static_cast<float>(friction));
    material->SetRestitution(static_cast<float>(restitution));
    material->SetCompliance(static_cast<float>(GetScopedEnvDouble(env_prefix, "COMPLIANCE", contact_compliance)));
    material->SetComplianceT(static_cast<float>(GetScopedEnvDouble(env_prefix, "COMPLIANCE_T", contact_compliance_t)));
    material->SetDampingF(static_cast<float>(GetScopedEnvDouble(env_prefix, "DAMPING_F", contact_damping_f)));

    const float resolved_collision_envelope =
        static_cast<float>(GetScopedEnvDouble(env_prefix, "COLLISION_ENVELOPE", collision_envelope));

    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
    (void)sdf_build_tuning;
    (void)sample_tuning;
    (void)contact_regime;
    (void)env_prefix;
    (void)dynamics_substeps;
#endif

    const auto body1_pos_W = chrono::ChVector3d(body1_pos[0], body1_pos[1], body1_pos[2]);
    const auto body3_pos_W = chrono::ChVector3d(body3_pos[0], body3_pos[1], body3_pos[2]);
    const auto body2_cm_ref = chrono::ChVector3d(body2_cm_offset[0], body2_cm_offset[1], body2_cm_offset[2]);

#if defined(SPCC_ENABLE_VDB)
    std::shared_ptr<chrono::ChTriangleMeshConnected> body1_mesh;
    std::shared_ptr<chrono::ChTriangleMeshConnected> body3_mesh;
    if (use_vdb) {
        body1_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        body1_mesh->LoadWavefrontMesh(body1_obj, true, false);
        body3_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        body3_mesh->LoadWavefrontMesh(body3_obj, true, false);
    }
#endif

    std::shared_ptr<chrono::ChBodyEasyMesh> body1;
    if (use_vdb) {
        body1 = std::make_shared<chrono::ChBodyEasyMesh>(body1_obj, 1000.0, true, true, false, material, 1.0e-6);
        auto dummy_shape = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        body1->AddCollisionShape(dummy_shape);
        body1->EnableCollision(true);
    } else {
        body1 = std::make_shared<chrono::ChBodyEasyMesh>(body1_obj, 1000.0, true, true, true, material, 1.0e-6);
    }
    body1->SetFixed(true);
    auto body1_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body1);
    if (body1_aux) {
        body1_aux->SetFrameRefToAbs(chrono::ChFrame<>(body1_pos_W, chrono::QUNIT));
    } else {
        body1->SetPos(body1_pos_W);
    }
    if (body1->GetCollisionModel()) {
        body1->GetCollisionModel()->SetEnvelope(resolved_collision_envelope);
    }
    m_system->AddBody(body1);
    m_clearance_body1 = body1;

    std::shared_ptr<chrono::ChBodyEasyMesh> body3;
    if (use_vdb) {
        body3 = std::make_shared<chrono::ChBodyEasyMesh>(body3_obj, 1000.0, true, true, false, material, 1.0e-6);
        auto dummy_shape = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        body3->AddCollisionShape(dummy_shape);
        body3->EnableCollision(true);
    } else {
        body3 = std::make_shared<chrono::ChBodyEasyMesh>(body3_obj, 1000.0, true, true, true, material, 1.0e-6);
    }
    auto body3_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body3);
    if (body3_aux) {
        body3_aux->SetFrameRefToAbs(chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    } else {
        body3->SetPos(body3_pos_W);
    }
    body3->SetMass(body3_mass);
    body3->SetInertiaXX(chrono::ChVector3d(body3_inertia_xx[0], body3_inertia_xx[1], body3_inertia_xx[2]));
    body3->SetInertiaXY(chrono::ChVector3d(body3_inertia_xy[0], body3_inertia_xy[1], body3_inertia_xy[2]));
    body3->SetSleepingAllowed(false);
    if (body3->GetCollisionModel()) {
        body3->GetCollisionModel()->SetEnvelope(resolved_collision_envelope);
    }
    m_system->AddBody(body3);
    m_clearance_body3 = body3;

    auto body2 = std::make_shared<chrono::ChBodyAuxRef>();
    body2->SetFixed(false);
    body2->EnableCollision(false);
    body2->SetMass(body2_mass);
    body2->SetInertiaXX(chrono::ChVector3d(body2_inertia_xx[0], body2_inertia_xx[1], body2_inertia_xx[2]));
    body2->SetInertiaXY(chrono::ChVector3d(body2_inertia_xy[0], body2_inertia_xy[1], body2_inertia_xy[2]));
    body2->SetFrameCOMToRef(chrono::ChFrame<>(body2_cm_ref, chrono::QUNIT));
    body2->SetFrameRefToAbs(chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    body2->SetSleepingAllowed(false);
    m_system->AddBody(body2);
    m_clearance_body2 = body2;

    auto fixed_link = std::make_shared<chrono::ChLinkLockLock>();
    fixed_link->Initialize(body2, body3, chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    m_system->AddLink(fixed_link);

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));
        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);

        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*body1_mesh, sdf_options);

        const double surface_res = GetScopedEnvDouble(env_prefix, "SURFACE_RES", sample_tuning.surface_res);
        const int max_samples = GetScopedEnvInt(env_prefix, "MAX_SAMPLES", sample_tuning.max_samples);
        m_contact_surface_samples = spcc::DenseSurfaceSampler::BuildFromMesh(
            *body3_mesh, static_cast<std::size_t>(std::max(0, max_samples)), surface_res);
        m_contact_pipeline.Configure(ResolveCompressedContactConfig(env_prefix, contact_regime));
        m_contact_pipeline.SetSlaveSurfaceSamples(m_contact_surface_samples);
        m_contact_mu_default = friction;

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            body1, body3, m_gear1_sdf.get(), &m_contact_pipeline, material, &m_reduced_contacts);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif
}

void ChronoRigidSystemNSC::InitializeHeadOnSphereCase(
    const std::string& sphere_a_obj,
    const std::string& sphere_b_obj,
    double sphere_radius,
    double sphere_a_density,
    double sphere_b_density,
    const double* sphere_a_pos,
    const double* sphere_b_pos,
    const double* sphere_a_vel,
    const double* sphere_b_vel,
    double friction,
    double restitution,
    double gravity_y,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::CompressedContactConfig& contact_regime) {
    m_sphere.reset();
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_reduced_contacts.clear();
    m_contact_surface_samples.clear();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));
    m_headon_sphere_radius = sphere_radius;
    m_headon_restitution = restitution;

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction(static_cast<float>(friction));
    material->SetRestitution(static_cast<float>(restitution));

    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    const auto sphere_a_pos_W = chrono::ChVector3d(sphere_a_pos[0], sphere_a_pos[1], sphere_a_pos[2]);
    const auto sphere_b_pos_W = chrono::ChVector3d(sphere_b_pos[0], sphere_b_pos[1], sphere_b_pos[2]);
    const auto sphere_a_vel_W = chrono::ChVector3d(sphere_a_vel[0], sphere_a_vel[1], sphere_a_vel[2]);
    const auto sphere_b_vel_W = chrono::ChVector3d(sphere_b_vel[0], sphere_b_vel[1], sphere_b_vel[2]);
    const auto x_axis_prismatic_rot = chrono::QuatFromAngleY(-chrono::CH_PI_2);
    const double sphere_volume =
        (4.0 / 3.0) * chrono::CH_PI * sphere_radius * sphere_radius * sphere_radius;
    const double sphere_a_mass = sphere_a_density * sphere_volume;
    const double sphere_b_mass = sphere_b_density * sphere_volume;
    const double sphere_a_inertia = 0.4 * sphere_a_mass * sphere_radius * sphere_radius;
    const double sphere_b_inertia = 0.4 * sphere_b_mass * sphere_radius * sphere_radius;

#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
    m_headon_manual_elastic_resolution =
        use_vdb && (GetScopedEnvDouble(env_prefix, "USE_EXACT_CORRECTION", 0.0) > 0.5);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
    (void)sdf_build_tuning;
    (void)sample_tuning;
    (void)contact_regime;
    (void)env_prefix;
    (void)dynamics_substeps;
#endif

    if (!use_vdb) {
        auto sphere_a =
            std::make_shared<chrono::ChBodyEasySphere>(sphere_radius, sphere_a_density, true, true, material);
        auto sphere_b =
            std::make_shared<chrono::ChBodyEasySphere>(sphere_radius, sphere_b_density, true, true, material);
        sphere_a->SetFixed(false);
        sphere_b->SetFixed(false);
        sphere_a->SetSleepingAllowed(false);
        sphere_b->SetSleepingAllowed(false);
        sphere_a->SetPos(sphere_a_pos_W);
        sphere_b->SetPos(sphere_b_pos_W);
        sphere_a->SetPosDt(sphere_a_vel_W);
        sphere_b->SetPosDt(sphere_b_vel_W);
        sphere_a->SetMass(sphere_a_mass);
        sphere_b->SetMass(sphere_b_mass);
        sphere_a->SetInertiaXX(chrono::ChVector3d(sphere_a_inertia, sphere_a_inertia, sphere_a_inertia));
        sphere_b->SetInertiaXX(chrono::ChVector3d(sphere_b_inertia, sphere_b_inertia, sphere_b_inertia));
        m_system->AddBody(sphere_a);
        m_system->AddBody(sphere_b);
        m_headon_a = sphere_a;
        m_headon_b = sphere_b;
    }

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        auto mesh_a = std::make_shared<chrono::ChTriangleMeshConnected>();
        auto mesh_b = std::make_shared<chrono::ChTriangleMeshConnected>();
        bool loaded_a = !sphere_a_obj.empty() && mesh_a->LoadWavefrontMesh(sphere_a_obj, true, false);
        bool loaded_b = !sphere_b_obj.empty() && mesh_b->LoadWavefrontMesh(sphere_b_obj, true, false);
        if (!loaded_a) {
            mesh_a = BuildUvSphereMesh(sphere_radius, 64, 32);
        }
        if (!loaded_b) {
            mesh_b = BuildUvSphereMesh(sphere_radius, 64, 32);
        }

        auto sphere_a = std::make_shared<chrono::ChBodyEasyMesh>(mesh_a, sphere_a_density, true, true, false, material,
                                                                 1.0e-6);
        auto sphere_b = std::make_shared<chrono::ChBodyEasyMesh>(mesh_b, sphere_b_density, true, true, false, material,
                                                                 1.0e-6);

        auto dummy_shape_a = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        auto dummy_shape_b = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        sphere_a->AddCollisionShape(dummy_shape_a);
        sphere_a->EnableCollision(true);
        sphere_b->AddCollisionShape(dummy_shape_b);
        sphere_b->EnableCollision(true);

        SetBodyFrameOrPos(sphere_a, sphere_a_pos_W);
        SetBodyFrameOrPos(sphere_b, sphere_b_pos_W);
        sphere_a->SetFixed(false);
        sphere_b->SetFixed(false);
        sphere_a->SetSleepingAllowed(false);
        sphere_b->SetSleepingAllowed(false);
        sphere_a->SetPosDt(sphere_a_vel_W);
        sphere_b->SetPosDt(sphere_b_vel_W);
        sphere_a->SetMass(sphere_a_mass);
        sphere_b->SetMass(sphere_b_mass);
        sphere_a->SetInertiaXX(chrono::ChVector3d(sphere_a_inertia, sphere_a_inertia, sphere_a_inertia));
        sphere_b->SetInertiaXX(chrono::ChVector3d(sphere_b_inertia, sphere_b_inertia, sphere_b_inertia));
        m_system->AddBody(sphere_a);
        m_system->AddBody(sphere_b);
        m_headon_a = sphere_a;
        m_headon_b = sphere_b;

        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));

        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);

        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*mesh_a, sdf_options);

        m_contact_surface_samples.clear();
        spcc::DenseSurfaceSample sample;
        sample.xi_slave_S = chrono::ChVector3d(-sphere_radius, 0.0, 0.0);
        sample.normal_slave_S = chrono::ChVector3d(-1.0, 0.0, 0.0);
        sample.area_weight = 4.0 * chrono::CH_PI * sphere_radius * sphere_radius;
        m_contact_surface_samples.push_back(sample);
        auto resolved_contact_regime = ResolveCompressedContactConfig(env_prefix, contact_regime);
        m_contact_pipeline.Configure(resolved_contact_regime);
        m_contact_pipeline.SetSlaveSurfaceSamples(m_contact_surface_samples);
        m_contact_mu_default = friction;

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            sphere_a, sphere_b, m_gear1_sdf.get(), &m_contact_pipeline, material, &m_reduced_contacts,
            true, resolved_contact_regime.analytic_sphere_toi_contact, sphere_radius, sphere_radius);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif

    auto prism_a = std::make_shared<chrono::ChLinkLockPrismatic>();
    prism_a->Initialize(m_headon_a, ground, chrono::ChFrame<>(m_headon_a->GetPos(), x_axis_prismatic_rot));
    m_system->AddLink(prism_a);

    auto prism_b = std::make_shared<chrono::ChLinkLockPrismatic>();
    prism_b->Initialize(m_headon_b, ground, chrono::ChFrame<>(m_headon_b->GetPos(), x_axis_prismatic_rot));
    m_system->AddLink(prism_b);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetGear1AngVelParent() const {
    if (m_gear1) return m_gear1->GetAngVelParent();
    return chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetGear2AngVelParent() const {
    if (m_gear2) return m_gear2->GetAngVelParent();
    return chrono::ChVector3d(0, 0, 0);
}

} // namespace backend
} // namespace platform


