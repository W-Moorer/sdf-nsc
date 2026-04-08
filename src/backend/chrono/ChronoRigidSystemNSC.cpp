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
#include "platform/backend/spcc/ContactActivation.h"
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

std::vector<chrono::ChVector3d> BuildDownsampledVertexSamples(const chrono::ChTriangleMeshConnected& mesh,
                                                               std::size_t max_samples,
                                                               double surface_res = 0.0005) {
    const auto& vertices = mesh.GetCoordsVertices();
    const auto& faces = mesh.GetIndicesVertexes();
    std::vector<chrono::ChVector3d> samples;
    if (vertices.empty()) {
        std::cout << "[DEBUG] BuildDownsampledVertexSamples returning 0 samples! Mesh is empty!" << std::endl;
        return samples;
    }

    // Add original vertices
    for (const auto& v : vertices) {
        if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
            continue;
        }
        samples.push_back(v);
    }
    
    // Add extra samples inside triangles to prevent tunneling/aliasing
    // We sample faces by barycentric coordinates to create a denser point cloud
    for (const auto& face : faces) {
        auto v0 = vertices[face.x()];
        auto v1 = vertices[face.y()];
        auto v2 = vertices[face.z()];
        
        // Edge lengths
        double l01 = (v1 - v0).Length();
        double l12 = (v2 - v1).Length();
        double l20 = (v0 - v2).Length();
        double area = chrono::Vcross(v1 - v0, v2 - v0).Length() * 0.5;
        
        // Target face-interior sample spacing for contact manifold coverage.
        const double res = std::max(1.0e-5, surface_res);
        int target_pts = std::min(500, std::max(0, (int)(area / (res * res * 0.5))));
        
        if (target_pts > 0) {
            // Systematic or random barycentric sampling. Using a simple grid approach on u,v
            int steps = std::max(2, (int)std::sqrt(target_pts * 2));
            for (int i = 1; i < steps; ++i) {
                for (int j = 1; j < steps - i; ++j) {
                    double u = (double)i / steps;
                    double v = (double)j / steps;
                    double w = 1.0 - u - v;
                    chrono::ChVector3d p = v0 * u + v1 * v + v2 * w;
                    samples.push_back(p);
                }
            }
        }
    }

    // Downsample if we still exceeded max_samples by a lot, but usually we want all of them.
    if (samples.size() > max_samples && max_samples > 0) {
        std::vector<chrono::ChVector3d> filtered;
        filtered.reserve(max_samples);
        // keep all original vertices
        std::size_t num_orig = vertices.size();
        for (std::size_t i = 0; i < num_orig; ++i) {
            filtered.push_back(samples[i]);
        }
        std::size_t remaining = max_samples > num_orig ? max_samples - num_orig : 0;
        if (remaining > 0) {
            std::size_t generated = samples.size() - num_orig;
            std::size_t stride = std::max((std::size_t)1, generated / remaining);
            for(std::size_t i = num_orig; i < samples.size(); i += stride) {
                filtered.push_back(samples[i]);
                if(filtered.size() >= max_samples) break;
            }
        }
        std::cout << "[DEBUG] BuildDownsampledVertexSamples returning " << filtered.size() << " samples. downsampled from " << samples.size() << std::endl;
        return filtered;
    }

    std::cout << "[DEBUG] BuildDownsampledVertexSamples returning " << samples.size() << " samples." << std::endl;
    return samples;
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

class SDFCollisionCallback : public chrono::ChSystem::CustomCollisionCallback {
public:
    std::shared_ptr<chrono::ChBody> m_master;
    std::shared_ptr<chrono::ChBody> m_slave;
    spcc::VDBSDFField* m_sdf;
    spcc::ContactActivation* m_activation;
    std::vector<chrono::ChVector3d> m_slave_vertices_local;
    std::shared_ptr<chrono::ChContactMaterial> m_material;
    std::vector<spcc::ActiveContactSample>* m_active_contacts_out;
    bool m_enable_curvature_term;

    SDFCollisionCallback(std::shared_ptr<chrono::ChBody> master,
                         std::shared_ptr<chrono::ChBody> slave,
                         spcc::VDBSDFField* sdf,
                         spcc::ContactActivation* activation,
                         const std::vector<chrono::ChVector3d>& slave_vertices,
                         std::shared_ptr<chrono::ChContactMaterial> mat,
                         std::vector<spcc::ActiveContactSample>* active_contacts_out,
                         bool enable_curvature_term)
        : m_master(master), m_slave(slave), m_sdf(sdf), m_activation(activation),
          m_slave_vertices_local(slave_vertices), m_material(mat),
          m_active_contacts_out(active_contacts_out),
          m_enable_curvature_term(enable_curvature_term) {}

    void OnCustomCollision(chrono::ChSystem* sys) override {
        if (!m_master || !m_slave || !m_sdf || !m_activation || m_slave_vertices_local.empty()) {
            if (m_active_contacts_out) {
                m_active_contacts_out->clear();
            }
            return;
        }

        // 1. Build world samples W
        auto world_samples_W = BuildWorldSamplesFromLocal(m_slave, m_slave_vertices_local);
        
        // 2. Build kinematic states
        auto master_state = MakeRigidBodyStateW(m_master, 1);
        auto slave_state = MakeRigidBodyStateW(m_slave, 2);

        // 3. Build active set with ContactActivation
        std::vector<spcc::ActiveContactSample> active_contacts;
        int penetration_count = 0;
        m_activation->BuildActiveSet(master_state, slave_state, *m_sdf,
                                     m_slave_vertices_local, world_samples_W,
                                     (float)std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetSlidingFriction(),
                                     sys->GetStep(),
                                     m_enable_curvature_term,
                                     active_contacts);
        for (auto& contact : active_contacts) {
            penetration_count++;
            
            chrono::ChCollisionInfo cinfo;
            cinfo.modelA = m_master->GetCollisionModel().get();
            cinfo.modelB = m_slave->GetCollisionModel().get();
            cinfo.shapeA = nullptr;
            cinfo.shapeB = nullptr;
            cinfo.vN = contact.n_W; 
            // n_W points from master to slave.
            
            // Relative velocity of the slave point relative to the master
            chrono::ChVector3d v_master = master_state.v_com_W + chrono::Vcross(master_state.w_W, contact.x_W - master_state.x_com_W);
            chrono::ChVector3d v_slave  = slave_state.v_com_W + chrono::Vcross(slave_state.w_W, contact.x_W - slave_state.x_com_W);
            chrono::ChVector3d v_rel = v_slave - v_master;

            double curvature_term = 0.0;
            if (m_enable_curvature_term) {
                const double dt = sys->GetStep();
                chrono::ChMatrix33<> H_eff = contact.hessian_W;
                if (contact.curvature_tangential_only) {
                    H_eff = contact.P_W * H_eff * contact.P_W;
                }
                curvature_term = 0.5 * dt * dt * chrono::Vdot(v_rel, H_eff * v_rel);
                curvature_term *= contact.curvature_gate;
                if (!std::isfinite(curvature_term)) {
                    curvature_term = 0.0;
                }
                if (contact.curvature_term_abs_max > 0.0 && std::isfinite(contact.curvature_term_abs_max)) {
                    curvature_term =
                        std::clamp(curvature_term, -contact.curvature_term_abs_max, contact.curvature_term_abs_max);
                }
                if (contact.curvature_term_ratio_max > 0.0 && std::isfinite(contact.curvature_term_ratio_max)) {
                    const double gap_scale =
                        std::max(std::abs(contact.phi), std::max(0.0, contact.curvature_gap_floor));
                    const double ratio_cap = contact.curvature_term_ratio_max * gap_scale;
                    curvature_term = std::clamp(curvature_term, -ratio_cap, ratio_cap);
                }
            }
            cinfo.vpB = contact.x_W; // point on slave
            
            // Inject the curvature compensation into the perceived distance
            cinfo.distance = contact.phi + curvature_term;

            contact.x_master_surface_W = cinfo.vpA;
            contact.v_rel_W = v_rel;
            contact.curvature_term = curvature_term;
            contact.phi_eff = cinfo.distance;
             
            sys->GetContactContainer()->AddContact(cinfo, m_material, m_material);
        }

        if (m_active_contacts_out) {
            *m_active_contacts_out = active_contacts;
        }
        
        static int step_id = 0;
        if (step_id % 10 == 0) {
             std::cout << "[SDF] Step " << step_id << " active contacts: " << penetration_count;
#if defined(SPCC_ENABLE_VDB)
             if (GetEnvDouble("SPCC_DEBUG_CONTACT_STATS", 0.0) > 0.5) {
                 const auto& stats = m_activation->GetStats();
                 std::cout << " queried=" << stats.queried << " after_cap=" << stats.accepted_after_cap
                           << " fit_attempted=" << stats.local_fit_attempted
                           << " fit_applied=" << stats.local_fit_applied
                           << " fit_reject=" << stats.local_fit_rejected_positive_gap;
             }
#endif
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
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
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
    const int substeps = std::max(1, m_dynamics_substeps);
    const double sub_dt = step_size / static_cast<double>(substeps);
    for (int i = 0; i < substeps; ++i) {
        m_system->DoStepDynamics(sub_dt);
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
    const platform::backend::spcc::ContactRegimeConfig& contact_regime,
    const platform::models::FollowerPreloadConfig& follower_preload
) {
    m_cam_body.reset();
    m_follower.reset();
    m_dynamics_substeps = 1;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
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

    const bool enable_curvature_term =
        (contact_algorithm == platform::common::ContactAlgorithm::SdfSecondOrder);
#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = true;
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
        sdf_options.direct_phi_hessian =
            (GetScopedEnvDouble(env_prefix, "DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) >
             0.5);
        
        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>(); // reusing the gear1_sdf for Cam
        std::cout << "[DEBUG] Building VDB SDF for Cam with wide band..." << std::endl;
        m_gear1_sdf->BuildFromTriangleMesh(*cam_mesh, sdf_options);
        
        m_contact_activation.SetEnvPrefix(env_prefix);
        auto resolved_contact_regime = contact_regime;
        resolved_contact_regime.activation.delta_on =
            GetScopedEnvDouble(env_prefix, "DELTA_ON", resolved_contact_regime.activation.delta_on);
        resolved_contact_regime.activation.delta_off =
            GetScopedEnvDouble(env_prefix, "DELTA_OFF", resolved_contact_regime.activation.delta_off);
        resolved_contact_regime.activation.hold_steps =
            GetScopedEnvInt(env_prefix, "HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
        resolved_contact_regime.activation.max_active_keep =
            GetScopedEnvInt(env_prefix, "MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
        resolved_contact_regime.curvature.enabled =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_ENABLED",
                                resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.tangential_only =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_TANGENTIAL_ONLY",
                                resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.normal_alignment_cos_min =
            GetScopedEnvDouble(env_prefix, "CURVATURE_ALIGNMENT_COS_MIN",
                               resolved_contact_regime.curvature.normal_alignment_cos_min);
        resolved_contact_regime.curvature.max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.max_hessian_frobenius);
        resolved_contact_regime.curvature.max_curvature_term_abs =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_ABS",
                               resolved_contact_regime.curvature.max_curvature_term_abs);
        resolved_contact_regime.curvature.max_curvature_term_ratio =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_RATIO",
                               resolved_contact_regime.curvature.max_curvature_term_ratio);
        resolved_contact_regime.curvature.gap_floor =
            GetScopedEnvDouble(env_prefix, "CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
        m_contact_activation.SetPolicy(resolved_contact_regime);

        const double cam_surface_res = GetScopedEnvDouble(env_prefix, "SURFACE_RES", sample_tuning.surface_res);
        const int cam_max_samples = GetScopedEnvInt(env_prefix, "MAX_SAMPLES", sample_tuning.max_samples);
        std::vector<chrono::ChVector3d> slave_samples =
            BuildDownsampledVertexSamples(*follower_mesh, static_cast<std::size_t>(std::max(0, cam_max_samples)),
                                          cam_surface_res); 

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            cam_body, m_follower, m_gear1_sdf.get(), &m_contact_activation, slave_samples, material,
            &m_active_contacts, enable_curvature_term);
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
const std::vector<platform::backend::spcc::ActiveContactSample>& ChronoRigidSystemNSC::GetActiveContacts() const {
    return m_active_contacts;
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
    bool enable_curvature_term,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::ContactRegimeConfig& contact_regime
) {
    m_cam_body.reset();
    m_follower.reset();
    m_dynamics_substeps = 1;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
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
    std::cout << "[DEBUG] Using VDB SDF instead of NSC Mesh Collision" << std::endl;
    // Set create_collision flag (5th param) to FALSE to avoid HACD OOM
    // Params: mesh, density, compute_mass, create_visualization, create_collision, material, sphere_swept
    m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh1, density, true, true, false, material, 0.000001);
    m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh2, density, true, true, false, material, 0.000001);

    // Provide a dummy collision shape so the bodies have a valid ChCollisionModel
    // otherwise the custom SDF contact injection will crash due to null models.
    auto dummy_shape1 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
    m_gear1->AddCollisionShape(dummy_shape1);
    m_gear1->EnableCollision(true);

    auto dummy_shape2 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
    m_gear2->AddCollisionShape(dummy_shape2);
    m_gear2->EnableCollision(true);

    // Build SDF for gear 1
    spcc::VDBSDFField::BuildOptions sdf_options;
    // Use an adequate voxel size based on mesh_scale or global size.
    // Make sure it's smaller than the gear mesh features!
    sdf_options.voxel_size = GetEnvDouble("SPCC_GEAR_VOXEL_SIZE", sdf_build_tuning.voxel_size);
    sdf_options.half_band_width_voxels =
        GetEnvDouble("SPCC_GEAR_HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
    sdf_options.direct_phi_hessian =
        (GetEnvDouble("SPCC_GEAR_DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) > 0.5);
    m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
    std::cout << "[DEBUG] Building VDB SDF for gear 1..." << std::endl;
    m_gear1_sdf->BuildFromTriangleMesh(*mesh1, sdf_options);
    
    // Configure contact activation
    m_contact_activation.SetEnvPrefix("SPCC_GEAR");
    auto resolved_contact_regime = contact_regime;
    resolved_contact_regime.activation.delta_on =
        GetEnvDouble("SPCC_GEAR_DELTA_ON", resolved_contact_regime.activation.delta_on);
    resolved_contact_regime.activation.delta_off =
        GetEnvDouble("SPCC_GEAR_DELTA_OFF", resolved_contact_regime.activation.delta_off);
    resolved_contact_regime.activation.hold_steps =
        GetEnvInt("SPCC_GEAR_HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
    resolved_contact_regime.activation.max_active_keep =
        GetEnvInt("SPCC_GEAR_MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
    resolved_contact_regime.curvature.enabled =
        (GetEnvDouble("SPCC_GEAR_CURVATURE_ENABLED", resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
    resolved_contact_regime.curvature.tangential_only =
        (GetEnvDouble("SPCC_GEAR_CURVATURE_TANGENTIAL_ONLY",
                      resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
    resolved_contact_regime.curvature.normal_alignment_cos_min =
        GetEnvDouble("SPCC_GEAR_CURVATURE_ALIGNMENT_COS_MIN",
                     resolved_contact_regime.curvature.normal_alignment_cos_min);
    resolved_contact_regime.curvature.max_hessian_frobenius =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_HESSIAN_FROBENIUS",
                     resolved_contact_regime.curvature.max_hessian_frobenius);
    resolved_contact_regime.curvature.max_curvature_term_abs =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_ABS", resolved_contact_regime.curvature.max_curvature_term_abs);
    resolved_contact_regime.curvature.max_curvature_term_ratio =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_RATIO", resolved_contact_regime.curvature.max_curvature_term_ratio);
    resolved_contact_regime.curvature.gap_floor =
        GetEnvDouble("SPCC_GEAR_CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
    m_contact_activation.SetPolicy(resolved_contact_regime);

    // Register custom collision callback and sample slave vertices
    std::vector<chrono::ChVector3d> slave_samples;
    // Downsample slightly or use all vertices
    const double gear_surface_res = GetEnvDouble("SPCC_GEAR_SURFACE_RES", sample_tuning.surface_res);
    const int gear_max_samples = GetEnvInt("SPCC_GEAR_MAX_SAMPLES", sample_tuning.max_samples);
    slave_samples = BuildDownsampledVertexSamples(*mesh2, static_cast<std::size_t>(std::max(0, gear_max_samples)),
                                                  gear_surface_res); 

    auto sdf_callback = std::make_shared<SDFCollisionCallback>(
        m_gear1, m_gear2, m_gear1_sdf.get(), &m_contact_activation, slave_samples, material, &m_active_contacts,
        enable_curvature_term);
    m_system->RegisterCustomCollisionCallback(sdf_callback);
    std::cout << "[DEBUG] VDB SDF initialized successfully." << std::endl;
#else
    m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh1, density, true, true, true, material, 0.000001);
    std::cout << "[DEBUG] Created ChBodyEasyMesh gear1" << std::endl;
    m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh2, density, true, true, true, material, 0.000001);
    std::cout << "[DEBUG] Created ChBodyEasyMesh gear2" << std::endl;
#endif

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


