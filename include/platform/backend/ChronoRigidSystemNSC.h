#pragma once

#include "platform/backend/IRigidSystem.h"
#include "platform/common/ContactAlgorithm.h"
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/physics/ChBodyEasy.h>
#include <chrono/core/ChQuaternion.h>
#include <string>
#include <memory>
#include <vector>

#if defined(SPCC_ENABLE_VDB)
#include "platform/backend/spcc/ContactActivation.h"
#include "platform/backend/spcc/ContactTuning.h"
#include "platform/backend/spcc/SampleBVH.h"
#include "platform/backend/spcc/VDBSDFField.h"
#endif

namespace platform {
namespace models {
struct FollowerPreloadConfig;
}
namespace backend {

struct BodyDebugSnapshot {
    bool valid = false;
    chrono::ChVector3d x_com_W;
    chrono::ChVector3d x_ref_W;
    chrono::ChQuaternion<> q_WL = chrono::QUNIT;
    chrono::ChQuaternion<> q_WRef = chrono::QUNIT;
    chrono::ChVector3d v_com_W;
    chrono::ChVector3d w_W;
};

class ChronoRigidSystemNSC : public IRigidSystem {
public:
    ChronoRigidSystemNSC();
    virtual ~ChronoRigidSystemNSC() override;

    // Use this specific initialize for this basic drop case
    void Initialize(
        double ball_radius, double ball_density, double ball_height,
        double ground_x, double ground_y, double ground_z,
        double friction, double restitution, double gravity_y
    );

    // Initializer for the Cam case
    void InitializeCamCase(
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
    );

    // Initializer for the Simple Gear case
    void InitializeSimpleGearCase(
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
    );

    void InitializeHeadOnSphereCase(
        const std::string& sphere_a_obj, const std::string& sphere_b_obj,
        double sphere_radius, double sphere_a_density, double sphere_b_density,
        const double* sphere_a_pos, const double* sphere_b_pos,
        const double* sphere_a_vel, const double* sphere_b_vel,
        double friction, double restitution, double gravity_y,
        int dynamics_substeps,
        const std::string& env_prefix,
        platform::common::ContactAlgorithm contact_algorithm,
        const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
        const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
        const platform::backend::spcc::ContactRegimeConfig& contact_regime);

    void InitializeRevoluteClearanceCase(
        const std::string& body1_obj, const std::string& body3_obj,
        const double* body1_pos, const double* body3_pos, const double* body2_cm_offset,
        double body3_mass, const double* body3_inertia_xx, const double* body3_inertia_xy,
        double body2_mass, const double* body2_inertia_xx, const double* body2_inertia_xy,
        double friction, double restitution, double gravity_y,
        double contact_compliance, double contact_compliance_t, double contact_damping_f, double collision_envelope,
        int dynamics_substeps, const std::string& env_prefix,
        platform::common::ContactAlgorithm contact_algorithm,
        const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
        const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
        const platform::backend::spcc::ContactRegimeConfig& contact_regime,
        bool use_lcp_manifold_quadrature,
        int manifold_quadrature_contacts,
        double manifold_quadrature_span_scale,
        double manifold_quadrature_min_half_span);

    // Legacy initialize to satisfy IRigidSystem base if needed
    void Initialize() override;

    void StepDynamics(double step_size) override;
    double GetTime() const override;
    unsigned int GetNumContacts() const override;

    // Direct accessors for demonstration printing (Drop Case)
    double GetDynamicSpherePosY() const;
    double GetDynamicSphereVelY() const;
    double GetHeadOnSphereAPosX() const;
    double GetHeadOnSphereAVelX() const;
    double GetHeadOnSphereBPosX() const;
    double GetHeadOnSphereBVelX() const;
    chrono::ChVector3d GetClearanceBody2Pos() const;
    chrono::ChVector3d GetClearanceBody2Vel() const;
    chrono::ChVector3d GetClearanceBody3Pos() const;
    chrono::ChVector3d GetClearanceBody3Vel() const;
    chrono::ChVector3d GetClearanceBody3AngVel() const;

    // Direct accessors for Cam Case
    double GetFollowerPosY() const;
    double GetFollowerVelY() const;
    bool GetCamBodySnapshot(BodyDebugSnapshot& out) const;
    bool GetFollowerBodySnapshot(BodyDebugSnapshot& out) const;
#if defined(SPCC_ENABLE_VDB)
    const std::vector<platform::backend::spcc::ActiveContactSample>& GetActiveContacts() const;
#endif

    // Direct accessors for Simple Gear Case (angular velocity in parent/world frame)
    chrono::ChVector3d GetGear1AngVelParent() const;
    chrono::ChVector3d GetGear2AngVelParent() const;

private:
    std::unique_ptr<chrono::ChSystemNSC> m_system;
    std::shared_ptr<chrono::ChBody> m_sphere;      // used in drop case
    std::shared_ptr<chrono::ChBody> m_headon_a;    // used in head-on sphere case
    std::shared_ptr<chrono::ChBody> m_headon_b;    // used in head-on sphere case
    std::shared_ptr<chrono::ChBody> m_cam_body;    // used in cam case
    std::shared_ptr<chrono::ChBody> m_follower;    // used in cam case
    std::shared_ptr<chrono::ChBody> m_gear1;       // used in simple gear case
    std::shared_ptr<chrono::ChBody> m_gear2;       // used in simple gear case
    std::shared_ptr<chrono::ChBody> m_clearance_body1;
    std::shared_ptr<chrono::ChBody> m_clearance_body2;
    std::shared_ptr<chrono::ChBody> m_clearance_body3;

#if defined(SPCC_ENABLE_VDB)
    std::unique_ptr<platform::backend::spcc::VDBSDFField> m_gear1_sdf;
    std::shared_ptr<platform::backend::spcc::SampleBVH> m_sample_bvh;
    platform::backend::spcc::ContactActivation m_contact_activation;
    std::vector<chrono::ChVector3d> m_gear2_samples_S;
    std::vector<platform::backend::spcc::ActiveContactSample> m_active_contacts;
    double m_contact_mu_default = 0.0;
#endif
    int m_dynamics_substeps = 1;
    bool m_headon_manual_elastic_resolution = false;
    double m_headon_sphere_radius = 0.0;
    double m_headon_restitution = 1.0;
};

} // namespace backend
} // namespace platform
