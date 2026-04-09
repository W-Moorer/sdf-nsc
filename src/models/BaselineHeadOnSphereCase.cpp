#include "platform/models/BaselineHeadOnSphereCase.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace platform {
namespace models {

BaselineHeadOnSphereCase::BaselineHeadOnSphereCase(const HeadOnSphereCaseConfig& config) : m_config(config) {
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineHeadOnSphereCase::SetupSystem() {
    m_backend->InitializeHeadOnSphereCase(
        m_config.sphere_a_mesh_path, m_config.sphere_b_mesh_path, m_config.sphere_radius, m_config.sphere_a_density,
        m_config.sphere_b_density,
        m_config.sphere_a_init_pos, m_config.sphere_b_init_pos, m_config.sphere_a_init_vel, m_config.sphere_b_init_vel,
        m_config.friction, m_config.restitution, m_config.gravity_y, m_config.dynamics_substeps, m_config.env_prefix,
        m_config.contact_algorithm, m_config.sdf_build, m_config.sample_tuning, m_config.contact_regime);
}

void BaselineHeadOnSphereCase::SaveCSV(double t,
                                       double pos_a_x,
                                       double vel_a_x,
                                       double pos_b_x,
                                       double vel_b_x,
                                       unsigned int nc) {
    m_csv_buffer += std::to_string(t) + "," + std::to_string(pos_a_x) + "," + std::to_string(vel_a_x) + "," +
                    std::to_string(pos_b_x) + "," + std::to_string(vel_b_x) + "," + std::to_string(nc) + "\n";
}

void BaselineHeadOnSphereCase::Run() {
    SetupSystem();

    std::cout << "===========================================" << std::endl;
    std::cout << " Head-On Sphere Collision (NSC)" << std::endl;
    std::cout << " - Contact Algorithm: "
              << platform::common::ContactAlgorithmToCliName(m_config.contact_algorithm) << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << " s, End Time: " << m_config.total_time << " s"
              << std::endl;
    std::cout << " - Density A / B: " << m_config.sphere_a_density << " / " << m_config.sphere_b_density
              << " kg/m^3" << std::endl;
    std::cout << " - Initial Vel A: " << m_config.sphere_a_init_vel[0] << " m/s, Initial Vel B: "
              << m_config.sphere_b_init_vel[0] << " m/s" << std::endl;
    std::cout << "===========================================" << std::endl;

    m_csv_buffer = "Time,SphereA_Pos_X,SphereA_Vel_X,SphereB_Pos_X,SphereB_Vel_X,Num_Contacts\n";
    std::cout << std::fixed << std::setprecision(6);

    while (m_backend->GetTime() < m_config.total_time) {
        m_backend->StepDynamics(m_config.step_size);

        const double t = m_backend->GetTime();
        const double pos_a_x = m_backend->GetHeadOnSphereAPosX();
        const double vel_a_x = m_backend->GetHeadOnSphereAVelX();
        const double pos_b_x = m_backend->GetHeadOnSphereBPosX();
        const double vel_b_x = m_backend->GetHeadOnSphereBVelX();
        const unsigned int nc = m_backend->GetNumContacts();

        if (std::fmod(t, 0.05) < (m_config.step_size * 0.5) || t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << "Time: " << std::setw(8) << t << " s | A(vx): " << std::setw(10) << vel_a_x
                      << " | B(vx): " << std::setw(10) << vel_b_x << " | Contacts: " << nc << std::endl;
        }

        SaveCSV(t, pos_a_x, vel_a_x, pos_b_x, vel_b_x, nc);
    }

    const std::filesystem::path out_path(m_config.output_csv_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(out_path);
    if (out.is_open()) {
        out << m_csv_buffer;
        out.close();
        std::cout << "\nResults saved to: " << m_config.output_csv_path << std::endl;
    }

    std::cout << "===========================================" << std::endl;
    std::cout << " SUMMARY:" << std::endl;
    std::cout << " - Final Sphere A Velocity X : " << m_backend->GetHeadOnSphereAVelX() << " m/s" << std::endl;
    std::cout << " - Final Sphere B Velocity X : " << m_backend->GetHeadOnSphereBVelX() << " m/s" << std::endl;
    std::cout << " - Output CSV path           : " << m_config.output_csv_path << std::endl;
    std::cout << "===========================================" << std::endl;
}

}  // namespace models
}  // namespace platform
