#include "platform/models/BaselineDropCase.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace platform {
namespace models {

BaselineDropCase::BaselineDropCase(const DropCaseConfig& config)
    : m_config(config),
      m_min_y(1e9), m_final_y(0), m_final_vy(0),
      m_had_contact(false), m_first_contact_time(-1.0)
{
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineDropCase::SetupSystem() {
    m_backend->Initialize(
        m_config.ball_radius, m_config.ball_density, m_config.ball_initial_height,
        m_config.ground_size_x, m_config.ground_size_y, m_config.ground_size_z,
        m_config.friction, m_config.restitution, m_config.gravity
    );
}

void BaselineDropCase::SaveCSV(double t, double p_y, double v_y, unsigned int nc) {
    // Append to in-memory buffer
    m_csv_buffer += std::to_string(t) + "," +
                    std::to_string(p_y) + "," +
                    std::to_string(v_y) + "," +
                    std::to_string(nc) + "\n";
}

void BaselineDropCase::Run() {
    SetupSystem();

    std::cout << "===========================================" << std::endl;
    std::cout << " Baseline Drop Case (NSC)" << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << "s, End Time: " << m_config.total_time << "s" << std::endl;
    std::cout << " - Pos Y Init: " << m_config.ball_initial_height << "m" << std::endl;
    std::cout << "===========================================" << std::endl;

    m_csv_buffer = "Time,Pos_Y,Vel_Y,Num_Contacts\n";

    std::cout << std::fixed << std::setprecision(4);

    while (m_backend->GetTime() < m_config.total_time) {
        m_backend->StepDynamics(m_config.step_size);
        
        double t = m_backend->GetTime();
        double p_y = m_backend->GetDynamicSpherePosY();
        double v_y = m_backend->GetDynamicSphereVelY();
        unsigned int nc = m_backend->GetNumContacts();
        
        if (nc > 0 && !m_had_contact) {
            m_had_contact = true;
            m_first_contact_time = t;
        }
        
        m_min_y = std::min(m_min_y, p_y);
        m_final_y = p_y;
        m_final_vy = v_y;
        
        // Console output every 0.1s roughly
        if (std::abs(std::fmod(t, 0.1)) < (m_config.step_size * 0.5) || t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << "Time: " << std::setw(6) << t 
                      << " s | Pos Y: " << std::setw(8) << p_y 
                      << " m | Vel Y: " << std::setw(8) << v_y << " m/s" 
                      << " | Contacts: " << nc << std::endl;
        }

        SaveCSV(t, p_y, v_y, nc);
    }

    // Write to disk
    std::ofstream out(m_config.output_csv_path);
    if (out.is_open()) {
        out << m_csv_buffer;
        out.close();
        std::cout << "\nResults saved to: " << m_config.output_csv_path << std::endl;
    }

    // Console Summary
    std::cout << "===========================================" << std::endl;
    std::cout << " SUMMARY:" << std::endl;
    std::cout << " - Min Y Position   : " << m_min_y << " m" << std::endl;
    std::cout << " - Final Y Pos      : " << m_final_y << " m" << std::endl;
    std::cout << " - Final Velocity   : " << m_final_vy << " m/s" << std::endl;
    std::cout << " - Has Contact?     : " << (m_had_contact ? "YES" : "NO") << std::endl;
    if (m_had_contact) {
        std::cout << " - First Contact @  : " << m_first_contact_time << " s" << std::endl;
    }
    std::cout << " - Output CSV path  : " << m_config.output_csv_path << std::endl;
    std::cout << "===========================================" << std::endl;
}

} // namespace models
} // namespace platform
