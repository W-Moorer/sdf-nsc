#include "platform/models/BaselineSimpleGearCase.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace platform {
namespace models {

BaselineSimpleGearCase::BaselineSimpleGearCase(const SimpleGearCaseConfig& config)
    : m_config(config),
      m_min_wrx(1e9),
      m_max_wrx(-1e9),
      m_final_wrx(0.0),
      m_ratio_sum(0.0),
      m_ratio_count(0) {
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineSimpleGearCase::SetupSystem() {
    m_backend->InitializeSimpleGearCase(
        m_config.gear1_mesh_path,
        m_config.gear2_mesh_path,
        m_config.gear1_ref_pos,
        m_config.gear2_ref_pos,
        m_config.joint1_pos,
        m_config.joint2_pos,
        m_config.mesh_scale,
        m_config.density,
        m_config.friction,
        m_config.restitution,
        m_config.gravity_y,
        m_config.motor_speed,
        m_config.gear1_mass,
        m_config.gear1_inertia_xx,
        m_config.gear2_mass,
        m_config.gear2_inertia_xx,
        m_config.sdf_type == 2,
        m_config.sdf_build,
        m_config.sample_tuning,
        m_config.contact_regime);
}

void BaselineSimpleGearCase::SaveCSV(double t, double wrx, double wry, double wrz) {
    auto fmt = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(9) << v;
        return ss.str();
    };

    // Match the RecurDyn Gear22.csv style: first 3 columns are time vectors, last 3 are Y values.
    m_csv_buffer += fmt(t) + "," +
                    fmt(t) + "," +
                    fmt(t) + "," +
                    fmt(wrx) + "," +
                    fmt(wry) + "," +
                    fmt(wrz) + "\n";
}

void BaselineSimpleGearCase::Run() {
    std::cout << "Loading Simple Gear Meshes..." << std::endl;
    try {
        SetupSystem();
    } catch (const std::exception& e) {
        std::cerr << "Exception during setup: " << e.what() << std::endl;
        return;
    }

    std::cout << "===========================================" << std::endl;
    std::cout << " Baseline Simple Gear Case (NSC)" << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << "s, End Time: " << m_config.total_time << "s" << std::endl;
    std::cout << " - Drive Speed: " << m_config.motor_speed << " rad/s" << std::endl;
    std::cout << "===========================================" << std::endl;

    m_csv_buffer = "X:Vel_RX-GEAR22-chrono(rad/s),X:Vel_RY-GEAR22-chrono(rad/s),X:Vel_RZ-GEAR22-chrono(rad/s),"
                   "Y:Vel_RX-GEAR22-chrono(rad/s),Y:Vel_RY-GEAR22-chrono(rad/s),Y:Vel_RZ-GEAR22-chrono(rad/s)\n";

    std::cout << std::fixed << std::setprecision(4);

    auto w1_0 = m_backend->GetGear1AngVelParent();
    auto w2_0 = m_backend->GetGear2AngVelParent();
    SaveCSV(m_backend->GetTime(), w2_0.x(), w2_0.y(), w2_0.z());

    while (m_backend->GetTime() < m_config.total_time) {
        m_backend->StepDynamics(m_config.step_size);

        double t = m_backend->GetTime();
        auto w1 = m_backend->GetGear1AngVelParent();
        auto w2 = m_backend->GetGear2AngVelParent();

        m_min_wrx = std::min(m_min_wrx, w2.x());
        m_max_wrx = std::max(m_max_wrx, w2.x());
        m_final_wrx = w2.x();

        if (std::abs(w1.x()) > 1e-8) {
            m_ratio_sum += (w2.x() / w1.x());
            m_ratio_count += 1;
        }

        if (std::abs(std::fmod(t, 0.1)) < (m_config.step_size * 0.5) ||
            t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << "Time: " << std::setw(6) << t
                      << " s | Gear22 Wrx: " << std::setw(9) << w2.x()
                      << " rad/s | Gear22 Wry: " << std::setw(9) << w2.y()
                      << " rad/s | Gear22 Wrz: " << std::setw(9) << w2.z()
                      << " rad/s | Contacts: " << m_backend->GetNumContacts()
                      << std::endl;
        }

        SaveCSV(t, w2.x(), w2.y(), w2.z());
    }

    std::ofstream out(m_config.output_csv_path);
    if (out.is_open()) {
        out << m_csv_buffer;
        out.close();
        std::cout << "\nResults saved to: " << m_config.output_csv_path << std::endl;
    }

    double avg_ratio = (m_ratio_count > 0) ? (m_ratio_sum / m_ratio_count) : 0.0;
    std::cout << "===========================================" << std::endl;
    std::cout << " SUMMARY:" << std::endl;
    std::cout << " - Min Gear22 Wrx   : " << m_min_wrx << " rad/s" << std::endl;
    std::cout << " - Max Gear22 Wrx   : " << m_max_wrx << " rad/s" << std::endl;
    std::cout << " - Final Gear22 Wrx : " << m_final_wrx << " rad/s" << std::endl;
    std::cout << " - Avg Ratio W2/W1x : " << avg_ratio << std::endl;
    std::cout << " - Output CSV path  : " << m_config.output_csv_path << std::endl;
    std::cout << "===========================================" << std::endl;
}

} // namespace models
} // namespace platform
