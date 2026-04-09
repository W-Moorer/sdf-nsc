#include "platform/models/BaselineRevoluteClearanceCase.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace platform {
namespace models {

namespace {

double WrapToPi(double angle) {
    while (angle > chrono::CH_PI) {
        angle -= chrono::CH_2PI;
    }
    while (angle < -chrono::CH_PI) {
        angle += chrono::CH_2PI;
    }
    return angle;
}

}  // namespace

BaselineRevoluteClearanceCase::BaselineRevoluteClearanceCase(const RevoluteClearanceCaseConfig& config)
    : m_config(config) {
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineRevoluteClearanceCase::SetupSystem() {
    m_backend->InitializeRevoluteClearanceCase(
        m_config.body1_mesh_path, m_config.body3_mesh_path, m_config.body1_init_pos, m_config.body3_init_pos,
        m_config.body2_cm_offset, m_config.body3_mass, m_config.body3_inertia_xx, m_config.body3_inertia_xy,
        m_config.body2_mass, m_config.body2_inertia_xx, m_config.body2_inertia_xy, m_config.friction,
        m_config.restitution, m_config.gravity_y, m_config.contact_compliance, m_config.contact_compliance_t,
        m_config.contact_damping_f, m_config.collision_envelope, m_config.dynamics_substeps, m_config.env_prefix,
        m_config.contact_algorithm, m_config.sdf_build, m_config.sample_tuning, m_config.contact_regime,
        m_config.use_lcp_manifold_quadrature, m_config.manifold_quadrature_contacts,
        m_config.manifold_quadrature_span_scale, m_config.manifold_quadrature_min_half_span);
}

double BaselineRevoluteClearanceCase::ComputeBody3SwingAngleX(const chrono::ChVector3d& body2_pos,
                                                              const chrono::ChVector3d& body3_pos) const {
    const chrono::ChVector3d rel = body2_pos - body3_pos;
    const double ref_angle = std::atan2(m_config.body2_cm_offset[2], m_config.body2_cm_offset[1]);
    const double cur_angle = std::atan2(rel.z(), rel.y());
    return WrapToPi(cur_angle - ref_angle);
}

double BaselineRevoluteClearanceCase::UnwrapBody3SwingAngle(double raw_angle) {
    if (!m_has_prev_body3_angle) {
        m_has_prev_body3_angle = true;
        m_prev_body3_raw_angle = raw_angle;
        m_body3_angle_offset = 0.0;
        return raw_angle;
    }

    const double delta = raw_angle - m_prev_body3_raw_angle;
    if (delta > chrono::CH_PI) {
        m_body3_angle_offset -= chrono::CH_2PI;
    } else if (delta < -chrono::CH_PI) {
        m_body3_angle_offset += chrono::CH_2PI;
    }
    m_prev_body3_raw_angle = raw_angle;
    return raw_angle + m_body3_angle_offset;
}

void BaselineRevoluteClearanceCase::SaveCSV(double t,
                                            const chrono::ChVector3d& body2_pos,
                                            const chrono::ChVector3d& body2_vel,
                                            const chrono::ChVector3d& body3_pos,
                                            const chrono::ChVector3d& body3_vel,
                                            double body3_angle_x,
                                            double body3_angvel_x,
                                            unsigned int nc) {
    auto fmt = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(9) << v;
        return ss.str();
    };

    m_csv_buffer += fmt(t) + "," + fmt(body2_pos.x()) + "," + fmt(body2_pos.y()) + "," + fmt(body2_pos.z()) + "," +
                    fmt(body2_vel.x()) + "," + fmt(body2_vel.y()) + "," + fmt(body2_vel.z()) + "," +
                    fmt(body3_pos.x()) + "," + fmt(body3_pos.y()) + "," + fmt(body3_pos.z()) + "," +
                    fmt(body3_vel.x()) + "," + fmt(body3_vel.y()) + "," + fmt(body3_vel.z()) + "," +
                    fmt(body3_angle_x) + "," + fmt(body3_angvel_x) + "," +
                    std::to_string(nc) + "\n";
}

void BaselineRevoluteClearanceCase::Run() {
    SetupSystem();

    std::cout << "===========================================" << std::endl;
    std::cout << " Revolute Clearance Joint Case (NSC)" << std::endl;
    std::cout << " - Contact Algorithm: "
              << platform::common::ContactAlgorithmToCliName(m_config.contact_algorithm) << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << " s, End Time: " << m_config.total_time << " s"
              << std::endl;
    std::cout << "===========================================" << std::endl;

    m_has_prev_body3_angle = false;
    m_prev_body3_raw_angle = 0.0;
    m_body3_angle_offset = 0.0;

    m_csv_buffer =
        "Time,Body2_Pos_TX,Body2_Pos_TY,Body2_Pos_TZ,Body2_Vel_TX,Body2_Vel_TY,Body2_Vel_TZ,"
        "Body3_Pos_TX,Body3_Pos_TY,Body3_Pos_TZ,Body3_Vel_TX,Body3_Vel_TY,Body3_Vel_TZ,"
        "Body3_Ang_X,Body3_AngVel_X,Num_Contacts\n";

    while (m_backend->GetTime() < m_config.total_time - 1.0e-12) {
        m_backend->StepDynamics(m_config.step_size);

        const double t = m_backend->GetTime();
        const chrono::ChVector3d body2_pos = m_backend->GetClearanceBody2Pos();
        const chrono::ChVector3d body2_vel = m_backend->GetClearanceBody2Vel();
        const chrono::ChVector3d body3_pos = m_backend->GetClearanceBody3Pos();
        const chrono::ChVector3d body3_vel = m_backend->GetClearanceBody3Vel();
        const chrono::ChVector3d body3_ang_vel = m_backend->GetClearanceBody3AngVel();
        const unsigned int nc = m_backend->GetNumContacts();
        const double body3_angle_x = UnwrapBody3SwingAngle(ComputeBody3SwingAngleX(body2_pos, body3_pos));

        if ((static_cast<int>(t * 1000.0 + 0.5) % 50) == 0 || t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << std::fixed << std::setprecision(6) << "Time: " << std::setw(8) << t
                      << " s | Body2(TY,TZ): (" << std::setw(10) << body2_pos.y() << ", " << std::setw(10)
                      << body2_pos.z() << ") | Body3(TY,TZ): (" << std::setw(10) << body3_pos.y() << ", "
                      << std::setw(10) << body3_pos.z() << ") | AngleX: " << std::setw(10) << body3_angle_x
                      << " | Contacts: " << nc << std::endl;
        }

        SaveCSV(t, body2_pos, body2_vel, body3_pos, body3_vel, body3_angle_x, body3_ang_vel.x(), nc);
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
}

}  // namespace models
}  // namespace platform
