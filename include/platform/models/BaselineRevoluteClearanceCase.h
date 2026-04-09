#pragma once

#include "platform/backend/ChronoRigidSystemNSC.h"
#include "platform/models/RevoluteClearanceCaseConfig.h"

#include <memory>
#include <string>

namespace platform {
namespace models {

class BaselineRevoluteClearanceCase {
  public:
    explicit BaselineRevoluteClearanceCase(const RevoluteClearanceCaseConfig& config);
    ~BaselineRevoluteClearanceCase() = default;

    void Run();

  private:
    void SetupSystem();
    double ComputeBody3SwingAngleX(const chrono::ChVector3d& body2_pos,
                                   const chrono::ChVector3d& body3_pos) const;
    double UnwrapBody3SwingAngle(double raw_angle);
    void SaveCSV(double t,
                 const chrono::ChVector3d& body2_pos,
                 const chrono::ChVector3d& body2_vel,
                 const chrono::ChVector3d& body3_pos,
                 const chrono::ChVector3d& body3_vel,
                 double body3_angle_x,
                 double body3_angvel_x,
                 unsigned int nc);

    RevoluteClearanceCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;
    std::string m_csv_buffer;
    bool m_has_prev_body3_angle = false;
    double m_prev_body3_raw_angle = 0.0;
    double m_body3_angle_offset = 0.0;
};

}  // namespace models
}  // namespace platform
