#pragma once

#include "platform/backend/ChronoRigidSystemNSC.h"
#include "platform/models/HeadOnSphereCaseConfig.h"

#include <memory>
#include <string>

namespace platform {
namespace models {

class BaselineHeadOnSphereCase {
  public:
    explicit BaselineHeadOnSphereCase(const HeadOnSphereCaseConfig& config);
    ~BaselineHeadOnSphereCase() = default;

    void Run();

  private:
    void SetupSystem();
    void SaveCSV(double t, double pos_a_x, double vel_a_x, double pos_b_x, double vel_b_x, unsigned int nc);

    HeadOnSphereCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;
    std::string m_csv_buffer;
};

}  // namespace models
}  // namespace platform
