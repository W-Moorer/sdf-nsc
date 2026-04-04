#pragma once

#include "platform/models/SimpleGearCaseConfig.h"
#include "platform/backend/ChronoRigidSystemNSC.h"
#include <memory>
#include <string>

namespace platform {
namespace models {

class BaselineSimpleGearCase {
public:
    BaselineSimpleGearCase(const SimpleGearCaseConfig& config);

    void Run();

private:
    void SetupSystem();
    void SaveCSV(double t, double wrx, double wry, double wrz);

    SimpleGearCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;

    std::string m_csv_buffer;
    double m_min_wrx;
    double m_max_wrx;
    double m_final_wrx;
    double m_ratio_sum;
    int m_ratio_count;
};

} // namespace models
} // namespace platform
