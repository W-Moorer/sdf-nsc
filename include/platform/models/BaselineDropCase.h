#pragma once

#include "platform/models/DropCaseConfig.h"
#include "platform/backend/ChronoRigidSystemNSC.h"
#include <memory>
#include <string>

namespace platform {
namespace models {

class BaselineDropCase {
public:
    BaselineDropCase(const DropCaseConfig& config);
    ~BaselineDropCase() = default;

    // Run the case loop until total_time is reached
    void Run();

private:
    void SetupSystem();
    void SaveCSV(double t, double p_y, double v_y, unsigned int nc);

    DropCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;

    // Statistics
    double m_min_y;
    double m_final_y;
    double m_final_vy;
    bool m_had_contact;
    double m_first_contact_time;

    // CSV File pointer
    std::string m_csv_buffer; // Build up buffer for easy disk IO later
};

} // namespace models
} // namespace platform
