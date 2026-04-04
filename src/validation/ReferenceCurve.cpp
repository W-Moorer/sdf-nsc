#include "platform/validation/ReferenceCurve.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace platform {
namespace validation {

bool ReferenceCurve::LoadFromCSV(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[ReferenceCurve] Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;
    // Read header
    if (!std::getline(file, line)) return false;

    // Check columns
    bool has_contacts = (line.find("Num_Contacts") != std::string::npos);

    m_data.clear();
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        
        DataPoint pt;
        pt.num_contacts = -1;

        // Time
        if (!std::getline(ss, token, ',')) continue;
        pt.time = std::stod(token);

        // Pos_Y
        if (!std::getline(ss, token, ',')) continue;
        pt.pos_y = std::stod(token);

        // Vel_Y
        if (!std::getline(ss, token, ',')) continue;
        pt.vel_y = std::stod(token);

        // Num_Contacts
        if (has_contacts && std::getline(ss, token, ',')) {
            pt.num_contacts = std::stoi(token);
        }

        m_data.push_back(pt);
    }

    // Ensure sorted by time
    std::sort(m_data.begin(), m_data.end(), [](const DataPoint& a, const DataPoint& b) {
        return a.time < b.time;
    });

    return true;
}

bool ReferenceCurve::GetInterpolated(double time, DataPoint& out_point) const {
    if (m_data.empty()) return false;
    if (time < m_data.front().time || time > m_data.back().time) return false;

    auto it = std::lower_bound(m_data.begin(), m_data.end(), time, [](const DataPoint& pt, double t) {
        return pt.time < t;
    });

    if (it == m_data.end()) {
        out_point = m_data.back();
        return true;
    }
    if (it == m_data.begin()) {
        out_point = m_data.front();
        return true;
    }

    auto prev = it - 1;
    double t0 = prev->time;
    double t1 = it->time;
    double factor = (time - t0) / (t1 - t0);

    out_point.time = time;
    out_point.pos_y = prev->pos_y + factor * (it->pos_y - prev->pos_y);
    out_point.vel_y = prev->vel_y + factor * (it->vel_y - prev->vel_y);
    
    // Nearest integer for contacts
    out_point.num_contacts = (factor < 0.5) ? prev->num_contacts : it->num_contacts;

    return true;
}

} // namespace validation
} // namespace platform
