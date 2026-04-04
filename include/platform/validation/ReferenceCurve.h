#pragma once

#include <string>
#include <vector>
#include <map>

namespace platform {
namespace validation {

struct DataPoint {
    double time;
    double pos_y;
    double vel_y;
    int num_contacts; // -1 means missing
};

class ReferenceCurve {
public:
    ReferenceCurve() = default;

    // Load from CSV. Expected columns: Time, Pos_Y, Vel_Y, [Num_Contacts]
    bool LoadFromCSV(const std::string& filepath);

    // Get interpolated data point at specific time
    // Returns false if time is out of bounds
    bool GetInterpolated(double time, DataPoint& out_point) const;

    const std::vector<DataPoint>& GetData() const { return m_data; }

private:
    std::vector<DataPoint> m_data;
};

} // namespace validation
} // namespace platform
