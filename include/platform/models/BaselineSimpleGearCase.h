#pragma once

#include "platform/models/SimpleGearCaseConfig.h"
#include "platform/backend/ChronoRigidSystemNSC.h"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace platform {
namespace models {

class BaselineSimpleGearCase {
public:
    BaselineSimpleGearCase(const SimpleGearCaseConfig& config);

    void Run();

private:
    struct TriangleMesh {
        std::vector<chrono::ChVector3d> vertices;
        std::vector<std::array<int, 3>> faces;
    };

    struct ExportedVTKFrame {
        std::string filename;
        double time = 0.0;
    };

    void SetupSystem();
    void SaveCSV(double t, double wrx, double wry, double wrz);
    bool LoadTriangleMesh(const std::string& obj_path, double scale, TriangleMesh& out_mesh);
    void TryExportVTKFrame(double t);
    void FinalizeVTKSeriesOutput();

    SimpleGearCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;

    std::string m_csv_buffer;
    TriangleMesh m_gear1_mesh;
    TriangleMesh m_gear2_mesh;
    std::vector<ExportedVTKFrame> m_exported_vtk_frames;
    std::string m_vtk_frame_prefix;
    std::size_t m_sim_step_index = 0;
    std::size_t m_vtk_export_index = 0;
    double m_min_wrx;
    double m_max_wrx;
    double m_final_wrx;
    double m_ratio_sum;
    int m_ratio_count;
};

} // namespace models
} // namespace platform
