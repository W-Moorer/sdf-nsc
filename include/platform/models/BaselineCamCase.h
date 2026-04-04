#pragma once

#include "platform/models/CamCaseConfig.h"
#include "platform/backend/ChronoRigidSystemNSC.h"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace platform {
namespace models {

class BaselineCamCase {
public:
    BaselineCamCase(const CamCaseConfig& config);
    
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
    void SaveCSV(double t, double p_y, double v_y, double a_y, unsigned int nc);
    void TryCaptureSnapshots(double t, double p_y, double v_y, double a_y, unsigned int nc);
    void FinalizeSnapshotOutput();
    bool LoadTriangleMesh(const std::string& obj_path, TriangleMesh& out_mesh);
    void TryExportVTKFrame(double t);
    void FinalizeVTKSeriesOutput();

    CamCaseConfig m_config;
    std::unique_ptr<platform::backend::ChronoRigidSystemNSC> m_backend;
    
    std::string m_csv_buffer;
    std::vector<std::string> m_snapshot_frame_json;
    std::vector<double> m_snapshot_times;
    std::size_t m_next_snapshot_index = 0;
    TriangleMesh m_cam_mesh;
    TriangleMesh m_follower_mesh;
    std::vector<ExportedVTKFrame> m_exported_vtk_frames;
    std::string m_vtk_frame_prefix;
    std::size_t m_sim_step_index = 0;
    std::size_t m_vtk_export_index = 0;
    double m_min_y;
    double m_max_y;
    double m_final_y;
    bool m_had_contact;
};

} // namespace models
} // namespace platform
