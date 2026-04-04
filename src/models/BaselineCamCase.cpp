#include "platform/models/BaselineCamCase.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <cctype>

namespace platform {
namespace models {

namespace {

std::string JsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char ch : s) {
        switch (ch) {
            case '\\':
                out += "\\\\";
                break;
            case '"':
                out += "\\\"";
                break;
            case '\n':
                out += "\\n";
                break;
            case '\r':
                out += "\\r";
                break;
            case '\t':
                out += "\\t";
                break;
            default:
                out += ch;
                break;
        }
    }
    return out;
}

bool ParseObjIndexToken(const std::string& token, int& out_index_zero_based) {
    if (token.empty()) {
        return false;
    }
    const auto slash_pos = token.find('/');
    const std::string index_token = (slash_pos == std::string::npos) ? token : token.substr(0, slash_pos);
    if (index_token.empty()) {
        return false;
    }
    const int idx = std::stoi(index_token);
    if (idx <= 0) {
        return false;
    }
    out_index_zero_based = idx - 1;
    return true;
}

chrono::ChVector3d TransformRefVertex(const backend::BodyDebugSnapshot& body,
                                      const chrono::ChVector3d& x_ref_local) {
    return body.x_ref_W + body.q_WRef.Rotate(x_ref_local);
}

std::string SanitizeToken(const std::string& token) {
    std::string out;
    out.reserve(token.size());
    bool last_was_underscore = false;
    for (unsigned char ch : token) {
        if (std::isalnum(ch) != 0) {
            out.push_back(static_cast<char>(std::tolower(ch)));
            last_was_underscore = false;
        } else if (!last_was_underscore) {
            out.push_back('_');
            last_was_underscore = true;
        }
    }
    while (!out.empty() && out.front() == '_') {
        out.erase(out.begin());
    }
    while (!out.empty() && out.back() == '_') {
        out.pop_back();
    }
    return out.empty() ? std::string("frame") : out;
}

std::string BuildFramePrefixFromMeshes(const std::string& mesh_a_path, const std::string& mesh_b_path) {
    const std::string stem_a = SanitizeToken(std::filesystem::path(mesh_a_path).stem().string());
    const std::string stem_b = SanitizeToken(std::filesystem::path(mesh_b_path).stem().string());
    if (stem_a == stem_b) {
        return stem_a;
    }
    return stem_a + "_" + stem_b;
}

std::string JsonNumber(double v) {
    std::ostringstream ss;
    ss << std::setprecision(12) << (std::isfinite(v) ? v : 0.0);
    return ss.str();
}

std::string JsonVec3(const chrono::ChVector3d& v) {
    return "[" + JsonNumber(v.x()) + ", " + JsonNumber(v.y()) + ", " + JsonNumber(v.z()) + "]";
}

std::string JsonQuat(const chrono::ChQuaternion<>& q) {
    return "[" + JsonNumber(q.e0()) + ", " + JsonNumber(q.e1()) + ", " + JsonNumber(q.e2()) + ", " +
           JsonNumber(q.e3()) + "]";
}

std::string JsonMat33(const chrono::ChMatrix33<>& m) {
    std::ostringstream ss;
    ss << "[";
    for (int r = 0; r < 3; ++r) {
        if (r > 0) {
            ss << ", ";
        }
        ss << "[" << JsonNumber(m(r, 0)) << ", " << JsonNumber(m(r, 1)) << ", " << JsonNumber(m(r, 2)) << "]";
    }
    ss << "]";
    return ss.str();
}

std::string BodySnapshotToJson(const backend::BodyDebugSnapshot& body) {
    std::ostringstream ss;
    ss << "{"
       << "\"valid\": " << (body.valid ? "true" : "false") << ", "
       << "\"x_com_W\": " << JsonVec3(body.x_com_W) << ", "
       << "\"x_ref_W\": " << JsonVec3(body.x_ref_W) << ", "
       << "\"q_WL\": " << JsonQuat(body.q_WL) << ", "
       << "\"q_WRef\": " << JsonQuat(body.q_WRef) << ", "
       << "\"v_com_W\": " << JsonVec3(body.v_com_W) << ", "
       << "\"w_W\": " << JsonVec3(body.w_W)
       << "}";
    return ss.str();
}

#if defined(SPCC_ENABLE_VDB)
std::string ContactSnapshotToJson(const backend::spcc::ActiveContactSample& contact) {
    std::ostringstream ss;
    ss << "{"
       << "\"sample_id\": " << contact.sample_id << ", "
       << "\"active_age\": " << contact.active_age << ", "
       << "\"point_slave_world\": " << JsonVec3(contact.x_W) << ", "
       << "\"point_master_surface_world\": " << JsonVec3(contact.x_master_surface_W) << ", "
       << "\"point_master_local\": " << JsonVec3(contact.x_master_M) << ", "
       << "\"normal_world\": " << JsonVec3(contact.n_W) << ", "
       << "\"grad_master_local\": " << JsonVec3(contact.grad_M) << ", "
       << "\"u_pred_world\": " << JsonVec3(contact.u_pred_W) << ", "
       << "\"u_tau_pred_world\": " << JsonVec3(contact.u_tau_pred_W) << ", "
       << "\"v_rel_world\": " << JsonVec3(contact.v_rel_W) << ", "
       << "\"phi\": " << JsonNumber(contact.phi) << ", "
       << "\"phi_eff\": " << JsonNumber(contact.phi_eff) << ", "
       << "\"curvature_term\": " << JsonNumber(contact.curvature_term) << ", "
       << "\"u_n_pred\": " << JsonNumber(contact.u_n_pred) << ", "
       << "\"mu\": " << JsonNumber(contact.mu) << ", "
       << "\"hessian_world\": " << JsonMat33(contact.hessian_W)
       << "}";
    return ss.str();
}
#endif

}  // namespace

BaselineCamCase::BaselineCamCase(const CamCaseConfig& config)
    : m_config(config),
      m_min_y(1e9), m_max_y(-1e9), m_final_y(0),
      m_had_contact(false)
{
    m_snapshot_times = m_config.snapshot_times;
    std::sort(m_snapshot_times.begin(), m_snapshot_times.end());
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineCamCase::SetupSystem() {
    if (!m_config.vtk_output_dir.empty()) {
        m_vtk_frame_prefix = BuildFramePrefixFromMeshes(m_config.cam_mesh_path, m_config.follower_mesh_path);
        if (!LoadTriangleMesh(m_config.cam_mesh_path, m_cam_mesh)) {
            std::cerr << "Failed to load cam mesh for VTK export: " << m_config.cam_mesh_path << std::endl;
            m_config.vtk_output_dir.clear();
        }
        if (!LoadTriangleMesh(m_config.follower_mesh_path, m_follower_mesh)) {
            std::cerr << "Failed to load follower mesh for VTK export: " << m_config.follower_mesh_path << std::endl;
            m_config.vtk_output_dir.clear();
        }
        if (!m_config.vtk_output_dir.empty()) {
            std::filesystem::create_directories(m_config.vtk_output_dir);
        }
    }

    m_backend->InitializeCamCase(
        m_config.cam_mesh_path, m_config.follower_mesh_path,
        m_config.cam_init_pos, m_config.follower_init_pos,
        m_config.density, m_config.friction, m_config.restitution,
        m_config.gravity_y, m_config.motor_speed, m_config.contact_algorithm
    );
}

void BaselineCamCase::SaveCSV(double t, double p_y, double v_y, double a_y, unsigned int nc) {
    (void)nc;

    auto fmt = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(9) << v;
        return ss.str();
    };

    // Export rows in a reference-compatible structure to make field-by-field comparison unambiguous.
    m_csv_buffer += fmt(t) + "," +   // X:Pos_TX-Body2(m) -> time-like field in provided reference
                    fmt(t) + "," +   // X:Pos_TY-Body2(m)
                    fmt(t) + "," +   // X:Pos_TZ-Body2(m)
                    fmt(t) + "," +   // X:Vel_TX-Body2(m/s)
                    fmt(t) + "," +   // X:Vel_TY-Body2(m/s)
                    fmt(t) + "," +   // X:Vel_TZ-Body2(m/s)
                    fmt(t) + "," +   // X:Acc_TX-Body2(m/s^2)
                    fmt(t) + "," +   // X:Acc_TY-Body2(m/s^2)
                    fmt(t) + "," +   // X:Acc_TZ-Body2(m/s^2)
                    fmt(0.0) + "," + // Y:Pos_TX-Body2(m)
                    fmt(p_y) + "," + // Y:Pos_TY-Body2(m)
                    fmt(0.0) + "," + // Y:Pos_TZ-Body2(m)
                    fmt(0.0) + "," + // Y:Vel_TX-Body2(m/s)
                    fmt(v_y) + "," + // Y:Vel_TY-Body2(m/s)
                    fmt(0.0) + "," + // Y:Vel_TZ-Body2(m/s)
                    fmt(0.0) + "," + // Y:Acc_TX-Body2(m/s^2)
                    fmt(a_y) + "," + // Y:Acc_TY-Body2(m/s^2)
                    fmt(0.0) + "\n"; // Y:Acc_TZ-Body2(m/s^2)
}

void BaselineCamCase::TryCaptureSnapshots(double t, double p_y, double v_y, double a_y, unsigned int nc) {
    if (m_config.snapshot_output_path.empty() || m_snapshot_times.empty()) {
        return;
    }

    const double tol = std::max(1e-9, m_config.step_size * 0.51);
    while (m_next_snapshot_index < m_snapshot_times.size() &&
           t + tol >= m_snapshot_times[m_next_snapshot_index]) {
        const double requested_time = m_snapshot_times[m_next_snapshot_index];

        backend::BodyDebugSnapshot cam_body;
        backend::BodyDebugSnapshot follower_body;
        m_backend->GetCamBodySnapshot(cam_body);
        m_backend->GetFollowerBodySnapshot(follower_body);

        std::ostringstream ss;
        ss << "{"
           << "\"requested_time\": " << JsonNumber(requested_time) << ", "
           << "\"captured_time\": " << JsonNumber(t) << ", "
           << "\"follower_pos_y\": " << JsonNumber(p_y) << ", "
           << "\"follower_vel_y\": " << JsonNumber(v_y) << ", "
           << "\"follower_acc_y\": " << JsonNumber(a_y) << ", "
           << "\"system_num_contacts\": " << nc << ", "
           << "\"cam\": " << BodySnapshotToJson(cam_body) << ", "
           << "\"follower\": " << BodySnapshotToJson(follower_body) << ", "
           << "\"contacts\": [";

#if defined(SPCC_ENABLE_VDB)
        const auto& contacts = m_backend->GetActiveContacts();
        for (std::size_t i = 0; i < contacts.size(); ++i) {
            if (i > 0) {
                ss << ", ";
            }
            ss << ContactSnapshotToJson(contacts[i]);
        }
#endif
        ss << "]"
           << "}";

        m_snapshot_frame_json.push_back(ss.str());
        std::cout << "[SNAPSHOT] Captured frame near t=" << requested_time << " s at simulation time " << t
                  << " s" << std::endl;
        ++m_next_snapshot_index;
    }
}

void BaselineCamCase::FinalizeSnapshotOutput() {
    if (m_config.snapshot_output_path.empty()) {
        return;
    }

    const std::filesystem::path out_path(m_config.snapshot_output_path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open snapshot output path: " << m_config.snapshot_output_path << std::endl;
        return;
    }

    out << "{\n"
        << "  \"case\": \"cam\",\n"
        << "  \"contact_algorithm\": \""
        << platform::common::ContactAlgorithmToCliName(m_config.contact_algorithm) << "\",\n"
        << "  \"cam_mesh_path\": \"" << JsonEscape(m_config.cam_mesh_path) << "\",\n"
        << "  \"follower_mesh_path\": \"" << JsonEscape(m_config.follower_mesh_path) << "\",\n"
        << "  \"step_size\": " << JsonNumber(m_config.step_size) << ",\n"
        << "  \"total_time\": " << JsonNumber(m_config.total_time) << ",\n"
        << "  \"motor_speed\": " << JsonNumber(m_config.motor_speed) << ",\n"
        << "  \"view_plane\": \"xy\",\n"
        << "  \"frames\": [\n";

    for (std::size_t i = 0; i < m_snapshot_frame_json.size(); ++i) {
        out << "    " << m_snapshot_frame_json[i];
        if (i + 1 < m_snapshot_frame_json.size()) {
            out << ",";
        }
        out << "\n";
    }

    out << "  ]\n"
        << "}\n";

    std::cout << "Snapshot JSON saved to: " << m_config.snapshot_output_path << " ("
              << m_snapshot_frame_json.size() << " frames)" << std::endl;
}

bool BaselineCamCase::LoadTriangleMesh(const std::string& obj_path, TriangleMesh& out_mesh) {
    out_mesh = TriangleMesh{};

    std::ifstream in(obj_path);
    if (!in.is_open()) {
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2) {
            continue;
        }

        std::istringstream iss(line);
        std::string tag;
        iss >> tag;

        if (tag == "v") {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            if (iss >> x >> y >> z) {
                out_mesh.vertices.emplace_back(x, y, z);
            }
        } else if (tag == "f") {
            std::vector<int> face_indices;
            std::string token;
            while (iss >> token) {
                int idx = -1;
                if (ParseObjIndexToken(token, idx)) {
                    face_indices.push_back(idx);
                }
            }
            if (face_indices.size() >= 3) {
                for (std::size_t i = 1; i + 1 < face_indices.size(); ++i) {
                    out_mesh.faces.push_back({face_indices[0], face_indices[i], face_indices[i + 1]});
                }
            }
        }
    }

    return !out_mesh.vertices.empty() && !out_mesh.faces.empty();
}

void BaselineCamCase::TryExportVTKFrame(double t) {
    if (m_config.vtk_output_dir.empty()) {
        return;
    }

    if ((m_sim_step_index % static_cast<std::size_t>(std::max(1, m_config.vtk_stride))) != 0) {
        return;
    }

    backend::BodyDebugSnapshot cam_body;
    backend::BodyDebugSnapshot follower_body;
    if (!m_backend->GetCamBodySnapshot(cam_body) || !m_backend->GetFollowerBodySnapshot(follower_body)) {
        return;
    }

    std::vector<chrono::ChVector3d> world_points;
    world_points.reserve(m_cam_mesh.vertices.size() + m_follower_mesh.vertices.size());

    for (const auto& v : m_cam_mesh.vertices) {
        world_points.push_back(TransformRefVertex(cam_body, v));
    }
    const int follower_vertex_offset = static_cast<int>(world_points.size());
    for (const auto& v : m_follower_mesh.vertices) {
        world_points.push_back(TransformRefVertex(follower_body, v));
    }

    const std::size_t num_polys = m_cam_mesh.faces.size() + m_follower_mesh.faces.size();
    const std::size_t polygon_int_count = num_polys * 4;

    std::ostringstream filename;
    filename << m_vtk_frame_prefix << "_" << std::setw(4) << std::setfill('0') << m_vtk_export_index << ".vtk";
    const auto out_path = std::filesystem::path(m_config.vtk_output_dir) / filename.str();

    std::ofstream out(out_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open VTK output path: " << out_path.string() << std::endl;
        return;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "cam frame t=" << std::setprecision(12) << t << "\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << world_points.size() << " float\n";
    out << std::fixed << std::setprecision(9);
    for (const auto& p : world_points) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    out << "POLYGONS " << num_polys << " " << polygon_int_count << "\n";
    for (const auto& f : m_cam_mesh.faces) {
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }
    for (const auto& f : m_follower_mesh.faces) {
        out << "3 " << (follower_vertex_offset + f[0]) << " " << (follower_vertex_offset + f[1]) << " "
            << (follower_vertex_offset + f[2]) << "\n";
    }

    out << "CELL_DATA " << num_polys << "\n";
    out << "SCALARS part_id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < m_cam_mesh.faces.size(); ++i) {
        out << "0\n";
    }
    for (std::size_t i = 0; i < m_follower_mesh.faces.size(); ++i) {
        out << "1\n";
    }

    std::cout << "[VTK] Exported frame " << filename.str() << " at t=" << t << " s" << std::endl;
    m_exported_vtk_frames.push_back(ExportedVTKFrame{filename.str(), t});
    ++m_vtk_export_index;
}

void BaselineCamCase::FinalizeVTKSeriesOutput() {
    if (m_config.vtk_output_dir.empty() || m_exported_vtk_frames.empty()) {
        return;
    }

    const auto series_path =
        std::filesystem::path(m_config.vtk_output_dir) / (m_vtk_frame_prefix + ".vtk.series");
    std::ofstream out(series_path);
    if (!out.is_open()) {
        std::cerr << "Failed to open VTK series output path: " << series_path.string() << std::endl;
        return;
    }

    out << "{\n";
    out << "  \"file-series-version\": \"1.0\",\n";
    out << "  \"files\": [\n";
    for (std::size_t i = 0; i < m_exported_vtk_frames.size(); ++i) {
        const auto& frame = m_exported_vtk_frames[i];
        out << "    {\"name\": \"" << JsonEscape(frame.filename) << "\", \"time\": " << JsonNumber(frame.time) << "}";
        if (i + 1 < m_exported_vtk_frames.size()) {
            out << ",";
        }
        out << "\n";
    }
    out << "  ]\n";
    out << "}\n";

    std::cout << "[VTK] Series manifest saved to: " << series_path.string() << std::endl;
}

void BaselineCamCase::Run() {
    std::cout << "Loading Cam Meshes..." << std::endl;
    try {
        SetupSystem();
    } catch(const std::exception& e) {
        std::cerr << "Exception during setup: " << e.what() << std::endl;
        return;
    }

    std::cout << "===========================================" << std::endl;
    std::cout << " Baseline Cam Case (NSC)" << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << "s, End Time: " << m_config.total_time << "s" << std::endl;
    std::cout << " - Motor Speed: " << m_config.motor_speed << " rad/s" << std::endl;
    std::cout << " - Contact Algorithm: " << platform::common::ContactAlgorithmToCliName(m_config.contact_algorithm)
              << std::endl;
    std::cout << "===========================================" << std::endl;

    m_csv_buffer = "X:Pos_TX-Body2(m),X:Pos_TY-Body2(m),X:Pos_TZ-Body2(m),"
                   "X:Vel_TX-Body2(m/s),X:Vel_TY-Body2(m/s),X:Vel_TZ-Body2(m/s),"
                   "X:Acc_TX-Body2(m/s^2),X:Acc_TY-Body2(m/s^2),X:Acc_TZ-Body2(m/s^2),"
                   "Y:Pos_TX-Body2(m),Y:Pos_TY-Body2(m),Y:Pos_TZ-Body2(m),"
                   "Y:Vel_TX-Body2(m/s),Y:Vel_TY-Body2(m/s),Y:Vel_TZ-Body2(m/s),"
                   "Y:Acc_TX-Body2(m/s^2),Y:Acc_TY-Body2(m/s^2),Y:Acc_TZ-Body2(m/s^2)\n";

    std::cout << std::fixed << std::setprecision(4);

    double t0 = m_backend->GetTime();
    double p_y0 = m_backend->GetFollowerPosY();
    double v_y0 = m_backend->GetFollowerVelY();
    unsigned int nc0 = m_backend->GetNumContacts();
    double a_y0 = m_config.gravity_y;

    if (nc0 > 0) m_had_contact = true;
    m_min_y = std::min(m_min_y, p_y0);
    m_max_y = std::max(m_max_y, p_y0);
    m_final_y = p_y0;
    SaveCSV(t0, p_y0, v_y0, a_y0, nc0);
    TryCaptureSnapshots(t0, p_y0, v_y0, a_y0, nc0);
    TryExportVTKFrame(t0);

    double prev_v_y = v_y0;

    const double step_tol = std::max(1e-12, m_config.step_size * 1e-9);
    const std::size_t max_steps =
        static_cast<std::size_t>(std::floor((m_config.total_time + step_tol) / m_config.step_size));

    for (std::size_t step = 0; step < max_steps; ++step) {
        m_backend->StepDynamics(m_config.step_size);
        ++m_sim_step_index;
        
        double t = m_backend->GetTime();
        double p_y = m_backend->GetFollowerPosY();
        double v_y = m_backend->GetFollowerVelY();
        unsigned int nc = m_backend->GetNumContacts();
        double a_y = (v_y - prev_v_y) / m_config.step_size;
        prev_v_y = v_y;
        
        if (nc > 0) m_had_contact = true;
        
        m_min_y = std::min(m_min_y, p_y);
        m_max_y = std::max(m_max_y, p_y);
        m_final_y = p_y;
        
        // Print roughly every 0.1s
        if (std::abs(std::fmod(t, 0.1)) < (m_config.step_size * 0.5) || t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << "Time: " << std::setw(6) << t 
                      << " s | Pos Y: " << std::setw(8) << p_y 
                      << " m | Vel Y: " << std::setw(8) << v_y << " m/s" 
                      << " | Contacts: " << nc << std::endl;
        }

        SaveCSV(t, p_y, v_y, a_y, nc);
        TryCaptureSnapshots(t, p_y, v_y, a_y, nc);
        TryExportVTKFrame(t);
    }

    const std::filesystem::path csv_path(m_config.output_csv_path);
    if (csv_path.has_parent_path()) {
        std::filesystem::create_directories(csv_path.parent_path());
    }

    std::ofstream out(m_config.output_csv_path);
    if (out.is_open()) {
        out << m_csv_buffer;
        out.close();
        std::cout << "\nResults saved to: " << m_config.output_csv_path << std::endl;
    }

    FinalizeSnapshotOutput();
    FinalizeVTKSeriesOutput();

    std::cout << "===========================================" << std::endl;
    std::cout << " SUMMARY:" << std::endl;
    std::cout << " - Min Y Position   : " << m_min_y << " m" << std::endl;
    std::cout << " - Max Y Position   : " << m_max_y << " m" << std::endl;
    std::cout << " - Final Y Pos      : " << m_final_y << " m" << std::endl;
    std::cout << " - Has Contact?     : " << (m_had_contact ? "YES" : "NO") << std::endl;
    std::cout << " - Output CSV path  : " << m_config.output_csv_path << std::endl;
    std::cout << "===========================================" << std::endl;
}

} // namespace models
} // namespace platform
