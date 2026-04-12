#include "platform/models/BaselineSimpleGearCase.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace platform {
namespace models {

namespace {

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

std::string XmlEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 8);
    for (char ch : s) {
        switch (ch) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            case '\'':
                out += "&apos;";
                break;
            default:
                out += ch;
                break;
        }
    }
    return out;
}

std::string JsonNumber(double v) {
    std::ostringstream ss;
    ss << std::setprecision(12) << (std::isfinite(v) ? v : 0.0);
    return ss.str();
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

}  // namespace

BaselineSimpleGearCase::BaselineSimpleGearCase(const SimpleGearCaseConfig& config)
    : m_config(config),
      m_min_wrx(1e9),
      m_max_wrx(-1e9),
      m_final_wrx(0.0),
      m_ratio_sum(0.0),
      m_ratio_count(0) {
    m_backend = std::make_unique<platform::backend::ChronoRigidSystemNSC>();
}

void BaselineSimpleGearCase::SetupSystem() {
    if (!m_config.vtk_output_dir.empty()) {
        m_vtk_frame_prefix = BuildFramePrefixFromMeshes(m_config.gear1_mesh_path, m_config.gear2_mesh_path);
        if (!LoadTriangleMesh(m_config.gear1_mesh_path, m_config.mesh_scale, m_gear1_mesh)) {
            std::cerr << "Failed to load gear 1 mesh for VTK export: " << m_config.gear1_mesh_path << std::endl;
            m_config.vtk_output_dir.clear();
        }
        if (!LoadTriangleMesh(m_config.gear2_mesh_path, m_config.mesh_scale, m_gear2_mesh)) {
            std::cerr << "Failed to load gear 2 mesh for VTK export: " << m_config.gear2_mesh_path << std::endl;
            m_config.vtk_output_dir.clear();
        }
        if (!m_config.vtk_output_dir.empty()) {
            std::filesystem::create_directories(m_config.vtk_output_dir);
        }
    }

    m_backend->InitializeSimpleGearCase(
        m_config.gear1_mesh_path,
        m_config.gear2_mesh_path,
        m_config.gear1_ref_pos,
        m_config.gear2_ref_pos,
        m_config.joint1_pos,
        m_config.joint2_pos,
        m_config.mesh_scale,
        m_config.density,
        m_config.friction,
        m_config.restitution,
        m_config.gravity_y,
        m_config.motor_speed,
        m_config.gear1_mass,
        m_config.gear1_inertia_xx,
        m_config.gear2_mass,
        m_config.gear2_inertia_xx,
        m_config.sdf_type == 2,
        m_config.sdf_build,
        m_config.sample_tuning,
        m_config.contact_regime);
}

void BaselineSimpleGearCase::SaveCSV(double t, double wrx, double wry, double wrz) {
    auto fmt = [](double v) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(9) << v;
        return ss.str();
    };

    // Match the RecurDyn Gear22.csv style: first 3 columns are time vectors, last 3 are Y values.
    m_csv_buffer += fmt(t) + "," +
                    fmt(t) + "," +
                    fmt(t) + "," +
                    fmt(wrx) + "," +
                    fmt(wry) + "," +
                    fmt(wrz) + "\n";
}

bool BaselineSimpleGearCase::LoadTriangleMesh(const std::string& obj_path,
                                              double scale,
                                              TriangleMesh& out_mesh) {
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
                out_mesh.vertices.emplace_back(scale * x, scale * y, scale * z);
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

void BaselineSimpleGearCase::TryExportVTKFrame(double t) {
    if (m_config.vtk_output_dir.empty()) {
        return;
    }

    if ((m_sim_step_index % static_cast<std::size_t>(std::max(1, m_config.vtk_stride))) != 0) {
        return;
    }

    backend::BodyDebugSnapshot gear1_body;
    backend::BodyDebugSnapshot gear2_body;
    if (!m_backend->GetGear1BodySnapshot(gear1_body) || !m_backend->GetGear2BodySnapshot(gear2_body)) {
        return;
    }

    std::vector<chrono::ChVector3d> world_points;
    world_points.reserve(m_gear1_mesh.vertices.size() + m_gear2_mesh.vertices.size());

    for (const auto& v : m_gear1_mesh.vertices) {
        world_points.push_back(TransformRefVertex(gear1_body, v));
    }
    const int gear2_vertex_offset = static_cast<int>(world_points.size());
    for (const auto& v : m_gear2_mesh.vertices) {
        world_points.push_back(TransformRefVertex(gear2_body, v));
    }

    const std::size_t num_polys = m_gear1_mesh.faces.size() + m_gear2_mesh.faces.size();
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
    out << "simple gear frame t=" << std::setprecision(12) << t << "\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << world_points.size() << " float\n";
    out << std::fixed << std::setprecision(9);
    for (const auto& p : world_points) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    out << "POLYGONS " << num_polys << " " << polygon_int_count << "\n";
    for (const auto& f : m_gear1_mesh.faces) {
        out << "3 " << f[0] << " " << f[1] << " " << f[2] << "\n";
    }
    for (const auto& f : m_gear2_mesh.faces) {
        out << "3 " << (gear2_vertex_offset + f[0]) << " " << (gear2_vertex_offset + f[1]) << " "
            << (gear2_vertex_offset + f[2]) << "\n";
    }

    out << "CELL_DATA " << num_polys << "\n";
    out << "SCALARS part_id int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < m_gear1_mesh.faces.size(); ++i) {
        out << "0\n";
    }
    for (std::size_t i = 0; i < m_gear2_mesh.faces.size(); ++i) {
        out << "1\n";
    }

    std::cout << "[VTK] Exported frame " << filename.str() << " at t=" << t << " s" << std::endl;
    m_exported_vtk_frames.push_back(ExportedVTKFrame{filename.str(), t});
    ++m_vtk_export_index;
}

void BaselineSimpleGearCase::FinalizeVTKSeriesOutput() {
    if (m_config.vtk_output_dir.empty() || m_exported_vtk_frames.empty()) {
        return;
    }

    const auto output_dir = std::filesystem::path(m_config.vtk_output_dir);
    const auto series_path = output_dir / (m_vtk_frame_prefix + ".vtk.series");
    std::ofstream series_out(series_path);
    if (!series_out.is_open()) {
        std::cerr << "Failed to open VTK series output path: " << series_path.string() << std::endl;
        return;
    }

    series_out << "{\n";
    series_out << "  \"file-series-version\": \"1.0\",\n";
    series_out << "  \"files\": [\n";
    for (std::size_t i = 0; i < m_exported_vtk_frames.size(); ++i) {
        const auto& frame = m_exported_vtk_frames[i];
        series_out << "    {\"name\": \"" << JsonEscape(frame.filename) << "\", \"time\": "
                   << JsonNumber(frame.time) << "}";
        if (i + 1 < m_exported_vtk_frames.size()) {
            series_out << ",";
        }
        series_out << "\n";
    }
    series_out << "  ]\n";
    series_out << "}\n";

    std::cout << "[VTK] Series manifest saved to: " << series_path.string() << std::endl;

    const auto pvd_path = output_dir / (m_vtk_frame_prefix + ".pvd");
    std::ofstream pvd_out(pvd_path);
    if (!pvd_out.is_open()) {
        std::cerr << "Failed to open PVD output path: " << pvd_path.string() << std::endl;
        return;
    }

    pvd_out << "<?xml version=\"1.0\"?>\n";
    pvd_out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvd_out << "  <Collection>\n";
    for (const auto& frame : m_exported_vtk_frames) {
        pvd_out << "    <DataSet timestep=\"" << JsonNumber(frame.time)
                << "\" group=\"\" part=\"0\" file=\"" << XmlEscape(frame.filename) << "\"/>\n";
    }
    pvd_out << "  </Collection>\n";
    pvd_out << "</VTKFile>\n";

    std::cout << "[VTK] PVD manifest saved to: " << pvd_path.string() << std::endl;
}

void BaselineSimpleGearCase::Run() {
    std::cout << "Loading Simple Gear Meshes..." << std::endl;
    try {
        SetupSystem();
    } catch (const std::exception& e) {
        std::cerr << "Exception during setup: " << e.what() << std::endl;
        return;
    }

    std::cout << "===========================================" << std::endl;
    std::cout << " Baseline Simple Gear Case (NSC)" << std::endl;
    std::cout << " - Step Size: " << m_config.step_size << "s, End Time: " << m_config.total_time << "s" << std::endl;
    std::cout << " - Drive Speed: " << m_config.motor_speed << " rad/s" << std::endl;
    std::cout << "===========================================" << std::endl;

    m_csv_buffer = "X:Vel_RX-GEAR22-chrono(rad/s),X:Vel_RY-GEAR22-chrono(rad/s),X:Vel_RZ-GEAR22-chrono(rad/s),"
                   "Y:Vel_RX-GEAR22-chrono(rad/s),Y:Vel_RY-GEAR22-chrono(rad/s),Y:Vel_RZ-GEAR22-chrono(rad/s)\n";

    std::cout << std::fixed << std::setprecision(4);

    auto w1_0 = m_backend->GetGear1AngVelParent();
    auto w2_0 = m_backend->GetGear2AngVelParent();
    const double t0 = m_backend->GetTime();
    m_min_wrx = std::min(m_min_wrx, w2_0.x());
    m_max_wrx = std::max(m_max_wrx, w2_0.x());
    m_final_wrx = w2_0.x();
    SaveCSV(t0, w2_0.x(), w2_0.y(), w2_0.z());
    TryExportVTKFrame(t0);

    while (m_backend->GetTime() < m_config.total_time) {
        m_backend->StepDynamics(m_config.step_size);
        ++m_sim_step_index;

        double t = m_backend->GetTime();
        auto w1 = m_backend->GetGear1AngVelParent();
        auto w2 = m_backend->GetGear2AngVelParent();

        m_min_wrx = std::min(m_min_wrx, w2.x());
        m_max_wrx = std::max(m_max_wrx, w2.x());
        m_final_wrx = w2.x();

        if (std::abs(w1.x()) > 1e-8) {
            m_ratio_sum += (w2.x() / w1.x());
            m_ratio_count += 1;
        }

        if (std::abs(std::fmod(t, 0.1)) < (m_config.step_size * 0.5) ||
            t >= m_config.total_time - m_config.step_size * 0.5) {
            std::cout << "Time: " << std::setw(6) << t
                      << " s | Gear22 Wrx: " << std::setw(9) << w2.x()
                      << " rad/s | Gear22 Wry: " << std::setw(9) << w2.y()
                      << " rad/s | Gear22 Wrz: " << std::setw(9) << w2.z()
                      << " rad/s | Contacts: " << m_backend->GetNumContacts()
                      << std::endl;
        }

        SaveCSV(t, w2.x(), w2.y(), w2.z());
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

    FinalizeVTKSeriesOutput();

    double avg_ratio = (m_ratio_count > 0) ? (m_ratio_sum / m_ratio_count) : 0.0;
    std::cout << "===========================================" << std::endl;
    std::cout << " SUMMARY:" << std::endl;
    std::cout << " - Min Gear22 Wrx   : " << m_min_wrx << " rad/s" << std::endl;
    std::cout << " - Max Gear22 Wrx   : " << m_max_wrx << " rad/s" << std::endl;
    std::cout << " - Final Gear22 Wrx : " << m_final_wrx << " rad/s" << std::endl;
    std::cout << " - Avg Ratio W2/W1x : " << avg_ratio << std::endl;
    std::cout << " - Output CSV path  : " << m_config.output_csv_path << std::endl;
    std::cout << "===========================================" << std::endl;
}

} // namespace models
} // namespace platform
