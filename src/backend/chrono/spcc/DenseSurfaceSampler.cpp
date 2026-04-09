#include "platform/backend/spcc/DenseSurfaceSampler.h"

#include <algorithm>
#include <cmath>

namespace platform {
namespace backend {
namespace spcc {

namespace {

struct PendingSample {
    DenseSurfaceSample sample;
    double rank = 0.0;
};

}  // namespace

std::vector<DenseSurfaceSample> DenseSurfaceSampler::BuildFromMesh(const chrono::ChTriangleMeshConnected& mesh,
                                                                   std::size_t max_samples,
                                                                   double surface_res) {
    const auto& vertices = mesh.GetCoordsVertices();
    const auto& faces = mesh.GetIndicesVertexes();

    std::vector<PendingSample> pending;
    pending.reserve(faces.size() * 4);

    const double res = std::max(1.0e-6, surface_res);
    for (std::size_t face_index = 0; face_index < faces.size(); ++face_index) {
        const auto& face = faces[face_index];
        const chrono::ChVector3d v0 = vertices[face.x()];
        const chrono::ChVector3d v1 = vertices[face.y()];
        const chrono::ChVector3d v2 = vertices[face.z()];

        chrono::ChVector3d normal = chrono::Vcross(v1 - v0, v2 - v0);
        const double double_area = normal.Length();
        if (!(double_area > 1.0e-16) || !std::isfinite(double_area)) {
            continue;
        }
        const double area = 0.5 * double_area;
        normal *= (1.0 / double_area);

        const int nominal_points = std::max(1, static_cast<int>(std::ceil(area / (res * res))));
        const int steps = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(nominal_points)))));

        std::vector<chrono::ChVector3d> face_points;
        face_points.reserve(static_cast<std::size_t>(steps * steps));
        for (int i = 0; i <= steps; ++i) {
            for (int j = 0; j <= steps - i; ++j) {
                const double a = static_cast<double>(i) / static_cast<double>(steps);
                const double b = static_cast<double>(j) / static_cast<double>(steps);
                const double c = 1.0 - a - b;
                face_points.push_back(a * v0 + b * v1 + c * v2);
            }
        }

        const double area_weight =
            face_points.empty() ? 0.0 : (area / static_cast<double>(face_points.size()));
        for (const auto& point : face_points) {
            DenseSurfaceSample sample;
            sample.source_face_index = face_index;
            sample.xi_slave_S = point;
            sample.normal_slave_S = normal;
            sample.area_weight = area_weight;

            PendingSample ranked;
            ranked.sample = sample;
            ranked.rank = area_weight;
            pending.push_back(ranked);
        }
    }

    if (max_samples > 0 && pending.size() > max_samples) {
        std::stable_sort(pending.begin(), pending.end(), [](const PendingSample& a, const PendingSample& b) {
            return a.rank > b.rank;
        });
        pending.resize(max_samples);
    }

    std::vector<DenseSurfaceSample> out;
    out.reserve(pending.size());
    for (const auto& ranked : pending) {
        out.push_back(ranked.sample);
    }
    return out;
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
