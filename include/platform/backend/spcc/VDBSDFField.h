#pragma once

#include <memory>
#include <string>

#include <chrono/core/ChVector3.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>

#include "platform/backend/spcc/ContactActivation.h"

namespace platform {
namespace backend {
namespace spcc {

class VDBSDFField final : public SDFField {
  public:
    struct BuildOptions {
        double voxel_size = 5e-4;
        double half_band_width_voxels = 3.0;
        std::string grid_name = "master_sdf";
        bool direct_phi_hessian = false;
    };

    VDBSDFField();
    ~VDBSDFField() override;

    bool BuildFromTriangleMesh(const chrono::ChTriangleMeshConnected& mesh_M);

    bool BuildFromTriangleMesh(const chrono::ChTriangleMeshConnected& mesh_M,
                   const BuildOptions& options);

    // Optional path if a prebuilt .vdb asset already exists.
    bool LoadFromVDBFile(const std::string& vdb_path,
                         const std::string& grid_name = std::string());

    bool IsReady() const;
    const std::string& LastError() const;

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override;

    bool QueryPhiGradHessianM(const chrono::ChVector3d& x_M,
                              double& phi,
                              chrono::ChVector3d& grad_M,
                              chrono::ChMatrix33<>& hessian_M) const override;

  private:
    bool BuildNanoFromOpenGrid();

    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace spcc
}  // namespace backend
}  // namespace platform
