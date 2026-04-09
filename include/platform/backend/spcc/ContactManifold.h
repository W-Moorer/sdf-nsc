#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <chrono/core/ChMatrix33.h>
#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/ContactTuning.h"

namespace platform {
namespace backend {
namespace spcc {

struct ContactManifoldHeader {
    ContactManifoldKind kind = ContactManifoldKind::CompactPoint;
    std::size_t manifold_id = 0;
    std::size_t source_contact_index = 0;
    int age = 0;
    bool matched = false;

    chrono::ChVector3d x_W;
    chrono::ChVector3d n_W;
    chrono::ChMatrix33<> P_W = chrono::ChMatrix33<>(1);

    chrono::ChVector3d x_master_M;
    double phi = 0.0;
    double phi_eff = 0.0;

    chrono::ChVector3d v_rel_W;
    double u_n_pred = 0.0;
    chrono::ChVector3d u_tau_pred_W;

    chrono::ChMatrix33<> hessian_W = chrono::ChMatrix33<>(0);
    double curvature_term = 0.0;
    double curvature_gate = 1.0;
    bool curvature_tangential_only = false;
    double curvature_term_abs_max = 0.0;
    double curvature_term_ratio_max = 0.0;
    double curvature_gap_floor = 0.0;

    double mu = 0.0;
    double weight = 1.0;
};

// Final solver-ready contact constraint emitted from manifold quadrature.
struct QuadratureContact {
    std::size_t manifold_id = 0;
    std::size_t source_contact_index = 0;
    int quadrature_id = 0;

    chrono::ChVector3d x_W;
    chrono::ChVector3d vpA_W;
    chrono::ChVector3d vpB_W;
    chrono::ChVector3d n_W;

    double phi = 0.0;
    double phi_eff = 0.0;
    double weight = 1.0;
    bool requery_geometry = false;

    chrono::ChMatrix33<> hessian_W = chrono::ChMatrix33<>(0);
    bool curvature_tangential_only = false;
    double mu = 0.0;
};

struct ImpactPairManifoldState {
    chrono::ChVector3d impact_normal_W;
    double toi = 0.0;
};

struct CompactPointManifoldState {
    int cluster_size = 1;
    chrono::ChVector3d deepest_x_W;
};

struct SlidingPatchManifoldState {
    chrono::ChVector3d tangent1_W;
    chrono::ChVector3d tangent2_W;
    chrono::ChVector3d coverage_center_W;
    double coverage_radius = 0.0;
    int cluster_size = 1;
};

struct ArcSlidingManifoldState {
    chrono::ChVector3d axis_W;
    chrono::ChVector3d radial_center_W;
    chrono::ChVector3d tangent_arc_W;

    double theta_center = 0.0;
    double theta_span = 0.0;
    double axial_center = 0.0;
    double axial_span = 0.0;

    double cylinder_radius = 0.0;
    int cluster_size = 1;
};

struct ContactManifold {
    ContactManifoldHeader header;
    std::vector<std::size_t> member_contact_indices;
    std::vector<QuadratureContact> seed_contacts;
    ImpactPairManifoldState impact_pair;
    CompactPointManifoldState compact_point;
    SlidingPatchManifoldState sliding_patch;
    ArcSlidingManifoldState arc_sliding;
};

class SDFField;
struct ActiveContactSample;
struct RigidBodyStateW;

class IContactManifoldBuilder {
  public:
    virtual ~IContactManifoldBuilder() = default;

    virtual void BuildManifolds(const RigidBodyStateW& master_pred,
                                const RigidBodyStateW& slave_pred,
                                const SDFField& sdf,
                                const std::vector<ActiveContactSample>& candidates,
                                double step_size,
                                std::vector<ContactManifold>& out_manifolds) = 0;
};

class IManifoldQuadrature {
  public:
    virtual ~IManifoldQuadrature() = default;

    virtual void Discretize(const ContactManifold& manifold,
                            std::vector<QuadratureContact>& out_contacts) const = 0;
};

std::unique_ptr<IContactManifoldBuilder> MakeContactManifoldBuilder(ContactManifoldKind kind);
std::unique_ptr<IManifoldQuadrature> MakeManifoldQuadrature(ContactManifoldKind kind,
                                                            bool enable_multi_contact,
                                                            int target_contacts,
                                                            double span_scale,
                                                            double min_half_span);

}  // namespace spcc
}  // namespace backend
}  // namespace platform
