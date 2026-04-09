#include "platform/backend/ChronoRigidSystemNSC.h"
#include "platform/models/CamCaseConfig.h"
#include <chrono/physics/ChContactMaterialNSC.h>
#include <chrono/physics/ChLinkMotorRotationSpeed.h>
#include <chrono/physics/ChLinkLock.h>
#include <chrono/physics/ChLinkTSDA.h>
#include <chrono/geometry/ChTriangleMeshConnected.h>
#include <chrono/collision/ChCollisionShapeSphere.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#if defined(SPCC_ENABLE_VDB)
#include "platform/backend/spcc/ContactActivation.h"
#include "platform/backend/spcc/ContactManifold.h"
#include "platform/backend/spcc/VDBSDFField.h"
#endif

namespace platform {
namespace backend {

#if defined(SPCC_ENABLE_VDB)
namespace {

double GetEnvDouble(const char* name, double fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || !std::isfinite(parsed)) {
        return fallback;
    }
    return parsed;
}

int GetEnvInt(const char* name, int fallback) {
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return fallback;
    }
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value) {
        return fallback;
    }
    return static_cast<int>(parsed);
}

std::string MakeScopedEnvName(const std::string& prefix, const char* suffix) {
    if (prefix.empty()) {
        return std::string(suffix);
    }
    return prefix + "_" + suffix;
}

double GetScopedEnvDouble(const std::string& prefix, const char* suffix, double fallback) {
    return GetEnvDouble(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

int GetScopedEnvInt(const std::string& prefix, const char* suffix, int fallback) {
    return GetEnvInt(MakeScopedEnvName(prefix, suffix).c_str(), fallback);
}

std::vector<chrono::ChVector3d> BuildDownsampledVertexSamples(const chrono::ChTriangleMeshConnected& mesh,
                                                               std::size_t max_samples,
                                                               double surface_res = 0.0005) {
    const auto& vertices = mesh.GetCoordsVertices();
    const auto& faces = mesh.GetIndicesVertexes();
    std::vector<chrono::ChVector3d> samples;
    if (vertices.empty()) {
        std::cout << "[DEBUG] BuildDownsampledVertexSamples returning 0 samples! Mesh is empty!" << std::endl;
        return samples;
    }

    // Add original vertices
    for (const auto& v : vertices) {
        if (!std::isfinite(v.x()) || !std::isfinite(v.y()) || !std::isfinite(v.z())) {
            continue;
        }
        samples.push_back(v);
    }
    
    // Add extra samples inside triangles to prevent tunneling/aliasing
    // We sample faces by barycentric coordinates to create a denser point cloud
    for (const auto& face : faces) {
        auto v0 = vertices[face.x()];
        auto v1 = vertices[face.y()];
        auto v2 = vertices[face.z()];
        
        // Edge lengths
        double l01 = (v1 - v0).Length();
        double l12 = (v2 - v1).Length();
        double l20 = (v0 - v2).Length();
        double area = chrono::Vcross(v1 - v0, v2 - v0).Length() * 0.5;
        
        // Target face-interior sample spacing for contact manifold coverage.
        const double res = std::max(1.0e-5, surface_res);
        int target_pts = std::min(500, std::max(0, (int)(area / (res * res * 0.5))));
        
        if (target_pts > 0) {
            // Systematic or random barycentric sampling. Using a simple grid approach on u,v
            int steps = std::max(2, (int)std::sqrt(target_pts * 2));
            for (int i = 1; i < steps; ++i) {
                for (int j = 1; j < steps - i; ++j) {
                    double u = (double)i / steps;
                    double v = (double)j / steps;
                    double w = 1.0 - u - v;
                    chrono::ChVector3d p = v0 * u + v1 * v + v2 * w;
                    samples.push_back(p);
                }
            }
        }
    }

    // Downsample if we still exceeded max_samples by a lot, but usually we want all of them.
    if (samples.size() > max_samples && max_samples > 0) {
        std::vector<chrono::ChVector3d> filtered;
        filtered.reserve(max_samples);
        // keep all original vertices
        std::size_t num_orig = vertices.size();
        for (std::size_t i = 0; i < num_orig; ++i) {
            filtered.push_back(samples[i]);
        }
        std::size_t remaining = max_samples > num_orig ? max_samples - num_orig : 0;
        if (remaining > 0) {
            std::size_t generated = samples.size() - num_orig;
            std::size_t stride = std::max((std::size_t)1, generated / remaining);
            for(std::size_t i = num_orig; i < samples.size(); i += stride) {
                filtered.push_back(samples[i]);
                if(filtered.size() >= max_samples) break;
            }
        }
        std::cout << "[DEBUG] BuildDownsampledVertexSamples returning " << filtered.size() << " samples. downsampled from " << samples.size() << std::endl;
        return filtered;
    }

    std::cout << "[DEBUG] BuildDownsampledVertexSamples returning " << samples.size() << " samples." << std::endl;
    return samples;
}

spcc::RigidBodyStateW MakeRigidBodyStateW(const std::shared_ptr<chrono::ChBody>& body,
                                          int body_id) {
    spcc::RigidBodyStateW state;
    state.body_id = body_id;
    if (!body) {
        return state;
    }

    state.x_com_W = body->GetPos();
    state.R_WL = body->GetRotMat();
    state.v_com_W = body->GetPosDt();
    state.w_W = body->GetAngVelParent();
    
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        state.x_ref_W = aux_ref->GetFrameRefToAbs().GetPos();
        state.R_WRef  = aux_ref->GetFrameRefToAbs().GetRotMat();
    } else {
        state.x_ref_W = state.x_com_W;
        state.R_WRef  = state.R_WL;
    }

    state.mass = body->GetMass();
    state.inv_mass = (state.mass > 0.0) ? (1.0 / state.mass) : 0.0;

    state.I_inv_L = body->GetInvInertia();
    state.I_inv_W = state.R_WL * state.I_inv_L * state.R_WL.transpose();
    return state;
}

BodyDebugSnapshot MakeBodyDebugSnapshot(const std::shared_ptr<chrono::ChBody>& body) {
    BodyDebugSnapshot out;
    if (!body) {
        return out;
    }

    out.valid = true;
    out.x_com_W = body->GetPos();
    out.q_WL = body->GetRot();
    out.v_com_W = body->GetPosDt();
    out.w_W = body->GetAngVelParent();

    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        const auto frame_ref = aux_ref->GetFrameRefToAbs();
        out.x_ref_W = frame_ref.GetPos();
        out.q_WRef = frame_ref.GetRot();
    } else {
        out.x_ref_W = out.x_com_W;
        out.q_WRef = out.q_WL;
    }

    return out;
}

std::vector<chrono::ChVector3d> BuildWorldSamplesFromLocal(const std::shared_ptr<chrono::ChBody>& slave_body,
                                                            const std::vector<chrono::ChVector3d>& local_samples_S) {
    std::vector<chrono::ChVector3d> world_samples_W;
    if (!slave_body || local_samples_S.empty()) {
        return world_samples_W;
    }

    world_samples_W.reserve(local_samples_S.size());
    // Crucial fix: ChBodyEasyMesh places the actual mesh at the REF frame, not COG frame.
    // If the body has a REF frame offset, we must apply it.
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(slave_body);
    chrono::ChFrame<> frame_W = aux_ref ? aux_ref->GetFrameRefToAbs() : slave_body->GetFrameCOMToAbs();

    for (const auto& xi_S : local_samples_S) {
        world_samples_W.push_back(frame_W.TransformPointLocalToParent(xi_S));
    }
    return world_samples_W;
}

std::shared_ptr<chrono::ChTriangleMeshConnected> BuildUvSphereMesh(double radius, int slices, int stacks) {
    auto mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
    auto& vertices = mesh->GetCoordsVertices();
    auto& faces = mesh->GetIndicesVertexes();

    slices = std::max(8, slices);
    stacks = std::max(4, stacks);
    vertices.reserve(static_cast<std::size_t>((stacks + 1) * slices));
    faces.reserve(static_cast<std::size_t>(2 * slices * (stacks - 1)));

    for (int stack = 0; stack <= stacks; ++stack) {
        const double v = static_cast<double>(stack) / static_cast<double>(stacks);
        const double phi = chrono::CH_PI * v;
        const double sin_phi = std::sin(phi);
        const double cos_phi = std::cos(phi);
        for (int slice = 0; slice < slices; ++slice) {
            const double u = static_cast<double>(slice) / static_cast<double>(slices);
            const double theta = chrono::CH_2PI * u;
            vertices.emplace_back(radius * sin_phi * std::cos(theta), radius * cos_phi,
                                  radius * sin_phi * std::sin(theta));
        }
    }

    auto vid = [slices](int stack, int slice) { return stack * slices + (slice % slices); };
    for (int stack = 0; stack < stacks; ++stack) {
        for (int slice = 0; slice < slices; ++slice) {
            const int a = vid(stack, slice);
            const int b = vid(stack, slice + 1);
            const int c = vid(stack + 1, slice);
            const int d = vid(stack + 1, slice + 1);
            if (stack != 0) {
                faces.emplace_back(a, c, b);
            }
            if (stack != stacks - 1) {
                faces.emplace_back(b, c, d);
            }
        }
    }

    return mesh;
}

void SetBodyFrameOrPos(const std::shared_ptr<chrono::ChBody>& body, const chrono::ChVector3d& pos_W) {
    if (!body) {
        return;
    }
    auto aux_ref = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body);
    if (aux_ref) {
        aux_ref->SetFrameRefToAbs(chrono::ChFrame<>(pos_W, chrono::QUNIT));
    } else {
        body->SetPos(pos_W);
    }
}

class SDFCollisionCallback : public chrono::ChSystem::CustomCollisionCallback {
public:
    std::shared_ptr<chrono::ChBody> m_master;
    std::shared_ptr<chrono::ChBody> m_slave;
    spcc::VDBSDFField* m_sdf;
    spcc::ContactActivation* m_activation;
    std::vector<chrono::ChVector3d> m_slave_vertices_local;
    std::shared_ptr<const spcc::SampleBVH> m_sample_bvh;
    std::shared_ptr<chrono::ChContactMaterial> m_material;
    std::vector<spcc::ActiveContactSample>* m_active_contacts_out;
    bool m_enable_curvature_term;
    bool m_project_master_point_from_phi;
    bool m_flip_contact_normal;
    bool m_use_predictive_gap;
    bool m_use_analytic_sphere_toi_contact;
    double m_master_sphere_radius;
    double m_slave_sphere_radius;
    spcc::ContactManifoldKind m_manifold_kind;
    std::unique_ptr<spcc::IContactManifoldBuilder> m_manifold_builder;
    std::unique_ptr<spcc::IManifoldQuadrature> m_manifold_quadrature_builder;
    bool m_use_manifold_quadrature;
    int m_manifold_quadrature_contacts;
    double m_manifold_quadrature_span_scale;
    double m_manifold_quadrature_min_half_span;

    SDFCollisionCallback(std::shared_ptr<chrono::ChBody> master,
                         std::shared_ptr<chrono::ChBody> slave,
                         spcc::VDBSDFField* sdf,
                         spcc::ContactActivation* activation,
                         const std::vector<chrono::ChVector3d>& slave_vertices,
                         std::shared_ptr<const spcc::SampleBVH> sample_bvh,
                         std::shared_ptr<chrono::ChContactMaterial> mat,
                         std::vector<spcc::ActiveContactSample>* active_contacts_out,
                         bool enable_curvature_term,
                         bool project_master_point_from_phi = false,
                         bool flip_contact_normal = false,
                         bool use_predictive_gap = false,
                         bool use_analytic_sphere_toi_contact = false,
                         double master_sphere_radius = 0.0,
                         double slave_sphere_radius = 0.0,
                         spcc::ContactManifoldKind manifold_kind = spcc::ContactManifoldKind::CompactPoint,
                         bool use_manifold_quadrature = false,
                         int manifold_quadrature_contacts = 1,
                         double manifold_quadrature_span_scale = 1.0,
                         double manifold_quadrature_min_half_span = 0.0)
        : m_master(master), m_slave(slave), m_sdf(sdf), m_activation(activation),
          m_slave_vertices_local(slave_vertices), m_sample_bvh(std::move(sample_bvh)), m_material(mat),
          m_active_contacts_out(active_contacts_out),
          m_enable_curvature_term(enable_curvature_term),
          m_project_master_point_from_phi(project_master_point_from_phi),
          m_flip_contact_normal(flip_contact_normal),
          m_use_predictive_gap(use_predictive_gap),
          m_use_analytic_sphere_toi_contact(use_analytic_sphere_toi_contact),
          m_master_sphere_radius(master_sphere_radius),
          m_slave_sphere_radius(slave_sphere_radius),
          m_manifold_kind(manifold_kind),
          m_manifold_builder(spcc::MakeContactManifoldBuilder(manifold_kind)),
          m_manifold_quadrature_builder(spcc::MakeManifoldQuadrature(manifold_kind,
                                                                    use_manifold_quadrature,
                                                                    manifold_quadrature_contacts,
                                                                    manifold_quadrature_span_scale,
                                                                    manifold_quadrature_min_half_span)),
          m_use_manifold_quadrature(use_manifold_quadrature),
          m_manifold_quadrature_contacts(std::max(1, manifold_quadrature_contacts)),
          m_manifold_quadrature_span_scale(manifold_quadrature_span_scale),
          m_manifold_quadrature_min_half_span(manifold_quadrature_min_half_span) {}

    void BuildSolverContacts(const spcc::RigidBodyStateW& master_state,
                             const spcc::RigidBodyStateW& slave_state,
                             double step_size,
                             const std::vector<spcc::ActiveContactSample>& active_contacts,
                             std::vector<spcc::ContactManifold>& manifolds,
                             std::vector<spcc::QuadratureContact>& quadrature_contacts) {
        manifolds.clear();
        quadrature_contacts.clear();
        if (!m_manifold_builder || !m_manifold_quadrature_builder || active_contacts.empty()) {
            return;
        }

        m_manifold_builder->BuildManifolds(master_state, slave_state, *m_sdf, active_contacts, step_size, manifolds);
        for (const auto& manifold : manifolds) {
            m_manifold_quadrature_builder->Discretize(manifold, quadrature_contacts);
        }
        std::stable_sort(quadrature_contacts.begin(), quadrature_contacts.end(),
                         [](const spcc::QuadratureContact& a, const spcc::QuadratureContact& b) {
                             if (a.source_contact_index != b.source_contact_index) {
                                 return a.source_contact_index < b.source_contact_index;
                             }
                             return a.quadrature_id < b.quadrature_id;
                         });
    }

    void AddAnalyticSphereTOIContact(chrono::ChSystem* sys,
                                     const spcc::RigidBodyStateW& master_state,
                                     const spcc::RigidBodyStateW& slave_state) {
        if (m_active_contacts_out) {
            m_active_contacts_out->clear();
        }

        const double radius_sum = m_master_sphere_radius + m_slave_sphere_radius;
        if (!(radius_sum > 0.0) || !std::isfinite(radius_sum)) {
            return;
        }

        chrono::ChVector3d delta_W = slave_state.x_com_W - master_state.x_com_W;
        double center_dist = delta_W.Length();
        chrono::ChVector3d n_W;
        if (center_dist > 1.0e-12 && std::isfinite(center_dist)) {
            n_W = delta_W * (1.0 / center_dist);
        } else {
            n_W = slave_state.v_com_W - master_state.v_com_W;
            const double n_len = n_W.Length();
            if (!(n_len > 1.0e-12) || !std::isfinite(n_len)) {
                n_W = chrono::ChVector3d(1.0, 0.0, 0.0);
            } else {
                n_W *= (1.0 / n_len);
            }
            center_dist = 0.0;
        }

        const chrono::ChVector3d v_rel = slave_state.v_com_W - master_state.v_com_W;
        const double vn = chrono::Vdot(n_W, v_rel);
        const double gap = center_dist - radius_sum;
        const double dt = sys->GetStep();

        const bool overlapping = gap <= 0.0;
        const bool crosses_this_step = (gap > 0.0) && (vn < -1.0e-12) && ((gap + dt * vn) <= 0.0);
        if ((!overlapping && !crosses_this_step) || !(vn < -1.0e-12)) {
            return;
        }

        double tau = 0.0;
        if (crosses_this_step) {
            tau = std::clamp(gap / (-vn), 0.0, dt);
        }

        const chrono::ChVector3d cA_toi_W = master_state.x_com_W + master_state.v_com_W * tau;
        const chrono::ChVector3d cB_toi_W = slave_state.x_com_W + slave_state.v_com_W * tau;
        const double restitution = static_cast<double>(
            std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetRestitution());
        const double denom = master_state.inv_mass + slave_state.inv_mass;
        if (!(denom > 0.0) || !std::isfinite(denom)) {
            return;
        }

        const double j = -(1.0 + restitution) * vn / denom;
        const chrono::ChVector3d vA_after_W = master_state.v_com_W - j * master_state.inv_mass * n_W;
        const chrono::ChVector3d vB_after_W = slave_state.v_com_W + j * slave_state.inv_mass * n_W;
        const double remain = std::max(0.0, dt - tau);
        const chrono::ChVector3d xA_after_W = cA_toi_W + vA_after_W * remain;
        const chrono::ChVector3d xB_after_W = cB_toi_W + vB_after_W * remain;

        SetBodyFrameOrPos(m_master, xA_after_W);
        SetBodyFrameOrPos(m_slave, xB_after_W);
        m_master->SetPosDt(vA_after_W);
        m_slave->SetPosDt(vB_after_W);
        m_master->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_slave->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_master->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
        m_slave->SetAngVelParent(chrono::ChVector3d(0, 0, 0));

        const chrono::ChVector3d vpA_W = cA_toi_W + m_master_sphere_radius * n_W;
        const chrono::ChVector3d vpB_W = cB_toi_W - m_slave_sphere_radius * n_W;
        const chrono::ChVector3d contact_normal_W = m_flip_contact_normal ? (-n_W) : n_W;

        if (m_active_contacts_out) {
            spcc::ActiveContactSample sample;
            sample.sample_id = 0;
            sample.x_W = 0.5 * (vpA_W + vpB_W);
            sample.x_master_surface_W = vpA_W;
            sample.phi = 0.0;
            sample.phi_eff = 0.0;
            sample.n_W = contact_normal_W;
            sample.v_rel_W = vB_after_W - vA_after_W;
            sample.mu = static_cast<double>(
                std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetSlidingFriction());
            sample.rA_W = sample.x_W - xA_after_W;
            sample.rB_W = sample.x_W - xB_after_W;
            sample.u_pred_W = sample.v_rel_W;
            sample.u_n_pred = chrono::Vdot(contact_normal_W, sample.u_pred_W);
            sample.active_age = 1;
            sample.cluster_size = 1;
            m_active_contacts_out->push_back(sample);
        }
    }

    void OnCustomCollision(chrono::ChSystem* sys) override {
        if (!m_master || !m_slave || !m_sdf || !m_activation || m_slave_vertices_local.empty()) {
            if (m_active_contacts_out) {
                m_active_contacts_out->clear();
            }
            return;
        }

        // 1. Build kinematic states
        auto master_state = MakeRigidBodyStateW(m_master, 1);
        auto slave_state = MakeRigidBodyStateW(m_slave, 2);

        if (m_use_analytic_sphere_toi_contact) {
            AddAnalyticSphereTOIContact(sys, master_state, slave_state);
            return;
        }

        // 2. Build active set with ContactActivation
        std::vector<spcc::ActiveContactSample> active_contacts;
        int penetration_count = 0;
        const auto friction =
            (float)std::static_pointer_cast<chrono::ChContactMaterialNSC>(m_material)->GetSlidingFriction();
        if (m_sample_bvh && !m_sample_bvh->Empty()) {
            m_activation->BuildActiveSetWithBVH(master_state, slave_state, *m_sdf, m_slave_vertices_local,
                                                *m_sample_bvh, friction, sys->GetStep(), m_enable_curvature_term,
                                                active_contacts);
        } else {
            auto world_samples_W = BuildWorldSamplesFromLocal(m_slave, m_slave_vertices_local);
            m_activation->BuildActiveSet(master_state, slave_state, *m_sdf, m_slave_vertices_local, world_samples_W,
                                         friction, sys->GetStep(), m_enable_curvature_term, active_contacts);
        }
        auto emit_contact = [&](spcc::ActiveContactSample& contact,
                                const chrono::ChVector3d& sample_W,
                                bool requery_geometry) {
            chrono::ChVector3d contact_normal_W = contact.n_W;
            double contact_phi = contact.phi;

            if (requery_geometry) {
                const chrono::ChVector3d x_master_M =
                    master_state.R_WRef.transpose() * (sample_W - master_state.x_ref_W);
                double phi_q = contact_phi;
                chrono::ChVector3d grad_q;
                if (!m_sdf->QueryPhiGradM(x_master_M, phi_q, grad_q)) {
                    return false;
                }
                chrono::ChVector3d n_q_W = master_state.R_WRef * grad_q;
                const double n_len = n_q_W.Length();
                if (!(n_len > 1.0e-12) || !std::isfinite(n_len)) {
                    return false;
                }
                n_q_W *= (1.0 / n_len);
                if (!std::isfinite(phi_q)) {
                    return false;
                }
                contact_phi = phi_q;
                contact_normal_W = n_q_W;
            }

            chrono::ChCollisionInfo cinfo;
            cinfo.modelA = m_master->GetCollisionModel().get();
            cinfo.modelB = m_slave->GetCollisionModel().get();
            cinfo.shapeA = nullptr;
            cinfo.shapeB = nullptr;
            cinfo.vN = m_flip_contact_normal ? (-contact_normal_W) : contact_normal_W;

            const chrono::ChVector3d rA_W = sample_W - master_state.x_com_W;
            const chrono::ChVector3d rB_W = sample_W - slave_state.x_com_W;
            const chrono::ChVector3d v_master =
                master_state.v_com_W + chrono::Vcross(master_state.w_W, rA_W);
            const chrono::ChVector3d v_slave =
                slave_state.v_com_W + chrono::Vcross(slave_state.w_W, rB_W);
            const chrono::ChVector3d v_rel = v_slave - v_master;

            chrono::ChMatrix33<> P_W = contact.P_W;
            if (requery_geometry) {
                P_W = chrono::ChMatrix33<>(1);
                P_W(0, 0) -= contact_normal_W.x() * contact_normal_W.x();
                P_W(0, 1) -= contact_normal_W.x() * contact_normal_W.y();
                P_W(0, 2) -= contact_normal_W.x() * contact_normal_W.z();
                P_W(1, 0) -= contact_normal_W.y() * contact_normal_W.x();
                P_W(1, 1) -= contact_normal_W.y() * contact_normal_W.y();
                P_W(1, 2) -= contact_normal_W.y() * contact_normal_W.z();
                P_W(2, 0) -= contact_normal_W.z() * contact_normal_W.x();
                P_W(2, 1) -= contact_normal_W.z() * contact_normal_W.y();
                P_W(2, 2) -= contact_normal_W.z() * contact_normal_W.z();
            }

            double curvature_term = 0.0;
            if (m_enable_curvature_term) {
                const double dt = sys->GetStep();
                chrono::ChMatrix33<> H_eff = contact.hessian_W;
                if (contact.curvature_tangential_only) {
                    H_eff = P_W * H_eff * P_W;
                }
                curvature_term = 0.5 * dt * dt * chrono::Vdot(v_rel, H_eff * v_rel);
                curvature_term *= contact.curvature_gate;
                if (!std::isfinite(curvature_term)) {
                    curvature_term = 0.0;
                }
                if (contact.curvature_term_abs_max > 0.0 && std::isfinite(contact.curvature_term_abs_max)) {
                    curvature_term =
                        std::clamp(curvature_term, -contact.curvature_term_abs_max, contact.curvature_term_abs_max);
                }
                if (contact.curvature_term_ratio_max > 0.0 && std::isfinite(contact.curvature_term_ratio_max)) {
                    const double gap_scale =
                        std::max(std::abs(contact_phi), std::max(0.0, contact.curvature_gap_floor));
                    const double ratio_cap = contact.curvature_term_ratio_max * gap_scale;
                    curvature_term = std::clamp(curvature_term, -ratio_cap, ratio_cap);
                }
            }

            cinfo.vpB = sample_W;
            cinfo.vpA = sample_W - contact_phi * contact_normal_W;

            double predicted_gap = contact_phi;
            if (m_use_predictive_gap) {
                const double dt = sys->GetStep();
                const double closing_speed = std::min(0.0, chrono::Vdot(contact_normal_W, v_rel));
                predicted_gap += dt * closing_speed;
            }
            cinfo.distance = predicted_gap + curvature_term;

            contact.x_master_surface_W = cinfo.vpA;
            contact.v_rel_W = v_rel;
            contact.curvature_term = curvature_term;
            contact.phi_eff = cinfo.distance;

            sys->GetContactContainer()->AddContact(cinfo, m_material, m_material);
            penetration_count++;
            return true;
        };

        std::vector<spcc::ContactManifold> manifolds;
        std::vector<spcc::QuadratureContact> quadrature_contacts;
        BuildSolverContacts(master_state, slave_state, sys->GetStep(), active_contacts, manifolds,
                            quadrature_contacts);

        for (const auto& quadrature_contact : quadrature_contacts) {
            if (quadrature_contact.source_contact_index >= active_contacts.size()) {
                continue;
            }
            auto& contact = active_contacts[quadrature_contact.source_contact_index];
            emit_contact(contact, quadrature_contact.x_W, quadrature_contact.requery_geometry);
        }

        if (m_active_contacts_out) {
            *m_active_contacts_out = active_contacts;
        }
        
        static int step_id = 0;
        if (step_id % 10 == 0) {
             std::cout << "[SDF] Step " << step_id << " active contacts: " << penetration_count;
#if defined(SPCC_ENABLE_VDB)
             if (GetEnvDouble("SPCC_DEBUG_CONTACT_STATS", 0.0) > 0.5) {
                 const auto& stats = m_activation->GetStats();
                 std::cout << " queried=" << stats.queried << " after_cap=" << stats.accepted_after_cap
                           << " fit_attempted=" << stats.local_fit_attempted
                           << " fit_applied=" << stats.local_fit_applied
                           << " fit_reject=" << stats.local_fit_rejected_positive_gap
                           << " bvh_nodes=" << stats.bvh_nodes_tested
                           << " bvh_pruned=" << stats.bvh_nodes_pruned
                           << " bvh_leaf=" << stats.bvh_leaf_nodes
                           << " bvh_marked=" << stats.bvh_samples_marked;
             }
             if (GetEnvDouble("SPCC_DEBUG_MANIFOLD_QUADRATURE", 0.0) > 0.5 && !active_contacts.empty()) {
                 double avg_span = 0.0;
                 double max_span = 0.0;
                 double avg_cluster = 0.0;
                 int quad_eligible = 0;
                 for (const auto& contact : active_contacts) {
                     avg_span += contact.coverage_half_span;
                     max_span = std::max(max_span, contact.coverage_half_span);
                     avg_cluster += static_cast<double>(contact.cluster_size);
                     if (contact.coverage_half_span >= m_manifold_quadrature_min_half_span &&
                         contact.coverage_tangent_W.Length2() > 1.0e-12) {
                         quad_eligible += 1;
                     }
                 }
                 avg_span /= static_cast<double>(active_contacts.size());
                 avg_cluster /= static_cast<double>(active_contacts.size());
                 std::cout << " quad_eligible=" << quad_eligible << "/" << active_contacts.size()
                           << " avg_span=" << avg_span << " max_span=" << max_span
                           << " avg_cluster=" << avg_cluster;
             }
#endif
             std::cout << std::endl;
        }
        step_id++;
    }
};

}  // namespace
#endif

ChronoRigidSystemNSC::ChronoRigidSystemNSC() {
    m_system = std::make_unique<chrono::ChSystemNSC>();
}

ChronoRigidSystemNSC::~ChronoRigidSystemNSC() = default;

void ChronoRigidSystemNSC::Initialize() {
    // Default call with baseline hardcoded params if Interface is strictly used
    Initialize(0.5, 1000.0, 5.0, 10.0, 1.0, 10.0, 0.5, 0.0, -9.81);
}

void ChronoRigidSystemNSC::Initialize(
    double ball_radius, double ball_density, double ball_height,
    double ground_x, double ground_y, double ground_z,
    double friction, double restitution, double gravity_y
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
    m_sample_bvh.reset();
#endif

    // 0. Initialize Collision System (Crucial since Chrono 9+)
    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);

    // 1. Setup Gravity
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    // 2. Setup Base Contact Material
    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution); 

    // 3. Fixed Ground Body
    auto ground = std::make_shared<chrono::ChBodyEasyBox>(ground_x, ground_y, ground_z, 1000.0, true, true, material);
    // Align ground top surface to Y=0
    ground->SetPos(chrono::ChVector3d(0, -ground_y / 2.0, 0)); 
    ground->SetFixed(true);
    m_system->AddBody(ground);

    // 4. Dynamic Sphere
    m_sphere = std::make_shared<chrono::ChBodyEasySphere>(ball_radius, ball_density, true, true, material);
    m_sphere->SetPos(chrono::ChVector3d(0, ball_height, 0));
    m_sphere->SetFixed(false);
    m_system->AddBody(m_sphere);
}

void ChronoRigidSystemNSC::StepDynamics(double step_size) {
    double headon_xA0 = 0.0;
    double headon_xB0 = 0.0;
    double headon_vA0 = 0.0;
    double headon_vB0 = 0.0;
    bool use_headon_exact_step = m_headon_manual_elastic_resolution && m_headon_a && m_headon_b;
    if (use_headon_exact_step) {
        headon_xA0 = m_headon_a->GetPos().x();
        headon_xB0 = m_headon_b->GetPos().x();
        headon_vA0 = m_headon_a->GetPosDt().x();
        headon_vB0 = m_headon_b->GetPosDt().x();
    }

    const int substeps = std::max(1, m_dynamics_substeps);
    const double sub_dt = step_size / static_cast<double>(substeps);
    for (int i = 0; i < substeps; ++i) {
        m_system->DoStepDynamics(sub_dt);
    }

    if (use_headon_exact_step) {
        const double mA = std::max(0.0, m_headon_a->GetMass());
        const double mB = std::max(0.0, m_headon_b->GetMass());
        const double radius_sum = 2.0 * m_headon_sphere_radius;
        const double gap0 = (headon_xB0 - headon_xA0) - radius_sum;
        const double closing_speed = headon_vA0 - headon_vB0;

        double xA1 = headon_xA0 + headon_vA0 * step_size;
        double xB1 = headon_xB0 + headon_vB0 * step_size;
        double vA1 = headon_vA0;
        double vB1 = headon_vB0;

        if (mA > 0.0 && mB > 0.0 && closing_speed > 1.0e-12) {
            const double gap1 = gap0 - closing_speed * step_size;
            if (gap0 <= 0.0 || gap1 <= 0.0) {
                double tau = 0.0;
                if (gap0 > 0.0) {
                    tau = std::clamp(gap0 / closing_speed, 0.0, step_size);
                }
                const double xA_toi = headon_xA0 + headon_vA0 * tau;
                const double xB_toi = headon_xB0 + headon_vB0 * tau;
                const double e = std::clamp(m_headon_restitution, 0.0, 1.0);
                vA1 = ((mA - e * mB) * headon_vA0 + (1.0 + e) * mB * headon_vB0) / (mA + mB);
                vB1 = ((mB - e * mA) * headon_vB0 + (1.0 + e) * mA * headon_vA0) / (mA + mB);
                const double remain = std::max(0.0, step_size - tau);
                xA1 = xA_toi + vA1 * remain;
                xB1 = xB_toi + vB1 * remain;
            }
        }

        auto posA = m_headon_a->GetPos();
        auto posB = m_headon_b->GetPos();
        posA.x() = xA1;
        posB.x() = xB1;
        m_headon_a->SetPos(posA);
        m_headon_b->SetPos(posB);

        auto velA = m_headon_a->GetPosDt();
        auto velB = m_headon_b->GetPosDt();
        velA.x() = vA1;
        velB.x() = vB1;
        velA.y() = 0.0;
        velA.z() = 0.0;
        velB.y() = 0.0;
        velB.z() = 0.0;
        m_headon_a->SetPosDt(velA);
        m_headon_b->SetPosDt(velB);
        m_headon_a->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_headon_b->SetPosDt2(chrono::ChVector3d(0, 0, 0));
        m_headon_a->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
        m_headon_b->SetAngVelParent(chrono::ChVector3d(0, 0, 0));
    }
}

double ChronoRigidSystemNSC::GetTime() const {
    return m_system->GetChTime();
}

unsigned int ChronoRigidSystemNSC::GetNumContacts() const {
    return m_system->GetNumContacts();
}

double ChronoRigidSystemNSC::GetDynamicSpherePosY() const {
    return m_sphere->GetPos().y();
}

double ChronoRigidSystemNSC::GetDynamicSphereVelY() const {
    return m_sphere->GetPosDt().y();
}

double ChronoRigidSystemNSC::GetHeadOnSphereAPosX() const {
    return m_headon_a ? m_headon_a->GetPos().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereAVelX() const {
    return m_headon_a ? m_headon_a->GetPosDt().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereBPosX() const {
    return m_headon_b ? m_headon_b->GetPos().x() : 0.0;
}

double ChronoRigidSystemNSC::GetHeadOnSphereBVelX() const {
    return m_headon_b ? m_headon_b->GetPosDt().x() : 0.0;
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody2Pos() const {
    return m_clearance_body2 ? m_clearance_body2->GetPos() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody2Vel() const {
    return m_clearance_body2 ? m_clearance_body2->GetPosDt() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3Pos() const {
    return m_clearance_body3 ? m_clearance_body3->GetPos() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3Vel() const {
    return m_clearance_body3 ? m_clearance_body3->GetPosDt() : chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetClearanceBody3AngVel() const {
    return m_clearance_body3 ? m_clearance_body3->GetAngVelParent() : chrono::ChVector3d(0, 0, 0);
}

void ChronoRigidSystemNSC::InitializeCamCase(
    const std::string& cam_obj, const std::string& follower_obj,
    const double* cam_pos, const double* follower_pos, const double* motor_joint_pos,
    double density, double friction, double restitution,
    double gravity_y, double motor_speed,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::ContactRegimeConfig& contact_regime,
    const platform::models::FollowerPreloadConfig& follower_preload
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
    m_sample_bvh.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution);

    // Ground anchor for joints
    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    const bool enable_curvature_term =
        (contact_algorithm == platform::common::ContactAlgorithm::SdfSecondOrder);
#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = true;
    (void)contact_algorithm;
#endif

#if defined(SPCC_ENABLE_VDB)
    std::shared_ptr<chrono::ChTriangleMeshConnected> cam_mesh;
    std::shared_ptr<chrono::ChTriangleMeshConnected> follower_mesh;
    if (use_vdb) {
        cam_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        cam_mesh->LoadWavefrontMesh(cam_obj, true, false);
        follower_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        follower_mesh->LoadWavefrontMesh(follower_obj, true, false);
    }
#endif

    // Load Cam Body
    std::shared_ptr<chrono::ChBodyEasyMesh> cam_body;
    if (use_vdb) {
        cam_body = std::make_shared<chrono::ChBodyEasyMesh>(
            cam_obj, density, true, true, false, material, 0.000001);
        auto dummy_sphere_c = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        cam_body->AddCollisionShape(dummy_sphere_c);
        cam_body->EnableCollision(true);
    } else {
        cam_body = std::make_shared<chrono::ChBodyEasyMesh>(
            cam_obj, density, true, true, true, material, 0.000001);
    }

    // Move frame to align with provided global starting pos
    cam_body->SetFrameRefToAbs(chrono::ChFrame<>(chrono::ChVector3d(cam_pos[0], cam_pos[1], cam_pos[2]), chrono::QUNIT));
    m_system->AddBody(cam_body);
    m_cam_body = cam_body;

    // Motor for Cam to Ground (Revolute around Z axis approx)
    auto motor = std::make_shared<chrono::ChLinkMotorRotationSpeed>();      
    // The revolute motor frame
    motor->Initialize(cam_body, ground,
                      chrono::ChFrame<>(chrono::ChVector3d(motor_joint_pos[0], motor_joint_pos[1], motor_joint_pos[2]),
                                        chrono::QuatFromAngleX(0)));
    auto mfun = std::make_shared<chrono::ChFunctionConst>(motor_speed);     
    motor->SetSpeedFunction(mfun);
    m_system->AddLink(motor);

    // Load Follower Body
    if (use_vdb) {
        m_follower = std::make_shared<chrono::ChBodyEasyMesh>(
            follower_obj, density, true, true, false, material, 0.000001);
        auto dummy_sphere_f = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
        m_follower->AddCollisionShape(dummy_sphere_f);
        m_follower->EnableCollision(true);
    } else {
        m_follower = std::make_shared<chrono::ChBodyEasyMesh>(
            follower_obj, density, true, true, true, material, 0.000001);       
    }

    auto follower_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_follower);
    if(follower_aux) {
        follower_aux->SetFrameRefToAbs(chrono::ChFrame<>(chrono::ChVector3d(follower_pos[0], follower_pos[1], follower_pos[2]), chrono::QUNIT));
    }
    m_system->AddBody(m_follower);

    // Prismatic joint: Follower to Ground (along Y axis)
    auto prismatic = std::make_shared<chrono::ChLinkLockPrismatic>();       
    // Initialize passing the center of the joint at follower's current position to allow sliding along Y
    prismatic->Initialize(m_follower, ground, chrono::ChFrame<>(m_follower->GetPos(), chrono::QuatFromAngleX(chrono::CH_PI_2)));
    m_system->AddLink(prismatic);

    if (follower_preload.enabled) {
        auto preload = std::make_shared<chrono::ChLinkTSDA>();
        preload->Initialize(m_follower, ground, false, m_follower->GetPos(),
                            chrono::ChVector3d(follower_preload.anchor_pos[0], follower_preload.anchor_pos[1],
                                               follower_preload.anchor_pos[2]));
        preload->SetRestLength(follower_preload.rest_length);
        preload->SetSpringCoefficient(follower_preload.stiffness);
        preload->SetDampingCoefficient(follower_preload.damping);
        preload->SetActuatorForce(0.0);
        m_system->AddLink(preload);
    }

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));
        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
        sdf_options.direct_phi_hessian =
            (GetScopedEnvDouble(env_prefix, "DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) >
             0.5);
        
        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>(); // reusing the gear1_sdf for Cam
        std::cout << "[DEBUG] Building VDB SDF for Cam with wide band..." << std::endl;
        m_gear1_sdf->BuildFromTriangleMesh(*cam_mesh, sdf_options);
        
        m_contact_activation.SetEnvPrefix(env_prefix);
        auto resolved_contact_regime = contact_regime;
        resolved_contact_regime.activation.delta_on =
            GetScopedEnvDouble(env_prefix, "DELTA_ON", resolved_contact_regime.activation.delta_on);
        resolved_contact_regime.activation.delta_off =
            GetScopedEnvDouble(env_prefix, "DELTA_OFF", resolved_contact_regime.activation.delta_off);
        resolved_contact_regime.activation.hold_steps =
            GetScopedEnvInt(env_prefix, "HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
        resolved_contact_regime.activation.max_active_keep =
            GetScopedEnvInt(env_prefix, "MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
        resolved_contact_regime.activation.use_sample_bvh =
            (GetScopedEnvDouble(env_prefix, "USE_SAMPLE_BVH",
                                resolved_contact_regime.activation.use_sample_bvh ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.activation.sample_bvh_leaf_size =
            GetScopedEnvInt(env_prefix, "SAMPLE_BVH_LEAF_SIZE",
                            resolved_contact_regime.activation.sample_bvh_leaf_size);
        resolved_contact_regime.activation.sample_bvh_margin_scale =
            GetScopedEnvDouble(env_prefix, "SAMPLE_BVH_MARGIN_SCALE",
                               resolved_contact_regime.activation.sample_bvh_margin_scale);
        resolved_contact_regime.curvature.enabled =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_ENABLED",
                                resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.tangential_only =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_TANGENTIAL_ONLY",
                                resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.normal_alignment_cos_min =
            GetScopedEnvDouble(env_prefix, "CURVATURE_ALIGNMENT_COS_MIN",
                               resolved_contact_regime.curvature.normal_alignment_cos_min);
        resolved_contact_regime.curvature.max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.max_hessian_frobenius);
        resolved_contact_regime.curvature.small_step_dt_threshold =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_DT_THRESHOLD",
                               resolved_contact_regime.curvature.small_step_dt_threshold);
        resolved_contact_regime.curvature.small_step_max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.small_step_max_hessian_frobenius);
        resolved_contact_regime.curvature.max_curvature_term_abs =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_ABS",
                               resolved_contact_regime.curvature.max_curvature_term_abs);
        resolved_contact_regime.curvature.max_curvature_term_ratio =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_RATIO",
                               resolved_contact_regime.curvature.max_curvature_term_ratio);
        resolved_contact_regime.curvature.gap_floor =
            GetScopedEnvDouble(env_prefix, "CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
        m_contact_activation.SetPolicy(resolved_contact_regime);

        const double cam_surface_res = GetScopedEnvDouble(env_prefix, "SURFACE_RES", sample_tuning.surface_res);
        const int cam_max_samples = GetScopedEnvInt(env_prefix, "MAX_SAMPLES", sample_tuning.max_samples);
        std::vector<chrono::ChVector3d> slave_samples =
            BuildDownsampledVertexSamples(*follower_mesh, static_cast<std::size_t>(std::max(0, cam_max_samples)),
                                          cam_surface_res); 
        m_sample_bvh.reset();
        if (resolved_contact_regime.activation.use_sample_bvh && !slave_samples.empty()) {
            m_sample_bvh = std::make_shared<spcc::SampleBVH>();
            m_sample_bvh->Build(slave_samples, resolved_contact_regime.activation.sample_bvh_leaf_size);
        }

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            cam_body, m_follower, m_gear1_sdf.get(), &m_contact_activation, slave_samples, m_sample_bvh, material,
            &m_active_contacts, enable_curvature_term, false, false, false, false, 0.0, 0.0,
            resolved_contact_regime.quadrature.manifold_kind);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif
}

double ChronoRigidSystemNSC::GetFollowerPosY() const {
    if (m_follower) return m_follower->GetPos().y();
    return 0.0;
}

double ChronoRigidSystemNSC::GetFollowerVelY() const {
    if (m_follower) return m_follower->GetPosDt().y();
    return 0.0;
}

bool ChronoRigidSystemNSC::GetCamBodySnapshot(BodyDebugSnapshot& out) const {
    out = MakeBodyDebugSnapshot(m_cam_body);
    return out.valid;
}

bool ChronoRigidSystemNSC::GetFollowerBodySnapshot(BodyDebugSnapshot& out) const {
    out = MakeBodyDebugSnapshot(m_follower);
    return out.valid;
}

#if defined(SPCC_ENABLE_VDB)
const std::vector<platform::backend::spcc::ActiveContactSample>& ChronoRigidSystemNSC::GetActiveContacts() const {
    return m_active_contacts;
}
#endif

void ChronoRigidSystemNSC::InitializeSimpleGearCase(
    const std::string& gear1_obj, const std::string& gear2_obj,
    const double* gear1_ref_pos, const double* gear2_ref_pos,
    const double* joint1_pos, const double* joint2_pos,
    double mesh_scale, double density,
    double friction, double restitution,
    double gravity_y, double motor_speed,
    double gear1_mass, const double* gear1_inertia_xx,
    double gear2_mass, const double* gear2_inertia_xx,
    bool enable_curvature_term,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::ContactRegimeConfig& contact_regime
) {
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
    m_sample_bvh.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction((float)friction);
    material->SetRestitution((float)restitution);

    auto load_scaled_mesh = [mesh_scale](const std::string& mesh_file) {
        auto mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        mesh->LoadWavefrontMesh(mesh_file, true, false);

        chrono::ChMatrix33<> scale_matrix(1);
        scale_matrix(0, 0) = mesh_scale;
        scale_matrix(1, 1) = mesh_scale;
        scale_matrix(2, 2) = mesh_scale;
        mesh->Transform(chrono::ChVector3d(0, 0, 0), scale_matrix);
        return mesh;
    };

    // Ground anchor for joints
    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    std::cout << "[DEBUG] Loading gear meshes..." << std::endl;
    auto mesh1 = load_scaled_mesh(gear1_obj);
    auto mesh2 = load_scaled_mesh(gear2_obj);
    std::cout << "[DEBUG] Loaded gear meshes" << std::endl;

#if defined(SPCC_ENABLE_VDB)
    std::cout << "[DEBUG] Using VDB SDF instead of NSC Mesh Collision" << std::endl;
    // Set create_collision flag (5th param) to FALSE to avoid HACD OOM
    // Params: mesh, density, compute_mass, create_visualization, create_collision, material, sphere_swept
    m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh1, density, true, true, false, material, 0.000001);
    m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh2, density, true, true, false, material, 0.000001);

    // Provide a dummy collision shape so the bodies have a valid ChCollisionModel
    // otherwise the custom SDF contact injection will crash due to null models.
    auto dummy_shape1 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
    m_gear1->AddCollisionShape(dummy_shape1);
    m_gear1->EnableCollision(true);

    auto dummy_shape2 = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1e-6);
    m_gear2->AddCollisionShape(dummy_shape2);
    m_gear2->EnableCollision(true);

    // Build SDF for gear 1
    spcc::VDBSDFField::BuildOptions sdf_options;
    // Use an adequate voxel size based on mesh_scale or global size.
    // Make sure it's smaller than the gear mesh features!
    sdf_options.voxel_size = GetEnvDouble("SPCC_GEAR_VOXEL_SIZE", sdf_build_tuning.voxel_size);
    sdf_options.half_band_width_voxels =
        GetEnvDouble("SPCC_GEAR_HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
    sdf_options.direct_phi_hessian =
        (GetEnvDouble("SPCC_GEAR_DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) > 0.5);
    m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
    std::cout << "[DEBUG] Building VDB SDF for gear 1..." << std::endl;
    m_gear1_sdf->BuildFromTriangleMesh(*mesh1, sdf_options);
    
    // Configure contact activation
    m_contact_activation.SetEnvPrefix("SPCC_GEAR");
    auto resolved_contact_regime = contact_regime;
    resolved_contact_regime.activation.delta_on =
        GetEnvDouble("SPCC_GEAR_DELTA_ON", resolved_contact_regime.activation.delta_on);
    resolved_contact_regime.activation.delta_off =
        GetEnvDouble("SPCC_GEAR_DELTA_OFF", resolved_contact_regime.activation.delta_off);
    resolved_contact_regime.activation.hold_steps =
        GetEnvInt("SPCC_GEAR_HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
    resolved_contact_regime.activation.max_active_keep =
        GetEnvInt("SPCC_GEAR_MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
    resolved_contact_regime.activation.use_sample_bvh =
        (GetEnvDouble("SPCC_GEAR_USE_SAMPLE_BVH",
                      resolved_contact_regime.activation.use_sample_bvh ? 1.0 : 0.0) > 0.5);
    resolved_contact_regime.activation.sample_bvh_leaf_size =
        GetEnvInt("SPCC_GEAR_SAMPLE_BVH_LEAF_SIZE", resolved_contact_regime.activation.sample_bvh_leaf_size);
    resolved_contact_regime.activation.sample_bvh_margin_scale =
        GetEnvDouble("SPCC_GEAR_SAMPLE_BVH_MARGIN_SCALE",
                     resolved_contact_regime.activation.sample_bvh_margin_scale);
    resolved_contact_regime.curvature.enabled =
        (GetEnvDouble("SPCC_GEAR_CURVATURE_ENABLED", resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
    resolved_contact_regime.curvature.tangential_only =
        (GetEnvDouble("SPCC_GEAR_CURVATURE_TANGENTIAL_ONLY",
                      resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
    resolved_contact_regime.curvature.normal_alignment_cos_min =
        GetEnvDouble("SPCC_GEAR_CURVATURE_ALIGNMENT_COS_MIN",
                     resolved_contact_regime.curvature.normal_alignment_cos_min);
    resolved_contact_regime.curvature.max_hessian_frobenius =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_HESSIAN_FROBENIUS",
                     resolved_contact_regime.curvature.max_hessian_frobenius);
    resolved_contact_regime.curvature.small_step_dt_threshold =
        GetEnvDouble("SPCC_GEAR_CURVATURE_SMALL_STEP_DT_THRESHOLD",
                     resolved_contact_regime.curvature.small_step_dt_threshold);
    resolved_contact_regime.curvature.small_step_max_hessian_frobenius =
        GetEnvDouble("SPCC_GEAR_CURVATURE_SMALL_STEP_MAX_HESSIAN_FROBENIUS",
                     resolved_contact_regime.curvature.small_step_max_hessian_frobenius);
    resolved_contact_regime.curvature.max_curvature_term_abs =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_ABS", resolved_contact_regime.curvature.max_curvature_term_abs);
    resolved_contact_regime.curvature.max_curvature_term_ratio =
        GetEnvDouble("SPCC_GEAR_CURVATURE_MAX_RATIO", resolved_contact_regime.curvature.max_curvature_term_ratio);
    resolved_contact_regime.curvature.gap_floor =
        GetEnvDouble("SPCC_GEAR_CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
    m_contact_activation.SetPolicy(resolved_contact_regime);

    // Register custom collision callback and sample slave vertices
    std::vector<chrono::ChVector3d> slave_samples;
    // Downsample slightly or use all vertices
    const double gear_surface_res = GetEnvDouble("SPCC_GEAR_SURFACE_RES", sample_tuning.surface_res);
    const int gear_max_samples = GetEnvInt("SPCC_GEAR_MAX_SAMPLES", sample_tuning.max_samples);
    slave_samples = BuildDownsampledVertexSamples(*mesh2, static_cast<std::size_t>(std::max(0, gear_max_samples)),
                                                  gear_surface_res); 
    m_sample_bvh.reset();
    if (resolved_contact_regime.activation.use_sample_bvh && !slave_samples.empty()) {
        m_sample_bvh = std::make_shared<spcc::SampleBVH>();
        m_sample_bvh->Build(slave_samples, resolved_contact_regime.activation.sample_bvh_leaf_size);
    }

    auto sdf_callback = std::make_shared<SDFCollisionCallback>(
        m_gear1, m_gear2, m_gear1_sdf.get(), &m_contact_activation, slave_samples, m_sample_bvh, material,
        &m_active_contacts, enable_curvature_term, false, false, false, false, 0.0, 0.0,
        resolved_contact_regime.quadrature.manifold_kind);
    m_system->RegisterCustomCollisionCallback(sdf_callback);
    std::cout << "[DEBUG] VDB SDF initialized successfully." << std::endl;
#else
    m_gear1 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh1, density, true, true, true, material, 0.000001);
    std::cout << "[DEBUG] Created ChBodyEasyMesh gear1" << std::endl;
    m_gear2 = std::make_shared<chrono::ChBodyEasyMesh>(
        mesh2, density, true, true, true, material, 0.000001);
    std::cout << "[DEBUG] Created ChBodyEasyMesh gear2" << std::endl;
#endif

    auto gear1_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_gear1);
    if (gear1_aux) {
        gear1_aux->SetFrameRefToAbs(
            chrono::ChFrame<>(
                chrono::ChVector3d(gear1_ref_pos[0], gear1_ref_pos[1], gear1_ref_pos[2]),
                chrono::QUNIT));
    }

    auto gear2_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(m_gear2);
    if (gear2_aux) {
        gear2_aux->SetFrameRefToAbs(
            chrono::ChFrame<>(
                chrono::ChVector3d(gear2_ref_pos[0], gear2_ref_pos[1], gear2_ref_pos[2]),
                chrono::QUNIT));
    }

    // Overwrite mass/inertia with values parsed from RecurDyn model.
    m_gear1->SetMass(gear1_mass);
    m_gear1->SetInertiaXX(
        chrono::ChVector3d(gear1_inertia_xx[0], gear1_inertia_xx[1], gear1_inertia_xx[2]));
    m_gear2->SetMass(gear2_mass);
    m_gear2->SetInertiaXX(
        chrono::ChVector3d(gear2_inertia_xx[0], gear2_inertia_xx[1], gear2_inertia_xx[2]));

    m_gear1->SetSleepingAllowed(false);
    m_gear2->SetSleepingAllowed(false);

    m_system->AddBody(m_gear1);
    m_system->AddBody(m_gear2);

    // RecurDyn results indicate RX is dominant, so use joint axis along global X.
    auto x_axis_revolute_rot = chrono::QuatFromAngleY(-chrono::CH_PI_2);

    auto motor = std::make_shared<chrono::ChLinkMotorRotationSpeed>();
    motor->Initialize(
        m_gear1,
        ground,
        chrono::ChFrame<>(
            chrono::ChVector3d(joint1_pos[0], joint1_pos[1], joint1_pos[2]),
            x_axis_revolute_rot));
    auto mfun = std::make_shared<chrono::ChFunctionConst>(motor_speed);
    motor->SetSpeedFunction(mfun);
    m_system->AddLink(motor);

    auto rev2 = std::make_shared<chrono::ChLinkLockRevolute>();
    rev2->Initialize(
        m_gear2,
        ground,
        chrono::ChFrame<>(
            chrono::ChVector3d(joint2_pos[0], joint2_pos[1], joint2_pos[2]),
            x_axis_revolute_rot));
    m_system->AddLink(rev2);
}

void ChronoRigidSystemNSC::InitializeRevoluteClearanceCase(
    const std::string& body1_obj,
    const std::string& body3_obj,
    const double* body1_pos,
    const double* body3_pos,
    const double* body2_cm_offset,
    double body3_mass,
    const double* body3_inertia_xx,
    const double* body3_inertia_xy,
    double body2_mass,
    const double* body2_inertia_xx,
    const double* body2_inertia_xy,
    double friction,
    double restitution,
    double gravity_y,
    double contact_compliance,
    double contact_compliance_t,
    double contact_damping_f,
    double collision_envelope,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::ContactRegimeConfig& contact_regime,
    bool use_lcp_manifold_quadrature,
    int manifold_quadrature_contacts,
    double manifold_quadrature_span_scale,
    double manifold_quadrature_min_half_span) {
    m_sphere.reset();
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
    m_sample_bvh.reset();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction(static_cast<float>(friction));
    material->SetRestitution(static_cast<float>(restitution));
    material->SetCompliance(static_cast<float>(GetScopedEnvDouble(env_prefix, "COMPLIANCE", contact_compliance)));
    material->SetComplianceT(static_cast<float>(GetScopedEnvDouble(env_prefix, "COMPLIANCE_T", contact_compliance_t)));
    material->SetDampingF(static_cast<float>(GetScopedEnvDouble(env_prefix, "DAMPING_F", contact_damping_f)));

    const float resolved_collision_envelope =
        static_cast<float>(GetScopedEnvDouble(env_prefix, "COLLISION_ENVELOPE", collision_envelope));

    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    const bool enable_curvature_term =
        (contact_algorithm == platform::common::ContactAlgorithm::SdfSecondOrder);
#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
    (void)sdf_build_tuning;
    (void)sample_tuning;
    (void)contact_regime;
    (void)env_prefix;
    (void)dynamics_substeps;
#endif

    const auto body1_pos_W = chrono::ChVector3d(body1_pos[0], body1_pos[1], body1_pos[2]);
    const auto body3_pos_W = chrono::ChVector3d(body3_pos[0], body3_pos[1], body3_pos[2]);
    const auto body2_cm_ref = chrono::ChVector3d(body2_cm_offset[0], body2_cm_offset[1], body2_cm_offset[2]);

#if defined(SPCC_ENABLE_VDB)
    std::shared_ptr<chrono::ChTriangleMeshConnected> body1_mesh;
    std::shared_ptr<chrono::ChTriangleMeshConnected> body3_mesh;
    if (use_vdb) {
        body1_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        body1_mesh->LoadWavefrontMesh(body1_obj, true, false);
        body3_mesh = std::make_shared<chrono::ChTriangleMeshConnected>();
        body3_mesh->LoadWavefrontMesh(body3_obj, true, false);
    }
#endif

    std::shared_ptr<chrono::ChBodyEasyMesh> body1;
    if (use_vdb) {
        body1 = std::make_shared<chrono::ChBodyEasyMesh>(body1_obj, 1000.0, true, true, false, material, 1.0e-6);
        auto dummy_shape = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        body1->AddCollisionShape(dummy_shape);
        body1->EnableCollision(true);
    } else {
        body1 = std::make_shared<chrono::ChBodyEasyMesh>(body1_obj, 1000.0, true, true, true, material, 1.0e-6);
    }
    body1->SetFixed(true);
    auto body1_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body1);
    if (body1_aux) {
        body1_aux->SetFrameRefToAbs(chrono::ChFrame<>(body1_pos_W, chrono::QUNIT));
    } else {
        body1->SetPos(body1_pos_W);
    }
    if (body1->GetCollisionModel()) {
        body1->GetCollisionModel()->SetEnvelope(resolved_collision_envelope);
    }
    m_system->AddBody(body1);
    m_clearance_body1 = body1;

    std::shared_ptr<chrono::ChBodyEasyMesh> body3;
    if (use_vdb) {
        body3 = std::make_shared<chrono::ChBodyEasyMesh>(body3_obj, 1000.0, true, true, false, material, 1.0e-6);
        auto dummy_shape = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        body3->AddCollisionShape(dummy_shape);
        body3->EnableCollision(true);
    } else {
        body3 = std::make_shared<chrono::ChBodyEasyMesh>(body3_obj, 1000.0, true, true, true, material, 1.0e-6);
    }
    auto body3_aux = std::dynamic_pointer_cast<chrono::ChBodyAuxRef>(body3);
    if (body3_aux) {
        body3_aux->SetFrameRefToAbs(chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    } else {
        body3->SetPos(body3_pos_W);
    }
    body3->SetMass(body3_mass);
    body3->SetInertiaXX(chrono::ChVector3d(body3_inertia_xx[0], body3_inertia_xx[1], body3_inertia_xx[2]));
    body3->SetInertiaXY(chrono::ChVector3d(body3_inertia_xy[0], body3_inertia_xy[1], body3_inertia_xy[2]));
    body3->SetSleepingAllowed(false);
    if (body3->GetCollisionModel()) {
        body3->GetCollisionModel()->SetEnvelope(resolved_collision_envelope);
    }
    m_system->AddBody(body3);
    m_clearance_body3 = body3;

    auto body2 = std::make_shared<chrono::ChBodyAuxRef>();
    body2->SetFixed(false);
    body2->EnableCollision(false);
    body2->SetMass(body2_mass);
    body2->SetInertiaXX(chrono::ChVector3d(body2_inertia_xx[0], body2_inertia_xx[1], body2_inertia_xx[2]));
    body2->SetInertiaXY(chrono::ChVector3d(body2_inertia_xy[0], body2_inertia_xy[1], body2_inertia_xy[2]));
    body2->SetFrameCOMToRef(chrono::ChFrame<>(body2_cm_ref, chrono::QUNIT));
    body2->SetFrameRefToAbs(chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    body2->SetSleepingAllowed(false);
    m_system->AddBody(body2);
    m_clearance_body2 = body2;

    auto fixed_link = std::make_shared<chrono::ChLinkLockLock>();
    fixed_link->Initialize(body2, body3, chrono::ChFrame<>(body3_pos_W, chrono::QUNIT));
    m_system->AddLink(fixed_link);

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));
        const bool resolved_use_lcp_manifold_quadrature =
            (GetScopedEnvDouble(env_prefix, "LCP_MANIFOLD_QUADRATURE",
                                use_lcp_manifold_quadrature ? 1.0 : 0.0) > 0.5);
        const int resolved_manifold_quadrature_contacts =
            std::max(1, GetScopedEnvInt(env_prefix, "LCP_MANIFOLD_QUADRATURE_CONTACTS",
                                        manifold_quadrature_contacts));
        const double resolved_manifold_quadrature_span_scale =
            GetScopedEnvDouble(env_prefix, "LCP_MANIFOLD_QUADRATURE_SPAN_SCALE",
                               manifold_quadrature_span_scale);
        const double resolved_manifold_quadrature_min_half_span =
            GetScopedEnvDouble(env_prefix, "LCP_MANIFOLD_QUADRATURE_MIN_HALF_SPAN",
                               manifold_quadrature_min_half_span);

        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
        sdf_options.direct_phi_hessian =
            (GetScopedEnvDouble(env_prefix, "DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) >
             0.5);

        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*body1_mesh, sdf_options);

        m_contact_activation.SetEnvPrefix(env_prefix);
        auto resolved_contact_regime = contact_regime;
        resolved_contact_regime.activation.delta_on =
            GetScopedEnvDouble(env_prefix, "DELTA_ON", resolved_contact_regime.activation.delta_on);
        resolved_contact_regime.activation.delta_off =
            GetScopedEnvDouble(env_prefix, "DELTA_OFF", resolved_contact_regime.activation.delta_off);
        resolved_contact_regime.activation.hold_steps =
            GetScopedEnvInt(env_prefix, "HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
        resolved_contact_regime.activation.max_active_keep =
            GetScopedEnvInt(env_prefix, "MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
        resolved_contact_regime.activation.use_sample_bvh =
            (GetScopedEnvDouble(env_prefix, "USE_SAMPLE_BVH",
                                resolved_contact_regime.activation.use_sample_bvh ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.activation.sample_bvh_leaf_size =
            GetScopedEnvInt(env_prefix, "SAMPLE_BVH_LEAF_SIZE",
                            resolved_contact_regime.activation.sample_bvh_leaf_size);
        resolved_contact_regime.activation.sample_bvh_margin_scale =
            GetScopedEnvDouble(env_prefix, "SAMPLE_BVH_MARGIN_SCALE",
                               resolved_contact_regime.activation.sample_bvh_margin_scale);
        resolved_contact_regime.activation.sample_bvh_use_persistent_seeds_only =
            (GetScopedEnvDouble(env_prefix, "SAMPLE_BVH_PERSISTENT_ONLY",
                                resolved_contact_regime.activation.sample_bvh_use_persistent_seeds_only ? 1.0 : 0.0) >
             0.5);
        resolved_contact_regime.curvature.enabled =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_ENABLED",
                                resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.tangential_only =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_TANGENTIAL_ONLY",
                                resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.normal_alignment_cos_min =
            GetScopedEnvDouble(env_prefix, "CURVATURE_ALIGNMENT_COS_MIN",
                               resolved_contact_regime.curvature.normal_alignment_cos_min);
        resolved_contact_regime.curvature.max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.max_hessian_frobenius);
        resolved_contact_regime.curvature.small_step_dt_threshold =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_DT_THRESHOLD",
                               resolved_contact_regime.curvature.small_step_dt_threshold);
        resolved_contact_regime.curvature.small_step_max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.small_step_max_hessian_frobenius);
        resolved_contact_regime.curvature.max_curvature_term_abs =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_ABS",
                               resolved_contact_regime.curvature.max_curvature_term_abs);
        resolved_contact_regime.curvature.max_curvature_term_ratio =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_RATIO",
                               resolved_contact_regime.curvature.max_curvature_term_ratio);
        resolved_contact_regime.curvature.gap_floor =
            GetScopedEnvDouble(env_prefix, "CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
        m_contact_activation.SetPolicy(resolved_contact_regime);

        const double surface_res = GetScopedEnvDouble(env_prefix, "SURFACE_RES", sample_tuning.surface_res);
        const int max_samples = GetScopedEnvInt(env_prefix, "MAX_SAMPLES", sample_tuning.max_samples);
        std::vector<chrono::ChVector3d> slave_samples =
            BuildDownsampledVertexSamples(*body3_mesh, static_cast<std::size_t>(std::max(0, max_samples)), surface_res);
        m_sample_bvh.reset();
        if (resolved_contact_regime.activation.use_sample_bvh && !slave_samples.empty()) {
            m_sample_bvh = std::make_shared<spcc::SampleBVH>();
            m_sample_bvh->Build(slave_samples, resolved_contact_regime.activation.sample_bvh_leaf_size);
        }

        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            body1, body3, m_gear1_sdf.get(), &m_contact_activation, slave_samples, m_sample_bvh, material,
            &m_active_contacts, enable_curvature_term, false, false, false, false, 0.0, 0.0,
            resolved_contact_regime.quadrature.manifold_kind,
            resolved_use_lcp_manifold_quadrature, resolved_manifold_quadrature_contacts,
            resolved_manifold_quadrature_span_scale, resolved_manifold_quadrature_min_half_span);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif
}

void ChronoRigidSystemNSC::InitializeHeadOnSphereCase(
    const std::string& sphere_a_obj,
    const std::string& sphere_b_obj,
    double sphere_radius,
    double sphere_a_density,
    double sphere_b_density,
    const double* sphere_a_pos,
    const double* sphere_b_pos,
    const double* sphere_a_vel,
    const double* sphere_b_vel,
    double friction,
    double restitution,
    double gravity_y,
    int dynamics_substeps,
    const std::string& env_prefix,
    platform::common::ContactAlgorithm contact_algorithm,
    const platform::backend::spcc::SdfBuildTuning& sdf_build_tuning,
    const platform::backend::spcc::SurfaceSampleTuning& sample_tuning,
    const platform::backend::spcc::ContactRegimeConfig& contact_regime) {
    m_sphere.reset();
    m_headon_a.reset();
    m_headon_b.reset();
    m_cam_body.reset();
    m_follower.reset();
    m_gear1.reset();
    m_gear2.reset();
    m_clearance_body1.reset();
    m_clearance_body2.reset();
    m_clearance_body3.reset();
    m_dynamics_substeps = 1;
    m_headon_manual_elastic_resolution = false;
    m_headon_sphere_radius = 0.0;
    m_headon_restitution = 1.0;
#if defined(SPCC_ENABLE_VDB)
    m_active_contacts.clear();
    m_sample_bvh.reset();
    m_gear1_sdf.reset();
#endif

    m_system->SetCollisionSystemType(chrono::ChCollisionSystem::Type::BULLET);
    m_system->SetGravitationalAcceleration(chrono::ChVector3d(0, gravity_y, 0));
    m_headon_sphere_radius = sphere_radius;
    m_headon_restitution = restitution;

    auto material = std::make_shared<chrono::ChContactMaterialNSC>();
    material->SetFriction(static_cast<float>(friction));
    material->SetRestitution(static_cast<float>(restitution));

    auto ground = std::make_shared<chrono::ChBody>();
    ground->SetFixed(true);
    m_system->AddBody(ground);

    const auto sphere_a_pos_W = chrono::ChVector3d(sphere_a_pos[0], sphere_a_pos[1], sphere_a_pos[2]);
    const auto sphere_b_pos_W = chrono::ChVector3d(sphere_b_pos[0], sphere_b_pos[1], sphere_b_pos[2]);
    const auto sphere_a_vel_W = chrono::ChVector3d(sphere_a_vel[0], sphere_a_vel[1], sphere_a_vel[2]);
    const auto sphere_b_vel_W = chrono::ChVector3d(sphere_b_vel[0], sphere_b_vel[1], sphere_b_vel[2]);
    const auto x_axis_prismatic_rot = chrono::QuatFromAngleY(-chrono::CH_PI_2);
    const double sphere_volume =
        (4.0 / 3.0) * chrono::CH_PI * sphere_radius * sphere_radius * sphere_radius;
    const double sphere_a_mass = sphere_a_density * sphere_volume;
    const double sphere_b_mass = sphere_b_density * sphere_volume;
    const double sphere_a_inertia = 0.4 * sphere_a_mass * sphere_radius * sphere_radius;
    const double sphere_b_inertia = 0.4 * sphere_b_mass * sphere_radius * sphere_radius;

#if defined(SPCC_ENABLE_VDB)
    const bool use_vdb = (contact_algorithm != platform::common::ContactAlgorithm::NativeMesh);
    m_headon_manual_elastic_resolution =
        use_vdb && (GetScopedEnvDouble(env_prefix, "USE_EXACT_CORRECTION", 0.0) > 0.5);
#else
    const bool use_vdb = false;
    (void)contact_algorithm;
    (void)sdf_build_tuning;
    (void)sample_tuning;
    (void)contact_regime;
    (void)env_prefix;
    (void)dynamics_substeps;
#endif

    if (!use_vdb) {
        auto sphere_a =
            std::make_shared<chrono::ChBodyEasySphere>(sphere_radius, sphere_a_density, true, true, material);
        auto sphere_b =
            std::make_shared<chrono::ChBodyEasySphere>(sphere_radius, sphere_b_density, true, true, material);
        sphere_a->SetFixed(false);
        sphere_b->SetFixed(false);
        sphere_a->SetSleepingAllowed(false);
        sphere_b->SetSleepingAllowed(false);
        sphere_a->SetPos(sphere_a_pos_W);
        sphere_b->SetPos(sphere_b_pos_W);
        sphere_a->SetPosDt(sphere_a_vel_W);
        sphere_b->SetPosDt(sphere_b_vel_W);
        sphere_a->SetMass(sphere_a_mass);
        sphere_b->SetMass(sphere_b_mass);
        sphere_a->SetInertiaXX(chrono::ChVector3d(sphere_a_inertia, sphere_a_inertia, sphere_a_inertia));
        sphere_b->SetInertiaXX(chrono::ChVector3d(sphere_b_inertia, sphere_b_inertia, sphere_b_inertia));
        m_system->AddBody(sphere_a);
        m_system->AddBody(sphere_b);
        m_headon_a = sphere_a;
        m_headon_b = sphere_b;
    }

#if defined(SPCC_ENABLE_VDB)
    if (use_vdb) {
        auto mesh_a = std::make_shared<chrono::ChTriangleMeshConnected>();
        auto mesh_b = std::make_shared<chrono::ChTriangleMeshConnected>();
        bool loaded_a = !sphere_a_obj.empty() && mesh_a->LoadWavefrontMesh(sphere_a_obj, true, false);
        bool loaded_b = !sphere_b_obj.empty() && mesh_b->LoadWavefrontMesh(sphere_b_obj, true, false);
        if (!loaded_a) {
            mesh_a = BuildUvSphereMesh(sphere_radius, 64, 32);
        }
        if (!loaded_b) {
            mesh_b = BuildUvSphereMesh(sphere_radius, 64, 32);
        }

        auto sphere_a = std::make_shared<chrono::ChBodyEasyMesh>(mesh_a, sphere_a_density, true, true, false, material,
                                                                 1.0e-6);
        auto sphere_b = std::make_shared<chrono::ChBodyEasyMesh>(mesh_b, sphere_b_density, true, true, false, material,
                                                                 1.0e-6);

        auto dummy_shape_a = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        auto dummy_shape_b = chrono_types::make_shared<chrono::ChCollisionShapeSphere>(material, 1.0e-6);
        sphere_a->AddCollisionShape(dummy_shape_a);
        sphere_a->EnableCollision(true);
        sphere_b->AddCollisionShape(dummy_shape_b);
        sphere_b->EnableCollision(true);

        SetBodyFrameOrPos(sphere_a, sphere_a_pos_W);
        SetBodyFrameOrPos(sphere_b, sphere_b_pos_W);
        sphere_a->SetFixed(false);
        sphere_b->SetFixed(false);
        sphere_a->SetSleepingAllowed(false);
        sphere_b->SetSleepingAllowed(false);
        sphere_a->SetPosDt(sphere_a_vel_W);
        sphere_b->SetPosDt(sphere_b_vel_W);
        sphere_a->SetMass(sphere_a_mass);
        sphere_b->SetMass(sphere_b_mass);
        sphere_a->SetInertiaXX(chrono::ChVector3d(sphere_a_inertia, sphere_a_inertia, sphere_a_inertia));
        sphere_b->SetInertiaXX(chrono::ChVector3d(sphere_b_inertia, sphere_b_inertia, sphere_b_inertia));
        m_system->AddBody(sphere_a);
        m_system->AddBody(sphere_b);
        m_headon_a = sphere_a;
        m_headon_b = sphere_b;

        m_dynamics_substeps = std::max(1, GetScopedEnvInt(env_prefix, "DYNAMICS_SUBSTEPS", dynamics_substeps));

        spcc::VDBSDFField::BuildOptions sdf_options;
        sdf_options.voxel_size = GetScopedEnvDouble(env_prefix, "VOXEL_SIZE", sdf_build_tuning.voxel_size);
        sdf_options.half_band_width_voxels =
            GetScopedEnvDouble(env_prefix, "HALF_BAND_VOXELS", sdf_build_tuning.half_band_width_voxels);
        sdf_options.direct_phi_hessian =
            (GetScopedEnvDouble(env_prefix, "DIRECT_PHI_HESSIAN", sdf_build_tuning.direct_phi_hessian ? 1.0 : 0.0) >
             0.5);

        m_gear1_sdf = std::make_unique<spcc::VDBSDFField>();
        m_gear1_sdf->BuildFromTriangleMesh(*mesh_a, sdf_options);

        m_contact_activation.SetEnvPrefix(env_prefix);
        auto resolved_contact_regime = contact_regime;
        resolved_contact_regime.activation.delta_on =
            GetScopedEnvDouble(env_prefix, "DELTA_ON", resolved_contact_regime.activation.delta_on);
        resolved_contact_regime.activation.delta_off =
            GetScopedEnvDouble(env_prefix, "DELTA_OFF", resolved_contact_regime.activation.delta_off);
        resolved_contact_regime.activation.hold_steps =
            GetScopedEnvInt(env_prefix, "HOLD_STEPS", resolved_contact_regime.activation.hold_steps);
        resolved_contact_regime.activation.max_active_keep =
            GetScopedEnvInt(env_prefix, "MAX_CONTACTS", resolved_contact_regime.activation.max_active_keep);
        resolved_contact_regime.activation.use_sample_bvh =
            (GetScopedEnvDouble(env_prefix, "USE_SAMPLE_BVH",
                                resolved_contact_regime.activation.use_sample_bvh ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.activation.sample_bvh_leaf_size =
            GetScopedEnvInt(env_prefix, "SAMPLE_BVH_LEAF_SIZE",
                            resolved_contact_regime.activation.sample_bvh_leaf_size);
        resolved_contact_regime.activation.sample_bvh_margin_scale =
            GetScopedEnvDouble(env_prefix, "SAMPLE_BVH_MARGIN_SCALE",
                               resolved_contact_regime.activation.sample_bvh_margin_scale);
        resolved_contact_regime.activation.sample_bvh_use_persistent_seeds_only =
            (GetScopedEnvDouble(env_prefix, "SAMPLE_BVH_PERSISTENT_ONLY",
                                resolved_contact_regime.activation.sample_bvh_use_persistent_seeds_only ? 1.0 : 0.0) >
             0.5);
        resolved_contact_regime.curvature.enabled =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_ENABLED",
                                resolved_contact_regime.curvature.enabled ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.tangential_only =
            (GetScopedEnvDouble(env_prefix, "CURVATURE_TANGENTIAL_ONLY",
                                resolved_contact_regime.curvature.tangential_only ? 1.0 : 0.0) > 0.5);
        resolved_contact_regime.curvature.normal_alignment_cos_min =
            GetScopedEnvDouble(env_prefix, "CURVATURE_ALIGNMENT_COS_MIN",
                               resolved_contact_regime.curvature.normal_alignment_cos_min);
        resolved_contact_regime.curvature.max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.max_hessian_frobenius);
        resolved_contact_regime.curvature.small_step_dt_threshold =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_DT_THRESHOLD",
                               resolved_contact_regime.curvature.small_step_dt_threshold);
        resolved_contact_regime.curvature.small_step_max_hessian_frobenius =
            GetScopedEnvDouble(env_prefix, "CURVATURE_SMALL_STEP_MAX_HESSIAN_FROBENIUS",
                               resolved_contact_regime.curvature.small_step_max_hessian_frobenius);
        resolved_contact_regime.curvature.max_curvature_term_abs =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_ABS",
                               resolved_contact_regime.curvature.max_curvature_term_abs);
        resolved_contact_regime.curvature.max_curvature_term_ratio =
            GetScopedEnvDouble(env_prefix, "CURVATURE_MAX_RATIO",
                               resolved_contact_regime.curvature.max_curvature_term_ratio);
        resolved_contact_regime.curvature.gap_floor =
            GetScopedEnvDouble(env_prefix, "CURVATURE_GAP_FLOOR", resolved_contact_regime.curvature.gap_floor);
        m_contact_activation.SetPolicy(resolved_contact_regime);

        std::vector<chrono::ChVector3d> slave_samples;
        slave_samples.emplace_back(-sphere_radius, 0.0, 0.0);
        m_sample_bvh.reset();
        if (resolved_contact_regime.activation.use_sample_bvh && !slave_samples.empty()) {
            m_sample_bvh = std::make_shared<spcc::SampleBVH>();
            m_sample_bvh->Build(slave_samples, resolved_contact_regime.activation.sample_bvh_leaf_size);
        }

        const bool use_analytic_toi_contact =
            (GetScopedEnvDouble(env_prefix, "USE_ANALYTIC_TOI_CONTACT", 0.0) > 0.5);
        auto sdf_callback = std::make_shared<SDFCollisionCallback>(
            sphere_a, sphere_b, m_gear1_sdf.get(), &m_contact_activation, slave_samples, m_sample_bvh, material,
            &m_active_contacts, contact_algorithm == platform::common::ContactAlgorithm::SdfSecondOrder,
            true, false, true, use_analytic_toi_contact, sphere_radius, sphere_radius,
            resolved_contact_regime.quadrature.manifold_kind);
        m_system->RegisterCustomCollisionCallback(sdf_callback);
    }
#endif

    auto prism_a = std::make_shared<chrono::ChLinkLockPrismatic>();
    prism_a->Initialize(m_headon_a, ground, chrono::ChFrame<>(m_headon_a->GetPos(), x_axis_prismatic_rot));
    m_system->AddLink(prism_a);

    auto prism_b = std::make_shared<chrono::ChLinkLockPrismatic>();
    prism_b->Initialize(m_headon_b, ground, chrono::ChFrame<>(m_headon_b->GetPos(), x_axis_prismatic_rot));
    m_system->AddLink(prism_b);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetGear1AngVelParent() const {
    if (m_gear1) return m_gear1->GetAngVelParent();
    return chrono::ChVector3d(0, 0, 0);
}

chrono::ChVector3d ChronoRigidSystemNSC::GetGear2AngVelParent() const {
    if (m_gear2) return m_gear2->GetAngVelParent();
    return chrono::ChVector3d(0, 0, 0);
}

} // namespace backend
} // namespace platform


