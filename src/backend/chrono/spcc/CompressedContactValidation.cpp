#include "platform/backend/spcc/CompressedContactValidation.h"

#include <algorithm>
#include <cmath>

namespace platform {
namespace backend {
namespace spcc {

namespace {

chrono::ChVector3d SafeNormalized(const chrono::ChVector3d& v, const chrono::ChVector3d& fallback) {
    const double len = v.Length();
    if (!(len > 1.0e-12) || !std::isfinite(len)) {
        return fallback;
    }
    return v * (1.0 / len);
}

double ProxyLoad(const DenseValidationContact& point) {
    return std::max(0.0, -point.phi_eff) * std::max(0.0, point.area_weight);
}

double ProxyLoad(const ReducedContactPoint& point) {
    if (point.allocated_load > 0.0) {
        return point.allocated_load;
    }
    return std::max(0.0, -point.phi_eff) * std::max(0.0, point.support_weight);
}

chrono::ChVector3d ProxyForce(const std::vector<DenseValidationContact>& points, double& total_load) {
    chrono::ChVector3d force(0.0, 0.0, 0.0);
    total_load = 0.0;
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        force += load * point.n_W;
        total_load += load;
    }
    return force;
}

chrono::ChVector3d ProxyForce(const std::vector<ReducedContactPoint>& points, double& total_load) {
    chrono::ChVector3d force(0.0, 0.0, 0.0);
    total_load = 0.0;
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        force += load * point.n_W;
        total_load += load;
    }
    return force;
}

chrono::ChVector3d ProxyCoP(const std::vector<DenseValidationContact>& points) {
    chrono::ChVector3d cop(0.0, 0.0, 0.0);
    double total_load = 0.0;
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        cop += load * point.x_W;
        total_load += load;
    }
    if (!(total_load > 1.0e-12)) {
        return cop;
    }
    return cop * (1.0 / total_load);
}

chrono::ChVector3d ProxyCoP(const std::vector<ReducedContactPoint>& points) {
    chrono::ChVector3d cop(0.0, 0.0, 0.0);
    double total_load = 0.0;
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        cop += load * point.x_W;
        total_load += load;
    }
    if (!(total_load > 1.0e-12)) {
        return cop;
    }
    return cop * (1.0 / total_load);
}

chrono::ChVector3d ProxyMomentAtOrigin(const std::vector<DenseValidationContact>& points,
                                       const chrono::ChVector3d& origin_W) {
    chrono::ChVector3d moment(0.0, 0.0, 0.0);
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        moment += chrono::Vcross(point.x_W - origin_W, load * point.n_W);
    }
    return moment;
}

chrono::ChVector3d ProxyMomentAtOrigin(const std::vector<ReducedContactPoint>& points,
                                       const chrono::ChVector3d& origin_W) {
    chrono::ChVector3d moment(0.0, 0.0, 0.0);
    for (const auto& point : points) {
        const double load = ProxyLoad(point);
        moment += chrono::Vcross(point.x_W - origin_W, load * point.n_W);
    }
    return moment;
}

double WorstGap(const std::vector<DenseValidationContact>& points) {
    double worst_gap = 0.0;
    for (const auto& point : points) {
        worst_gap = std::min(worst_gap, point.phi_eff);
    }
    return worst_gap;
}

double WorstGap(const std::vector<ReducedContactPoint>& points) {
    double worst_gap = 0.0;
    for (const auto& point : points) {
        worst_gap = std::min(worst_gap, point.phi_eff);
    }
    return worst_gap;
}

void BuildDenseContactsImpl(const CompressedContactConfig& cfg,
                            const std::vector<DenseSurfaceSample>& slave_surface_samples,
                            const RigidBodyStateW& master_state,
                            const RigidBodyStateW& slave_state,
                            const FirstOrderSDF& sdf,
                            double step_size,
                            std::vector<DenseValidationContact>& out_dense_contacts) {
    out_dense_contacts.clear();
    out_dense_contacts.reserve(slave_surface_samples.size());

    for (std::size_t sample_id = 0; sample_id < slave_surface_samples.size(); ++sample_id) {
        const auto& sample = slave_surface_samples[sample_id];
        const chrono::ChVector3d x_W = slave_state.x_ref_W + slave_state.R_WRef * sample.xi_slave_S;
        const chrono::ChVector3d x_master_M = master_state.R_WRef.transpose() * (x_W - master_state.x_ref_W);

        double phi = 0.0;
        chrono::ChVector3d grad_M;
        if (!sdf.QueryPhiGradM(x_master_M, phi, grad_M)) {
            continue;
        }

        chrono::ChVector3d n_W = master_state.R_WRef * grad_M;
        n_W = SafeNormalized(n_W, chrono::ChVector3d(0.0, 1.0, 0.0));

        const chrono::ChVector3d rA_W = x_W - master_state.x_com_W;
        const chrono::ChVector3d rB_W = x_W - slave_state.x_com_W;
        const chrono::ChVector3d v_master_W = master_state.v_com_W + chrono::Vcross(master_state.w_W, rA_W);
        const chrono::ChVector3d v_slave_W = slave_state.v_com_W + chrono::Vcross(slave_state.w_W, rB_W);
        const chrono::ChVector3d v_rel_W = v_slave_W - v_master_W;

        double phi_eff = phi;
        if (cfg.predictive_gap) {
            phi_eff += step_size * std::min(0.0, chrono::Vdot(n_W, v_rel_W));
        }

        if (!(phi_eff <= cfg.delta_on || phi <= cfg.delta_on)) {
            continue;
        }

        DenseValidationContact point;
        point.sample_id = sample_id;
        point.x_W = x_W;
        point.x_master_M = x_master_M;
        point.x_master_surface_W = x_W - phi * n_W;
        point.n_W = n_W;
        point.v_rel_W = v_rel_W;
        point.phi = phi;
        point.phi_eff = phi_eff;
        point.area_weight = sample.area_weight;
        out_dense_contacts.push_back(point);
    }

    if (cfg.max_active_dense > 0 && static_cast<int>(out_dense_contacts.size()) > cfg.max_active_dense) {
        std::stable_sort(out_dense_contacts.begin(), out_dense_contacts.end(),
                         [](const DenseValidationContact& a, const DenseValidationContact& b) {
                             return a.phi_eff < b.phi_eff;
                         });
        out_dense_contacts.resize(static_cast<std::size_t>(cfg.max_active_dense));
    }
}

CompressionWrenchSummary SummarizeDense(const std::vector<DenseValidationContact>& points) {
    CompressionWrenchSummary summary;
    summary.force_W = ProxyForce(points, summary.total_proxy_load);
    summary.cop_W = ProxyCoP(points);
    summary.moment_at_origin_W = ProxyMomentAtOrigin(points, chrono::ChVector3d(0.0, 0.0, 0.0));
    summary.worst_gap = WorstGap(points);
    return summary;
}

CompressionWrenchSummary SummarizeReduced(const std::vector<ReducedContactPoint>& points) {
    CompressionWrenchSummary summary;
    summary.force_W = ProxyForce(points, summary.total_proxy_load);
    summary.cop_W = ProxyCoP(points);
    summary.moment_at_origin_W = ProxyMomentAtOrigin(points, chrono::ChVector3d(0.0, 0.0, 0.0));
    summary.worst_gap = WorstGap(points);
    return summary;
}

}  // namespace

void CompressedContactValidation::BuildDenseContacts(const CompressedContactConfig& cfg,
                                                     const std::vector<DenseSurfaceSample>& slave_surface_samples,
                                                     const RigidBodyStateW& master_state,
                                                     const RigidBodyStateW& slave_state,
                                                     const FirstOrderSDF& sdf,
                                                     double step_size,
                                                     std::vector<DenseValidationContact>& out_dense_contacts) {
    BuildDenseContactsImpl(cfg, slave_surface_samples, master_state, slave_state, sdf, step_size,
                           out_dense_contacts);
}

void CompressedContactValidation::Validate(const CompressedContactConfig& cfg,
                                           const std::vector<DenseSurfaceSample>& slave_surface_samples,
                                           const RigidBodyStateW& master_state,
                                           const RigidBodyStateW& slave_state,
                                           const FirstOrderSDF& sdf,
                                           double mu_default,
                                           double step_size,
                                           CompressionValidationReport& out_report) {
    out_report = CompressionValidationReport{};

    BuildDenseContactsImpl(cfg, slave_surface_samples, master_state, slave_state, sdf, step_size,
                           out_report.dense_contacts);

    CompressedContactPipeline pipeline;
    pipeline.Configure(cfg);
    pipeline.SetSlaveSurfaceSamples(slave_surface_samples);
    pipeline.BuildReducedContacts(master_state, slave_state, sdf, mu_default, step_size,
                                  out_report.reduced_contacts, &out_report.stats);

    out_report.dense_wrench = SummarizeDense(out_report.dense_contacts);
    out_report.reduced_wrench = SummarizeReduced(out_report.reduced_contacts);
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
