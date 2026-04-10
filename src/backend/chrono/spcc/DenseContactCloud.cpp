#include "platform/backend/spcc/DenseContactCloud.h"

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

struct RankedCandidate {
    std::size_t sample_id = 0;
    double phi = 0.0;
};

}  // namespace

void DenseContactCloudBuilder::Build(const CompressedContactConfig& cfg,
                                     const std::vector<DenseSurfaceSample>& slave_surface_samples,
                                     const DenseSampleBVH& dense_sample_bvh,
                                     const RigidBodyStateW& master_state,
                                     const RigidBodyStateW& slave_state,
                                     const FirstOrderSDF& sdf,
                                     double step_size,
                                     std::vector<DenseContactPoint>& out_dense_points,
                                     DenseContactCloudStats* out_stats) {
    out_dense_points.clear();

    DenseContactCloudStats stats;
    stats.total_samples = slave_surface_samples.size();

    std::vector<std::size_t> candidate_indices;
    if (!dense_sample_bvh.Empty()) {
        dense_sample_bvh.CollectCandidateSampleIndices(master_state, slave_state, sdf, cfg, step_size,
                                                       candidate_indices, &stats.bvh);
    } else {
        candidate_indices.resize(slave_surface_samples.size());
        for (std::size_t i = 0; i < slave_surface_samples.size(); ++i) {
            candidate_indices[i] = i;
        }
    }

    stats.candidate_samples = candidate_indices.size();
    if (cfg.max_exact_candidates > 0 &&
        static_cast<int>(candidate_indices.size()) > cfg.max_exact_candidates) {
        std::vector<RankedCandidate> ranked_candidates;
        ranked_candidates.reserve(candidate_indices.size());
        for (const auto sample_id : candidate_indices) {
            const auto& sample = slave_surface_samples[sample_id];
            const chrono::ChVector3d x_W = slave_state.x_ref_W + slave_state.R_WRef * sample.xi_slave_S;
            const chrono::ChVector3d x_master_M = master_state.R_WRef.transpose() * (x_W - master_state.x_ref_W);

            double phi = 0.0;
            if (!sdf.QueryPhiM(x_master_M, phi)) {
                continue;
            }

            RankedCandidate ranked;
            ranked.sample_id = sample_id;
            ranked.phi = phi;
            ranked_candidates.push_back(ranked);
        }

        stats.phi_prefilter_samples = ranked_candidates.size();
        const std::size_t exact_limit = static_cast<std::size_t>(std::max(cfg.max_exact_candidates, cfg.max_active_dense));
        if (ranked_candidates.size() > exact_limit) {
            std::stable_sort(ranked_candidates.begin(), ranked_candidates.end(),
                             [](const RankedCandidate& a, const RankedCandidate& b) { return a.phi < b.phi; });
            ranked_candidates.resize(exact_limit);
        }

        candidate_indices.clear();
        candidate_indices.reserve(ranked_candidates.size());
        for (const auto& ranked : ranked_candidates) {
            candidate_indices.push_back(ranked.sample_id);
        }
    }

    stats.exact_samples = candidate_indices.size();
    out_dense_points.reserve(candidate_indices.size());

    for (const auto sample_id : candidate_indices) {
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

        DenseContactPoint point;
        point.sample_id = sample_id;
        point.x_W = x_W;
        point.x_master_M = x_master_M;
        point.x_master_surface_W = x_W - phi * n_W;
        point.n_W = n_W;
        point.v_rel_W = v_rel_W;
        point.phi = phi;
        point.phi_eff = phi_eff;
        point.area_weight = sample.area_weight;
        out_dense_points.push_back(point);
    }

    if (cfg.max_active_dense > 0 && static_cast<int>(out_dense_points.size()) > cfg.max_active_dense) {
        std::stable_sort(out_dense_points.begin(), out_dense_points.end(), [](const DenseContactPoint& a,
                                                                              const DenseContactPoint& b) {
            return a.phi_eff < b.phi_eff;
        });
        out_dense_points.resize(static_cast<std::size_t>(cfg.max_active_dense));
    }

    stats.active_samples = out_dense_points.size();
    if (out_stats) {
        *out_stats = stats;
    }
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
