#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <chrono/core/ChVector3.h>

#include "platform/backend/spcc/CompressedContactValidation.h"

namespace fs = std::filesystem;

namespace {

using platform::backend::spcc::CompressedContactConfig;
using platform::backend::spcc::CompressedContactValidation;
using platform::backend::spcc::CompressionValidationReport;
using platform::backend::spcc::DenseSurfaceSample;
using platform::backend::spcc::FirstOrderSDF;
using platform::backend::spcc::RigidBodyStateW;

chrono::ChVector3d Normalized(const chrono::ChVector3d& v) {
    const double len = v.Length();
    if (!(len > 1.0e-12)) {
        return chrono::ChVector3d(1.0, 0.0, 0.0);
    }
    return v * (1.0 / len);
}

class PlaneSDF final : public FirstOrderSDF {
  public:
    PlaneSDF(const chrono::ChVector3d& origin_W, const chrono::ChVector3d& normal_W)
        : origin_W_(origin_W), normal_W_(Normalized(normal_W)) {}

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        grad_M = normal_W_;
        phi = chrono::Vdot(normal_W_, x_M - origin_W_);
        return true;
    }

  private:
    chrono::ChVector3d origin_W_;
    chrono::ChVector3d normal_W_;
};

class CylinderSDF final : public FirstOrderSDF {
  public:
    explicit CylinderSDF(double radius) : radius_(radius) {}

    bool QueryPhiGradM(const chrono::ChVector3d& x_M,
                       double& phi,
                       chrono::ChVector3d& grad_M) const override {
        const double radial = std::sqrt(x_M.x() * x_M.x() + x_M.y() * x_M.y());
        phi = radial - radius_;
        if (radial > 1.0e-12) {
            grad_M = chrono::ChVector3d(x_M.x() / radial, x_M.y() / radial, 0.0);
        } else {
            grad_M = chrono::ChVector3d(1.0, 0.0, 0.0);
        }
        return true;
    }

  private:
    double radius_ = 0.0;
};

struct ScenarioDefinition {
    std::string name;
    std::string description;
    CompressedContactConfig cfg;
    RigidBodyStateW master_state;
    RigidBodyStateW slave_state;
    std::vector<DenseSurfaceSample> samples;
    std::unique_ptr<FirstOrderSDF> sdf;
    double step_size = 1.0e-3;
    double mu_default = 0.2;
    std::size_t expected_patch_count = 1;
};

struct ScenarioResult {
    ScenarioDefinition definition;
    CompressionValidationReport report;
    bool pass = false;
};

void BuildBasis(const chrono::ChVector3d& n_W,
                chrono::ChVector3d& t1_W,
                chrono::ChVector3d& t2_W) {
    const chrono::ChVector3d seed =
        (std::abs(n_W.z()) < 0.9) ? chrono::ChVector3d(0.0, 0.0, 1.0) : chrono::ChVector3d(1.0, 0.0, 0.0);
    t1_W = Normalized(chrono::Vcross(seed, n_W));
    t2_W = Normalized(chrono::Vcross(n_W, t1_W));
}

std::vector<DenseSurfaceSample> BuildRectSamples(const chrono::ChVector3d& center_W,
                                                 const chrono::ChVector3d& axis_u_W,
                                                 const chrono::ChVector3d& axis_v_W,
                                                 double size_u,
                                                 double size_v,
                                                 int nu,
                                                 int nv) {
    std::vector<DenseSurfaceSample> samples;
    if (nu <= 0 || nv <= 0) {
        return samples;
    }

    const chrono::ChVector3d u_W = Normalized(axis_u_W);
    const chrono::ChVector3d v_W = Normalized(axis_v_W);
    const chrono::ChVector3d n_W = Normalized(chrono::Vcross(u_W, v_W));
    const double area_weight = (size_u * size_v) / static_cast<double>(nu * nv);

    samples.reserve(static_cast<std::size_t>(nu * nv));
    for (int iu = 0; iu < nu; ++iu) {
        for (int iv = 0; iv < nv; ++iv) {
            const double u = ((static_cast<double>(iu) + 0.5) / static_cast<double>(nu) - 0.5) * size_u;
            const double v = ((static_cast<double>(iv) + 0.5) / static_cast<double>(nv) - 0.5) * size_v;

            DenseSurfaceSample sample;
            sample.xi_slave_S = center_W + u * u_W + v * v_W;
            sample.normal_slave_S = n_W;
            sample.area_weight = area_weight;
            samples.push_back(sample);
        }
    }
    return samples;
}

CompressedContactConfig MakeValidationConfig() {
    CompressedContactConfig cfg;
    cfg.delta_on = 5.0e-3;
    cfg.delta_off = 6.0e-3;
    cfg.max_active_dense = 0;
    cfg.bvh_leaf_size = 24;
    cfg.bvh_query_margin = 3.0e-3;
    cfg.bvh_velocity_bound_scale = 1.0;
    cfg.bvh_enable_sdf_node_bound = true;
    cfg.patch_radius = 0.05;
    cfg.normal_cos_min = 0.95;
    cfg.max_patch_diameter = 0.14;
    cfg.max_subpatch_diameter = 0.0;
    cfg.max_plane_error = 0.0;
    cfg.sentinel_spacing = 6.0e-3;
    cfg.sentinel_margin = 1.5e-3;
    cfg.max_subpatch_depth = 4;
    cfg.min_dense_points_per_subpatch = 16;
    cfg.max_reduced_points_per_patch = 4;
    cfg.warm_start_match_radius = 2.0e-3;
    cfg.temporal_load_regularization = 1.0e-10;
    cfg.temporal_reference_blend = 0.0;
    cfg.max_wrench_error = 0.08;
    cfg.max_cop_error = 2.5e-3;
    cfg.max_gap_error = 0.25;
    cfg.predictive_gap = true;
    return cfg;
}

ScenarioDefinition MakeFlatSquareScenario() {
    ScenarioDefinition s;
    s.name = "flat_square";
    s.description = "single square patch against a plane";
    s.cfg = MakeValidationConfig();
    s.cfg.patch_radius = 0.06;
    s.cfg.max_patch_diameter = 0.08;
    s.cfg.max_subpatch_diameter = 0.05;
    s.cfg.max_plane_error = 1.0e-4;
    s.samples = BuildRectSamples(chrono::ChVector3d(0.0, -2.0e-3, 0.0),
                                 chrono::ChVector3d(1.0, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 4.0e-2, 4.0e-2, 24, 24);
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    return s;
}

ScenarioDefinition MakeTiltedSquareScenario() {
    ScenarioDefinition s;
    s.name = "tilted_square";
    s.description = "tilted plane with predictive gap active";
    s.cfg = MakeValidationConfig();
    s.cfg.patch_radius = 0.07;
    s.cfg.max_patch_diameter = 0.09;
    s.cfg.max_subpatch_diameter = 0.06;
    s.cfg.max_plane_error = 3.0e-4;
    const chrono::ChVector3d n_W = Normalized(chrono::ChVector3d(0.35, 1.0, 0.15));
    chrono::ChVector3d t1_W;
    chrono::ChVector3d t2_W;
    BuildBasis(n_W, t1_W, t2_W);
    s.samples = BuildRectSamples(-1.5e-3 * n_W, t1_W, t2_W, 5.0e-2, 5.0e-2, 24, 24);
    s.slave_state.v_com_W = -0.2 * n_W;
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), n_W);
    return s;
}

ScenarioDefinition MakeLongStripScenario() {
    ScenarioDefinition s;
    s.name = "long_strip";
    s.description = "elongated strip patch against a plane";
    s.cfg = MakeValidationConfig();
    s.cfg.patch_radius = 0.08;
    s.cfg.max_patch_diameter = 0.14;
    s.cfg.max_subpatch_diameter = 3.5e-2;
    s.cfg.max_plane_error = 1.0e-4;
    s.samples = BuildRectSamples(chrono::ChVector3d(0.0, -1.5e-3, 0.0),
                                 chrono::ChVector3d(1.0, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 1.2e-1, 1.2e-2, 60, 8);
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    return s;
}

ScenarioDefinition MakeDualPatchScenario() {
    ScenarioDefinition s;
    s.name = "dual_patch";
    s.description = "two disconnected square patches";
    s.cfg = MakeValidationConfig();
    s.cfg.patch_radius = 0.025;
    s.cfg.max_patch_diameter = 0.05;
    s.cfg.max_subpatch_diameter = 0.03;
    s.cfg.max_plane_error = 1.0e-4;
    s.expected_patch_count = 2;

    auto left = BuildRectSamples(chrono::ChVector3d(-4.0e-2, -1.5e-3, 0.0),
                                 chrono::ChVector3d(1.0, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 2.5e-2, 2.5e-2, 18, 18);
    auto right = BuildRectSamples(chrono::ChVector3d(4.0e-2, -1.5e-3, 0.0),
                                  chrono::ChVector3d(1.0, 0.0, 0.0),
                                  chrono::ChVector3d(0.0, 0.0, 1.0),
                                  2.5e-2, 2.5e-2, 18, 18);
    s.samples = left;
    s.samples.insert(s.samples.end(), right.begin(), right.end());
    s.sdf = std::make_unique<PlaneSDF>(chrono::ChVector3d(0.0, 0.0, 0.0), chrono::ChVector3d(0.0, 1.0, 0.0));
    return s;
}

ScenarioDefinition MakeCylinderBandScenario() {
    ScenarioDefinition s;
    s.name = "cylinder_band";
    s.description = "curved contact band against a cylinder";
    s.cfg = MakeValidationConfig();
    s.cfg.patch_radius = 0.06;
    s.cfg.normal_cos_min = 0.92;
    s.cfg.max_patch_diameter = 0.09;
    s.cfg.max_subpatch_diameter = 2.5e-2;
    s.cfg.max_plane_error = 1.0e-3;
    s.cfg.max_wrench_error = 0.12;
    s.cfg.max_cop_error = 3.5e-3;
    s.samples = BuildRectSamples(chrono::ChVector3d(4.8e-2, 0.0, 0.0),
                                 chrono::ChVector3d(0.0, 1.0, 0.0),
                                 chrono::ChVector3d(0.0, 0.0, 1.0),
                                 2.0e-2, 8.0e-2, 24, 40);
    s.sdf = std::make_unique<CylinderSDF>(5.0e-2);
    return s;
}

std::vector<ScenarioDefinition> BuildScenarios() {
    std::vector<ScenarioDefinition> scenarios;
    scenarios.push_back(MakeFlatSquareScenario());
    scenarios.push_back(MakeTiltedSquareScenario());
    scenarios.push_back(MakeLongStripScenario());
    scenarios.push_back(MakeDualPatchScenario());
    scenarios.push_back(MakeCylinderBandScenario());
    return scenarios;
}

bool IsPass(const ScenarioDefinition& scenario, const CompressionValidationReport& report) {
    return !report.dense_contacts.empty() && !report.reduced_contacts.empty() &&
           report.stats.patch_count == scenario.expected_patch_count &&
           report.stats.epsilon_F <= scenario.cfg.max_wrench_error &&
           report.stats.epsilon_M <= scenario.cfg.max_wrench_error &&
           report.stats.epsilon_CoP <= scenario.cfg.max_cop_error &&
           report.stats.epsilon_gap <= scenario.cfg.max_gap_error;
}

void WriteSummaryCsv(const std::string& path, const std::vector<ScenarioResult>& results) {
    if (path.empty()) {
        return;
    }

    const fs::path out_path(path);
    if (out_path.has_parent_path()) {
        fs::create_directories(out_path.parent_path());
    }

    std::ofstream out(path);
    out << "scenario,description,total_samples,candidate_count,dense_count,reduced_count,compression_ratio,"
           "patch_count,subpatch_count,expected_patch_count,bvh_nodes_visited,bvh_nodes_pruned_obb,bvh_nodes_pruned_sdf,"
           "bvh_leaf_samples_tested,epsilon_F,epsilon_M,epsilon_CoP,epsilon_gap,max_subpatch_plane_error,"
           "max_subpatch_gap_error,dense_worst_gap,reduced_worst_gap,max_subpatch_force_residual,max_subpatch_moment_residual,"
           "dense_force_norm,reduced_force_norm,dense_moment_norm,reduced_moment_norm,pass\n";
    for (const auto& result : results) {
        const double dense_count = static_cast<double>(result.report.stats.dense_count);
        const double reduced_count = static_cast<double>(result.report.stats.reduced_count);
        const double ratio = (dense_count > 0.0) ? (reduced_count / dense_count) : 0.0;
        out << result.definition.name << ','
            << '"' << result.definition.description << '"' << ','
            << result.report.stats.total_samples << ','
            << result.report.stats.candidate_count << ','
            << result.report.stats.dense_count << ','
            << result.report.stats.reduced_count << ','
            << ratio << ','
            << result.report.stats.patch_count << ','
            << result.report.stats.subpatch_count << ','
            << result.definition.expected_patch_count << ','
            << result.report.stats.bvh_nodes_visited << ','
            << result.report.stats.bvh_nodes_pruned_obb << ','
            << result.report.stats.bvh_nodes_pruned_sdf << ','
            << result.report.stats.bvh_leaf_samples_tested << ','
            << result.report.stats.epsilon_F << ','
            << result.report.stats.epsilon_M << ','
            << result.report.stats.epsilon_CoP << ','
            << result.report.stats.epsilon_gap << ','
            << result.report.stats.max_subpatch_plane_error << ','
            << result.report.stats.max_subpatch_gap_error << ','
            << result.report.stats.dense_worst_gap << ','
            << result.report.stats.reduced_worst_gap << ','
            << result.report.stats.max_subpatch_force_residual << ','
            << result.report.stats.max_subpatch_moment_residual << ','
            << result.report.dense_wrench.force_W.Length() << ','
            << result.report.reduced_wrench.force_W.Length() << ','
            << result.report.dense_wrench.moment_at_origin_W.Length() << ','
            << result.report.reduced_wrench.moment_at_origin_W.Length() << ','
            << (result.pass ? 1 : 0) << '\n';
    }
}

void PrintScenarioResult(const ScenarioResult& result) {
    const double dense_count = static_cast<double>(result.report.stats.dense_count);
    const double reduced_count = static_cast<double>(result.report.stats.reduced_count);
    const double ratio = (dense_count > 0.0) ? (reduced_count / dense_count) : 0.0;

    std::cout << result.definition.name << " : " << (result.pass ? "PASS" : "FAIL") << '\n';
    std::cout << "  desc            : " << result.definition.description << '\n';
    std::cout << "  dense query     : total=" << result.report.stats.total_samples
              << " candidates=" << result.report.stats.candidate_count
              << " bvhVisited=" << result.report.stats.bvh_nodes_visited
              << " prunedObb=" << result.report.stats.bvh_nodes_pruned_obb
              << " prunedSdf=" << result.report.stats.bvh_nodes_pruned_sdf << '\n';
    std::cout << "  contacts        : dense=" << result.report.stats.dense_count
              << " reduced=" << result.report.stats.reduced_count
              << " ratio=" << std::fixed << std::setprecision(4) << ratio << '\n';
    std::cout << "  patches         : actual=" << result.report.stats.patch_count
              << " expected=" << result.definition.expected_patch_count
              << " subpatches=" << result.report.stats.subpatch_count << '\n';
    std::cout << "  errors          : epsF=" << result.report.stats.epsilon_F
              << " epsM=" << result.report.stats.epsilon_M
              << " epsCoP=" << result.report.stats.epsilon_CoP
              << " epsGap=" << result.report.stats.epsilon_gap
              << " plane=" << result.report.stats.max_subpatch_plane_error
              << " sentGap=" << result.report.stats.max_subpatch_gap_error
              << " denseGap=" << result.report.stats.dense_worst_gap
              << " redGap=" << result.report.stats.reduced_worst_gap
              << " subF=" << result.report.stats.max_subpatch_force_residual
              << " subM=" << result.report.stats.max_subpatch_moment_residual << '\n';
    std::cout << "  dense wrench    : |F|=" << result.report.dense_wrench.force_W.Length()
              << " |M0|=" << result.report.dense_wrench.moment_at_origin_W.Length()
              << " CoP=(" << result.report.dense_wrench.cop_W.x() << ", "
              << result.report.dense_wrench.cop_W.y() << ", "
              << result.report.dense_wrench.cop_W.z() << ")\n";
    std::cout << "  reduced wrench  : |F|=" << result.report.reduced_wrench.force_W.Length()
              << " |M0|=" << result.report.reduced_wrench.moment_at_origin_W.Length()
              << " CoP=(" << result.report.reduced_wrench.cop_W.x() << ", "
              << result.report.reduced_wrench.cop_W.y() << ", "
              << result.report.reduced_wrench.cop_W.z() << ")\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    std::string csv_path = "data/generated/compressed_contact_validation.csv";
    std::string scenario_filter = "all";

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--csv" && i + 1 < argc) {
            csv_path = argv[++i];
        } else if (arg == "--scenario" && i + 1 < argc) {
            scenario_filter = argv[++i];
        } else if (arg == "--help") {
            std::cout << "Usage: " << argv[0]
                      << " [--scenario all|flat_square|tilted_square|long_strip|dual_patch|cylinder_band]"
                      << " [--csv output.csv]\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << '\n';
            return 1;
        }
    }

    auto scenarios = BuildScenarios();
    std::vector<ScenarioResult> results;

    for (auto& scenario : scenarios) {
        if (scenario_filter != "all" && scenario_filter != scenario.name) {
            continue;
        }

        ScenarioResult result;
        result.definition = std::move(scenario);
        CompressedContactValidation::Validate(result.definition.cfg,
                                              result.definition.samples,
                                              result.definition.master_state,
                                              result.definition.slave_state,
                                              *result.definition.sdf,
                                              result.definition.mu_default,
                                              result.definition.step_size,
                                              result.report);
        result.pass = IsPass(result.definition, result.report);
        PrintScenarioResult(result);
        results.push_back(std::move(result));
    }

    if (results.empty()) {
        std::cerr << "No scenarios matched filter: " << scenario_filter << '\n';
        return 1;
    }

    WriteSummaryCsv(csv_path, results);
    std::cout << "summary_csv : " << csv_path << '\n';

    for (const auto& result : results) {
        if (!result.pass) {
            return 2;
        }
    }
    return 0;
}
