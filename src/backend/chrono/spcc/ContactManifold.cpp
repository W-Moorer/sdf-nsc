#include "platform/backend/spcc/ContactManifold.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>

#include "platform/backend/spcc/SPCCProblem.h"

namespace platform {
namespace backend {
namespace spcc {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

chrono::ChVector3d NormalizeOrFallback(const chrono::ChVector3d& v_W,
                                       const chrono::ChVector3d& fallback_W) {
    const double len = v_W.Length();
    if (len > 1.0e-12 && std::isfinite(len)) {
        return v_W * (1.0 / len);
    }
    const double fallback_len = fallback_W.Length();
    if (fallback_len > 1.0e-12 && std::isfinite(fallback_len)) {
        return fallback_W * (1.0 / fallback_len);
    }
    return chrono::ChVector3d(1.0, 0.0, 0.0);
}

chrono::ChVector3d BuildFallbackTangent(const chrono::ChVector3d& n_W) {
    chrono::ChVector3d tangent = chrono::Vcross(n_W, chrono::ChVector3d(1.0, 0.0, 0.0));
    if (tangent.Length2() <= 1.0e-12) {
        tangent = chrono::Vcross(n_W, chrono::ChVector3d(0.0, 1.0, 0.0));
    }
    return NormalizeOrFallback(tangent, chrono::ChVector3d(0.0, 0.0, 1.0));
}

double WrapAngle(double angle) {
    while (angle <= -kPi) {
        angle += kTwoPi;
    }
    while (angle > kPi) {
        angle -= kTwoPi;
    }
    return angle;
}

double PolarAngleYZ(const chrono::ChVector3d& radial_S) {
    return std::atan2(radial_S.z(), radial_S.y());
}

double CircularMean(const std::vector<double>& angles, const std::vector<double>& weights) {
    double s = 0.0;
    double c = 0.0;
    for (std::size_t i = 0; i < angles.size(); ++i) {
        const double w = (i < weights.size()) ? weights[i] : 1.0;
        s += w * std::sin(angles[i]);
        c += w * std::cos(angles[i]);
    }
    if (std::abs(s) <= 1.0e-12 && std::abs(c) <= 1.0e-12) {
        return angles.empty() ? 0.0 : angles.front();
    }
    return std::atan2(s, c);
}

struct ArcContactParams {
    std::size_t source_index = 0;
    double theta = 0.0;
    double axial = 0.0;
    double radius = 0.0;
    double half_span = 0.0;
};

QuadratureContact BuildSeedContact(const ActiveContactSample& contact,
                                   std::size_t manifold_id,
                                   std::size_t source_contact_index);
ContactManifoldHeader BuildHeader(const ActiveContactSample& contact,
                                  ContactManifoldKind kind,
                                  std::size_t source_contact_index);
std::size_t GroupKeyForContact(const ActiveContactSample& contact, std::size_t source_contact_index);
std::size_t SelectRepresentativeIndex(const std::vector<std::size_t>& member_indices,
                                      const std::vector<ActiveContactSample>& candidates);

ContactManifold BuildArcManifoldFromGroup(const RigidBodyStateW& slave_pred,
                                          const std::vector<ActiveContactSample>& candidates,
                                          const std::vector<std::size_t>& member_indices) {
    ContactManifold manifold;
    const std::size_t rep = SelectRepresentativeIndex(member_indices, candidates);
    manifold.header = BuildHeader(candidates[rep], ContactManifoldKind::ArcSliding, rep);
    manifold.header.manifold_id = GroupKeyForContact(candidates[rep], rep);
    manifold.header.weight = 0.0;
    manifold.header.age = 0;
    manifold.header.matched = false;
    manifold.member_contact_indices = member_indices;

    chrono::ChVector3d weighted_center(0.0, 0.0, 0.0);
    chrono::ChVector3d weighted_normal(0.0, 0.0, 0.0);
    chrono::ChVector3d weighted_tangent(0.0, 0.0, 0.0);
    chrono::ChVector3d weighted_radial_W(0.0, 0.0, 0.0);
    double weight_sum = 0.0;
    double axial_sum = 0.0;
    double radius_sum = 0.0;
    double coverage_arc_span = 0.0;
    std::vector<double> theta_values;
    std::vector<double> theta_weights;
    theta_values.reserve(member_indices.size());
    theta_weights.reserve(member_indices.size());

    for (std::size_t idx : member_indices) {
        const auto& contact = candidates[idx];
        const double w = std::max(1.0, static_cast<double>(contact.cluster_size));
        weighted_center += w * contact.x_W;
        weighted_normal += w * contact.n_W;
        weighted_tangent += w * contact.coverage_tangent_W;
        weight_sum += w;
        manifold.header.weight += std::max(0.0, contact.patch_weight_sum);
        manifold.header.age = std::max(manifold.header.age, contact.active_age);
        manifold.header.matched = manifold.header.matched || contact.manifold_matched;
        manifold.seed_contacts.push_back(BuildSeedContact(contact, manifold.header.manifold_id, idx));

        chrono::ChVector3d radial_slave_S = contact.xi_slave_S;
        radial_slave_S.x() = 0.0;
        weighted_radial_W += w * (slave_pred.R_WRef * radial_slave_S);
        axial_sum += contact.xi_slave_S.x();
        radius_sum += radial_slave_S.Length();
        theta_values.push_back(PolarAngleYZ(radial_slave_S));
        theta_weights.push_back(w);

        if (radial_slave_S.Length() > 1.0e-12 && std::isfinite(contact.coverage_half_span) && contact.coverage_half_span > 0.0) {
            coverage_arc_span = std::max(coverage_arc_span, 2.0 * contact.coverage_half_span / radial_slave_S.Length());
        }
    }

    if (!(manifold.header.weight > 0.0) || !std::isfinite(manifold.header.weight)) {
        manifold.header.weight = static_cast<double>(member_indices.size());
    }

    const chrono::ChVector3d center_W =
        (weight_sum > 0.0) ? (weighted_center * (1.0 / weight_sum)) : manifold.header.x_W;
    const chrono::ChVector3d normal_W =
        NormalizeOrFallback(weighted_normal, manifold.header.n_W);
    const chrono::ChVector3d tangent_hint_W =
        NormalizeOrFallback(weighted_tangent, BuildFallbackTangent(normal_W));
    const chrono::ChVector3d axis_W =
        NormalizeOrFallback(slave_pred.R_WRef * chrono::ChVector3d(1.0, 0.0, 0.0),
                            chrono::ChVector3d(1.0, 0.0, 0.0));
    const chrono::ChVector3d radial_center_W =
        NormalizeOrFallback(weighted_radial_W, normal_W);
    const chrono::ChVector3d tangent_arc_W =
        NormalizeOrFallback(chrono::Vcross(axis_W, radial_center_W), tangent_hint_W);
    const chrono::ChVector3d tangent2_W =
        NormalizeOrFallback(chrono::Vcross(normal_W, tangent_arc_W), BuildFallbackTangent(normal_W));

    manifold.header.x_W = center_W;
    manifold.header.n_W = normal_W;
    manifold.sliding_patch.cluster_size = static_cast<int>(member_indices.size());
    manifold.sliding_patch.coverage_center_W = center_W;
    manifold.sliding_patch.coverage_radius =
        member_indices.empty() ? 0.0 : radius_sum / static_cast<double>(member_indices.size());
    manifold.sliding_patch.tangent1_W = tangent_arc_W;
    manifold.sliding_patch.tangent2_W = tangent2_W;

    const double theta_center = CircularMean(theta_values, theta_weights);
    double theta_half_span = 0.0;
    for (double theta : theta_values) {
        theta_half_span = std::max(theta_half_span, std::abs(WrapAngle(theta - theta_center)));
    }
    theta_half_span = std::max(theta_half_span, 0.5 * coverage_arc_span);

    double axial_center = 0.0;
    if (!member_indices.empty()) {
        axial_center = axial_sum / static_cast<double>(member_indices.size());
    }
    double axial_half_span = 0.0;
    for (std::size_t idx : member_indices) {
        axial_half_span = std::max(axial_half_span, std::abs(candidates[idx].xi_slave_S.x() - axial_center));
    }

    manifold.arc_sliding.cluster_size = static_cast<int>(member_indices.size());
    manifold.arc_sliding.axis_W = axis_W;
    manifold.arc_sliding.radial_center_W = radial_center_W;
    manifold.arc_sliding.tangent_arc_W = tangent_arc_W;
    manifold.arc_sliding.theta_center = theta_center;
    manifold.arc_sliding.theta_span = 2.0 * theta_half_span;
    manifold.arc_sliding.axial_center = axial_center;
    manifold.arc_sliding.axial_span = 2.0 * axial_half_span;
    manifold.arc_sliding.cylinder_radius =
        member_indices.empty() ? 0.0 : radius_sum / static_cast<double>(member_indices.size());

    return manifold;
}

QuadratureContact BuildSeedContact(const ActiveContactSample& contact,
                                   std::size_t manifold_id,
                                   std::size_t source_contact_index) {
    QuadratureContact out;
    out.manifold_id = manifold_id;
    out.source_contact_index = source_contact_index;
    out.quadrature_id = 0;
    out.x_W = contact.x_W;
    out.n_W = contact.n_W;
    out.phi = contact.phi;
    out.phi_eff = contact.phi_eff;
    out.weight = contact.patch_weight_sum;
    out.hessian_W = contact.hessian_W;
    out.curvature_tangential_only = contact.curvature_tangential_only;
    out.mu = contact.mu;
    return out;
}

ContactManifoldHeader BuildHeader(const ActiveContactSample& contact,
                                  ContactManifoldKind kind,
                                  std::size_t source_contact_index) {
    ContactManifoldHeader header;
    header.kind = kind;
    header.manifold_id = contact.manifold_id != 0 ? contact.manifold_id : (source_contact_index + 1);
    header.source_contact_index = source_contact_index;
    header.age = contact.active_age;
    header.matched = contact.manifold_matched;
    header.x_W = contact.x_W;
    header.n_W = contact.n_W;
    header.P_W = contact.P_W;
    header.x_master_M = contact.x_master_M;
    header.phi = contact.phi;
    header.phi_eff = contact.phi_eff;
    header.v_rel_W = contact.v_rel_W;
    header.u_n_pred = contact.u_n_pred;
    header.u_tau_pred_W = contact.u_tau_pred_W;
    header.hessian_W = contact.hessian_W;
    header.curvature_term = contact.curvature_term;
    header.curvature_gate = contact.curvature_gate;
    header.curvature_tangential_only = contact.curvature_tangential_only;
    header.curvature_term_abs_max = contact.curvature_term_abs_max;
    header.curvature_term_ratio_max = contact.curvature_term_ratio_max;
    header.curvature_gap_floor = contact.curvature_gap_floor;
    header.mu = contact.mu;
    header.weight = contact.patch_weight_sum;
    return header;
}

std::size_t GroupKeyForContact(const ActiveContactSample& contact, std::size_t source_contact_index) {
    return contact.manifold_id != 0 ? contact.manifold_id : (source_contact_index + 1);
}

std::size_t SelectRepresentativeIndex(const std::vector<std::size_t>& member_indices,
                                      const std::vector<ActiveContactSample>& candidates) {
    std::size_t best = member_indices.front();
    auto best_key = std::make_pair(candidates[best].phi_eff, candidates[best].phi);
    for (std::size_t idx : member_indices) {
        const auto key = std::make_pair(candidates[idx].phi_eff, candidates[idx].phi);
        if (key < best_key) {
            best = idx;
            best_key = key;
        }
    }
    return best;
}

class GroupedBuilderBase : public IContactManifoldBuilder {
  public:
    explicit GroupedBuilderBase(ContactManifoldKind kind) : kind_(kind) {}

    void BuildManifolds(const RigidBodyStateW& master_pred,
                        const RigidBodyStateW& slave_pred,
                        const SDFField& sdf,
                        const std::vector<ActiveContactSample>& candidates,
                        double step_size,
                        std::vector<ContactManifold>& out_manifolds) override {
        (void)master_pred;
        (void)sdf;
        (void)step_size;
        out_manifolds.clear();
        if (candidates.empty()) {
            return;
        }

        std::unordered_map<std::size_t, std::size_t> manifold_lookup;
        manifold_lookup.reserve(candidates.size());

        for (std::size_t i = 0; i < candidates.size(); ++i) {
            const auto key = GroupKeyForContact(candidates[i], i);
            auto [it, inserted] = manifold_lookup.emplace(key, out_manifolds.size());
            if (inserted) {
                ContactManifold manifold;
                manifold.header = BuildHeader(candidates[i], kind_, i);
                manifold.header.manifold_id = key;
                out_manifolds.push_back(manifold);
            }

            auto& manifold = out_manifolds[it->second];
            manifold.member_contact_indices.push_back(i);
            manifold.seed_contacts.push_back(BuildSeedContact(candidates[i], key, i));
            AbsorbContactIntoHeader(candidates[i], manifold.header);
        }

        for (auto& manifold : out_manifolds) {
            const std::size_t rep = SelectRepresentativeIndex(manifold.member_contact_indices, candidates);
            manifold.header = BuildHeader(candidates[rep], kind_, rep);
            manifold.header.manifold_id = GroupKeyForContact(candidates[rep], rep);
            manifold.header.weight = 0.0;
            manifold.header.age = 0;
            manifold.header.matched = false;
            for (std::size_t idx : manifold.member_contact_indices) {
                manifold.header.weight += std::max(0.0, candidates[idx].patch_weight_sum);
                manifold.header.age = std::max(manifold.header.age, candidates[idx].active_age);
                manifold.header.matched = manifold.header.matched || candidates[idx].manifold_matched;
            }
            if (!(manifold.header.weight > 0.0) || !std::isfinite(manifold.header.weight)) {
                manifold.header.weight = static_cast<double>(manifold.member_contact_indices.size());
            }
            BuildSpecificState(slave_pred, candidates, manifold);
        }
    }

  protected:
    virtual void BuildSpecificState(const RigidBodyStateW& slave_pred,
                                    const std::vector<ActiveContactSample>& candidates,
                                    ContactManifold& manifold) const = 0;

  private:
    static void AbsorbContactIntoHeader(const ActiveContactSample& contact, ContactManifoldHeader& header) {
        header.age = std::max(header.age, contact.active_age);
        header.matched = header.matched || contact.manifold_matched;
    }

    ContactManifoldKind kind_;
};

class ImpactPairBuilder final : public GroupedBuilderBase {
  public:
    ImpactPairBuilder() : GroupedBuilderBase(ContactManifoldKind::ImpactPair) {}

  protected:
    void BuildSpecificState(const RigidBodyStateW&,
                            const std::vector<ActiveContactSample>&,
                            ContactManifold& manifold) const override {
        manifold.impact_pair.impact_normal_W = manifold.header.n_W;
    }
};

class CompactPointBuilder final : public GroupedBuilderBase {
  public:
    CompactPointBuilder() : GroupedBuilderBase(ContactManifoldKind::CompactPoint) {}

  protected:
    void BuildSpecificState(const RigidBodyStateW&,
                            const std::vector<ActiveContactSample>& candidates,
                            ContactManifold& manifold) const override {
        manifold.compact_point.cluster_size = static_cast<int>(manifold.member_contact_indices.size());
        const auto rep = manifold.header.source_contact_index;
        manifold.compact_point.deepest_x_W = candidates[rep].x_W;
    }
};

class SlidingPatchBuilder final : public GroupedBuilderBase {
  public:
    explicit SlidingPatchBuilder(ContactManifoldKind kind) : GroupedBuilderBase(kind) {}

  protected:
    void BuildSpecificState(const RigidBodyStateW& slave_pred,
                            const std::vector<ActiveContactSample>& candidates,
                            ContactManifold& manifold) const override {
        const bool is_arc = manifold.header.kind == ContactManifoldKind::ArcSliding;
        chrono::ChVector3d weighted_center(0.0, 0.0, 0.0);
        chrono::ChVector3d weighted_normal(0.0, 0.0, 0.0);
        chrono::ChVector3d weighted_tangent(0.0, 0.0, 0.0);
        double weight_sum = 0.0;
        double max_radius = 0.0;
        chrono::ChVector3d radial_sum(0.0, 0.0, 0.0);
        double radial_radius_sum = 0.0;
        double axial_sum = 0.0;

        for (std::size_t idx : manifold.member_contact_indices) {
            const auto& contact = candidates[idx];
            const double w = std::max(1.0, static_cast<double>(contact.cluster_size));
            weighted_center += w * contact.x_W;
            weighted_normal += w * contact.n_W;
            weighted_tangent += w * contact.coverage_tangent_W;
            weight_sum += w;
            max_radius = std::max(max_radius, contact.coverage_half_span);

            if (is_arc) {
                chrono::ChVector3d radial_slave_S = contact.xi_slave_S;
                radial_slave_S.x() = 0.0;
                radial_sum += slave_pred.R_WRef * radial_slave_S;
                radial_radius_sum += radial_slave_S.Length();
                axial_sum += contact.xi_slave_S.x();
            }
        }

        const chrono::ChVector3d center_W =
            (weight_sum > 0.0) ? (weighted_center * (1.0 / weight_sum)) : manifold.header.x_W;
        const chrono::ChVector3d normal_W =
            NormalizeOrFallback(weighted_normal, manifold.header.n_W);
        const chrono::ChVector3d tangent1_W =
            NormalizeOrFallback(weighted_tangent, BuildFallbackTangent(normal_W));
        const chrono::ChVector3d tangent2_W =
            NormalizeOrFallback(chrono::Vcross(normal_W, tangent1_W), BuildFallbackTangent(normal_W));

        manifold.header.x_W = center_W;
        manifold.header.n_W = normal_W;

        manifold.sliding_patch.cluster_size = static_cast<int>(manifold.member_contact_indices.size());
        manifold.sliding_patch.coverage_center_W = center_W;
        manifold.sliding_patch.coverage_radius = max_radius;
        manifold.sliding_patch.tangent1_W = tangent1_W;
        manifold.sliding_patch.tangent2_W = tangent2_W;

        if (is_arc) {
            manifold.arc_sliding.cluster_size = manifold.sliding_patch.cluster_size;
            manifold.arc_sliding.axis_W =
                NormalizeOrFallback(slave_pred.R_WRef * chrono::ChVector3d(1.0, 0.0, 0.0),
                                    chrono::ChVector3d(1.0, 0.0, 0.0));
            manifold.arc_sliding.radial_center_W =
                NormalizeOrFallback(radial_sum, normal_W);
            manifold.arc_sliding.tangent_arc_W =
                NormalizeOrFallback(chrono::Vcross(manifold.arc_sliding.axis_W, manifold.arc_sliding.radial_center_W),
                                    tangent1_W);
            manifold.arc_sliding.axial_center =
                manifold.member_contact_indices.empty() ? 0.0 : axial_sum / static_cast<double>(manifold.member_contact_indices.size());
            manifold.arc_sliding.axial_span = 0.0;
            manifold.arc_sliding.theta_center = 0.0;
            manifold.arc_sliding.theta_span = 2.0 * max_radius;
            manifold.arc_sliding.cylinder_radius =
                manifold.member_contact_indices.empty() ? 0.0
                                                        : radial_radius_sum / static_cast<double>(manifold.member_contact_indices.size());
        }
    }
};

class ArcSlidingBuilder final : public IContactManifoldBuilder {
  public:
    void BuildManifolds(const RigidBodyStateW&,
                        const RigidBodyStateW& slave_pred,
                        const SDFField&,
                        const std::vector<ActiveContactSample>& candidates,
                        double,
                        std::vector<ContactManifold>& out_manifolds) override {
        out_manifolds.clear();
        if (candidates.empty()) {
            return;
        }

        std::vector<ArcContactParams> params;
        params.reserve(candidates.size());
        double max_local_span = 0.0;
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            ArcContactParams p;
            p.source_index = i;
            chrono::ChVector3d radial_slave_S = candidates[i].xi_slave_S;
            radial_slave_S.x() = 0.0;
            p.theta = PolarAngleYZ(radial_slave_S);
            p.axial = candidates[i].xi_slave_S.x();
            p.radius = radial_slave_S.Length();
            p.half_span = candidates[i].coverage_half_span;
            params.push_back(p);
            if (p.radius > 1.0e-12 && std::isfinite(p.half_span) && p.half_span > 0.0) {
                max_local_span = std::max(max_local_span, 2.0 * p.half_span / p.radius);
            }
        }

        std::sort(params.begin(), params.end(), [](const ArcContactParams& a, const ArcContactParams& b) {
            return a.theta < b.theta;
        });

        const double theta_gap_tol = std::max(0.35, max_local_span * 1.5);
        std::vector<std::vector<std::size_t>> groups;
        std::vector<std::size_t> current_group;
        current_group.push_back(params.front().source_index);
        for (std::size_t i = 1; i < params.size(); ++i) {
            const double gap = WrapAngle(params[i].theta - params[i - 1].theta);
            if (std::abs(gap) <= theta_gap_tol) {
                current_group.push_back(params[i].source_index);
            } else {
                groups.push_back(current_group);
                current_group.clear();
                current_group.push_back(params[i].source_index);
            }
        }
        groups.push_back(current_group);

        if (groups.size() > 1) {
            const double wrap_gap = WrapAngle((params.front().theta + kTwoPi) - params.back().theta);
            if (std::abs(wrap_gap) <= theta_gap_tol) {
                groups.front().insert(groups.front().begin(), groups.back().begin(), groups.back().end());
                groups.pop_back();
            }
        }

        out_manifolds.reserve(groups.size());
        for (const auto& group : groups) {
            out_manifolds.push_back(BuildArcManifoldFromGroup(slave_pred, candidates, group));
        }
    }
};

class SeedContactQuadrature final : public IManifoldQuadrature {
  public:
    void Discretize(const ContactManifold& manifold,
                    std::vector<QuadratureContact>& out_contacts) const override {
        out_contacts.insert(out_contacts.end(), manifold.seed_contacts.begin(), manifold.seed_contacts.end());
    }
};

class ArcSlidingQuadrature final : public IManifoldQuadrature {
  public:
    ArcSlidingQuadrature(int target_contacts, double span_scale, double min_half_span)
        : target_contacts_(std::max(1, target_contacts)),
          span_scale_(std::clamp(span_scale, 0.0, 1.0)),
          min_half_span_(std::max(0.0, min_half_span)) {}

    void Discretize(const ContactManifold& manifold,
                    std::vector<QuadratureContact>& out_contacts) const override {
        if (target_contacts_ <= 1 || manifold.seed_contacts.empty() ||
            manifold.header.kind != ContactManifoldKind::ArcSliding) {
            out_contacts.insert(out_contacts.end(), manifold.seed_contacts.begin(), manifold.seed_contacts.end());
            return;
        }

        if (static_cast<int>(manifold.seed_contacts.size()) <= target_contacts_) {
            out_contacts.insert(out_contacts.end(), manifold.seed_contacts.begin(), manifold.seed_contacts.end());
            return;
        }

        const auto tangent_W =
            NormalizeOrFallback(manifold.arc_sliding.tangent_arc_W, manifold.sliding_patch.tangent1_W);
        if (tangent_W.Length2() <= 1.0e-12) {
            out_contacts.insert(out_contacts.end(), manifold.seed_contacts.begin(), manifold.seed_contacts.end());
            return;
        }

        std::vector<std::pair<double, std::size_t>> order;
        order.reserve(manifold.seed_contacts.size());
        for (std::size_t i = 0; i < manifold.seed_contacts.size(); ++i) {
            const double coord = chrono::Vdot(manifold.seed_contacts[i].x_W - manifold.header.x_W, tangent_W);
            order.emplace_back(coord, i);
        }
        std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

        const int denom = std::max(1, target_contacts_ - 1);
        for (int qi = 0; qi < target_contacts_; ++qi) {
            const double alpha = (target_contacts_ == 1) ? 0.5 : static_cast<double>(qi) / static_cast<double>(denom);
            const std::size_t pick =
                static_cast<std::size_t>(std::round(alpha * static_cast<double>(order.size() - 1)));
            QuadratureContact contact = manifold.seed_contacts[order[pick].second];
            contact.quadrature_id = qi;
            contact.weight = manifold.header.weight / static_cast<double>(target_contacts_);
            contact.requery_geometry = false;
            out_contacts.push_back(contact);
        }
    }

  private:
    int target_contacts_ = 1;
    double span_scale_ = 1.0;
    double min_half_span_ = 0.0;
};

}  // namespace

std::unique_ptr<IContactManifoldBuilder> MakeContactManifoldBuilder(ContactManifoldKind kind) {
    switch (kind) {
        case ContactManifoldKind::ImpactPair:
            return std::make_unique<ImpactPairBuilder>();
        case ContactManifoldKind::CompactPoint:
            return std::make_unique<CompactPointBuilder>();
        case ContactManifoldKind::SlidingPatch:
            return std::make_unique<SlidingPatchBuilder>(ContactManifoldKind::SlidingPatch);
        case ContactManifoldKind::ArcSliding:
            return std::make_unique<ArcSlidingBuilder>();
    }
    return std::make_unique<CompactPointBuilder>();
}

std::unique_ptr<IManifoldQuadrature> MakeManifoldQuadrature(ContactManifoldKind kind,
                                                            bool enable_multi_contact,
                                                            int target_contacts,
                                                            double span_scale,
                                                            double min_half_span) {
    if (enable_multi_contact && kind == ContactManifoldKind::ArcSliding && target_contacts > 1) {
        return std::make_unique<ArcSlidingQuadrature>(target_contacts, span_scale, min_half_span);
    }
    return std::make_unique<SeedContactQuadrature>();
}

}  // namespace spcc
}  // namespace backend
}  // namespace platform
