#pragma once

#include <algorithm>
#include <cmath>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace gqh {

struct Parameters {
    double gmax = 1.0;
    double tau_max = 1.0;
    double theta1 = 0.0;
    double theta2 = 0.0;
    double theta3 = 1.0;
    double theta4 = 1.0;
    double theta5 = 1.0;
    double mrdf_p1 = 1.0;
    double mrdf_p2 = 0.0;
    double mrdf_p3 = 0.0;
};

inline double clamp01(const double value) {
    return std::max(0.0, std::min(1.0, value));
}

inline double gqh_reduction(const double gamma_abs, const Parameters &p) {
    if (gamma_abs <= 0.0)
        return 1.0;

    const double gamma_ref = p.tau_max / p.gmax;
    const double r = gamma_abs / gamma_ref;
    const double r_pow = std::pow(r, p.theta5); // theta5 == 0 is valid: r^0 == 1.
    const double theta_tau = p.theta1 + p.theta2 *
        (p.theta4 * r_pow / (std::pow(p.theta3, p.theta5) + p.theta4 * r_pow));
    const double radicand = std::max(
        std::pow(1.0 + r, 2.0) - 4.0 * theta_tau * r, 0.0);
    const double ratio = std::abs(theta_tau) < 1.0e-10
        ? r / (1.0 + r)
        : ((1.0 + r) - std::sqrt(radicand)) / (2.0 * theta_tau);
    const double tau_abs = p.tau_max * clamp01(ratio);
    const double denominator = std::max(p.gmax * gamma_abs, 1.0e-12);
    return tau_abs / denominator;
}

inline double gqh_backbone(const double gamma, const Parameters &p) {
    if (gamma == 0.0)
        return 0.0;
    return std::copysign(
        p.tau_max * clamp01(gqh_reduction(std::abs(gamma), p) *
                            p.gmax * std::abs(gamma) / p.tau_max),
        gamma);
}

inline double mrdf_factor(const double gamma_m, const Parameters &p) {
    if (gamma_m <= 0.0)
        return 1.0;
    const double g_ratio = clamp01(
        gqh_backbone(gamma_m, p) / std::max(p.gmax * gamma_m, 1.0e-30));
    const double factor = p.mrdf_p1 - p.mrdf_p2 *
        std::pow(1.0 - g_ratio, p.mrdf_p3);
    return std::max(0.0, std::min(1.5, factor));
}

struct ReversalPoint {
    double gamma = 0.0;
    double tau = 0.0;
    double gamma_m = 1.0e-12;
};

// Event-resolved translated-branch Masing state with UIUC MRDF scaling.
class MasingState {
public:
    explicit MasingState(const Parameters &params) : params_(params) {}

    double update(const double gamma_target) {
        MasingState next = *this;
        const double tau = next.preview(gamma_target);
        *this = next;
        return tau;
    }

    double currentGamma() const { return gamma_prev_; }
    double currentTau() const { return tau_prev_; }
    int direction() const { return direction_; }
    std::size_t branchDepth() const { return stack_.size(); }

    void save(std::ostream &out) const {
        out.write(reinterpret_cast<const char *>(&gamma_prev_), sizeof(gamma_prev_));
        out.write(reinterpret_cast<const char *>(&tau_prev_), sizeof(tau_prev_));
        out.write(reinterpret_cast<const char *>(&gamma_max_seen_), sizeof(gamma_max_seen_));
        out.write(reinterpret_cast<const char *>(&gamma_min_seen_), sizeof(gamma_min_seen_));
        out.write(reinterpret_cast<const char *>(&last_extreme_gamma_), sizeof(last_extreme_gamma_));
        out.write(reinterpret_cast<const char *>(&direction_), sizeof(direction_));
        out.write(reinterpret_cast<const char *>(&initialized_), sizeof(initialized_));
        const std::size_t count = stack_.size();
        out.write(reinterpret_cast<const char *>(&count), sizeof(count));
        if (count)
            out.write(reinterpret_cast<const char *>(stack_.data()),
                      static_cast<std::streamsize>(count * sizeof(ReversalPoint)));
    }

    void restore(std::istream &in) {
        in.read(reinterpret_cast<char *>(&gamma_prev_), sizeof(gamma_prev_));
        in.read(reinterpret_cast<char *>(&tau_prev_), sizeof(tau_prev_));
        in.read(reinterpret_cast<char *>(&gamma_max_seen_), sizeof(gamma_max_seen_));
        in.read(reinterpret_cast<char *>(&gamma_min_seen_), sizeof(gamma_min_seen_));
        in.read(reinterpret_cast<char *>(&last_extreme_gamma_), sizeof(last_extreme_gamma_));
        in.read(reinterpret_cast<char *>(&direction_), sizeof(direction_));
        in.read(reinterpret_cast<char *>(&initialized_), sizeof(initialized_));
        std::size_t count = 0;
        in.read(reinterpret_cast<char *>(&count), sizeof(count));
        if (count > 100000)
            throw std::runtime_error("GQ/H-MRDF restore: reversal stack is invalid.");
        stack_.resize(count);
        if (count)
            in.read(reinterpret_cast<char *>(stack_.data()),
                    static_cast<std::streamsize>(count * sizeof(ReversalPoint)));
    }

private:
    double branchTau(const double gamma) const {
        const ReversalPoint &tip = stack_.back();
        const double delta = (gamma - tip.gamma) / 2.0;
        const double factor = mrdf_factor(tip.gamma_m, params_);
        if (factor == 1.0)
            return tip.tau + 2.0 * gqh_backbone(delta, params_);
        const double tau_b_m = gqh_backbone(tip.gamma_m, params_);
        const double secant_term = tau_b_m * (gamma - tip.gamma) /
                                   std::max(tip.gamma_m, 1.0e-12);
        return tip.tau + factor *
                   (2.0 * gqh_backbone(delta, params_) - secant_term) +
               secant_term;
    }

    double tauAt(const double gamma) const {
        return stack_.empty() ? gqh_backbone(gamma, params_) : branchTau(gamma);
    }

    double branchGammaM() const {
        return std::max(gamma_max_seen_, std::max(std::abs(gamma_prev_), 1.0e-12));
    }

    void refreshBranchGammaM(const double gamma) {
        if (stack_.empty())
            return;
        ReversalPoint &tip = stack_.back();
        tip.gamma_m = std::max(tip.gamma_m, std::max(gamma_max_seen_, std::abs(gamma)));
        tip.gamma_m = std::max(tip.gamma_m, 1.0e-12);
    }

    double preview(const double gamma_target) {
        const double gamma0 = gamma_prev_;
        if (!initialized_) {
            direction_ = gamma_target >= gamma0 ? 1 : -1;
            initialized_ = true;
        }

        int new_direction = 0;
        if (gamma_target > gamma_prev_ + 1.0e-12)
            new_direction = 1;
        else if (gamma_target < gamma_prev_ - 1.0e-12)
            new_direction = -1;
        if (new_direction == 0)
            return tauAt(gamma_prev_);

        if (direction_ != 0 && new_direction != direction_) {
            const double gamma_m = branchGammaM();
            stack_.push_back({gamma_prev_, tau_prev_, gamma_m});
            last_extreme_gamma_ = gamma_prev_;
        }
        direction_ = new_direction;

        while (stack_.size() >= 2) {
            const double previous_turn = stack_[stack_.size() - 2].gamma;
            const double guard = static_cast<double>(new_direction) *
                                     (gamma_target - previous_turn) - 1.0e-12;
            if (guard <= 0.0)
                break;
            stack_.pop_back();
            stack_.pop_back();
        }

        if (stack_.size() == 1) {
            const double tip_gamma = stack_.back().gamma;
            const double guard = static_cast<double>(new_direction) *
                                     (gamma_target + tip_gamma) - 1.0e-12;
            if (guard > 0.0)
                stack_.pop_back();
        }

        refreshBranchGammaM(gamma_target);
        const double tau = tauAt(gamma_target);
        gamma_prev_ = gamma_target;
        tau_prev_ = tau;
        gamma_max_seen_ = std::max(gamma_max_seen_, std::abs(gamma_target));
        gamma_min_seen_ = std::min(gamma_min_seen_, gamma_target);
        return tau;
    }

    Parameters params_;
    double gamma_prev_ = 0.0;
    double tau_prev_ = 0.0;
    double gamma_max_seen_ = 0.0;
    double gamma_min_seen_ = 0.0;
    double last_extreme_gamma_ = 0.0;
    int direction_ = 0;
    bool initialized_ = false;
    std::vector<ReversalPoint> stack_;
};

inline double masing_damping(const Parameters &p, const double gamma_a,
                             const int n = 512) {
    if (gamma_a <= 0.0)
        return 0.0;
    const double tau_a = gqh_backbone(gamma_a, p);
    const double denominator = tau_a * gamma_a;
    if (denominator <= 0.0)
        return 0.0;

    const int nn = std::max(n, 32);
    double previous_gamma = 0.0;
    double previous_tau = 0.0;
    double area = 0.0;
    for (int i = 0; i <= nn; ++i) {
        const double gamma = i == 0
            ? 0.0
            : gamma_a * std::pow(1.0e-9, 1.0 -
                                 static_cast<double>(i - 1) / (nn - 1));
        const double tau = gqh_backbone(gamma, p);
        if (i > 0)
            area += 0.5 * (previous_tau + tau) * (gamma - previous_gamma);
        previous_gamma = gamma;
        previous_tau = tau;
    }
    return std::max(0.0, (2.0 / 3.14159265358979323846) *
                           ((2.0 * area / denominator) - 1.0));
}

} // namespace gqh
