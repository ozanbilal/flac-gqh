#include "refcurve.h"

#include "gqh_math.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gqh::refcurve {

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kAtmosphericPressureKpa = 101.3;
constexpr double kFitStrainMin = 1.0e-6;
constexpr double kFitStrainMax = 5.0e-4;
constexpr double kTargetStrengthRatio = 0.95;
constexpr double kTargetStrengthStrain = 1.0e-1;
constexpr int kMaxIterations = 200;
constexpr double kLogFloor = 1.0e-6;

double clamp(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(hi, value));
}

double positiveOrThrow(const char *name, const double value) {
    if (!std::isfinite(value) || value <= 0.0)
        throw std::runtime_error(std::string("GQ/H refcurve: ") + name + " must be > 0.");
    return value;
}

double nonnegativeOrThrow(const char *name, const double value) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::runtime_error(std::string("GQ/H refcurve: ") + name + " must be >= 0.");
    return value;
}

double logSpace(const int index, const double lo, const double hi) {
    const double fraction = static_cast<double>(index) /
                            static_cast<double>(kReferencePointCount - 1);
    return std::exp(std::log(lo) + fraction * (std::log(hi) - std::log(lo)));
}

double interpolateLog(const double x, const double *anchor_x, const double *anchor_y,
                      const int count) {
    const double lx = std::log10(x);
    if (lx <= std::log10(anchor_x[0]))
        return anchor_y[0];
    if (lx >= std::log10(anchor_x[count - 1]))
        return anchor_y[count - 1];
    for (int i = 1; i < count; ++i) {
        const double x0 = std::log10(anchor_x[i - 1]);
        const double x1 = std::log10(anchor_x[i]);
        if (lx <= x1) {
            const double fraction = (lx - x0) / (x1 - x0);
            return anchor_y[i - 1] + fraction * (anchor_y[i] - anchor_y[i - 1]);
        }
    }
    return anchor_y[count - 1];
}

void fillLogGrid(ReferenceCurves *curves) {
    for (int i = 0; i < kReferencePointCount; ++i)
        curves->strain[static_cast<std::size_t>(i)] = logSpace(i, 1.0e-6, 1.0e-1);
}

double modifiedHyperbolicReduction(const double gamma, const double strain_ref,
                                   const double curvature) {
    const double safe_gamma = std::max(std::abs(gamma), 1.0e-12);
    return 1.0 / (1.0 + std::pow(safe_gamma / strain_ref, curvature));
}

double modifiedHyperbolicMasingDamping(const double gamma, const double strain_ref,
                                       const double curvature, const double scaling) {
    const double strain_pct = std::max(std::abs(gamma), 1.0e-12) * 100.0;
    const double strain_ref_pct = strain_ref * 100.0;
    const double damping_masing_a1 = (100.0 / kPi) *
        (4.0 * (strain_pct - strain_ref_pct *
                       std::log((strain_pct + strain_ref_pct) / strain_ref_pct)) /
             (strain_pct * strain_pct / (strain_pct + strain_ref_pct)) - 2.0);
    const double c1 = (-1.1143 * curvature * curvature) + (1.8618 * curvature) + 0.2523;
    const double c2 = (0.0805 * curvature * curvature) - (0.0710 * curvature) - 0.0095;
    const double c3 = (-0.0005 * curvature * curvature) + (0.0002 * curvature) + 0.0003;
    const double damping_masing = (c1 * damping_masing_a1) +
        (c2 * damping_masing_a1 * damping_masing_a1) +
        (c3 * damping_masing_a1 * damping_masing_a1 * damping_masing_a1);
    const double reduction = modifiedHyperbolicReduction(gamma, strain_ref, curvature);
    const double damping = (scaling * damping_masing * std::pow(reduction, 0.1)) / 100.0;
    return std::isfinite(damping) ? clamp(damping, 0.0, 0.5) : 0.0;
}

void setAnchoredCurve(ReferenceCurves *curves, const double *anchor_gg,
                      const double *anchor_damping, const int count) {
    static constexpr double anchor_strain[] = {
        1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4,
        1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1,
    };
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        curves->modulus_reduction[index] = clamp(
            interpolateLog(curves->strain[index], anchor_strain, anchor_gg, count), 0.0, 1.0);
        curves->damping_ratio[index] = clamp(
            interpolateLog(curves->strain[index], anchor_strain, anchor_damping, count), 0.0, 0.5);
    }
}

void setAnchoredCurve17(ReferenceCurves *curves, const double *anchor_gg,
                        const double *anchor_damping) {
    static constexpr double anchor_strain[] = {
        1.0e-7, 1.0e-6, 2.0e-6, 5.0e-6, 1.0e-5, 2.0e-5,
        5.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3,
        5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 1.0e-1,
    };
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        curves->modulus_reduction[index] = clamp(
            interpolateLog(curves->strain[index], anchor_strain, anchor_gg, 17), 0.0, 1.0);
        curves->damping_ratio[index] = clamp(
            interpolateLog(curves->strain[index], anchor_strain, anchor_damping, 17), 0.0, 0.5);
    }
}

void makeDarendeli(ReferenceCurves *curves, const CalibrationInput &input) {
    const double pi = nonnegativeOrThrow("pi", input.plasticity_index);
    const double ocr = positiveOrThrow("ocr", input.ocr);
    const double stress = positiveOrThrow("sigma_m_eff", input.mean_effective_stress_kpa);
    const double strain_ref = ((0.0352 + 0.0010 * pi * std::pow(ocr, 0.3246)) *
                               std::pow(stress / kAtmosphericPressureKpa, 0.3483)) / 100.0;
    curves->damping_min = clamp(
        ((0.8005 + 0.0129 * pi * std::pow(ocr, -0.1069)) *
         std::pow(stress / kAtmosphericPressureKpa, -0.2889)) / 100.0,
        0.0, 0.5);
    const double scaling = 0.6329 - 0.00566 * std::log(10.0);
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        curves->modulus_reduction[index] = modifiedHyperbolicReduction(
            curves->strain[index], strain_ref, 0.9190);
        curves->damping_ratio[index] = clamp(
            curves->damping_min + modifiedHyperbolicMasingDamping(
                curves->strain[index], strain_ref, 0.9190, scaling),
            curves->damping_min, 0.5);
    }
    curves->name = "Darendeli";
}

void makeVuceticDobry(ReferenceCurves *curves, const double plasticity_index) {
    const double pi = std::max(0.0, plasticity_index);
    static constexpr double gg_table[4][17] = {
        {1.000, 1.000, 1.000, 0.990, 0.984, 0.916, 0.818, 0.711, 0.578, 0.381, 0.256, 0.160, 0.067, 0.027, 0.008, 0.004, 0.004},
        {1.000, 1.000, 1.000, 1.000, 0.992, 0.965, 0.898, 0.818, 0.719, 0.549, 0.408, 0.287, 0.158, 0.096, 0.055, 0.028, 0.028},
        {1.000, 1.000, 1.000, 1.000, 1.000, 0.992, 0.953, 0.898, 0.816, 0.664, 0.537, 0.416, 0.266, 0.162, 0.090, 0.045, 0.045},
        {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.982, 0.953, 0.898, 0.781, 0.676, 0.535, 0.377, 0.246, 0.135, 0.068, 0.068},
    };
    static constexpr double damping_table[4][17] = {
        {0.010, 0.012, 0.012, 0.015, 0.018, 0.025, 0.038, 0.054, 0.078, 0.120, 0.152, 0.184, 0.218, 0.239, 0.254, 0.267, 0.267},
        {0.010, 0.011, 0.012, 0.013, 0.016, 0.021, 0.032, 0.046, 0.063, 0.091, 0.116, 0.142, 0.177, 0.200, 0.221, 0.243, 0.243},
        {0.010, 0.010, 0.011, 0.012, 0.014, 0.017, 0.027, 0.037, 0.050, 0.069, 0.086, 0.108, 0.141, 0.169, 0.199, 0.226, 0.226},
        {0.009, 0.010, 0.010, 0.011, 0.013, 0.016, 0.023, 0.029, 0.037, 0.049, 0.061, 0.078, 0.109, 0.134, 0.163, 0.192, 0.192},
    };
    const double pi_keys[] = {0.0, 15.0, 30.0, 50.0};
    int lo = 0;
    int hi = 0;
    double fraction = 0.0;
    if (pi <= pi_keys[0]) {
        lo = hi = 0;
    } else if (pi >= pi_keys[3]) {
        lo = hi = 3;
    } else {
        for (int i = 1; i < 4; ++i) {
            if (pi <= pi_keys[i]) {
                lo = i - 1;
                hi = i;
                fraction = (pi - pi_keys[lo]) / (pi_keys[hi] - pi_keys[lo]);
                break;
            }
        }
    }
    double gg[17] = {};
    double damping[17] = {};
    for (int i = 0; i < 17; ++i) {
        gg[i] = gg_table[lo][i] + fraction * (gg_table[hi][i] - gg_table[lo][i]);
        damping[i] = damping_table[lo][i] + fraction *
            (damping_table[hi][i] - damping_table[lo][i]);
    }
    setAnchoredCurve17(curves, gg, damping);
    curves->damping_min = *std::min_element(curves->damping_ratio.begin(), curves->damping_ratio.end());
    curves->name = "Vucetic-Dobry";
}

void makeSeedIdriss(ReferenceCurves *curves, const char *which) {
    static constexpr double upper_gg[] = {1.000, 1.000, 0.990, 0.960, 0.850, 0.655, 0.370, 0.190, 0.085, 0.085, 0.085};
    static constexpr double upper_damping[] = {0.007, 0.012, 0.027, 0.055, 0.099, 0.148, 0.210, 0.255, 0.279, 0.279, 0.279};
    static constexpr double mean_gg[] = {1.000, 0.980, 0.950, 0.890, 0.730, 0.520, 0.290, 0.140, 0.060, 0.060, 0.060};
    static constexpr double mean_damping[] = {0.009, 0.021, 0.028, 0.031, 0.065, 0.107, 0.165, 0.219, 0.257, 0.257, 0.257};
    static constexpr double lower_gg[] = {1.000, 1.000, 0.990, 0.930, 0.770, 0.550, 0.330, 0.180, 0.080, 0.080, 0.080};
    static constexpr double lower_damping[] = {0.0024, 0.0042, 0.0080, 0.0150, 0.0310, 0.0550, 0.0960, 0.1510, 0.2100, 0.2100, 0.2100};
    const bool upper = std::string(which) == "upper";
    const bool mean = std::string(which) == "mean";
    setAnchoredCurve(curves, upper ? upper_gg : (mean ? mean_gg : lower_gg),
                     upper ? upper_damping : (mean ? mean_damping : lower_damping), 11);
    curves->damping_min = *std::min_element(curves->damping_ratio.begin(), curves->damping_ratio.end());
    curves->name = upper ? "Seed-Idriss Sand (Upper)" :
        (mean ? "Seed-Idriss Sand (Mean)" : "Seed-Idriss Sand (Lower)");
}

void makeIdrissRock(ReferenceCurves *curves) {
    static constexpr double gg[] = {1.000, 1.000, 0.980, 0.940, 0.840, 0.710, 0.530, 0.350, 0.180, 0.180, 0.180};
    static constexpr double damping[] = {0.0040, 0.0050, 0.0080, 0.0120, 0.0210, 0.0340, 0.0590, 0.0920, 0.1250, 0.1250, 0.1250};
    setAnchoredCurve(curves, gg, damping, 11);
    curves->damping_min = *std::min_element(curves->damping_ratio.begin(), curves->damping_ratio.end());
    curves->name = "Idriss (1990) Rock";
}

void makeMenq(ReferenceCurves *curves, const CalibrationInput &input) {
    const double cu = positiveOrThrow("cu", input.coef_unif);
    const double stress = positiveOrThrow("sigma_m_eff", input.mean_effective_stress_kpa);
    const double stress_atm = stress / kAtmosphericPressureKpa;
    const double strain_ref = (0.12 * std::pow(cu, -0.6) *
                               std::pow(stress_atm, 0.5 * std::pow(cu, -0.15))) / 100.0;
    const double curvature = 0.86 + 0.1 * std::log10(stress_atm);
    const double scaling = 0.6329 - 0.00566 * std::log(10.0);
    curves->damping_min = clamp(
        (0.55 * std::pow(cu, 0.1) * std::pow(0.3, -0.3) *
         std::pow(stress_atm, -0.08)) / 100.0, 0.0, 0.5);
    double running_damping = 0.0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        curves->modulus_reduction[index] = modifiedHyperbolicReduction(
            curves->strain[index], strain_ref, curvature);
        const double correction = modifiedHyperbolicMasingDamping(
            curves->strain[index], strain_ref, curvature, scaling);
        running_damping = std::max(running_damping, correction);
        curves->damping_ratio[index] = clamp(curves->damping_min + running_damping, 0.0, 0.5);
    }
    curves->name = "Menq";
}

void makeIshibashiZhang(ReferenceCurves *curves, const CalibrationInput &input) {
    const double pi = nonnegativeOrThrow("pi", input.plasticity_index);
    const double stress = positiveOrThrow("sigma_m_eff", input.mean_effective_stress_kpa);
    double n_pi = 0.0;
    if (pi <= 15.0)
        n_pi = pi <= 0.0 ? 0.0 : 3.37e-6 * std::pow(pi, 1.404);
    else if (pi <= 70.0)
        n_pi = 7.0e-7 * std::pow(pi, 1.976);
    else
        n_pi = 2.7e-5 * std::pow(pi, 1.115);
    const double pi_power = std::pow(pi, 1.3);
    const double pi_factor = (1.0 + std::exp(-0.0145 * pi_power)) / 2.0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const double gamma = std::max(std::abs(curves->strain[index]), 1.0e-12);
        const double k_gamma = 0.5 * (1.0 + std::tanh(
            std::log(std::pow((0.000102 + n_pi) / gamma, 0.492))));
        const double m_minus_m0 = 0.272 *
            (1.0 - std::tanh(std::log(std::pow(0.000556 / gamma, 0.4)))) *
            std::exp(-0.0145 * pi_power);
        const double gg = clamp(k_gamma * std::pow(stress, m_minus_m0), 0.0, 1.0);
        curves->modulus_reduction[index] = gg;
        curves->damping_ratio[index] = clamp(
            0.333 * pi_factor * (0.586 * gg * gg - 1.547 * gg + 1.0), 0.0, 0.5);
    }
    curves->damping_min = *std::min_element(curves->damping_ratio.begin(), curves->damping_ratio.end());
    curves->name = "Ishibashi-Zhang";
}

struct FitBounds {
    std::array<double, 5> lower = {-25.0, -25.0, 0.05, 0.05, 0.0};
    std::array<double, 5> upper = {1.0, 25.0, 150.0, 10.0, 4.0};
};

struct FitContext {
    const CalibrationInput &input;
    const ReferenceCurves &target;
    FitBounds bounds;
    double damping_min = 0.0;
    std::string procedure;
    std::string mrdf_mode;
};

gqh::Parameters makeGqhParameters(const std::array<double, 5> &theta,
                                   const std::array<double, 3> &mrdf,
                                   const FitContext &context) {
    gqh::Parameters p;
    p.gmax = context.input.gmax;
    p.tau_max = context.input.tau_max;
    p.theta1 = theta[0];
    p.theta2 = theta[1];
    p.theta3 = theta[2];
    p.theta4 = theta[3];
    p.theta5 = theta[4];
    p.mrdf_p1 = mrdf[0];
    p.mrdf_p2 = mrdf[1];
    p.mrdf_p3 = mrdf[2];
    return p;
}

double dampingFor(const double gamma, const gqh::Parameters &p, const double damping_min,
                  const std::string &mrdf_mode) {
    const double reduction = clamp(gqh::gqh_reduction(std::abs(gamma), p), 0.0, 1.0);
    const double masing = gqh::masing_damping(p, std::abs(gamma), 512);
    double factor = 1.0;
    if (mrdf_mode == "darendeli")
        factor = p.mrdf_p1 * std::pow(std::max(reduction, 0.0), p.mrdf_p2);
    else
        factor = p.mrdf_p1 - p.mrdf_p2 * std::pow(std::max(1.0 - reduction, 0.0), p.mrdf_p3);
    factor = clamp(factor, 0.0, 1.5);
    return clamp(damping_min + masing * factor, 0.0, 0.4);
}

double dampingMaxFit(const std::array<double, kReferencePointCount> &reduction,
                     const ReferenceCurves &target, const double damping_min) {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const double x = clamp(1.0 - reduction[static_cast<std::size_t>(i)], 0.0, 1.0);
        const double y = std::max(target.damping_ratio[static_cast<std::size_t>(i)] - damping_min, 0.0);
        numerator += x * y;
        denominator += x * x;
    }
    if (denominator <= 1.0e-14)
        return clamp(damping_min, 0.0, 0.5);
    return clamp(damping_min + std::max(numerator / denominator, 0.0), damping_min, 0.5);
}

double logResidual(const double predicted, const double target) {
    return std::log(clamp(predicted, kLogFloor, 1.0)) -
           std::log(clamp(target, kLogFloor, 1.0));
}

double thetaCost(const std::array<double, 5> &theta, const FitContext &context,
                 const std::array<double, 3> &mrdf, const bool actual_damping) {
    const gqh::Parameters p = makeGqhParameters(theta, mrdf, context);
    std::array<double, kReferencePointCount> reduction = {};
    for (int i = 0; i < kReferencePointCount; ++i)
        reduction[static_cast<std::size_t>(i)] = clamp(gqh::gqh_reduction(
            context.target.strain[static_cast<std::size_t>(i)], p), 0.0, 1.0);

    double cost = 0.0;
    const double damping_max = dampingMaxFit(reduction, context.target, context.damping_min);
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const double gamma = context.target.strain[index];
        const bool in_modulus_window = gamma >= kFitStrainMin && gamma <= kFitStrainMax;
        const double damping = actual_damping
            ? dampingFor(gamma, p, context.damping_min, context.mrdf_mode)
            : context.damping_min + (1.0 - reduction[index]) * (damping_max - context.damping_min);
        if (in_modulus_window) {
            const double residual = logResidual(reduction[index], context.target.modulus_reduction[index]);
            cost += residual * residual;
        }
        const double damping_residual = logResidual(damping, context.target.damping_ratio[index]);
        const double damping_weight = context.procedure == "MR" ? 0.35 :
            (context.procedure == "DC" ? 1.0 : 1.0);
        cost += (damping_weight * damping_residual) * (damping_weight * damping_residual);
    }
    const double strength = gqh::gqh_backbone(kTargetStrengthStrain, p) /
                            std::max(context.input.tau_max, 1.0e-30);
    const double strength_residual = 6.0 * (strength - kTargetStrengthRatio);
    cost += strength_residual * strength_residual;
    cost += 0.02 * (theta[3] - 1.0) * (theta[3] - 1.0);
    cost += 0.02 * (theta[4] - 1.0) * (theta[4] - 1.0);
    cost += 50.0 * std::max(theta[0] + theta[1] - 1.0, 0.0) *
            std::max(theta[0] + theta[1] - 1.0, 0.0);
    cost += 0.01 * std::max(theta[2] - 80.0, 0.0) * std::max(theta[2] - 80.0, 0.0);
    return std::isfinite(cost) ? cost : std::numeric_limits<double>::max() / 4.0;
}

double mrdfCost(const std::array<double, 3> &mrdf, const std::array<double, 5> &theta,
                const FitContext &context) {
    const gqh::Parameters p = makeGqhParameters(theta, mrdf, context);
    double cost = 0.0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const double reduction = clamp(gqh::gqh_reduction(context.target.strain[index], p), 0.0, 1.0);
        double factor = 1.0;
        if (context.mrdf_mode == "darendeli")
            factor = mrdf[0] * std::pow(std::max(reduction, 0.0), mrdf[1]);
        else
            factor = mrdf[0] - mrdf[1] * std::pow(std::max(1.0 - reduction, 0.0), mrdf[2]);
        factor = clamp(factor, 0.02, 1.0);
        const double predicted = clamp(context.damping_min +
            gqh::masing_damping(p, context.target.strain[index], 512) * factor, 0.0, 0.5);
        const double residual = logResidual(predicted, context.target.damping_ratio[index]);
        cost += residual * residual;
    }
    if (context.mrdf_mode == "darendeli") {
        cost += 0.02 * (mrdf[0] - 1.0) * (mrdf[0] - 1.0);
        cost += 0.02 * mrdf[1] * mrdf[1];
    } else {
        cost += 0.01 * (mrdf[0] - 1.0) * (mrdf[0] - 1.0);
        cost += 0.01 * mrdf[2] * mrdf[2] / 100.0;
    }
    return std::isfinite(cost) ? cost : std::numeric_limits<double>::max() / 4.0;
}

template <std::size_t N>
struct SimplexResult {
    std::array<double, N> point = {};
    double cost = std::numeric_limits<double>::max();
    int iterations = 0;
    bool converged = false;
};

template <std::size_t N, typename CostFunction>
SimplexResult<N> nelderMead(const std::array<double, N> &initial,
                            const std::array<double, N> &lower,
                            const std::array<double, N> &upper,
                            CostFunction cost_function) {
    std::array<std::array<double, N>, N + 1> simplex = {};
    std::array<double, N + 1> costs = {};
    simplex[0] = initial;
    for (std::size_t i = 0; i < N; ++i) {
        simplex[i + 1] = initial;
        const double scale = std::max(std::abs(initial[i]), upper[i] - lower[i]);
        const double step = std::max(1.0e-4 * (upper[i] - lower[i]), 0.05 * scale);
        simplex[i + 1][i] = clamp(initial[i] + step, lower[i], upper[i]);
        if (simplex[i + 1][i] == initial[i])
            simplex[i + 1][i] = clamp(initial[i] - step, lower[i], upper[i]);
    }
    for (std::size_t i = 0; i < N + 1; ++i)
        costs[i] = cost_function(simplex[i]);

    SimplexResult<N> result;
    for (int iteration = 0; iteration < kMaxIterations; ++iteration) {
        std::array<std::size_t, N + 1> order = {};
        for (std::size_t i = 0; i < N + 1; ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(), [&](const std::size_t left, const std::size_t right) {
            return costs[left] < costs[right];
        });
        std::array<std::array<double, N>, N + 1> ordered_simplex = {};
        std::array<double, N + 1> ordered_costs = {};
        for (std::size_t i = 0; i < N + 1; ++i) {
            ordered_simplex[i] = simplex[order[i]];
            ordered_costs[i] = costs[order[i]];
        }
        simplex = ordered_simplex;
        costs = ordered_costs;
        double diameter = 0.0;
        for (std::size_t i = 1; i < N + 1; ++i)
            for (std::size_t j = 0; j < N; ++j)
                diameter = std::max(diameter, std::abs(simplex[i][j] - simplex[0][j]) /
                    std::max(1.0, upper[j] - lower[j]));
        if (diameter < 1.0e-6 || std::abs(costs[N] - costs[0]) < 1.0e-10) {
            result.point = simplex[0];
            result.cost = costs[0];
            result.iterations = iteration + 1;
            result.converged = true;
            return result;
        }

        std::array<double, N> centroid = {};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                centroid[j] += simplex[i][j] / static_cast<double>(N);
        std::array<double, N> reflected = {};
        for (std::size_t j = 0; j < N; ++j)
            reflected[j] = clamp(centroid[j] + (centroid[j] - simplex[N][j]), lower[j], upper[j]);
        const double reflected_cost = cost_function(reflected);
        if (reflected_cost < costs[0]) {
            std::array<double, N> expanded = {};
            for (std::size_t j = 0; j < N; ++j)
                expanded[j] = clamp(centroid[j] + 2.0 * (reflected[j] - centroid[j]), lower[j], upper[j]);
            const double expanded_cost = cost_function(expanded);
            simplex[N] = expanded_cost < reflected_cost ? expanded : reflected;
            costs[N] = std::min(expanded_cost, reflected_cost);
        } else if (reflected_cost < costs[N - 1]) {
            simplex[N] = reflected;
            costs[N] = reflected_cost;
        } else {
            std::array<double, N> contracted = {};
            for (std::size_t j = 0; j < N; ++j)
                contracted[j] = clamp(centroid[j] + 0.5 * (simplex[N][j] - centroid[j]), lower[j], upper[j]);
            const double contracted_cost = cost_function(contracted);
            if (contracted_cost < costs[N]) {
                simplex[N] = contracted;
                costs[N] = contracted_cost;
            } else {
                for (std::size_t i = 1; i < N + 1; ++i) {
                    for (std::size_t j = 0; j < N; ++j)
                        simplex[i][j] = clamp(simplex[0][j] + 0.5 *
                            (simplex[i][j] - simplex[0][j]), lower[j], upper[j]);
                    costs[i] = cost_function(simplex[i]);
                }
            }
        }
        result.iterations = iteration + 1;
    }
    std::size_t best = 0;
    for (std::size_t i = 1; i < N + 1; ++i)
        if (costs[i] < costs[best])
            best = i;
    result.point = simplex[best];
    result.cost = costs[best];
    result.converged = false;
    return result;
}

std::array<double, 5> clampTheta(std::array<double, 5> theta, const FitContext &context) {
    for (int i = 0; i < 5; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        if (context.input.fixed_theta[index])
            theta[index] = context.input.theta[index];
        else
            theta[index] = clamp(theta[index], context.bounds.lower[index], context.bounds.upper[index]);
    }
    return theta;
}

std::array<double, 5> fitTheta(const FitContext &context, const std::array<double, 3> &mrdf,
                                const bool actual_damping, int *iterations, bool *converged) {
    static constexpr std::array<std::array<double, 5>, 4> seeds = {{
        {{-2.5, 3.0, 5.0, 1.0, 1.0}},
        {{-5.0, 5.0, 15.0, 1.0, 1.0}},
        {{-1.0, 1.0, 1.0, 1.0, 1.0}},
        {{-8.0, 2.0, 20.0, 1.0, 0.99}},
    }};
    std::array<double, 5> best = clampTheta(seeds[0], context);
    double best_cost = thetaCost(best, context, mrdf, actual_damping);
    bool all_converged = true;
    int total_iterations = 0;
    for (const auto &seed : seeds) {
        std::array<int, 5> map = {};
        std::array<double, 5> initial_free = {};
        std::array<double, 5> lower_free = {};
        std::array<double, 5> upper_free = {};
        std::size_t free_count = 0;
        for (int i = 0; i < 5; ++i) {
            if (!context.input.fixed_theta[static_cast<std::size_t>(i)]) {
                map[free_count] = i;
                initial_free[free_count] = clamp(seed[static_cast<std::size_t>(i)],
                    context.bounds.lower[static_cast<std::size_t>(i)], context.bounds.upper[static_cast<std::size_t>(i)]);
                lower_free[free_count] = context.bounds.lower[static_cast<std::size_t>(i)];
                upper_free[free_count] = context.bounds.upper[static_cast<std::size_t>(i)];
                ++free_count;
            }
        }
        if (free_count == 0)
            break;
        std::array<double, 5> full = clampTheta(seed, context);
        auto makeFull = [&](const std::array<double, 5> &free) {
            std::array<double, 5> value = full;
            for (std::size_t j = 0; j < free_count; ++j)
                value[static_cast<std::size_t>(map[j])] = free[j];
            return clampTheta(value, context);
        };
        if (free_count == 1) {
            const auto nm = nelderMead<1>(
                std::array<double, 1>{initial_free[0]}, std::array<double, 1>{lower_free[0]},
                std::array<double, 1>{upper_free[0]}, [&](const std::array<double, 1> &free) {
                    return thetaCost(makeFull(std::array<double, 5>{free[0], 0, 0, 0, 0}),
                        context, mrdf, actual_damping);
                });
            const auto candidate = makeFull(std::array<double, 5>{nm.point[0], 0, 0, 0, 0});
            if (nm.cost < best_cost) { best = candidate; best_cost = nm.cost; }
            all_converged = all_converged && nm.converged;
            total_iterations = std::max(total_iterations, nm.iterations);
        } else if (free_count == 2) {
            const auto nm = nelderMead<2>(
                std::array<double, 2>{initial_free[0], initial_free[1]},
                std::array<double, 2>{lower_free[0], lower_free[1]},
                std::array<double, 2>{upper_free[0], upper_free[1]}, [&](const std::array<double, 2> &free) {
                    return thetaCost(makeFull(std::array<double, 5>{free[0], free[1], 0, 0, 0}),
                        context, mrdf, actual_damping);
                });
            const auto candidate = makeFull(std::array<double, 5>{nm.point[0], nm.point[1], 0, 0, 0});
            if (nm.cost < best_cost) { best = candidate; best_cost = nm.cost; }
            all_converged = all_converged && nm.converged;
            total_iterations = std::max(total_iterations, nm.iterations);
        } else if (free_count == 3) {
            const auto nm = nelderMead<3>(
                std::array<double, 3>{initial_free[0], initial_free[1], initial_free[2]},
                std::array<double, 3>{lower_free[0], lower_free[1], lower_free[2]},
                std::array<double, 3>{upper_free[0], upper_free[1], upper_free[2]}, [&](const std::array<double, 3> &free) {
                    return thetaCost(makeFull(std::array<double, 5>{free[0], free[1], free[2], 0, 0}),
                        context, mrdf, actual_damping);
                });
            const auto candidate = makeFull(std::array<double, 5>{nm.point[0], nm.point[1], nm.point[2], 0, 0});
            if (nm.cost < best_cost) { best = candidate; best_cost = nm.cost; }
            all_converged = all_converged && nm.converged;
            total_iterations = std::max(total_iterations, nm.iterations);
        } else if (free_count == 4) {
            const auto nm = nelderMead<4>(
                std::array<double, 4>{initial_free[0], initial_free[1], initial_free[2], initial_free[3]},
                std::array<double, 4>{lower_free[0], lower_free[1], lower_free[2], lower_free[3]},
                std::array<double, 4>{upper_free[0], upper_free[1], upper_free[2], upper_free[3]}, [&](const std::array<double, 4> &free) {
                    return thetaCost(makeFull(std::array<double, 5>{free[0], free[1], free[2], free[3], 0}),
                        context, mrdf, actual_damping);
                });
            const auto candidate = makeFull(std::array<double, 5>{nm.point[0], nm.point[1], nm.point[2], nm.point[3], 0});
            if (nm.cost < best_cost) { best = candidate; best_cost = nm.cost; }
            all_converged = all_converged && nm.converged;
            total_iterations = std::max(total_iterations, nm.iterations);
        } else {
            const auto nm = nelderMead<5>(
                clampTheta(seed, context), context.bounds.lower, context.bounds.upper,
                [&](const std::array<double, 5> &value) {
                    return thetaCost(clampTheta(value, context), context, mrdf, actual_damping);
                });
            const auto candidate = clampTheta(nm.point, context);
            if (nm.cost < best_cost) { best = candidate; best_cost = nm.cost; }
            all_converged = all_converged && nm.converged;
            total_iterations = std::max(total_iterations, nm.iterations);
        }
    }
    if (iterations) *iterations = total_iterations;
    if (converged) *converged = all_converged;
    return clampTheta(best, context);
}

std::array<double, 3> fitMrdf(const FitContext &context, const std::array<double, 5> &theta,
                              int *iterations, bool *converged) {
    const std::array<double, 3> lower = context.mrdf_mode == "darendeli"
        ? std::array<double, 3>{0.05, 0.0, 0.0}
        : std::array<double, 3>{0.05, 0.0, 0.5};
    const std::array<double, 3> upper = context.mrdf_mode == "darendeli"
        ? std::array<double, 3>{1.0, 4.0, 0.0}
        : std::array<double, 3>{1.0, 1.0, 40.0};
    const std::array<std::array<double, 3>, 4> seeds = context.mrdf_mode == "darendeli"
        ? std::array<std::array<double, 3>, 4>{{{{1.0, 0.0, 0.0}}, {{0.9, 0.1, 0.0}}, {{0.8, 0.3, 0.0}}, {{0.6, 1.0, 0.0}}}}
        : std::array<std::array<double, 3>, 4>{{{{1.0, 0.0, 20.0}}, {{0.9, 0.2, 15.0}}, {{0.87, 0.65, 8.0}}, {{0.8, 0.5, 25.0}}}};
    std::array<double, 3> best = seeds[0];
    double best_cost = mrdfCost(best, theta, context);
    bool all_converged = true;
    int total_iterations = 0;
    for (const auto &seed : seeds) {
        const auto nm = nelderMead<3>(seed, lower, upper, [&](const std::array<double, 3> &value) {
            return mrdfCost(value, theta, context);
        });
        if (nm.cost < best_cost) {
            best = nm.point;
            best_cost = nm.cost;
        }
        all_converged = all_converged && nm.converged;
        total_iterations = std::max(total_iterations, nm.iterations);
    }
    if (iterations) *iterations = total_iterations;
    if (converged) *converged = all_converged;
    return best;
}

double modulusRmse(const std::array<double, 5> &theta, const FitContext &context) {
    const gqh::Parameters p = makeGqhParameters(theta, {1.0, 0.0, 0.0}, context);
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        if (context.target.strain[index] < kFitStrainMin || context.target.strain[index] > kFitStrainMax)
            continue;
        const double residual = logResidual(gqh::gqh_reduction(context.target.strain[index], p),
                                            context.target.modulus_reduction[index]);
        sum += residual * residual;
        ++count;
    }
    return count > 0 ? std::sqrt(sum / static_cast<double>(count)) : 0.0;
}

double dampingRmse(const std::array<double, 5> &theta, const std::array<double, 3> &mrdf,
                   const FitContext &context) {
    const gqh::Parameters p = makeGqhParameters(theta, mrdf, context);
    double sum = 0.0;
    for (int i = 0; i < kReferencePointCount; ++i) {
        const std::size_t index = static_cast<std::size_t>(i);
        const double residual = logResidual(dampingFor(context.target.strain[index], p,
                                                        context.damping_min, context.mrdf_mode),
                                            context.target.damping_ratio[index]);
        sum += residual * residual;
    }
    return std::sqrt(sum / static_cast<double>(kReferencePointCount));
}

} // namespace

std::string normalizeCurveType(const std::string &value) {
    std::string result;
    result.reserve(value.size());
    for (const char character : value)
        result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(character))));
    if (result == "vucetic-dobry") return "vucetic_dobry";
    if (result == "ishibashi-zhang") return "ishibashi_zhang";
    if (result == "rock") return "idriss_rock";
    return result;
}

std::string normalizeMrdfMode(const std::string &value) {
    std::string result;
    result.reserve(value.size());
    for (const char character : value)
        result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(character))));
    return result;
}

std::string normalizeFitProcedure(const std::string &value) {
    std::string result;
    result.reserve(value.size());
    for (const char character : value)
        result.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(character))));
    return result;
}

bool isSupportedCurveType(const std::string &value) {
    const std::string curve = normalizeCurveType(value);
    return curve == "darendeli" || curve == "vucetic_dobry" ||
           curve == "seed_idriss_upper" || curve == "seed_idriss_mean" ||
           curve == "seed_idriss_lower" || curve == "idriss_rock" ||
           curve == "ishibashi_zhang" || curve == "menq" || curve == "none";
}

ReferenceCurves makeReferenceCurves(const CalibrationInput &input) {
    const std::string curve = normalizeCurveType(input.curve_type);
    if (!isSupportedCurveType(curve) || curve == "none")
        throw std::runtime_error("GQ/H refcurve: unknown curve_type '" + input.curve_type + "'.");
    ReferenceCurves curves;
    fillLogGrid(&curves);
    if (curve == "darendeli")
        makeDarendeli(&curves, input);
    else if (curve == "vucetic_dobry")
        makeVuceticDobry(&curves, input.plasticity_index);
    else if (curve == "seed_idriss_upper")
        makeSeedIdriss(&curves, "upper");
    else if (curve == "seed_idriss_mean")
        makeSeedIdriss(&curves, "mean");
    else if (curve == "seed_idriss_lower")
        makeSeedIdriss(&curves, "lower");
    else if (curve == "idriss_rock")
        makeIdrissRock(&curves);
    else if (curve == "ishibashi_zhang")
        makeIshibashiZhang(&curves, input);
    else if (curve == "menq")
        makeMenq(&curves, input);
    curves.damping_min = clamp(*std::min_element(curves.damping_ratio.begin(), curves.damping_ratio.end()),
                               0.0, 0.5);
    return curves;
}

CalibrationResult calibrate(const CalibrationInput &input) {
    positiveOrThrow("gmax", input.gmax);
    positiveOrThrow("tau_max", input.tau_max);
    nonnegativeOrThrow("pi", input.plasticity_index);
    positiveOrThrow("ocr", input.ocr);
    positiveOrThrow("sigma_m_eff", input.mean_effective_stress_kpa);
    positiveOrThrow("cu", input.coef_unif);
    const std::string mrdf_mode = normalizeMrdfMode(input.mrdf_mode);
    if (mrdf_mode != "uiuc" && mrdf_mode != "darendeli")
        throw std::runtime_error("GQ/H refcurve: mrdf_mode must be uiuc or darendeli.");
    const std::string procedure = normalizeFitProcedure(input.fit_procedure);
    if (procedure != "MR" && procedure != "MRD" && procedure != "DC")
        throw std::runtime_error("GQ/H refcurve: fit_procedure must be MR, MRD, or DC.");

    const ReferenceCurves target = makeReferenceCurves(input);
    FitContext context{input, target, FitBounds{}, target.damping_min, procedure, mrdf_mode};
    std::array<double, 3> mrdf = mrdf_mode == "darendeli"
        ? std::array<double, 3>{1.0, 0.0, 0.0}
        : std::array<double, 3>{1.0, 0.0, 20.0};
    int theta_iterations = 0;
    int mrdf_iterations = 0;
    bool theta_converged = false;
    bool mrdf_converged = false;
    std::array<double, 5> theta = fitTheta(context, mrdf, false, &theta_iterations, &theta_converged);
    mrdf = fitMrdf(context, theta, &mrdf_iterations, &mrdf_converged);
    if (procedure == "MRD" || procedure == "DC") {
        theta = fitTheta(context, mrdf, true, &theta_iterations, &theta_converged);
        mrdf = fitMrdf(context, theta, &mrdf_iterations, &mrdf_converged);
    }

    CalibrationResult result;
    result.theta = theta;
    result.mrdf = mrdf;
    result.damping_min = target.damping_min;
    result.damping_max = *std::max_element(target.damping_ratio.begin(), target.damping_ratio.end());
    result.modulus_rmse = modulusRmse(theta, context);
    result.damping_rmse = dampingRmse(theta, mrdf, context);
    const gqh::Parameters p = makeGqhParameters(theta, mrdf, context);
    result.strength_ratio = gqh::gqh_backbone(kTargetStrengthStrain, p) /
                            std::max(input.tau_max, 1.0e-30);
    result.iterations = std::max(theta_iterations, mrdf_iterations);
    result.converged = theta_converged && mrdf_converged;
    result.valid = true;
    result.curve_name = target.name;
    return result;
}

} // namespace gqh::refcurve
