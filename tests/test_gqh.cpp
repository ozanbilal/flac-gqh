#include "gqh_engine.h"
#include "gqh_math.h"
#include "refcurve.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double got, double expected, double tolerance, const std::string &message) {
    if (!std::isfinite(got) || std::abs(got - expected) > tolerance) {
        std::ostringstream out;
        out.precision(17);
        out << message << ": got " << got << ", expected " << expected;
        throw std::runtime_error(out.str());
    }
}

gqh::Parameters hyperbola() {
    gqh::Parameters p;
    p.gmax = 100.0;
    p.tau_max = 10.0;
    p.theta1 = 0.0;
    p.theta2 = 0.0;
    p.mrdf_p1 = 1.0;
    p.mrdf_p2 = 0.0;
    return p;
}

void test_backbone() {
    // With theta=0, tau = Gmax*gamma/(1+|gamma|/gamma_ref).
    // These numerical values are closed-form examples, not fitted outputs.
    const auto p = hyperbola();
    near(gqh::gqh_backbone(0.0, p), 0.0, 1e-14, "zero strain");
    near(gqh::gqh_backbone(0.05, p), 10.0 / 3.0, 1e-12, "half reference strain");
    near(gqh::gqh_backbone(0.1, p), 5.0, 1e-12, "reference strain");
    near(gqh::gqh_backbone(-0.2, p), -20.0 / 3.0, 1e-12, "negative strain");
    near(gqh::gqh_reduction(0.1, p), 0.5, 1e-12, "modulus reduction");
    require(gqh::gqh_backbone(1000.0, p) < p.tau_max, "strength asymptote");
    auto q = p;
    q.theta1 = 0.6;
    const double r = 1.0;
    const double x = gqh::gqh_backbone(0.1, q) / q.tau_max;
    // Independent implicit quadratic: theta*x^2 - (1+r)*x + r = 0.
    near(q.theta1 * x * x - (1.0 + r) * x + r, 0.0, 1e-12, "quadratic residual");
    require(x > 0.5 && x < 1.0, "quadratic physical root");
}

void test_mrdf_and_damping() {
    auto p = hyperbola();
    p.mrdf_p1 = 0.9;
    p.mrdf_p2 = 0.2;
    p.mrdf_p3 = 2.0;
    near(gqh::mrdf_factor(0.1, p), 0.85, 1e-12, "MRDF factor");
    // Analytic integral of the hyperbolic backbone: tau_max*(a-rf*ln(1+a/rf)).
    p = hyperbola();
    const double a = 0.1;
    const double area = p.tau_max * (a - a * std::log(2.0));
    const double expected = (2.0 / std::acos(-1.0)) * (2.0 * area / (5.0 * a) - 1.0);
    near(gqh::masing_damping(p, a, 8192), expected, 2e-5, "Masing damping integral");
}

void test_scalar_reversal_and_restore() {
    const auto p = hyperbola();
    gqh::MasingState state(p);
    near(state.update(0.1), 5.0, 1e-12, "virgin loading");
    near(state.update(0.0), -5.0 / 3.0, 1e-12, "first unloading");
    require(state.branchDepth() == 1, "record first reversal");
    std::stringstream saved(std::ios::in | std::ios::out | std::ios::binary);
    state.save(saved);
    saved.seekg(0);
    gqh::MasingState restored(p);
    restored.restore(saved);
    require(!saved.fail(), "successful restore");
    for (const double gamma : {-0.04, 0.01, -0.02, -0.1, -0.11}) {
        near(restored.update(gamma), state.update(gamma), 1e-12, "continuation after restore");
    }
    near(state.currentTau(), gqh::gqh_backbone(-0.11, p), 1e-12, "return to backbone");
}

void test_tensor_uniaxial_mapping() {
    const auto p = hyperbola();
    gqh::engine::Subzone subzone;
    const double c = 2.0 * std::sqrt(2.0 / 3.0);
    const double zero[3] = {0.0, 0.0, 0.0};
    const double peak[3] = {0.1 / c, 0.0, 0.0};
    double stress[3] = {};
    double tangent = 0.0;
    near(gqh::engine::gammaOctBetween(peak, zero), 0.1, 1e-14, "octahedral strain");
    gqh::engine::advance(subzone, zero, peak, p, p.gmax, stress, &tangent);
    near(stress[0], std::sqrt(1.5) * 5.0, 1e-11, "tensor virgin loading");
    near(stress[1], 0.0, 1e-14, "unused shear component");
    require(std::isfinite(tangent) && tangent > 0.0 && tangent <= p.gmax, "bounded tangent");
    gqh::engine::advance(subzone, peak, zero, p, p.gmax, stress, &tangent);
    near(stress[0], -std::sqrt(1.5) * 5.0 / 3.0, 1e-11, "tensor unloading");
    require(subzone.stack.size() == 1, "tensor reversal recorded");
}

void test_reference_curve_inputs() {
    using namespace gqh::refcurve;
    require(normalizeCurveType("VUCETIC-DOBRY") == "vucetic_dobry", "curve alias");
    require(isSupportedCurveType("DARENDELI"), "supported curve");
    CalibrationInput input;
    input.gmax = 100.0;
    input.tau_max = 10.0;
    input.curve_type = "darendeli";
    const auto curve = makeReferenceCurves(input);
    for (std::size_t i = 0; i < curve.strain.size(); ++i) {
        require(std::isfinite(curve.modulus_reduction[i]) && curve.modulus_reduction[i] >= 0.0 &&
                curve.modulus_reduction[i] <= 1.0, "finite reference reduction");
        require(std::isfinite(curve.damping_ratio[i]) && curve.damping_ratio[i] >= 0.0 &&
                curve.damping_ratio[i] <= 0.5, "finite reference damping");
        if (i) require(curve.strain[i] > curve.strain[i - 1], "increasing strain grid");
    }
    input.mean_effective_stress_kpa = -1.0;
    bool threw = false;
    try { (void)calibrate(input); } catch (const std::runtime_error &) { threw = true; }
    require(threw, "reject invalid effective stress");
}
} // namespace

int main() {
    try {
        test_backbone();
        test_mrdf_and_damping();
        test_scalar_reversal_and_restore();
        test_tensor_uniaxial_mapping();
        test_reference_curve_inputs();
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
    std::cout << "All constitutive core checks passed.\n";
    return 0;
}
