#pragma once

#include <array>
#include <string>

namespace gqh::refcurve {

constexpr int kReferencePointCount = 60;

struct CalibrationInput {
    double gmax = 0.0;
    double tau_max = 0.0;
    double plasticity_index = 0.0;
    double ocr = 1.0;
    double mean_effective_stress_kpa = 10.0;
    double coef_unif = 2.0;
    std::string curve_type = "none";
    std::string mrdf_mode = "uiuc";
    std::string fit_procedure = "MR";
    std::array<bool, 5> fixed_theta = {};
    std::array<double, 5> theta = {};
};

struct ReferenceCurves {
    std::array<double, kReferencePointCount> strain = {};
    std::array<double, kReferencePointCount> modulus_reduction = {};
    std::array<double, kReferencePointCount> damping_ratio = {};
    double damping_min = 0.0;
    std::string name;
};

struct CalibrationResult {
    std::array<double, 5> theta = {};
    std::array<double, 3> mrdf = {};
    double damping_min = 0.0;
    double damping_max = 0.0;
    double modulus_rmse = 0.0;
    double damping_rmse = 0.0;
    double strength_ratio = 0.0;
    int iterations = 0;
    bool converged = false;
    bool valid = false;
    std::string curve_name;
};

std::string normalizeCurveType(const std::string &value);
std::string normalizeMrdfMode(const std::string &value);
std::string normalizeFitProcedure(const std::string &value);
bool isSupportedCurveType(const std::string &value);

ReferenceCurves makeReferenceCurves(const CalibrationInput &input);
CalibrationResult calibrate(const CalibrationInput &input);

} // namespace gqh::refcurve
