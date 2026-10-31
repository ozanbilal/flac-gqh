#pragma once

// GQ/H + MRDF-UIUC hysteresis engine, free of any FLAC3D dependency.
//
// The FLAC3D v9 and v6 plug-ins are thin adapters over this header: they
// translate the host's State/property API and call into the functions below.
// Keeping the physics here reduces the risk of the two builds drifting apart;
// changes still require validation for each host version.
//
// Model summary:
//  * State lives PER TETRAHEDRON, not per zone.  FLAC3D splits a brick into
//    two overlays of five tets that can have different shear strain histories;
//    a single zone-level history cannot represent these differences.
//  * The hysteresis is a TOTAL map in the shear strain sub-space (e12,e13,e23):
//    stress is evaluated from the octahedral distance between the current shear
//    strain and the active branch anchor.  No per-increment scalar direction,
//    so there is no incremental drift.
//  * Normal deviatoric components stay elastic at gmax; the modulus reduction
//    is a shear phenomenon.

#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>

#include "gqh_math.h"

namespace gqh {
namespace engine {

constexpr double kSqrt3 = 1.7320508075688772935274463415059;
// Matches the 1e-12 gamma tolerance of the scalar MasingState so a purely
// uniaxial strain path reproduces the 1-D rules bit for bit.
constexpr double kGammaTolerance = 1.0e-12;
// A retreat smaller than this fraction of the active branch amplitude would be
// treated as numerical wiggle rather than unloading. Zero disables this hook.
constexpr double kReversalFraction = 0.0;

// One reversal point of the tensorial Masing history, in the shear sub-space.
struct Reversal {
    double strain[3] = {};
    double stress[3] = {};
    double gamma_m = 1.0e-12;
};

// Per-tetrahedron state.
struct Subzone {
    double strain[3] = {};         // accumulated shear strain of this tet
    double origin_strain[3] = {};  // shear strain at the hysteresis origin
    double origin_stress[3] = {};  // mapped shear stress at the origin
    double map_stress[3] = {};     // mapping value at strain
    double gamma_max = 0.0;        // max octahedral distance from origin
    // Farthest point on the ACTIVE branch and its mapped stress.  A reversal is
    // committed here, not at wherever the path was one step ago; storing the
    // previous step's stress instead puts the anchor one increment behind the
    // turning point and changes loop area.
    double peak_strain[3] = {};
    double peak_stress[3] = {};
    double peak_distance = 0.0;
    bool initialized = false;
    std::vector<Reversal> stack;
};

// gamma_oct = (2/sqrt3)*||e_s||, with the tensor-norm weight of two on the
// off-diagonal terms.  For a pure gamma_xz path this equals the scalar
// coordinate the 1-D MasingState is driven with, so the backbone parameters
// keep their octahedral meaning.
inline double gammaOctBetween(const double *a, const double *b) {
    const double d0 = a[0] - b[0];
    const double d1 = a[1] - b[1];
    const double d2 = a[2] - b[2];
    return (2.0 / kSqrt3) * std::sqrt(2.0 * (d0 * d0 + d1 * d1 + d2 * d2));
}

inline double branchTauOct(const Reversal &tip, const double gamma_r,
                           const Parameters &p) {
    const double half = gamma_r / 2.0;
    const double factor = mrdf_factor(tip.gamma_m, p);
    if (factor == 1.0)
        return 2.0 * gqh_backbone(half, p);
    const double tau_b_m = gqh_backbone(tip.gamma_m, p);
    const double secant = tau_b_m * gamma_r / std::max(tip.gamma_m, 1.0e-12);
    return factor * (2.0 * gqh_backbone(half, p) - secant) + secant;
}

inline void evaluateMapStress(const Subzone &sz, const double *shear,
                              const Parameters &p, double *stress_out) {
    const double *anchor_strain = sz.stack.empty() ? sz.origin_strain
                                                   : sz.stack.back().strain;
    const double *anchor_stress = sz.stack.empty() ? sz.origin_stress
                                                   : sz.stack.back().stress;
    double delta[3];
    for (int i = 0; i < 3; ++i)
        delta[i] = shear[i] - anchor_strain[i];
    const double n = std::sqrt(2.0 * (delta[0] * delta[0] +
                                      delta[1] * delta[1] +
                                      delta[2] * delta[2]));
    if (n <= 1.0e-30) {
        std::copy(anchor_stress, anchor_stress + 3, stress_out);
        return;
    }
    const double gamma = (2.0 / kSqrt3) * n;
    const double tau_oct = sz.stack.empty()
        ? gqh_backbone(gamma, p)
        : branchTauOct(sz.stack.back(), gamma, p);
    // ||s - s_anchor|| = sqrt3 * tau_oct reproduces s_xz = sqrt(3/2)*tau_oct
    // for a pure shear path -- exactly the stress the 1-D model produced.
    const double scale = kSqrt3 * tau_oct / n;
    for (int i = 0; i < 3; ++i)
        stress_out[i] = anchor_stress[i] + scale * delta[i];
}

// Advance one tetrahedron's hysteresis from shear_old to shear_new.
// Writes the mapped shear stress and the branch tangent modulus.
inline void advance(Subzone &sz, const double *shear_old,
                    const double *shear_new, const Parameters &p,
                    const double gmax, double *map_stress_out,
                    double *tangent_out) {
    const double step = gammaOctBetween(shear_new, shear_old);
    if (step > kGammaTolerance) {
        // 1. Reversal, committed at the branch's farthest point.
        const double *anchor = sz.stack.empty() ? sz.origin_strain
                                                : sz.stack.back().strain;
        const double d_new = gammaOctBetween(shear_new, anchor);
        if (d_new < sz.peak_distance * (1.0 - kReversalFraction) &&
            sz.peak_distance > kGammaTolerance) {
            Reversal turn;
            std::copy(std::begin(sz.peak_strain), std::end(sz.peak_strain),
                      std::begin(turn.strain));
            std::copy(std::begin(sz.peak_stress), std::end(sz.peak_stress),
                      std::begin(turn.stress));
            turn.gamma_m = std::max(
                {sz.gamma_max,
                 gammaOctBetween(sz.peak_strain, sz.origin_strain), 1.0e-12});
            sz.stack.push_back(turn);
            sz.peak_distance = 0.0;
        }

        // 2. Loop closure: leaving a nested loop reconnects to the outer branch
        //    (pop in pairs); passing the mirror point of the first reversal
        //    rejoins the virgin backbone.  Distance form of the 1-D pops.
        bool popped = false;
        while (sz.stack.size() >= 2) {
            const Reversal &tip = sz.stack[sz.stack.size() - 1];
            const Reversal &outer = sz.stack[sz.stack.size() - 2];
            if (gammaOctBetween(shear_new, tip.strain) >
                gammaOctBetween(tip.strain, outer.strain) + kGammaTolerance) {
                sz.stack.pop_back();
                sz.stack.pop_back();
                popped = true;
            } else {
                break;
            }
        }
        if (sz.stack.size() == 1) {
            const Reversal &tip = sz.stack.back();
            if (gammaOctBetween(shear_new, tip.strain) >
                2.0 * gammaOctBetween(tip.strain, sz.origin_strain) +
                    kGammaTolerance) {
                sz.stack.pop_back();
                popped = true;
            }
        }
        if (popped)
            sz.peak_distance = 0.0;

        // 3. MRDF amplitude bookkeeping.
        const double gamma_abs = gammaOctBetween(shear_new, sz.origin_strain);
        sz.gamma_max = std::max(sz.gamma_max, gamma_abs);
        if (!sz.stack.empty()) {
            Reversal &tip = sz.stack.back();
            tip.gamma_m = std::max({tip.gamma_m, sz.gamma_max, 1.0e-12});
        }
    }

    // 4. Total map.
    evaluateMapStress(sz, shear_new, p, map_stress_out);

    // 4b. Track the branch peak WITH its own mapped stress (after the map is
    //     evaluated -- see Subzone).
    {
        const double *anch = sz.stack.empty() ? sz.origin_strain
                                              : sz.stack.back().strain;
        const double d_now = gammaOctBetween(shear_new, anch);
        if (d_now >= sz.peak_distance) {
            sz.peak_distance = d_now;
            std::copy(shear_new, shear_new + 3, std::begin(sz.peak_strain));
            std::copy(map_stress_out, map_stress_out + 3,
                      std::begin(sz.peak_stress));
        }
    }

    // 5. Branch tangent at the current distance (Mohr-Coulomb return, stable
    //    timestep, stiffness damping).
    const double *anchor = sz.stack.empty() ? sz.origin_strain
                                            : sz.stack.back().strain;
    const double gamma = gammaOctBetween(shear_new, anchor);
    const double h = std::max(1.0e-9, 1.0e-3 * gamma);
    double tau0, tau1;
    if (sz.stack.empty()) {
        tau0 = gqh_backbone(gamma, p);
        tau1 = gqh_backbone(gamma + h, p);
    } else {
        tau0 = branchTauOct(sz.stack.back(), gamma, p);
        tau1 = branchTauOct(sz.stack.back(), gamma + h, p);
    }
    *tangent_out = std::min(gmax, std::max(1.0e-12, (tau1 - tau0) / h));
}

// Move the hysteresis origin to the current point, keeping the strain/stress
// baseline.  Used when the dynamic solve starts so that static settling is not
// replayed as earthquake history.
inline void restartAtCurrentPoint(Subzone &sz) {
    std::copy(std::begin(sz.strain), std::end(sz.strain),
              std::begin(sz.origin_strain));
    std::copy(std::begin(sz.map_stress), std::end(sz.map_stress),
              std::begin(sz.origin_stress));
    sz.stack.clear();
    sz.gamma_max = 0.0;
    std::copy(std::begin(sz.strain), std::end(sz.strain),
              std::begin(sz.peak_strain));
    std::copy(std::begin(sz.map_stress), std::end(sz.map_stress),
              std::begin(sz.peak_stress));
    sz.peak_distance = 0.0;
}

} // namespace engine
} // namespace gqh
