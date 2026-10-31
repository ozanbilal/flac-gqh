#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "models/src/conmodel.h"
#include "gqh_math.h"
#include "gqh_engine.h"
#include "refcurve.h"

namespace models {

class ModelGqhMrdf : public ConstitutiveModel {
public:
    ModelGqhMrdf(unsigned short option = 0);

    string getName() const override { return "gqh-mrdf"; }
    string getFullName() const override { return "GQ/H + MRDF-UIUC"; }
    uint32 getMinorVersion() const override;
    string getProperties() const override;
    string getStates() const override { return "yielding,reversal"; }
    base::Property getProperty(uint32 index) const override;
    void setProperty(uint32 index, const base::Property &p, uint32 restoreVersion = 0) override;
    ModelGqhMrdf *clone() const override { return new ModelGqhMrdf(); }

    double getConfinedModulus() const override;
    double getShearModulus() const override { return tangent_shear_; }
    double getBulkModulus() const override { return bulk_; }
    void copy(const ConstitutiveModel *mod) override;
    void run(uint32 dim, State *s) override;
    void initialize(uint32 dim, State *s) override;
    void save(std::ostream &out) const override;
    void restore(std::istream &in, uint32 restoreVersion) override;

    bool supportsHystereticDamping() const override { return true; }
    bool supportsStrengthStressRatio() const override { return true; }
    bool supportsPropertyScaling() const override { return true; }
    double getStrengthStressRatio(const SymTensor &st) const override;
    void scaleProperties(const double &scale, const std::vector<uint32> &props) override;

private:
    void updateParameters();
    void resetHistory();
    void invalidateCalibration();
    void calibrateFromReferenceCurve(uint32 dim, State *s);
    void applyMohrCoulomb(State *s, double tangent_shear);
    gqh::Parameters hysteresisParameters() const;


    static double deviatoricNorm(const double *value);
    static void deviatoricPart(const SymTensor &tensor, double *value);
    static double meanStress(const SymTensor &stress);

    // FLAC3D properties.  The first two provide the volumetric elastic response;
    // gmax is the small-strain shear modulus used by the GQ/H backbone.
    double bulk_ = 0.0;
    double shear_ = 0.0;
    double gmax_ = 1.0;
    double tau_max_ = 1.0;
    double theta1_ = 0.0;
    double theta2_ = 0.0;
    double theta3_ = 1.0;
    double theta4_ = 1.0;
    double theta5_ = 1.0;
    double mrdf_p1_ = 1.0;
    double mrdf_p2_ = 0.0;
    double mrdf_p3_ = 0.0;
    double cohesion_ = 0.0;
    double friction_ = 0.0;

    // Reference-curve calibration properties. curve_type=none preserves the
    // original direct theta/MRDF path.
    double pi_ = 0.0;
    double ocr_ = 1.0;
    std::string curve_type_ = "none";
    double cu_ = 2.0;
    std::string mrdf_mode_ = "uiuc";
    std::string fit_procedure_ = "MR";
    double sigma_m_min_ = 10.0;
    std::array<bool, 5> theta_explicit_ = {};
    bool calibration_initialized_ = false;
    bool calibration_converged_ = true;
    double calibrated_sigma_m_eff_ = 0.0;

    // Mohr-Coulomb cached parameters.
    double nph_ = 1.0;
    double csn_ = 0.0;
    double tangent_shear_ = 1.0;

    // Per-tetrahedron Masing states over the shear components (e12,e13,e23).
    // The hysteresis is a TOTAL map: the shear stress is evaluated from the
    // octahedral distance between the current shear strain and the active
    // branch anchor (last reversal point, or the origin on the virgin
    // backbone).  There is no per-increment scalar direction; reversals are
    // detected as a decrease of that distance, and nested loops close by
    // popping reversal pairs, exactly generalizing the 1-D MasingState rules
    // in gqh_math.h.  Normal deviatoric components stay elastic at gmax.
    std::vector<gqh::engine::Subzone> sub_;
    double initial_dev_stress_[6] = {};
    bool history_initialized_ = false;
    // initialize() is invoked repeatedly while FLAC3D steps (timestep and
    // multi-step reclassification), not only at solve boundaries.  The
    // "restart hysteresis at the dynamic solve" reset therefore must fire
    // exactly once; without this latch it wipes the hysteresis continuously
    // and static settling keeps erasing the earthquake history.
    bool dynamic_reset_done_ = false;
};

} // namespace models
