#pragma once

// FLAC3D 6.0 adapter for the GQ/H + MRDF-UIUC model.
//
// The physics lives in ../src/gqh_engine.h and is SHARED with the v9 plug-in --
// do not re-implement it here.  This file only translates the v6 host API
// (String/Variant/UInt/UByte, getStressStrengthRatio, ...) onto that engine.
//
// v6 vs v9 differences handled here:
//   * String is std::basic_string<Char> (wide) -> model/property names use L"".
//   * Variant instead of base::Property      -> p.toDouble() / p.toString().
//   * run/initialize take UByte, not uint32.
//   * getStressStrengthRatio  (v9: getStrengthStressRatio).
//   * supportsStressStrengthRatio (v9: supportsStrengthStressRatio).
//   * clone() uses the NEW() allocation macro.
//
// NOT ported: the reference-curve auto-calibration (curve_type/pi/ocr/...).
// That path pulls std::string plumbing and refcurve.cpp; the v6 build exposes
// the direct theta/MRDF parameters only.  Supplying curve_type in v6 is an
// error rather than a silent no-op.

// v6 SDK'da baslik "cmodels/src/conmodel.h" altindadir; v9'un
// "models/src/..." yolu BURADA CALISMAZ. build_v6.ps1 cmodels\src'i
// include yoluna koyar, bu yuzden duz ad yeterlidir.
#include "conmodel.h"

#include <vector>

#include "../src/gqh_engine.h"

namespace models {

class ModelGqhMrdfV6 : public ConstitutiveModel {
public:
    ModelGqhMrdfV6();

    virtual String getName() const { return L"gqh-mrdf"; }
    virtual String getFullName() const { return L"GQ/H + MRDF-UIUC"; }
    // Override, do not inherit: the v6 SDK header declares getPluginName inline
    // but the shipped conmodel006_64.lib does not export a definition, so
    // relying on the base class leaves an unresolved external at link time.
    virtual String getPluginName() const { return getName(); }
    virtual UInt getMinorVersion() const;
    virtual String getProperties() const {
        return L"bulk,shear,gmax,tau_max,theta1,theta2,theta3,theta4,theta5,"
               L"mrdf_p1,mrdf_p2,mrdf_p3,cohesion,friction";
    }
    virtual String getStates() const { return L"yielding,reversal"; }
    virtual Variant getProperty(UInt index) const;
    virtual void setProperty(UInt index, const Variant &p, UInt restoreVersion = 0);
    virtual ModelGqhMrdfV6 *clone() const { return NEW("ModelGqhMrdfV6") ModelGqhMrdfV6(); }

    virtual Double getConfinedModulus() const;
    virtual Double getShearModulus() const { return tangent_shear_; }
    virtual Double getBulkModulus() const { return bulk_; }
    virtual void copy(const ConstitutiveModel *mod);
    virtual void run(UByte dim, State *s);
    virtual void initialize(UByte dim, State *s);
    virtual void save(std::ostream &o) const;
    virtual void restore(std::istream &i, UInt restoreVersion);

    virtual bool supportsHystereticDamping() const { return true; }
    virtual bool supportsStressStrengthRatio() const { return true; }
    virtual bool supportsPropertyScaling() const { return true; }
    virtual Double getStressStrengthRatio(const SymTensor &st) const;
    virtual void scaleProperties(const Double &scale, const std::vector<UInt> &props);

private:
    void updateParameters();
    void resetHistory();
    void applyMohrCoulomb(State *s, Double tangent_shear);
    gqh::Parameters hysteresisParameters() const;

    static Double deviatoricNorm(const Double *value);
    static void deviatoricPart(const SymTensor &tensor, Double *value);
    static Double meanStress(const SymTensor &stress);

    Double bulk_, shear_, gmax_, tau_max_;
    Double theta1_, theta2_, theta3_, theta4_, theta5_;
    Double mrdf_p1_, mrdf_p2_, mrdf_p3_;
    Double cohesion_, friction_;
    Double nph_, csn_, tangent_shear_;

    std::vector<gqh::engine::Subzone> sub_;
    Double initial_dev_stress_[6];
    bool history_initialized_;
    bool dynamic_reset_done_;
};

} // namespace models
