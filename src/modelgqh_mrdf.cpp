#include "modelgqh_mrdf.h"
#include "models/src/state.h"
#include "models/src/convert.h"
#include "version.txt"
#include "gqh_engine.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

int __stdcall DllMain(void *, unsigned, void *) {
    return 1;
}

extern "C" EXPORT_TAG const char *getName() {
#ifdef MODELDEBUG
    return "cmodel";
#else
    return "cmodel";
#endif
}

extern "C" EXPORT_TAG unsigned getMajorVersion() {
    return MAJOR_VERSION;
}

extern "C" EXPORT_TAG unsigned getMinorVersion() {
    return UPDATE_VERSION;
}

extern "C" EXPORT_TAG void *createInstance() {
    return static_cast<void *>(new models::ModelGqhMrdf());
}

namespace models {

namespace {
constexpr double kPi = 3.1415926535897932384626433832795;
// Physics constants live in gqh_engine.h; only the Mohr-Coulomb cap needs pi.
std::atomic_flag sigma_floor_warning = ATOMIC_FLAG_INIT;
std::atomic_flag explicit_theta_warning = ATOMIC_FLAG_INIT;
// FLAC3D 9 passes restoreVersion=1 to restore(), not the plug-in minor version,
// so the archive carries its own tag (Itasca's cmodels do the same).  Archives
// written before the tag start directly with initial_dev_stress_.
constexpr uint32 kSaveMagic = 0x31485147u;  // "GQH1"
constexpr uint32 kSaveFormat = 4;           // 4: + calibration state
}

ModelGqhMrdf::ModelGqhMrdf(unsigned short option)
    : ConstitutiveModel(option) {
}

uint32 ModelGqhMrdf::getMinorVersion() const {
    return UPDATE_VERSION;
}

string ModelGqhMrdf::getProperties() const {
    return "bulk,shear,gmax,tau_max,theta1,theta2,theta3,theta4,theta5,"
           "mrdf_p1,mrdf_p2,mrdf_p3,cohesion,friction,"
           "pi,ocr,curve_type,cu,mrdf_mode,fit_procedure,sigma_m_min,"
           "d_gamma_abs,d_stack,d_map_sxz,d_e13,d_tangent";
}

base::Property ModelGqhMrdf::getProperty(uint32 index) const {
    switch (index) {
    case 1: return bulk_;
    case 2: return shear_;
    case 3: return gmax_;
    case 4: return tau_max_;
    case 5: return theta1_;
    case 6: return theta2_;
    case 7: return theta3_;
    case 8: return theta4_;
    case 9: return theta5_;
    case 10: return mrdf_p1_;
    case 11: return mrdf_p2_;
    case 12: return mrdf_p3_;
    case 13: return cohesion_;
    case 14: return friction_;
    case 15: return pi_;
    case 16: return ocr_;
    case 17: return curve_type_;
    case 18: return cu_;
    case 19: return mrdf_mode_;
    case 20: return fit_procedure_;
    case 21: return sigma_m_min_;
    // Read-only diagnostics of the per-tetrahedron Masing states.
    case 22: {
        double gmaxseen = 0.0;
        for (const gqh::engine::Subzone &sz : sub_)
            gmaxseen = std::max(gmaxseen,
                                gqh::engine::gammaOctBetween(sz.strain, sz.origin_strain));
        return gmaxseen;
    }
    case 23: {
        std::size_t depth = 0;
        for (const gqh::engine::Subzone &sz : sub_)
            depth = std::max(depth, sz.stack.size());
        return static_cast<double>(depth);
    }
    case 24: return sub_.size() > 1 ? sub_[1].map_stress[1]
             : (sub_.empty() ? 0.0 : sub_[0].map_stress[1]);
    case 25: return sub_.size() > 1 ? sub_[1].strain[1]
             : (sub_.empty() ? 0.0 : sub_[0].strain[1]);
    case 26: return tangent_shear_;
    default: return 0.0;
    }
}

void ModelGqhMrdf::setProperty(uint32 index, const base::Property &p, uint32 restoreVersion) {
    ConstitutiveModel::setProperty(index, p, restoreVersion);
    switch (index) {
    case 1: bulk_ = p.to<double>(); break;
    case 2: shear_ = p.to<double>(); break;
    case 3: gmax_ = p.to<double>(); break;
    case 4: tau_max_ = p.to<double>(); break;
    case 5: theta1_ = p.to<double>(); break;
    case 6: theta2_ = p.to<double>(); break;
    case 7: theta3_ = p.to<double>(); break;
    case 8: theta4_ = p.to<double>(); break;
    case 9: theta5_ = p.to<double>(); break;
    case 10: mrdf_p1_ = p.to<double>(); break;
    case 11: mrdf_p2_ = p.to<double>(); break;
    case 12: mrdf_p3_ = p.to<double>(); break;
    case 13: cohesion_ = p.to<double>(); break;
    case 14: friction_ = p.to<double>(); break;
    case 15: pi_ = p.to<double>(); break;
    case 16: ocr_ = p.to<double>(); break;
    case 17: curve_type_ = gqh::refcurve::normalizeCurveType(p.to<string>()); break;
    case 18: cu_ = p.to<double>(); break;
    case 19: mrdf_mode_ = gqh::refcurve::normalizeMrdfMode(p.to<string>()); break;
    case 20: fit_procedure_ = gqh::refcurve::normalizeFitProcedure(p.to<string>()); break;
    case 21: sigma_m_min_ = p.to<double>(); break;
    default: break;
    }
    // `model restore` replays every property through here with restoreVersion=1.
    // The values are the saved ones (calibrated theta included); the hysteresis
    // and calibration state come from restore().  Treating the replay as a user
    // edit would wipe that history and mark every theta as explicit.
    if (restoreVersion != 0)
        return;
    if (index >= 5 && index <= 9)
        theta_explicit_[static_cast<std::size_t>(index - 5)] = true;
    if (index >= 1 && index <= 21) {
        invalidateCalibration();
        resetHistory();
    }
}

double ModelGqhMrdf::getConfinedModulus() const {
    return bulk_ + (4.0 / 3.0) * std::max(gmax_, tangent_shear_);
}

gqh::Parameters ModelGqhMrdf::hysteresisParameters() const {
    gqh::Parameters p;
    p.gmax = gmax_;
    p.tau_max = tau_max_;
    p.theta1 = theta1_;
    p.theta2 = theta2_;
    p.theta3 = theta3_;
    p.theta4 = theta4_;
    p.theta5 = theta5_;
    p.mrdf_p1 = mrdf_p1_;
    p.mrdf_p2 = mrdf_p2_;
    p.mrdf_p3 = mrdf_p3_;
    return p;
}

void ModelGqhMrdf::updateParameters() {
    if (bulk_ <= 0.0)
        throw std::runtime_error("GQ/H-MRDF: bulk must be > 0.");
    if (shear_ <= 0.0)
        shear_ = gmax_;
    if (gmax_ <= 0.0)
        gmax_ = shear_;
    if (gmax_ <= 0.0)
        throw std::runtime_error("GQ/H-MRDF: gmax must be > 0.");
    if (tau_max_ <= 0.0)
        throw std::runtime_error("GQ/H-MRDF: tau_max must be > 0.");
    if (theta3_ <= 0.0 || theta4_ <= 0.0)
        throw std::runtime_error("GQ/H-MRDF: theta3 and theta4 must be > 0.");
    if (theta5_ < 0.0)
        throw std::runtime_error("GQ/H-MRDF: theta5 must be >= 0.");
    if (cohesion_ < 0.0)
        throw std::runtime_error("GQ/H-MRDF: cohesion must be >= 0.");
    if (friction_ < 0.0 || friction_ >= 90.0)
        throw std::runtime_error("GQ/H-MRDF: friction must be in [0,90) degrees.");
    if (pi_ < 0.0 || !std::isfinite(pi_))
        throw std::runtime_error("GQ/H-MRDF: pi must be >= 0.");
    if (ocr_ <= 0.0 || !std::isfinite(ocr_))
        throw std::runtime_error("GQ/H-MRDF: ocr must be > 0.");
    if (cu_ <= 0.0 || !std::isfinite(cu_))
        throw std::runtime_error("GQ/H-MRDF: cu must be > 0.");
    if (sigma_m_min_ <= 0.0 || !std::isfinite(sigma_m_min_))
        throw std::runtime_error("GQ/H-MRDF: sigma_m_min must be > 0.");
    curve_type_ = gqh::refcurve::normalizeCurveType(curve_type_);
    mrdf_mode_ = gqh::refcurve::normalizeMrdfMode(mrdf_mode_);
    fit_procedure_ = gqh::refcurve::normalizeFitProcedure(fit_procedure_);
    if (!gqh::refcurve::isSupportedCurveType(curve_type_))
        throw std::runtime_error("GQ/H-MRDF: unsupported curve_type '" + curve_type_ + "'.");
    if (mrdf_mode_ != "uiuc" && mrdf_mode_ != "darendeli")
        throw std::runtime_error("GQ/H-MRDF: mrdf_mode must be uiuc or darendeli.");
    if (fit_procedure_ != "MR" && fit_procedure_ != "MRD" && fit_procedure_ != "DC")
        throw std::runtime_error("GQ/H-MRDF: fit_procedure must be MR, MRD, or DC.");

    const double sin_phi = std::sin(friction_ * kPi / 180.0);
    nph_ = (1.0 + sin_phi) / (1.0 - sin_phi);
    csn_ = 2.0 * cohesion_ * std::sqrt(nph_);
    tangent_shear_ = gmax_;
}

void ModelGqhMrdf::resetHistory() {
    sub_.clear();
    std::fill(std::begin(initial_dev_stress_), std::end(initial_dev_stress_), 0.0);
    history_initialized_ = false;
    dynamic_reset_done_ = false;
    tangent_shear_ = gmax_;
}

void ModelGqhMrdf::invalidateCalibration() {
    calibration_initialized_ = false;
    calibration_converged_ = true;
    calibrated_sigma_m_eff_ = 0.0;
}

void ModelGqhMrdf::calibrateFromReferenceCurve(uint32 dim, State *s) {
    if (curve_type_ == "none" || calibration_initialized_)
        return;

    const double sig_m = -(s->stnS_.s11() + s->stnS_.s22() + s->stnS_.s33()) / 3.0;
    const double pp = s->getAveragePP();
    double sig_m_eff = sig_m - pp;
    const bool used_floor = sig_m_eff < sigma_m_min_;
    if (used_floor) {
        sig_m_eff = sigma_m_min_;
        if (!sigma_floor_warning.test_and_set())
            std::cerr << "GQ/H refcurve: sigma_m_eff was below sigma_m_min; "
                      << "using the configured floor once for this model.\n";
    }

    gqh::refcurve::CalibrationInput input;
    input.gmax = gmax_;
    input.tau_max = tau_max_;
    input.plasticity_index = pi_;
    input.ocr = ocr_;
    input.mean_effective_stress_kpa = sig_m_eff;
    input.coef_unif = cu_;
    input.curve_type = curve_type_;
    input.mrdf_mode = mrdf_mode_;
    input.fit_procedure = fit_procedure_;
    input.theta = {theta1_, theta2_, theta3_, theta4_, theta5_};
    input.fixed_theta = theta_explicit_;

    const gqh::refcurve::CalibrationResult result = gqh::refcurve::calibrate(input);
    theta1_ = result.theta[0];
    theta2_ = result.theta[1];
    theta3_ = result.theta[2];
    theta4_ = result.theta[3];
    theta5_ = result.theta[4];
    mrdf_p1_ = result.mrdf[0];
    mrdf_p2_ = result.mrdf[1];
    mrdf_p3_ = result.mrdf[2];
    calibrated_sigma_m_eff_ = sig_m_eff;
    calibration_initialized_ = true;
    calibration_converged_ = result.converged;
    if (used_floor)
        s->iworking_[0] |= 2;
    if (!result.converged) {
        s->iworking_[0] |= 1;
        std::cerr << "GQ/H refcurve: bounded fit reached its iteration limit; "
                  << "best candidate retained (iterations=" << result.iterations << ").\n";
    }
    bool explicit_theta = false;
    for (const bool is_explicit : theta_explicit_)
        explicit_theta = explicit_theta || is_explicit;
    if (explicit_theta && !explicit_theta_warning.test_and_set())
        std::cerr << "GQ/H refcurve: explicitly supplied theta value(s) override "
                  << "the reference-curve fit.\n";
    std::cerr << "GQ/H refcurve: zone=" << s->getZoneID()
              << " curve=" << result.curve_name
              << " sig_m=" << sig_m
              << " pp=" << pp
              << " sig_m_eff=" << sig_m_eff
              << " floor=" << (used_floor ? 1 : 0)
              << " theta=[" << theta1_ << "," << theta2_ << "," << theta3_
              << "," << theta4_ << "," << theta5_ << "]"
              << " mrdf=[" << mrdf_p1_ << "," << mrdf_p2_ << "," << mrdf_p3_ << "]"
              << " fit_procedure=" << fit_procedure_
              << " mrdf_mode=" << mrdf_mode_
              << " modulus_rmse=" << result.modulus_rmse
              << " damping_rmse=" << result.damping_rmse
              << " strength_ratio=" << result.strength_ratio
              << " iterations=" << result.iterations
              << " converged=" << (result.converged ? 1 : 0) << "\n";
    (void)dim;
}

void ModelGqhMrdf::initialize(uint32 dim, State *s) {
    const bool history_exists = history_initialized_;
    const bool calibration_exists = calibration_initialized_;
    ConstitutiveModel::initialize(dim, s);
    updateParameters();
    if (!history_exists)
        resetHistory();
    calibrateFromReferenceCurve(dim, s);
    if (!calibration_exists && calibration_initialized_)
        resetHistory();
    if (history_exists && !dynamic_reset_done_ &&
        deviatoricNorm(initial_dev_stress_) > 1.0e-8) {
        // FLAC3D calls initialize() again when the dynamic solve starts.  The
        // gravity solve may have accumulated shear history; retain the strain
        // and stress baseline of every tetrahedron but restart the hysteresis
        // from the current point, so static settling is not treated as
        // earthquake history.  One-shot: initialize() also fires repeatedly
        // during stepping (timestep updates), which must not wipe the state.
        // (With a zero initial deviatoric stress -- the element tests -- this
        // block never fires and cycles keep accumulating across commands.)
        for (gqh::engine::Subzone &sz : sub_)
            gqh::engine::restartAtCurrentPoint(sz);
        tangent_shear_ = gmax_;
        dynamic_reset_done_ = true;
    }
}

double ModelGqhMrdf::deviatoricNorm(const double *value) {
    const double norm2 = value[0] * value[0] + value[1] * value[1] + value[2] * value[2] +
                         2.0 * (value[3] * value[3] + value[4] * value[4] + value[5] * value[5]);
    return std::sqrt(std::max(norm2, 0.0));
}

void ModelGqhMrdf::deviatoricPart(const SymTensor &tensor, double *value) {
    const double mean = (tensor.s11() + tensor.s22() + tensor.s33()) / 3.0;
    value[0] = tensor.s11() - mean;
    value[1] = tensor.s22() - mean;
    value[2] = tensor.s33() - mean;
    value[3] = tensor.s12();
    value[4] = tensor.s13();
    value[5] = tensor.s23();
}

double ModelGqhMrdf::meanStress(const SymTensor &stress) {
    return (stress.s11() + stress.s22() + stress.s33()) / 3.0;
}

void ModelGqhMrdf::run(uint32 dim, State *s) {
    ConstitutiveModel::run(dim, s);

    // FLAC3D 9 calls run() once per tetrahedral sub-zone.  The ten tetrahedra
    // of a brick follow different strain paths whenever the zone does more
    // than pure shear (measured 01.09.2026 with the d_gsub* diagnostics: peak
    // strains coincide within ~12 % in a pure-shear column but spread by up to
    // 3x, central/corner ~0.63, when the column also bends), so every
    // tetrahedron carries its own Masing state, exactly like built-in models
    // keep per-tet stress.  (An earlier "~1e-3" figure here was not
    // reproduced and has been withdrawn.)
    const uint32 nsub = std::max<uint32>(1, s->total_sub_zones_);
    if (sub_.size() != nsub)
        sub_.assign(nsub, gqh::engine::Subzone{});
    gqh::engine::Subzone &sz = sub_[std::min<uint32>(s->sub_zone_, nsub - 1)];

    double increment_dev[6];
    deviatoricPart(s->stnE_, increment_dev);
    double current_dev[6];
    deviatoricPart(s->stnS_, current_dev);

    if (!sz.initialized) {
        for (int i = 0; i < 3; ++i) {
            sz.origin_stress[i] = current_dev[3 + i];
            sz.map_stress[i] = current_dev[3 + i];
        }
        sz.initialized = true;
    }
    if (!history_initialized_) {
        std::copy(std::begin(current_dev), std::end(current_dev),
                  std::begin(initial_dev_stress_));
        history_initialized_ = true;
    }

    double shear_old[3];
    std::copy(std::begin(sz.strain), std::end(sz.strain), std::begin(shear_old));
    for (int i = 0; i < 3; ++i)
        sz.strain[i] += increment_dev[3 + i];

    double map_new[3];
    double tangent = tangent_shear_;
    gqh::engine::advance(sz, shear_old, sz.strain, hysteresisParameters(),
                             gmax_, map_new, &tangent);

    // Apply the map incrementally on FLAC3D's current stress so external
    // adjustments between constitutive calls (Maxwell damping, boundary
    // corrections) are preserved; with none present this equals the total
    // map because map_stress tracks the map exactly.  Normal deviatoric
    // components respond elastically at gmax: the modulus reduction is a
    // shear phenomenon, and softening the normal response lets vertical
    // settling of a softened zone hijack the map direction.
    const double mean = meanStress(s->stnS_) +
                        bulk_ * (s->stnE_.s11() + s->stnE_.s22() + s->stnE_.s33());
    s->stnS_.rs11() = mean + current_dev[0] + 2.0 * gmax_ * increment_dev[0];
    s->stnS_.rs22() = mean + current_dev[1] + 2.0 * gmax_ * increment_dev[1];
    s->stnS_.rs33() = mean + current_dev[2] + 2.0 * gmax_ * increment_dev[2];
    s->stnS_.rs12() = current_dev[3] + (map_new[0] - sz.map_stress[0]);
    s->stnS_.rs13() = current_dev[4] + (map_new[1] - sz.map_stress[1]);
    s->stnS_.rs23() = current_dev[5] + (map_new[2] - sz.map_stress[2]);
    std::copy(std::begin(map_new), std::end(map_new), std::begin(sz.map_stress));
    tangent_shear_ = tangent;

    s->viscous_ = true;
    s->modulus_changed_ = true;
    applyMohrCoulomb(s, tangent);
    (void)dim;
}

void ModelGqhMrdf::applyMohrCoulomb(State *s, double tangent_shear) {
    if (!canFail())
        return;

    SymTensorInfo info;
    DVect3 principal = s->stnS_.getEigenInfo(&info);
    const double failure = nph_ * principal.z() - principal.x() - csn_;
    if (failure <= 0.0)
        return;

    // Zero dilation Mohr-Coulomb return.  With e1 - e2 = 2G the correction
    // ratios reduce to 1/(1+nph) and -1/(1+nph), independent of the tangent.
    // Computing them from e1 - e2 gives 0/0 = NaN once GQ/H saturates: the
    // tangent floor (1e-12) vanishes against bulk_ in double precision.
    (void)tangent_shear;
    const double sc1 = 1.0 / (1.0 + nph_);
    const double sc3 = -sc1;
    principal.rx() += failure * sc1;
    principal.ry() += 0.0;
    principal.rz() += failure * sc3;
    s->stnS_ = info.resolve(principal);
    s->state_ |= shear_now;
    s->viscous_ = false;
}

void ModelGqhMrdf::copy(const ConstitutiveModel *mod) {
    ConstitutiveModel::copy(mod);
    const ModelGqhMrdf *other = dynamic_cast<const ModelGqhMrdf *>(mod);
    if (!other)
        throw std::runtime_error("GQ/H-MRDF: constitutive model cast failed.");

    bulk_ = other->bulk_;
    shear_ = other->shear_;
    gmax_ = other->gmax_;
    tau_max_ = other->tau_max_;
    theta1_ = other->theta1_;
    theta2_ = other->theta2_;
    theta3_ = other->theta3_;
    theta4_ = other->theta4_;
    theta5_ = other->theta5_;
    mrdf_p1_ = other->mrdf_p1_;
    mrdf_p2_ = other->mrdf_p2_;
    mrdf_p3_ = other->mrdf_p3_;
    cohesion_ = other->cohesion_;
    friction_ = other->friction_;
    pi_ = other->pi_;
    ocr_ = other->ocr_;
    curve_type_ = other->curve_type_;
    cu_ = other->cu_;
    mrdf_mode_ = other->mrdf_mode_;
    fit_procedure_ = other->fit_procedure_;
    sigma_m_min_ = other->sigma_m_min_;
    theta_explicit_ = other->theta_explicit_;
    calibration_initialized_ = other->calibration_initialized_;
    calibration_converged_ = other->calibration_converged_;
    calibrated_sigma_m_eff_ = other->calibrated_sigma_m_eff_;
    nph_ = other->nph_;
    csn_ = other->csn_;
    tangent_shear_ = other->tangent_shear_;
    sub_ = other->sub_;
    std::copy(std::begin(other->initial_dev_stress_), std::end(other->initial_dev_stress_),
              std::begin(initial_dev_stress_));
    history_initialized_ = other->history_initialized_;
    dynamic_reset_done_ = other->dynamic_reset_done_;
}

void ModelGqhMrdf::save(std::ostream &out) const {
    ConstitutiveModel::save(out);
    out.write(reinterpret_cast<const char *>(&kSaveMagic), sizeof(kSaveMagic));
    out.write(reinterpret_cast<const char *>(&kSaveFormat), sizeof(kSaveFormat));
    out.write(reinterpret_cast<const char *>(initial_dev_stress_), sizeof(initial_dev_stress_));
    out.write(reinterpret_cast<const char *>(&history_initialized_), sizeof(history_initialized_));
    out.write(reinterpret_cast<const char *>(&dynamic_reset_done_), sizeof(dynamic_reset_done_));
    const std::size_t nsub = sub_.size();
    out.write(reinterpret_cast<const char *>(&nsub), sizeof(nsub));
    for (const gqh::engine::Subzone &sz : sub_) {
        out.write(reinterpret_cast<const char *>(sz.strain), sizeof(sz.strain));
        out.write(reinterpret_cast<const char *>(sz.origin_strain), sizeof(sz.origin_strain));
        out.write(reinterpret_cast<const char *>(sz.origin_stress), sizeof(sz.origin_stress));
        out.write(reinterpret_cast<const char *>(sz.map_stress), sizeof(sz.map_stress));
        out.write(reinterpret_cast<const char *>(sz.peak_strain), sizeof(sz.peak_strain));
        out.write(reinterpret_cast<const char *>(sz.peak_stress), sizeof(sz.peak_stress));
        out.write(reinterpret_cast<const char *>(&sz.peak_distance), sizeof(sz.peak_distance));
        out.write(reinterpret_cast<const char *>(&sz.gamma_max), sizeof(sz.gamma_max));
        out.write(reinterpret_cast<const char *>(&sz.initialized), sizeof(sz.initialized));
        const std::size_t count = sz.stack.size();
        out.write(reinterpret_cast<const char *>(&count), sizeof(count));
        if (count)
            out.write(reinterpret_cast<const char *>(sz.stack.data()),
                      static_cast<std::streamsize>(count * sizeof(gqh::engine::Reversal)));
    }
    // Without this block restore() recalibrates at the post-shaking stress
    // and initialize() then resets the history.
    out.write(reinterpret_cast<const char *>(&calibration_initialized_),
              sizeof(calibration_initialized_));
    out.write(reinterpret_cast<const char *>(&calibration_converged_),
              sizeof(calibration_converged_));
    out.write(reinterpret_cast<const char *>(&calibrated_sigma_m_eff_),
              sizeof(calibrated_sigma_m_eff_));
    for (const bool is_explicit : theta_explicit_)
        out.write(reinterpret_cast<const char *>(&is_explicit), sizeof(is_explicit));
}

void ModelGqhMrdf::restore(std::istream &in, uint32 restoreVersion) {
    ConstitutiveModel::restore(in, restoreVersion);
    uint32 tag = 0;
    in.read(reinterpret_cast<char *>(&tag), sizeof(tag));
    const bool legacy = tag != kSaveMagic;
    if (legacy) {
        // Untagged minor-version-3 archive: the four bytes already read are the
        // start of initial_dev_stress_.
        char *buffer = reinterpret_cast<char *>(initial_dev_stress_);
        std::copy(reinterpret_cast<const char *>(&tag),
                  reinterpret_cast<const char *>(&tag) + sizeof(tag), buffer);
        in.read(buffer + sizeof(tag),
                static_cast<std::streamsize>(sizeof(initial_dev_stress_) - sizeof(tag)));
    } else {
        uint32 format = 0;
        in.read(reinterpret_cast<char *>(&format), sizeof(format));
        if (format != kSaveFormat)
            throw std::runtime_error(
                "GQ/H-MRDF restore: unsupported save format; use the DLL that wrote the file.");
        in.read(reinterpret_cast<char *>(initial_dev_stress_), sizeof(initial_dev_stress_));
    }
    in.read(reinterpret_cast<char *>(&history_initialized_), sizeof(history_initialized_));
    in.read(reinterpret_cast<char *>(&dynamic_reset_done_), sizeof(dynamic_reset_done_));
    std::size_t nsub = 0;
    in.read(reinterpret_cast<char *>(&nsub), sizeof(nsub));
    if (nsub > 100)
        throw std::runtime_error("GQ/H-MRDF restore: subzone count is invalid.");
    sub_.assign(nsub, gqh::engine::Subzone{});
    for (gqh::engine::Subzone &sz : sub_) {
        in.read(reinterpret_cast<char *>(sz.strain), sizeof(sz.strain));
        in.read(reinterpret_cast<char *>(sz.origin_strain), sizeof(sz.origin_strain));
        in.read(reinterpret_cast<char *>(sz.origin_stress), sizeof(sz.origin_stress));
        in.read(reinterpret_cast<char *>(sz.map_stress), sizeof(sz.map_stress));
        in.read(reinterpret_cast<char *>(sz.peak_strain), sizeof(sz.peak_strain));
        in.read(reinterpret_cast<char *>(sz.peak_stress), sizeof(sz.peak_stress));
        in.read(reinterpret_cast<char *>(&sz.peak_distance), sizeof(sz.peak_distance));
        in.read(reinterpret_cast<char *>(&sz.gamma_max), sizeof(sz.gamma_max));
        in.read(reinterpret_cast<char *>(&sz.initialized), sizeof(sz.initialized));
        std::size_t count = 0;
        in.read(reinterpret_cast<char *>(&count), sizeof(count));
        if (count > 100000)
            throw std::runtime_error("GQ/H-MRDF restore: reversal stack is invalid.");
        sz.stack.resize(count);
        if (count)
            in.read(reinterpret_cast<char *>(sz.stack.data()),
                    static_cast<std::streamsize>(count * sizeof(gqh::engine::Reversal)));
    }
    if (legacy) {
        // The saved theta/MRDF properties are already the calibrated values,
        // and a zone with history has been through initialize() (calibration
        // runs before the first run()).  Explicit-theta flags are unknown.
        calibration_initialized_ = history_initialized_;
        calibration_converged_ = true;
        calibrated_sigma_m_eff_ = 0.0;
        theta_explicit_.fill(false);
    } else {
        in.read(reinterpret_cast<char *>(&calibration_initialized_),
                sizeof(calibration_initialized_));
        in.read(reinterpret_cast<char *>(&calibration_converged_),
                sizeof(calibration_converged_));
        in.read(reinterpret_cast<char *>(&calibrated_sigma_m_eff_),
                sizeof(calibrated_sigma_m_eff_));
        for (bool &is_explicit : theta_explicit_)
            in.read(reinterpret_cast<char *>(&is_explicit), sizeof(is_explicit));
    }
    if (!in)
        throw std::runtime_error("GQ/H-MRDF restore: truncated model state.");
    updateParameters();
}

double ModelGqhMrdf::getStrengthStressRatio(const SymTensor &stress) const {
    const DVect3 principal = stress.getEigenInfo();
    const double denominator = principal.z() - principal.x();
    if (denominator <= 0.0)
        return 10.0;
    const double sigma1_failure = nph_ * principal.z() - csn_;
    return std::min(10.0, std::max(0.0,
        (principal.z() - sigma1_failure) / denominator));
}

void ModelGqhMrdf::scaleProperties(const double &scale, const std::vector<uint32> &props) {
    for (const uint32 property : props) {
        switch (property) {
        case 1: bulk_ *= scale; break;
        case 2: shear_ *= scale; break;
        case 3: gmax_ *= scale; break;
        case 4: tau_max_ *= scale; break;
        case 13: cohesion_ *= scale; break;
        case 14:
            friction_ = std::max(0.0, std::min(89.0,
                std::atan(std::tan(friction_ * kPi / 180.0) * scale) * 180.0 / kPi));
            break;
        default: break;
        }
    }
    setValid(0);
}

} // namespace models
