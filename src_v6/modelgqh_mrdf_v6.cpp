#include "modelgqh_mrdf_v6.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace models {

namespace {
const Double kPi = 3.1415926535897932384626433832795;
}

// FLAC3D 6 discovers plug-in models through this factory symbol.
extern "C" __declspec(dllexport) ConstitutiveModel *createInstance() {
    return new ModelGqhMrdfV6();
}

extern "C" __declspec(dllexport) unsigned getMajorVersion() { return 6; }
extern "C" __declspec(dllexport) unsigned getMinorVersion() { return 1; }

// v6 probes a C-level getName() when loading the library (absent -> "Failed to
// find function getName()").  It must return the PLUG-IN CATEGORY, and for a
// constitutive model that string is literally "cmodel" -- returning anything
// else fails with 'name (X) failed to match (cmodel)'.  The DLL is therefore
// named cmodelgqh-mrdf006.dll, mirroring the v9 build.  The model's own name
// ("gqh-mrdf", used on the command line) comes from the class, not from here.
extern "C" __declspec(dllexport) const char *getName() { return "cmodel"; }

ModelGqhMrdfV6::ModelGqhMrdfV6()
    : bulk_(0.0), shear_(0.0), gmax_(1.0), tau_max_(1.0),
      theta1_(0.0), theta2_(0.0), theta3_(1.0), theta4_(1.0), theta5_(1.0),
      mrdf_p1_(1.0), mrdf_p2_(0.0), mrdf_p3_(0.0),
      cohesion_(0.0), friction_(0.0),
      nph_(1.0), csn_(0.0), tangent_shear_(1.0),
      history_initialized_(false), dynamic_reset_done_(false) {
    std::fill(initial_dev_stress_, initial_dev_stress_ + 6, 0.0);
}

UInt ModelGqhMrdfV6::getMinorVersion() const { return 1; }

Variant ModelGqhMrdfV6::getProperty(UInt index) const {
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
    default: return 0.0;
    }
}

void ModelGqhMrdfV6::setProperty(UInt index, const Variant &p, UInt restoreVersion) {
    ConstitutiveModel::setProperty(index, p, restoreVersion);
    switch (index) {
    case 1: bulk_ = p.toDouble(); break;
    case 2: shear_ = p.toDouble(); break;
    case 3: gmax_ = p.toDouble(); break;
    case 4: tau_max_ = p.toDouble(); break;
    case 5: theta1_ = p.toDouble(); break;
    case 6: theta2_ = p.toDouble(); break;
    case 7: theta3_ = p.toDouble(); break;
    case 8: theta4_ = p.toDouble(); break;
    case 9: theta5_ = p.toDouble(); break;
    case 10: mrdf_p1_ = p.toDouble(); break;
    case 11: mrdf_p2_ = p.toDouble(); break;
    case 12: mrdf_p3_ = p.toDouble(); break;
    case 13: cohesion_ = p.toDouble(); break;
    case 14: friction_ = p.toDouble(); break;
    default: break;
    }
    // `model restore` replays every property through here with restoreVersion=1;
    // only an interactive change may reset the history (Itasca cmodels use the
    // same `if (!restoreVersion)` guard).  FLAC3D 6 replays them before restore().
    if (!restoreVersion && index >= 1 && index <= 14)
        resetHistory();
}

Double ModelGqhMrdfV6::getConfinedModulus() const {
    return bulk_ + (4.0 / 3.0) * std::max(gmax_, tangent_shear_);
}

gqh::Parameters ModelGqhMrdfV6::hysteresisParameters() const {
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

void ModelGqhMrdfV6::updateParameters() {
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

    const Double sin_phi = std::sin(friction_ * kPi / 180.0);
    nph_ = (1.0 + sin_phi) / (1.0 - sin_phi);
    csn_ = 2.0 * cohesion_ * std::sqrt(nph_);
    tangent_shear_ = gmax_;
}

void ModelGqhMrdfV6::resetHistory() {
    sub_.clear();
    std::fill(initial_dev_stress_, initial_dev_stress_ + 6, 0.0);
    history_initialized_ = false;
    dynamic_reset_done_ = false;
    tangent_shear_ = gmax_;
}

Double ModelGqhMrdfV6::deviatoricNorm(const Double *value) {
    const Double norm2 = value[0] * value[0] + value[1] * value[1] + value[2] * value[2] +
                         2.0 * (value[3] * value[3] + value[4] * value[4] + value[5] * value[5]);
    return std::sqrt(std::max(norm2, 0.0));
}

void ModelGqhMrdfV6::deviatoricPart(const SymTensor &tensor, Double *value) {
    const Double mean = (tensor.s11() + tensor.s22() + tensor.s33()) / 3.0;
    value[0] = tensor.s11() - mean;
    value[1] = tensor.s22() - mean;
    value[2] = tensor.s33() - mean;
    value[3] = tensor.s12();
    value[4] = tensor.s13();
    value[5] = tensor.s23();
}

Double ModelGqhMrdfV6::meanStress(const SymTensor &stress) {
    return (stress.s11() + stress.s22() + stress.s33()) / 3.0;
}

void ModelGqhMrdfV6::initialize(UByte dim, State *s) {
    const bool history_exists = history_initialized_;
    ConstitutiveModel::initialize(dim, s);
    updateParameters();
    if (!history_exists)
        resetHistory();
    // One-shot restart when the dynamic solve begins: keep the strain/stress
    // baseline of every tetrahedron but drop the static settling history.
    // initialize() fires repeatedly during stepping, hence the latch.
    if (history_exists && !dynamic_reset_done_ &&
        deviatoricNorm(initial_dev_stress_) > 1.0e-8) {
        for (std::size_t i = 0; i < sub_.size(); ++i)
            gqh::engine::restartAtCurrentPoint(sub_[i]);
        tangent_shear_ = gmax_;
        dynamic_reset_done_ = true;
    }
}

void ModelGqhMrdfV6::run(UByte dim, State *s) {
    ConstitutiveModel::run(dim, s);

    // One Masing state per tetrahedron -- see gqh_engine.h.  Driving a single
    // zone-level history from sub-zone 0 starves the model in a column while
    // passing every uniform-strain element test.
    const UInt nsub = std::max<UInt>(1, s->total_sub_zones_);
    if (sub_.size() != nsub)
        sub_.assign(nsub, gqh::engine::Subzone());
    gqh::engine::Subzone &sz = sub_[std::min<UInt>(s->sub_zone_, nsub - 1)];

    Double increment_dev[6];
    deviatoricPart(s->stnE_, increment_dev);
    Double current_dev[6];
    deviatoricPart(s->stnS_, current_dev);

    if (!sz.initialized) {
        for (int i = 0; i < 3; ++i) {
            sz.origin_stress[i] = current_dev[3 + i];
            sz.map_stress[i] = current_dev[3 + i];
        }
        sz.initialized = true;
    }
    if (!history_initialized_) {
        std::copy(current_dev, current_dev + 6, initial_dev_stress_);
        history_initialized_ = true;
    }

    Double shear_old[3];
    std::copy(sz.strain, sz.strain + 3, shear_old);
    for (int i = 0; i < 3; ++i)
        sz.strain[i] += increment_dev[3 + i];

    Double map_new[3];
    Double tangent = tangent_shear_;
    gqh::engine::advance(sz, shear_old, sz.strain, hysteresisParameters(),
                         gmax_, map_new, &tangent);

    // The map is applied as a DIFFERENCE on FLAC3D's current stress so that
    // corrections made between constitutive calls survive.  Normal deviatoric
    // components stay elastic at gmax.
    const Double mean = meanStress(s->stnS_) +
                        bulk_ * (s->stnE_.s11() + s->stnE_.s22() + s->stnE_.s33());
    s->stnS_.rs11() = mean + current_dev[0] + 2.0 * gmax_ * increment_dev[0];
    s->stnS_.rs22() = mean + current_dev[1] + 2.0 * gmax_ * increment_dev[1];
    s->stnS_.rs33() = mean + current_dev[2] + 2.0 * gmax_ * increment_dev[2];
    s->stnS_.rs12() = current_dev[3] + (map_new[0] - sz.map_stress[0]);
    s->stnS_.rs13() = current_dev[4] + (map_new[1] - sz.map_stress[1]);
    s->stnS_.rs23() = current_dev[5] + (map_new[2] - sz.map_stress[2]);
    std::copy(map_new, map_new + 3, sz.map_stress);
    tangent_shear_ = tangent;

    s->viscous_ = true;
    applyMohrCoulomb(s, tangent);
}

void ModelGqhMrdfV6::applyMohrCoulomb(State *s, Double tangent_shear) {
    if (!canFail())
        return;

    SymTensorInfo info;
    DVect3 principal = s->stnS_.getEigenInfo(&info);
    const Double failure = nph_ * principal.z() - principal.x() - csn_;
    if (failure <= 0.0)
        return;

    // Zero-dilation return: ra = e1 - e2 = 2*G_t, rb = -ra, rd = ra*(1 + nph), so
    // ra/rd and rb/rd are +-1/(1 + nph) for any tangent.  Forming e1 - e2 with
    // bulk ~1e7 and the 1e-12 tangent floor (GQ/H saturated) gives 0/0 = NaN.
    (void)tangent_shear;
    principal.rx() += failure / (1.0 + nph_);
    principal.rz() -= failure / (1.0 + nph_);
    s->stnS_ = info.resolve(principal);
    s->state_ |= 1;
    s->viscous_ = false;
}

void ModelGqhMrdfV6::copy(const ConstitutiveModel *mod) {
    ConstitutiveModel::copy(mod);
    const ModelGqhMrdfV6 *other = dynamic_cast<const ModelGqhMrdfV6 *>(mod);
    if (!other)
        throw std::runtime_error("GQ/H-MRDF: constitutive model cast failed.");
    bulk_ = other->bulk_;      shear_ = other->shear_;
    gmax_ = other->gmax_;      tau_max_ = other->tau_max_;
    theta1_ = other->theta1_;  theta2_ = other->theta2_;
    theta3_ = other->theta3_;  theta4_ = other->theta4_;  theta5_ = other->theta5_;
    mrdf_p1_ = other->mrdf_p1_; mrdf_p2_ = other->mrdf_p2_; mrdf_p3_ = other->mrdf_p3_;
    cohesion_ = other->cohesion_; friction_ = other->friction_;
    nph_ = other->nph_; csn_ = other->csn_; tangent_shear_ = other->tangent_shear_;
    sub_ = other->sub_;
    std::copy(other->initial_dev_stress_, other->initial_dev_stress_ + 6,
              initial_dev_stress_);
    history_initialized_ = other->history_initialized_;
    dynamic_reset_done_ = other->dynamic_reset_done_;
}

void ModelGqhMrdfV6::save(std::ostream &o) const {
    ConstitutiveModel::save(o);
    o.write(reinterpret_cast<const char *>(initial_dev_stress_), sizeof(initial_dev_stress_));
    o.write(reinterpret_cast<const char *>(&history_initialized_), sizeof(history_initialized_));
    o.write(reinterpret_cast<const char *>(&dynamic_reset_done_), sizeof(dynamic_reset_done_));
    const std::size_t nsub = sub_.size();
    o.write(reinterpret_cast<const char *>(&nsub), sizeof(nsub));
    for (std::size_t k = 0; k < sub_.size(); ++k) {
        const gqh::engine::Subzone &sz = sub_[k];
        o.write(reinterpret_cast<const char *>(sz.strain), sizeof(sz.strain));
        o.write(reinterpret_cast<const char *>(sz.origin_strain), sizeof(sz.origin_strain));
        o.write(reinterpret_cast<const char *>(sz.origin_stress), sizeof(sz.origin_stress));
        o.write(reinterpret_cast<const char *>(sz.map_stress), sizeof(sz.map_stress));
        o.write(reinterpret_cast<const char *>(sz.peak_strain), sizeof(sz.peak_strain));
        o.write(reinterpret_cast<const char *>(sz.peak_stress), sizeof(sz.peak_stress));
        o.write(reinterpret_cast<const char *>(&sz.peak_distance), sizeof(sz.peak_distance));
        o.write(reinterpret_cast<const char *>(&sz.gamma_max), sizeof(sz.gamma_max));
        o.write(reinterpret_cast<const char *>(&sz.initialized), sizeof(sz.initialized));
        const std::size_t count = sz.stack.size();
        o.write(reinterpret_cast<const char *>(&count), sizeof(count));
        if (count)
            o.write(reinterpret_cast<const char *>(&sz.stack[0]),
                    static_cast<std::streamsize>(count * sizeof(gqh::engine::Reversal)));
    }
}

void ModelGqhMrdfV6::restore(std::istream &i, UInt restoreVersion) {
    ConstitutiveModel::restore(i, restoreVersion);
    i.read(reinterpret_cast<char *>(initial_dev_stress_), sizeof(initial_dev_stress_));
    i.read(reinterpret_cast<char *>(&history_initialized_), sizeof(history_initialized_));
    i.read(reinterpret_cast<char *>(&dynamic_reset_done_), sizeof(dynamic_reset_done_));
    // FLAC3D 6.00.69 and earlier lose plug-in save data (fixed in 6.00.70,
    // verified on 6.00.74): restore() gets an empty stream and every zone
    // would silently restart on the virgin curve.  A sentinel detects that
    // without relying on stream-state inlines across the VS2010/VS2019 runtimes.
    const std::size_t kUnread = static_cast<std::size_t>(-1);
    std::size_t nsub = kUnread;
    i.read(reinterpret_cast<char *>(&nsub), sizeof(nsub));
    if (nsub == kUnread)
        throw std::runtime_error(
            "GQ/H-MRDF restore: hysteresis state missing from the save file. "
            "FLAC3D 6.00.69 and earlier do not keep plug-in state; update to "
            "6.00.74 and re-run the analysis up to the save point.");
    if (nsub > 100)
        throw std::runtime_error("GQ/H-MRDF restore: subzone count is invalid.");
    sub_.assign(nsub, gqh::engine::Subzone());
    for (std::size_t k = 0; k < sub_.size(); ++k) {
        gqh::engine::Subzone &sz = sub_[k];
        i.read(reinterpret_cast<char *>(sz.strain), sizeof(sz.strain));
        i.read(reinterpret_cast<char *>(sz.origin_strain), sizeof(sz.origin_strain));
        i.read(reinterpret_cast<char *>(sz.origin_stress), sizeof(sz.origin_stress));
        i.read(reinterpret_cast<char *>(sz.map_stress), sizeof(sz.map_stress));
        i.read(reinterpret_cast<char *>(sz.peak_strain), sizeof(sz.peak_strain));
        i.read(reinterpret_cast<char *>(sz.peak_stress), sizeof(sz.peak_stress));
        i.read(reinterpret_cast<char *>(&sz.peak_distance), sizeof(sz.peak_distance));
        i.read(reinterpret_cast<char *>(&sz.gamma_max), sizeof(sz.gamma_max));
        i.read(reinterpret_cast<char *>(&sz.initialized), sizeof(sz.initialized));
        std::size_t count = 0;
        i.read(reinterpret_cast<char *>(&count), sizeof(count));
        if (count > 100000)
            throw std::runtime_error("GQ/H-MRDF restore: reversal stack is invalid.");
        sz.stack.resize(count);
        if (count)
            i.read(reinterpret_cast<char *>(&sz.stack[0]),
                   static_cast<std::streamsize>(count * sizeof(gqh::engine::Reversal)));
    }
    updateParameters();
}

Double ModelGqhMrdfV6::getStressStrengthRatio(const SymTensor &stress) const {
    const DVect3 principal = stress.getEigenInfo();
    const Double denominator = principal.z() - principal.x();
    if (denominator <= 0.0)
        return 10.0;
    const Double sigma1_failure = nph_ * principal.z() - csn_;
    return std::min<Double>(10.0, std::max<Double>(0.0,
        (principal.z() - sigma1_failure) / denominator));
}

void ModelGqhMrdfV6::scaleProperties(const Double &scale, const std::vector<UInt> &props) {
    for (std::size_t k = 0; k < props.size(); ++k) {
        switch (props[k]) {
        case 1: bulk_ *= scale; break;
        case 2: shear_ *= scale; break;
        case 3: gmax_ *= scale; break;
        case 4: tau_max_ *= scale; break;
        case 13: cohesion_ *= scale; break;
        case 14:
            friction_ = std::max<Double>(0.0, std::min<Double>(89.0,
                std::atan(std::tan(friction_ * kPi / 180.0) * scale) * 180.0 / kPi));
            break;
        default: break;
        }
    }
    setValid(0);
}

} // namespace models
