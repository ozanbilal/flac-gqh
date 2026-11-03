# FLAC3D GQ/H + MRDF

A user-defined constitutive model for FLAC3D, combining a generalized quadratic/hyperbolic (GQ/H) backbone and an MRDF-UIUC unloading/reloading rule. The portable constitutive core is shared by separate FLAC3D 9 and FLAC3D 6 adapters.

**Status:** preparation and independent validation are in progress. The portable tests do not establish FLAC3D compatibility, engineering suitability, or calibration for a particular project. This project's source is provided under the [MIT License](LICENSE).

## What is here

| Path | Purpose |
| --- | --- |
| `src/gqh_math.h`, `src/gqh_engine.h` | Shared backbone, hysteresis and tensor mapping; no FLAC3D dependency |
| `src/modelgqh_mrdf.*`, `src/refcurve.*` | FLAC3D 9 adapter and reference-curve calibration |
| `src_v6/` | FLAC3D 6 adapter; direct theta/MRDF parameters, without the v9 automatic reference-curve calibration |
| `tests/` | Portable numerical checks against simple closed-form cases |
| `build-tools/` | Windows build entry points that require a locally licensed, version-matched Itasca SDK |

The solver, SDK headers and libraries, compiled DLLs, project inputs, earthquake records and client results are not included. FLAC3D 6 and 9 use different plugin ABIs: do not interchange their DLLs.

## Run the portable checks

Requirements: CMake 3.16 or newer and a C++17 compiler. These commands do not require FLAC3D:

```sh
cmake -S . -B out -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build out --config Release
ctest --test-dir out -C Release --output-on-failure
```

The tests cover a closed-form hyperbolic backbone and damping integral, scalar and tensor reversal paths, continuation after scalar state save/restore, and selected reference-curve input checks. They are regression checks, not independent acceptance of the complete constitutive model. See [verification scope](docs/verification.md).

## Build a FLAC3D plugin on Windows

Both scripts require the corresponding licensed FLAC3D installation, its plugin SDK, MSVC and the Windows SDK. They only build a DLL; no solver run is started. Running a model with a constitutive plugin also requires the applicable FLAC3D C++ plugin option and configuration; see [Itasca's FLAC3D 9 plugin documentation](https://docs.itascacg.com/itasca900/flac3d/docproject/source/options/cppplugin/implementation.html). Review the SDK version and compiler requirements against the installed Itasca release before use. From the repository root in PowerShell:

```powershell
# FLAC3D 9: optionally override the default PluginFiles location.
$env:FLAC3D9_PLUGINFILES = 'C:\path\to\FLAC3D9\PluginFiles'
.\build-tools\build-v9.ps1

# FLAC3D 6: optionally override its pluginfiles location.
$env:FLAC3D6_PLUGINFILES = 'C:\path\to\FLAC3D6\pluginfiles'
.\build-tools\build-v6.ps1
```

The examples above contain placeholders and must be changed to the locally installed SDK paths. Omit each environment variable to use the script's default installation path. The v9 script uses C++17 and searches for MSVC 2019/2022; the v6 script uses C++14 and selects MSVC 2019. Output DLLs are written to `build-tools/`. The supplied scripts and adapters have not yet been compiled or run in the current verification environment. Check compiler output, binary hash, solver version, and a controlled element test before using a build.

## Engineering scope

Select parameters from reviewed material data or an explicitly reviewed calibration. Keep stress/strain units, initial effective stress, input-motion definition, and boundary conventions consistent. The cyclic rule already produces hysteretic response; avoid adding a second hysteretic damping model without an explicit modeling basis. Separate small-strain viscous damping is a project decision. Run project-specific verification before engineering use.

The GQ/H backbone is described by Groholski et al. (2016), *Simplified Model for Small-Strain Nonlinearity and Strength in 1D Seismic Site Response Analysis*, DOI [10.1061/(ASCE)GT.1943-5606.0001496](https://doi.org/10.1061/(ASCE)GT.1943-5606.0001496). See [scientific references](docs/references.md) for the named reference-curve models. These citations provide background and do not assert that the publications' authors validated this implementation.

## Known limitation: Mohr-Coulomb cap with free vertical motion

The optional Mohr-Coulomb cap (`cohesion`, `friction`) is not yet verified for 3D models in which neighbouring zones can move vertically relative to each other, for example models with piles, walls, foundations or other structural elements.

Both adapters use the same return rule. They apply GQ/H to the shear components only; the normal deviatoric components stay elastic at `gmax`. When the cap is active, the principal-space return reduces both the shear stress and the normal deviator. In a FLAC3D free-field test (periodic 5 × 5 m plan, 1 m zones, gravity with K0 = 0.667, static Tresca ratio about 0.9 near the base), this produced a zone-scale alternation of peak shear strain with depth and a higher surface acceleration than the equivalent 1 × 1 m column. The 3D and 1D results agreed when vertical velocity was fixed during shaking (spectral band ratio 0.999). FLAC3D's built-in `mohr-coulomb` model under the same setup gave matching 3D and 1D results. The same behaviour appeared at a lower static ratio (about 0.67) once a stronger input made the cap engage, so it is not limited to near-failure initial states. Alternative return rules tried so far (a closest-point return in the energy norm of the predictor, and returns that also remove the plastic shear strain from the hysteresis state) matched 3D and 1D in some records but not in others; none is adopted.

Until this is resolved:

- With the cap off (`cohesion` set very large), this behaviour was not observed.
- In free-field site-response models, fixing vertical velocity during the dynamic stage avoids it.
- In models that need vertical freedom, check the result against a model with the cap off or with a built-in plasticity model before using it.

## Contributions and security

See [CONTRIBUTING.md](CONTRIBUTING.md), [SECURITY.md](SECURITY.md), and [AGENTS.md](AGENTS.md). The workflows include an isolated Linux core-test design, build-free CodeQL, secret scanning, dependency review and GitHub Actions update proposals. The core-test job is configured for pushes to `main` and manual runs; pull requests receive static scans but do not run the core-test job. See [CI isolation](docs/ci-isolation.md) for the fail-closed preflight, limits, and current verification gaps. These workflows have not yet been run on GitHub. CodeQL's adapter coverage may be incomplete because the proprietary SDK is not present.
