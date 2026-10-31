# Verification scope

The portable checks use dimensionless examples with `Gmax = 100`, `tau_max = 10` in consistent units and `gamma_ref = tau_max/Gmax = 0.1`. For zero GQ/H curvature parameters, the backbone reduces to `tau = Gmax*gamma/(1+abs(gamma)/gamma_ref)`. The tests check this special case, the implicit quadratic residual for one nonzero curvature, an analytic hyperbolic-backbone integral for Masing damping, MRDF scaling, simple reversal states, and pure-shear tensor mapping.

Some expected values derived from those equations are:

| Check | Expected value |
| --- | ---: |
| Backbone stress at `gamma = 0.05` | 3.333333333333 |
| Backbone stress at `gamma = -0.2` | -6.666666666667 |
| Masing damping at amplitude `0.1` | 0.144774515882 |
| Pure-shear tensor stress at virgin peak | 6.123724356958 |
| Pure-shear tensor stress after unloading to zero | -2.041241452319 |

The damping test compares the numerical quadrature with the analytic integral `tau_max * (a - gamma_ref * ln(1 + a/gamma_ref))`; it uses an integration tolerance. These expected values are not recorded FLAC3D results.

These synthetic checks do not validate all nonlinear parameters, non-proportional multi-axial paths, calibration accuracy, FLAC3D save/restore, the Mohr-Coulomb adapter, or a complete site response model. They are not benchmarks for a particular project's design.

Before engineering use, record the exact source revision and toolchain; build separately with the intended SDK; test monotonic and cyclic element behavior, reversals, damping and save/restore in the matching FLAC3D version; and compare a controlled model against an independently computed reference with matching units and boundaries. Retain input, output, solver version, logs and DLL hash. Select acceptance tolerances through an explicit validation plan.
