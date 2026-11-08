# Benchmark report: WIMP Earth shielding

**Status:** Implemented and reproducible *within the declared idealized model*. The package is not a validated experimental dark-matter prediction.

## What was completed in the nine-section task

The scientific scenario and initial task specification predate this AI evaluation; numerical implementation is `solution.py` plus `transport_core.py`. The separate original task specification is not included in this compact package. Atmospheric and mixed-crust transport, speed-dependent SI/Helm scattering, flux-normalized xenon isotope-mixture recoil rate, three model comparisons, standalone verification, seed checks, uncertainty accounting, outputs, and graphical diagnostics are included.

### Frozen benchmark settings

- Reference WIMP mass: **100 GeV/c²**; xenon detector scoring plane at **1 km** depth.
- Four hypothetical proton cross sections: **1e-36, 1e-35, 1e-34, 4.44e-34 cm²**.
- Incoming halo: truncated 3-D Maxwellian, fixed laboratory orientation chosen anti-parallel to Galactic lab velocity.
- **5,000 independent histories per mass/cross section** (2 seeds × 2,500 histories); 200,000 independent draws for the incident-flux normalization.
- Detector: illustrative natural xenon isotope mixture, 5-keV step threshold, 100% efficiency above threshold, no detector resolution/background model.
- All published rate units: **events/(kg day keV)** for differential bins, **events/(kg day)** when integrated over recoil energy.

## Numerical comparison (integrated above the 5-keV threshold)

| σp (cm²) | Collisionless (events/kg/day) | Ballistic (events/kg/day) | Full transport (events/kg/day) | Scattered arrivals | Full/collisionless | Full/ballistic |
|---:|---:|---:|---:|---:|---:|---:|
| 1e-36 | 4.505e+05 ± 9.4e+03 | 4.483e+05 ± 8.9e+03 | 4.504e+05 ± 9.4e+03 | 14/5000 | 1.000 | 1.005 |
| 1e-35 | 4.505e+06 ± 9.4e+04 | 4.312e+06 ± 6e+04 | 4.427e+06 ± 5.9e+04 | 93/5000 | 0.983 | 1.027 |
| 1e-34 | 4.505e+07 ± 9.4e+05 | 3.298e+07 ± 2.9e+05 | 4.751e+07 ± 6.8e+05 | 991/5000 | 1.055 | 1.440 |
| 4.44e-34 | 2e+08 ± 4.2e+06 | 6.432e+07 ± 4.3e+05 | 2.342e+08 ± 5.1e+06 | 2976/5000 | 1.171 | 3.641 |

The quoted errors include sampling uncertainty and the independently estimated incident-flux normalisation error, combined in quadrature. All cross sections are illustrative, **not asserted to be consistent with experimental constraints**.

### Physical interpretation

- **σp = 1e-36 cm²:** zero-scatter arrivals = 99.72%; scattered arrivals = 0.28%; full/collisionless integrated rate ratio = 1.000; full/ballistic ratio = 1.005.
- **σp = 1e-35 cm²:** zero-scatter arrivals = 98.06%; scattered arrivals = 1.86%; full/collisionless integrated rate ratio = 0.983; full/ballistic ratio = 1.027.
- **σp = 1e-34 cm²:** zero-scatter arrivals = 79.90%; scattered arrivals = 19.82%; full/collisionless integrated rate ratio = 1.055; full/ballistic ratio = 1.440.
- **σp = 4.44e-34 cm²:** zero-scatter arrivals = 39.08%; scattered arrivals = 59.52%; full/collisionless integrated rate ratio = 1.171; full/ballistic ratio = 3.641.

The ballistic survival probability on initial straight trajectories **is not** the total transmitted WIMP fraction. In this particular infinite-plane reference geometry, some scattered trajectories contribute significantly to detector-plane crossing density; therefore the full rate may exceed the collisionless comparison. This is not a general prediction of realistic underground detectors. Recoil spectra and angle-sensitive residence-time weights are highly model-dependent.

### Uncertainty, convergence and rare-event flags

- Independent file-level verification: **48/48 tests passed**.
- Dedicated physics tests: incident-input validation; zero-momentum Helm normalisation; CM kinematics and conservation; no-scatter survival in a pure silicon slab; conditional first-collision distance CDF; uniform recoil fraction in the F=1 isotropic limit.
- No numerical collision-cap truncations occurred for this frozen grid.
- Two independent random-seed batches were compared with four-sigma uncertainty-aware acceptance rules. All integrated collisionless, ballistic, and transport comparisons passed.
- The benchmark reports an effective sample count per observable and standard errors from independent **initial histories**, not from correlated individual collisions.
- Increasing the budget from one seed (2,500) to both seeds (5,000) generally decreases estimated statistical uncertainty; individual observed estimates are not required to decrease monotonically. The file `summary.json` includes both half-budget and full-budget errors.

- **Low-statistics warning:** at σp=1e-36 cm² only 14 scattered histories reached the detector and the scattered-component effective sample count is **9.8**, below the declared diagnostic reference of 25. Do not interpret its detailed scattered-only spectrum as a resolved signal.

## Completion assessment against the nine sections

| Section | Reference package evidence | Status |
|---|---|---|
| 1. Scientific objective | Collisionless, ballistic, full-transport predictions | Implemented |
| 2. Difficulty | Coupled scattering, layered media, detector response | Documented |
| 3. Input data | Six validated, version-hashed input fixtures | Implemented; illustrative fixtures |
| 4. Numerical method | Velocity-dependent Helm kernel, exponential optical depth, 3-D propagation | Implemented in idealized geometry |
| 5. Interpretation & validation | Report, uncertainty diagnostics and plots | Implemented; limited convergence statistics |
| 6. Deliverables | Executables, CSV, JSON, procedure, tolerances, plots | Included |
| 7. Independent verification | Analytic first-collision, halo-flux integral, rate-grid & bin checks | Implemented checks passed |
| 8. Agentic requirements | Explicit transformations and numerical design decisions | Documented by code + guide |
| 9. Expected outcome | Quantitative comparison, event rates, uncertainty and limitations | Implemented in stated model |

## What this result **does not** establish

1. **A realistic detector response:** infinite scoring plane, perfect energy threshold and no backgrounds, exposure or resolution.
2. **A realistic location-dependent halo/geometry:** laboratory vertical orientation chosen for illustration; no sidereal modulation or spherical geological model.
3. **Complete geo-chemistry:** atmospheric profiles and multicomponent rock are provisional; the `Other` category is an explicitly fictional effective target.
4. **Precision at all energies and all cross sections:** near-grazing 1/cos(angle) weights increase variance, rare scattered contributions may be unresolved; speed tables and energy quadrature require deeper refinement studies for research-grade predictions.
5. **An experimental observation, exclusion or detection:** high benchmark cross sections are hypothetical parameter points.

## How to reproduce

```bash
pip install -r requirements.txt
python run_all.py
```

This command runs the frozen seed experiment, checks physical/schema/normalisation constraints and creates `recoil_spectra.png` and `transport_diagnostics.png`. See `procedure.txt` for the complete algorithm and `tolerance.txt` for acceptance decisions.
