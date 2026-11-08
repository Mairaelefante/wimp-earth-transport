![WIMP Earth Transport](cover.png)

# WIMP Earth Transport

**Monte Carlo simulation of dark-matter scattering through Earth's atmosphere and crust, with nuclear recoil predictions for an idealized underground xenon detector.**

Computational Astroparticle Physics | Monte Carlo Methods | Python | Dark Matter Direct Detection

## Overview

Weakly Interacting Massive Particles (WIMPs) are hypothetical dark-matter candidates that may interact with ordinary matter through rare nuclear scattering events.

When studying underground dark-matter detection, it is important to consider not only the incoming galactic dark-matter velocity distribution, but also how particles may interact with the material they traverse before reaching a detector.

This project investigates that problem through a numerical Monte Carlo framework, modeling WIMP propagation across atmospheric and crustal layers with different material compositions.

The central scientific question is:

**How do stochastic nuclear interactions in Earth's atmosphere and crust modify the WIMP population reaching an underground detector and its predicted nuclear recoil spectrum?**

## Scientific Method

The simulation follows six main stages.

**1. Galactic velocity distribution**

Incident WIMP velocities are sampled from a truncated Maxwellian halo model and transformed into the laboratory reference frame. The resulting velocity vectors determine the incoming particle energies and propagation directions.

**2. Atmospheric and crustal composition**

The propagation medium is modeled using atmospheric layers and a mixed crustal composition. Material densities and nuclear species determine the probability of scattering along a particle trajectory.

**3. Nuclear scattering**

Elastic WIMP–nucleus scattering is modeled using the interaction parameters, nuclear kinematics and a Helm form factor. Each scattering event may change the WIMP's kinetic energy and direction.

**4. Monte Carlo propagation**

Particle histories are followed stochastically through the propagation medium. Collision distances and scattering outcomes are sampled, and particle states are updated until they reach the scoring plane or satisfy a termination condition.

**5. Underground detector response**

Particles reaching the idealized detector scoring plane contribute to the calculation of nuclear recoil spectra for xenon.

**6. Numerical verification**

The calculation includes checks of physical consistency, numerical behavior, sampling uncertainty and the reproducibility of the simulation.

## Three Propagation Scenarios

The project compares three distinct physical descriptions:

| Scenario | Physical meaning |
|---|---|
| **Collisionless** | Reference model in which WIMPs propagate without interactions in the medium |
| **Ballistic** | Contribution from WIMPs that experience zero collisions in a medium where scattering is possible |
| **Full transport** | Combined contribution from unscattered and scattered trajectories reaching the scoring plane |

In the full-transport calculation, scattered trajectories may still contribute to the detector response after changing their energy or direction.

This distinction is important when studying how the propagation medium affects the predicted recoil spectrum.

## Scientific Results

The project generates recoil spectra and transport diagnostics across a set of illustrative WIMP–nucleon cross sections.

### Nuclear recoil spectra

![Nuclear recoil spectra](recoil_spectra.png)

The recoil spectra compare predicted detector responses under the different propagation scenarios.

### Transport diagnostics

![Transport diagnostics](transport_diagnostics.png)

The transport diagnostics illustrate how scattering modifies the trajectories contributing to the underground scoring plane.

The simulation is based on an idealized detector geometry and stated physical assumptions. The results are computational predictions within that model, not experimental measurements or dark-matter detection claims.

## Repository Structure

The main components are:

| File or folder | Purpose |
|---|---|
| `solution.py` | Main computational workflow |
| `transport_core.py` | Particle propagation, scattering and transport calculations |
| `run_all.py` | Runs the simulation and associated checks |
| `inputs/` | Halo, interaction, detector and material parameters |
| `verification.py` | Numerical and physical consistency checks |
| `test_physics.py` | Physics-specific tests |
| `results.csv` | Numerical recoil-spectrum results |
| `summary.json` | Aggregated simulation results and statistics |
| `recoil_spectra.png` | Recoil-spectrum visualization |
| `transport_diagnostics.png` | Transport-analysis visualization |

## Running the Simulation

The project uses Python with scientific computing libraries including NumPy, SciPy and Matplotlib.

Install the dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Run the computational workflow:

```bash
python3 run_all.py
```

The implementation also includes verification routines and physics tests.

The numerical configuration and model inputs are stored separately to make the assumptions inspectable and the calculations reproducible.

## Model Assumptions and Limitations

This is a computational study based on an idealized physical model.

Important limitations include the simplified terrestrial geometry, chosen atmospheric and crustal compositions, idealized detector response, finite Monte Carlo statistics and the use of illustrative interaction parameters.

Results should therefore be interpreted within the stated model rather than extrapolated directly to experimental detection rates or exclusion limits.

## Further Development

Possible extensions include more realistic terrestrial geometry and material profiles, refined detector response, larger Monte Carlo samples, improved treatment of rare scattering histories, and broader numerical convergence studies.

## Project Focus

This project brings together astroparticle physics, stochastic processes, nuclear scattering theory, numerical methods and scientific software development.

Its main contribution is a computational framework for exploring how interactions with terrestrial matter can affect dark-matter propagation and predicted underground recoil signals.
