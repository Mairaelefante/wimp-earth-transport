# Scientific foundation: WIMP transport and nuclear recoil simulation

This folder contains selected reproducible **original computational-physics work** from an earlier WIMP transport study. It predates the AI evaluation and exists independently of it.

The code explores a truncated Galactic halo velocity distribution transformed to the laboratory frame, propagation through stratified atmospheric and crustal matter, stochastic SI nuclear scattering including a Helm form factor, and a simplified xenon recoil scoring model. It compares collisionless, zero-collision (ballistic) and full-transport model outputs.

Start with `BENCHMARK_REPORT.md` for assumptions, parameter settings, interpretation and limitations. The original scientific plots are `transport_diagnostics.png` and `recoil_spectra.png`; they are model-generated simulated figures, not observed event data. Inputs are in `inputs/`. Main scripts: `transport_core.py`, `solution.py`, `run_all.py`, `make_plots.py`, `verification.py`, and `test_physics.py`.

Reproduce using an isolated Python environment with `python -m pip install -r requirements.txt`, then `python run_all.py`. This recomputes outputs and may be computationally expensive; run on a copy if you need to preserve the frozen published result files. Model assumptions and tested limits are documented in the included scientific report. The selected inputs are illustrative.

This scientific foundation is independent from the four synthetic evaluation prompts. No trial in the submitted AI benchmark measures full autonomous reproduction of this original simulation.
