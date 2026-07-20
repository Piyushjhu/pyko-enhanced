# pyKO

**Program Name:** pyKO  
**Author:** Sarah T. Stewart  
**Version:** 0.6.0  
**License:** GNU General Public License v3.0

pyKO is a 1-D elastic-plastic, 2nd order, Lagrangian hydrocode.  
Based on Mark Wilkins, *Computer Simulation of Dynamic Phenomena*, Springer-Verlag, 1999.  
Adapted from John Borg's Fortran KO code v11: https://www.eng.mu.edu/shockphysics/KO/

## Citation

[![DOI](https://zenodo.org/badge/602649996.svg)](https://zenodo.org/badge/latestdoi/602649996)

Stewart, S. T. pyKO code v0.6.1, doi:10.5281/zenodo.8092348, 2023.

## User manual

https://impactswiki.github.io/pyko/

KO code is available in Fortran, C, Python and Matlab. See the manual for links.

## Report bugs

https://github.com/ImpactsWiki/pyko/issues

## User mailing list

There is a user mailing list for KO code users (in any programming language).  
To sign up, send an email to **sympa@ucdavis.edu** from the address you want to subscribe:

- **Subject:** subscribe ko-code-users  
- **Body:** leave empty  

You will receive an email with a link to confirm your subscription.

## Authors and contributors

- **pyKO (original):** Sarah T. Stewart  
- **Interface separation:** Jake Diamond (contributor)  
- **Test 17, Test 18, Test 19:** Piyush Wanchoo (contributor)

---

## WARNING (interface separation)

This version includes an experimental interface separation criterion. The method has not been rigorously validated. Interface separation was added to help plan experiments; interpret results as approximations only.

---

## Example inputs: how to use the test cases

Example input files and scripts live in **`example-inputs/`**. Tabular EOS models are in **`eos/`**. Below is how to run each test.

**Step-by-step example runs and expected results:** see **[example-inputs/TUTORIAL.md](example-inputs/TUTORIAL.md)** for a short tutorial with commands, terminal output, and result files for Test 17, 18a, and 18b.

### Directory layout

```
example-inputs/
├── pko-test17-spall-interface.py      # Test 17: single impact, spall + interface
├── pko-test17-xt-animation.py         # Test 17: x-t diagram animation
├── pko-test18_a_multivelocity.py      # Test 18a: velocity sweep
├── pko-test18_b_multithickness.py     # Test 18b: thickness sweep
├── test17-spall-interface/            # Config and material files for Test 17
├── test18_a_multivelocity/           # Config and material DB for Test 18a
└── test18_b_multithickness/           # Config for Test 18b
```

---

### Test 17: Hybrid spall and interface analysis

Single impact with optional spall, interface tracking, FSV, and stress analysis.

**Run (from repo root):**

```bash
cd example-inputs
python pko-test17-spall-interface.py
```

**Configure:** Edit the script header to turn analyses on/off:

```python
ENABLE_SPALL_ANALYSIS = True      # Density + pressure-based spall
ENABLE_INTERFACE_ANALYSIS = True  # Material interface tracking
ENABLE_FSV_ANALYSIS = True        # Free surface velocity
ENABLE_STRESS_ANALYSIS = True     # Max stress in target
```

**Config files (in `test17-spall-interface/`):**

- `test17-with-interface-separation.yml` — with interface separation
- `test17-without-interface-separation.yml` — spall only
- `test17-material-config.yml` — material definitions

The script picks the YAML based on `ENABLE_INTERFACE_ANALYSIS`.  
**Outputs:** Pressure x–t plots (Eulerian/Lagrangian), FSV vs time, stress evolution, interface position, density ratio. Binary output and summary in the run directory.

---

### Test 18a: Multi-velocity analysis

Sweeps impact velocity at fixed target thickness and analyzes free surface velocity (FSV).

**Run (from repo root):**

```bash
cd example-inputs
python pko-test18_a_multivelocity.py
```

The script reads **`example-inputs/test18_a_multivelocity/test18.yml`**. Set the velocity range and number of steps there:

```yaml
velocity_sweep:
  min_velocity: 400.0   # m/s
  max_velocity: 800.0  # m/s
  velocity_steps: 5
```

**Optional – generate a new config with the material database:**

```bash
cd example-inputs/test18_a_multivelocity
python material_database.py --list-materials
python material_database.py --impactor Aluminum --target Copper \
  --min-velocity 100 --max-velocity 500 --velocity-steps 5
```

Then run `pko-test18_a_multivelocity.py` as above (from `example-inputs/`).  
**Outputs:** FSV vs time for each velocity, peak FSV vs impact velocity, zoomed FSV plot, and binary `.dat` files per velocity.

---

### Test 18b: Multi-thickness analysis

Sweeps target thickness at fixed impact velocity and analyzes FSV.

**Run (from repo root):**

```bash
cd example-inputs
python pko-test18_b_multithickness.py
```

The script reads **`example-inputs/test18_b_multithickness/test18b.yml`**. Set the thickness sweep there:

```yaml
thickness_sweep:
  constant_velocity: 300.0  # m/s
  min_thickness: 50.0      # μm
  max_thickness: 500.0     # μm
  thickness_steps: 5
```

**Optional – generate config with the interactive generator:**

```bash
cd example-inputs/test18_b_multithickness
python yml_generator_18b.py
```

Then run `pko-test18_b_multithickness.py` from `example-inputs/`.  
**Outputs:** FSV vs time for each thickness, peak FSV vs thickness, combined plots, and binary `.dat` files per thickness.

---

### Summary: which test to use

| Test   | What it does                          | Main output                          |
|--------|----------------------------------------|--------------------------------------|
| **17** | Single impact; spall, interface, FSV, stress | x–t diagrams, FSV, stress, interface |
| **18a**| Many velocities, one thickness         | FSV vs time; peak FSV vs velocity     |
| **18b**| One velocity, many thicknesses         | FSV vs time; peak FSV vs thickness   |

---

### Configuration tips

- **Mesh:** For shock/spall, ~1–5 μm per cell is typical; finer for spall, coarser for speed.
- **Spall (Test 17):** In the YAML, set `frac.pfrac` (e.g. Al ~2–5×10⁸ Pa, Cu ~8–12×10⁸ Pa) and `frac.nrhomin` (e.g. 0.9).
- **Paths:** Run the main scripts from **`example-inputs/`** so paths to config and output directories are correct.
- **Impact velocity:** A top-level `impact_velocity` key in the YAML sets the mat1 flyer speed directly (preferred over `mat1.init.up0`). `material_database.py` generates configs using this key. Test 18b reads its sweep velocity from `thickness_sweep.constant_velocity` instead.

---

## Beta release

v0.6.x is the first public release. The code needs user testing and development.
