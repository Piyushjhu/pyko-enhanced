# pyKO Example Runs and Results — Tutorial

This tutorial walks you through example runs for each test case and shows what to expect: commands, terminal output, and result files. Run commands from the **repository root** unless noted.

**Prerequisites:** Python 3.8+, pyKO and its dependencies installed, and the repo cloned. From the repo root, ensure you can `import pyko` (e.g. run from the directory that contains `pyko.py` or have it on `PYTHONPATH`).

---

## Overview of the three tests

| Test | What it does | Typical run time (approx.) |
|------|----------------|----------------------------|
| **Test 17** | Single impact: Al flyer → Cu target; spall, interface, FSV, stress | 1–3 min |
| **Test 18a** | Velocity sweep: several impact velocities, one target thickness | 2–10 min (depends on steps) |
| **Test 18b** | Thickness sweep: one velocity, several target thicknesses | 2–10 min (depends on steps) |

---

## Example 1: Test 17 — Single impact with spall and interface analysis

### What you will run

One planar impact with optional spall detection, interface tracking, free-surface velocity (FSV), and stress analysis. Good first run to see all analysis modules.

### Step 1: Go to example-inputs

```bash
cd example-inputs
```

### Step 2: Run Test 17

```bash
python pko-test17-spall-interface.py
```

The script uses config from **`test17-spall-interface/`**. It picks a YAML automatically: with interface separation if `ENABLE_INTERFACE_ANALYSIS = True`, or without if `False` (see script header).

### What you see in the terminal

- **Analysis configuration:** Which modules are ON/OFF (Spall, FSV, Stress, Interface).
- **Input parameters summary:** Timing (tstop, dtstart, dtoutput), material names, thicknesses, densities, EOS (c0, s1, gamma0), strength (model, shear modulus, yield), fracture (pfrac, nrhomin).
- **Physics:** Message indicating whether interface separation and spall are enabled or disabled.
- **Run progress:** PyKO runs; then the script loads the binary output and runs the requested analyses.
- **Spall analysis (if ON):** Density- and pressure-based spall checks; whether spall was detected; max tensile pressure.
- **FSV analysis (if ON):** Max free-surface velocity and time; optional FSV-based spall strength (½ρ₀c₀Δu) and comparison with YAML threshold.
- **Plots:** Matplotlib windows for Eulerian and Lagrangian x–t diagrams (pressure, particle velocity, density ratio, material ID), and any FSV/stress/interface plots. Close each window to continue.

### Result files

- **Binary output:** Path comes from the YAML `outputfilename` (e.g. `./pyko-test17-aluminum-copper-single-velocity-bin.dat`). Written in `example-inputs/` when the script runs from there.
- **Figures:** Test 17 uses `plt.show()`, so figures are displayed interactively. They are not saved to disk unless you add `plt.savefig(...)` in the script.

### Example summary lines (conceptual)

You might see lines similar to:

```
Maximum free surface velocity:  xxx m/s at x.xx μs
Peak FSV: xxx m/s at x.xxx μs
FSV-measured spall strength: x.xxx GPa
YAML Cu spall threshold:     x.xxx GPa
```

Exact values depend on impact velocity, materials, and mesh. The script also prints whether spall was detected (density or pressure criterion).

### If something goes wrong

- **“No module named 'pyko'”:** Run from the repo root or set `PYTHONPATH` so that the directory containing `pyko.py` is on the path. The script itself may do `sys.path.insert(1, ...)` to the repo; adjust if your layout differs.
- **YAML / config errors:** Check that `test17-spall-interface/test17-with-interface-separation.yml` (or `-without-`) and `test17-material-config.yml` exist and are valid.
- **Simulation fails:** Reduce impact velocity or increase `tstop` in the YAML.

---

## Example 2: Test 18a — Multi-velocity sweep

### What you will run

Several impacts at different velocities (e.g. 400–800 m/s in 5 steps), same target. The script computes FSV for each run and plots FSV vs time and peak FSV vs impact velocity.

### Step 1: Check or edit the velocity sweep

Open **`example-inputs/test18_a_multivelocity/test18.yml`**. For a short tutorial run you can use:

```yaml
velocity_sweep:
  min_velocity: 400.0   # m/s
  max_velocity: 800.0   # m/s
  velocity_steps: 5
```

So you get 5 velocities (400, 500, 600, 700, 800 m/s). Save and close.

### Step 2: Run Test 18a

From **repo root**:

```bash
cd example-inputs
python pko-test18_a_multivelocity.py
```

The script reads `test18.yml` from `test18_a_multivelocity/` and runs pyKO once per velocity.

### What you see in the terminal

- **Velocity sweep parameters:** Min/max velocity, number of steps, list of velocities.
- For each velocity:
  - Temporary YAML path, validation, “Starting PyKO simulation”, “PyKO simulation completed successfully”.
  - “Peak FSV: xxx m/s at x.xx μs”.
- At the end:
  - “Velocity sweep completed! Successfully simulated N velocities.”
  - “Three plots generated: fsv_vs_time.png, fsv_vs_time_zoomed.png, peak_fsv_vs_impact_velocity.png.”
  - **VELOCITY SWEEP SUMMARY:** Velocity range, peak FSV range, FSV amplification factor.
  - **MULTIVELOCITY ANALYSIS TEST 18a SUMMARY:** N/N simulations successful, velocity range.

### Result files (in `example-inputs/test18_a_multivelocity/`)

| File | Description |
|------|-------------|
| `fsv_vs_time.png` | FSV vs time for all velocities (overlaid). |
| `fsv_vs_time_zoomed.png` | Same, zoomed in time (e.g. 0.2–0.3 μs). |
| `peak_fsv_vs_impact_velocity.png` | Peak FSV vs impact velocity (correlation). |
| `pyko-test18-velocity-400-bin.dat`, `...-500-...`, ... | Binary output for each velocity. |

### Example summary (conceptual)

```
Velocity range: 400 - 800 m/s
Peak FSV range: ~xxx - xxx m/s
FSV amplification factor: ~x.xx
Velocity sweep completed: 5/5 simulations successful
```

Peak FSV should increase with impact velocity; the correlation plot should show a clear trend (often near-linear for this regime).

### Optional: New config from material database

To try a different material pair and velocity range:

```bash
cd example-inputs/test18_a_multivelocity
python material_database.py --list-materials
python material_database.py --impactor Aluminum --target Copper \
  --min-velocity 100 --max-velocity 500 --velocity-steps 5
```

Then run `pko-test18_a_multivelocity.py` again from `example-inputs/`; it will use the updated config in that directory.

---

## Example 3: Test 18b — Multi-thickness sweep

### What you will run

One impact velocity, several target thicknesses (e.g. 50–500 μm in 5 steps). The script computes FSV for each thickness and plots FSV vs time and peak FSV vs thickness.

### Step 1: Check or create config

Config file: **`example-inputs/test18_b_multithickness/test18b.yml`**. It must define `thickness_sweep` and material/mesh/timing. Example snippet:

```yaml
thickness_sweep:
  constant_velocity: 300.0  # m/s
  min_thickness: 50.0       # μm
  max_thickness: 500.0     # μm
  thickness_steps: 5
```

If the file does not exist, generate it:

```bash
cd example-inputs/test18_b_multithickness
python yml_generator_18b.py
```

Answer the prompts (impactor/target materials, velocity, thickness range, steps, mesh), then run the main script from `example-inputs/`.

### Step 2: Run Test 18b

From **repo root**:

```bash
cd example-inputs
python pko-test18_b_multithickness.py
```

### What you see in the terminal

- **Thickness sweep parameters:** Constant impact velocity, thickness range, number of steps, list of thicknesses (μm).
- For each thickness:
  - “Simulating thickness i/N: xxx μm”, timing (tstop, dtoutput), “PyKO simulation completed successfully”, “Peak FSV: xxx m/s at x.xx μs”.
- At the end:
  - “Thickness sweep completed! Successfully simulated N thicknesses.”
  - “Saved combined plot: thickness_sweep_fsv.png”, “Saved FSV vs Time plot: fsv_vs_time.png”, “Saved Peak FSV vs Thickness plot: peak_fsv_vs_thickness.png.”
  - **THICKNESS SWEEP SUMMARY:** Thickness range, peak FSV range, optional correlation coefficient.
  - **MULTITHICKNESS ANALYSIS TEST 18b SUMMARY:** N/N simulations successful, constant velocity, thickness range.

### Result files (in `example-inputs/test18_b_multithickness/`)

| File | Description |
|------|-------------|
| `fsv_vs_time.png` | FSV vs time for all thicknesses (overlaid). |
| `thickness_sweep_fsv.png` | Combined thickness-sweep visualization. |
| `peak_fsv_vs_thickness.png` | Peak FSV vs target thickness (correlation). |
| `pyko-test18b-thickness-50.0-bin.dat`, ... | Binary output for each thickness. |

### Example summary (conceptual)

```
Thickness range: 50 - 500 μm
Peak FSV range: ~xxx - xxx m/s
Thickness sweep completed: 5/5 simulations successful
```

Shock arrival at the free surface is later for thicker targets; peak FSV can vary with thickness and wave interactions.

---

## Quick reference: where things are

- **Test 17:** Run `python pko-test17-spall-interface.py` from `example-inputs/`. Config: `test17-spall-interface/*.yml`, `test17-material-config.yml`. Binary output: path from YAML (e.g. `./pyko-test17-...-bin.dat` in `example-inputs/`).
- **Test 18a:** Run `python pko-test18_a_multivelocity.py` from `example-inputs/`. Config: `test18_a_multivelocity/test18.yml`. Plots and `.dat` files: `test18_a_multivelocity/`.
- **Test 18b:** Run `python pko-test18_b_multithickness.py` from `example-inputs/`. Config: `test18_b_multithickness/test18b.yml`. Plots and `.dat` files: `test18_b_multithickness/`.

---

## Suggested order for a first-time tutorial

1. **Test 17** — One run, see full parameter summary and (if enabled) spall, FSV, stress, interface. Close plot windows to finish.
2. **Test 18a** — Small velocity sweep (e.g. 5 steps). Check the three PNGs in `test18_a_multivelocity/` and the printed peak FSV summary.
3. **Test 18b** — Small thickness sweep (e.g. 5 steps). Check the three PNGs in `test18_b_multithickness/` and the printed summary.

After that you can change materials (Test 17/18a material configs, or 18b generator), velocity ranges, thickness ranges, or mesh resolution and compare results.
