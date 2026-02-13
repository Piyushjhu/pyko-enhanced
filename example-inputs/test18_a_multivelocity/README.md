# pyKO Test 18a: Multi-Velocity Analysis

## Overview

Test 18a performs a **velocity sweep analysis** to study the effect of impact velocity on free surface velocity (FSV) response at constant target thickness. This test simulates a material flyer impacting a target across a range of velocities and analyzes how FSV scales with impact velocity.

## Key Features

- **Variable Impact Velocity**: Sweep through multiple impact velocities (typically 100-1000 m/s)
- **Constant Target Thickness**: Fixed target geometry across all simulations
- **FSV Analysis**: Track free surface velocity response for each velocity
- **Comparative Plotting**: Plot FSV vs time for all velocities on same graph
- **Correlation Analysis**: Analyze peak FSV vs impact velocity relationship
- **Material Database**: Easy material selection system with 8 common materials

## Physics Motivation

Understanding how impact velocity affects shock response is crucial for:
- **Shock strength characterization** as a function of impact velocity
- **Free surface velocity scaling** with impact conditions
- **Spall fracture thresholds** at different impact velocities
- **Impact protection design** across velocity ranges
- **Experimental planning** for shock physics experiments

## Files

### Core Scripts
- `pko-test18_a_multivelocity.py` - Main analysis script
- `material_database.py` - Material database and YAML configuration generator
- `example_material_usage.py` - Example usage demonstrations

### Configuration
- `test18.yml` - YAML configuration file with velocity sweep parameters
- `test18-with-interface-separation.yml` - Generated config with spall physics
- `test18-without-interface-separation.yml` - Generated config without spall physics

### Documentation
- `README.md` - This file
- `README_MATERIAL_SELECTION.md` - Detailed material selection guide

### Output Files
- `fsv_vs_time.png` - FSV vs time for all velocities
- `peak_fsv_vs_impact_velocity.png` - Peak FSV vs impact velocity correlation
- `fsv_vs_time_zoomed.png` - Zoomed FSV vs time plot
- `pyko-test18a-*-velocity-*.dat` - Binary output files for each velocity

## Usage

### Step 1: Generate Configuration (Optional)

If you want to use the material database system:

```bash
python material_database.py --impactor Aluminum --target Copper --min-velocity 100 --max-velocity 500 --velocity-steps 5
```

The generator will create YAML files with:
- Material properties (EOS, strength, fracture)
- Velocity sweep parameters
- Mesh configuration
- Timing parameters

### Step 2: Configure Velocity Sweep

Edit `test18.yml` to set velocity sweep parameters:

```yaml
velocity_sweep:
  min_velocity: 400.0  # Minimum impact velocity (m/s)
  max_velocity: 800.0  # Maximum impact velocity (m/s)
  velocity_steps: 5    # Number of velocity steps to simulate
```

### Step 3: Configure Analysis Modules

Edit `pko-test18_a_multivelocity.py` to enable/disable analysis:

```python
ENABLE_INTERFACE_ANALYSIS = True   # Interface analysis (usually enabled)
ENABLE_FSV_ANALYSIS = True         # Free surface velocity analysis (required)
ENABLE_STRESS_ANALYSIS = False     # Stress analysis (optional)
ENABLE_SPALL_ANALYSIS = False      # Spall analysis (optional for velocity sweep)
```

### Step 4: Run Analysis

```bash
python pko-test18_a_multivelocity.py
```

The script will:
1. Load the YAML configuration
2. Extract velocity sweep parameters
3. Run pyKO simulations for each impact velocity
4. Calculate FSV for each simulation
5. Generate comparative plots
6. Provide summary statistics

## Configuration Parameters

### Velocity Sweep
```yaml
velocity_sweep:
  min_velocity: 400.0  # m/s - minimum impact velocity
  max_velocity: 800.0  # m/s - maximum impact velocity
  velocity_steps: 5    # Number of velocity values to simulate
```

### Timing
- **tstop**: Total simulation time (typically 0.5 μs)
- **dtstart**: Initial time step (0.0001 μs)
- **dtoutput**: Output frequency (0.001 μs)

### Material Configuration
Each material requires:
- **EOS parameters**: Mie-Gruneisen equation of state
- **Strength model**: HYDRO, VM (Von Mises), or SG (Steinberg-Guinan)
- **Fracture parameters**: Spall strength (pfrac) and density ratio (nrhomin)

## Expected Results

### FSV vs Time Plot
- Multiple curves showing FSV evolution for different velocities
- Higher velocities show larger peak FSV
- Shock arrival time is similar (same target thickness)
- FSV rise time may vary with velocity

### Peak FSV vs Impact Velocity Plot
- Strong correlation between impact velocity and peak FSV
- Typically linear or power-law relationship
- Statistical correlation coefficient displayed
- Trend line showing scaling behavior

## Physics Insights

### Shock Strength Scaling
- **Higher velocities**: Stronger shocks, larger FSV
- **Lower velocities**: Weaker shocks, smaller FSV
- **Threshold effects**: Minimum velocity for spall or interface separation

### FSV Amplification
- **Velocity dependence**: Peak FSV scales with impact velocity
- **Material effects**: Different materials show different scaling
- **Strength model effects**: Elastic-plastic response affects FSV

### Spall Behavior
- **Velocity threshold**: Minimum velocity for spall initiation
- **Spall strength**: Material-dependent fracture threshold
- **Multiple spall**: Possible at high velocities

## Available Materials

The material database includes 8 common materials:

| Material | Density (kg/m³) | Sound Speed (m/s) | Spall Strength (MPa) |
|----------|----------------|-------------------|---------------------|
| **Aluminum** | 2700 | 5200 | 276 |
| **Copper** | 8930 | 3900 | 1500 |
| **Steel** | 7900 | 4570 | 2000 |
| **Iron** | 7870 | 4600 | 1800 |
| **Titanium** | 4500 | 4780 | 1200 |
| **Tungsten** | 19300 | 4020 | 3000 |
| **Gold** | 19320 | 3240 | 800 |
| **Lead** | 11340 | 2050 | 300 |

## Comparison with Test 18b

| Aspect | Test 18a (Multi-Velocity) | Test 18b (Multi-Thickness) |
|--------|---------------------------|----------------------------|
| **Variable Parameter** | Impact velocity | Target thickness |
| **Constant Parameter** | Target thickness | Impact velocity |
| **Analysis Focus** | Velocity effects on FSV | Thickness effects on FSV |
| **Physics Insight** | Shock strength dependence | Wave propagation effects |
| **Applications** | Impact characterization | Target design optimization |

## Troubleshooting

### Common Issues
1. **Simulation time too short**: Increase `tstop` in YAML
2. **Mesh resolution**: Ensure sufficient cells for accuracy (1-2 μm/cell recommended)
3. **Strength model stability**: VM is most stable, SG may have issues at high velocities
4. **File paths**: Ensure scripts run from correct directory
5. **Velocity range**: Very high velocities may require smaller time steps

### Performance Tips
- Start with fewer velocity steps for testing (3-5 steps)
- Use VM strength model for stability
- Monitor simulation completion for each velocity
- Check output files are generated before plotting
- Use appropriate mesh resolution (balance accuracy vs speed)

## Example Results

Typical output shows:
- **Velocity range**: 100-1000 m/s
- **Peak FSV range**: Scales with velocity (typically 0.5-2x impact velocity)
- **Correlation**: Strong positive correlation (R² > 0.9 typical)
- **Timing**: Shock arrival time constant (same target thickness)

## Next Steps

After running Test 18a, consider:
1. **Parameter sensitivity**: Vary target thickness or materials
2. **Strength model comparison**: Test different models (VM vs SG)
3. **Mesh convergence**: Study effect of cell count
4. **Combined analysis**: Use both 18a and 18b results together
5. **Material combinations**: Test different impactor-target pairs

## References

- **pyKO Manual**: [https://impactswiki.github.io/pyko/](https://impactswiki.github.io/pyko/)
- **Original KO Code**: [https://www.eng.mu.edu/shockphysics/KO/](https://www.eng.mu.edu/shockphysics/KO/)
- **Shock Physics**: Standard Hugoniot relations and shock wave theory
