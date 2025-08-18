# pyKO Test 18b: Multi-Thickness Analysis

## Overview

Test 18b performs a **thickness sweep analysis** to study the effect of target thickness on free surface velocity (FSV) response at constant impact velocity. This complements Test 18a (multi-velocity analysis) by keeping impact velocity constant while varying target thickness.

## Key Features

- **Constant Impact Velocity**: Fixed impact velocity across all simulations
- **Variable Target Thickness**: Sweep through multiple target thicknesses
- **FSV Analysis**: Track free surface velocity response for each thickness
- **Comparative Plotting**: Plot FSV vs time for all thicknesses on same graph
- **Correlation Analysis**: Analyze peak FSV vs target thickness relationship

## Physics Motivation

Understanding how target thickness affects shock response is crucial for:
- **Shock wave propagation** through different material thicknesses
- **Free surface velocity amplification** as a function of thickness
- **Spall fracture behavior** in thin vs thick targets
- **Design optimization** for impact protection systems

## Files

### Core Scripts
- `pko-test18_b_multithickness.py` - Main analysis script
- `yml_generator_18b.py` - Interactive YAML configuration generator

### Configuration
- `test18b.yml` - Generated YAML configuration file (created by generator)

### Output Files
- `fsv_vs_time.png` - FSV vs time for all thicknesses
- `peak_fsv_vs_thickness.png` - Peak FSV vs target thickness correlation
- `thickness_sweep_fsv.png` - Combined plots
- `pyko-test18b-thickness-*.dat` - Binary output files for each thickness

## Usage

### Step 1: Generate Configuration
```bash
python yml_generator_18b.py
```

The generator will prompt you for:
- **Impactor material** (e.g., Aluminum, Steel)
- **Target material** (e.g., Copper, Iron)
- **Strength models** for both materials (HYDRO, VM, SG)
- **Thickness parameters**:
  - Impactor thickness
  - Target thickness range (min, max, steps)
- **Impact velocity** (constant for all simulations)
- **Mesh resolution** (cells per material)

### Step 2: Run Analysis
```bash
python pko-test18_b_multithickness.py
```

The script will:
1. Load the YAML configuration
2. Extract thickness sweep parameters
3. Run pyKO simulations for each target thickness
4. Calculate FSV for each simulation
5. Generate comparative plots
6. Provide summary statistics

## Configuration Parameters

### Thickness Sweep
```yaml
thickness_sweep:
  constant_velocity: 300.0  # m/s - same for all simulations
  min_thickness: 50.0       # μm - minimum target thickness
  max_thickness: 500.0      # μm - maximum target thickness
  thickness_steps: 5        # Number of thickness values to simulate
```

### Timing
- **tstop**: Automatically calculated based on thickest target
- **dtstart**: Initial time step (0.1 ns)
- **dtoutput**: Output frequency (1 ns)

## Expected Results

### FSV vs Time Plot
- Multiple curves showing FSV evolution for different thicknesses
- Thicker targets show delayed shock arrival
- Peak FSV may vary with thickness

### Peak FSV vs Thickness Plot
- Correlation between target thickness and peak FSV
- Trend line showing relationship
- Statistical correlation coefficient

## Physics Insights

### Shock Arrival Time
- **Thicker targets**: Longer time for shock to reach free surface
- **Thinner targets**: Faster shock arrival, potential for multiple reflections

### FSV Amplification
- **Thickness effects**: May show optimal thickness for maximum FSV
- **Wave interactions**: Multiple reflections in thin targets
- **Material response**: Strength model effects on FSV

### Spall Behavior
- **Thin targets**: May experience complete spall
- **Thick targets**: Partial spall or no spall depending on conditions

## Comparison with Test 18a

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
2. **Mesh resolution**: Ensure sufficient cells for accuracy
3. **Strength model stability**: VM is most stable, SG may have issues
4. **File paths**: Ensure scripts run from correct directory

### Performance Tips
- Start with fewer thickness steps for testing
- Use VM strength model for stability
- Monitor simulation completion for each thickness
- Check output files are generated before plotting

## Example Results

Typical output shows:
- **Thickness range**: 50-500 μm
- **Peak FSV range**: Varies with thickness
- **Correlation**: May show positive or negative trend
- **Timing**: Shock arrival scales with thickness

## Next Steps

After running Test 18b, consider:
1. **Parameter sensitivity**: Vary impact velocity or materials
2. **Strength model comparison**: Test different models
3. **Mesh convergence**: Study effect of cell count
4. **Combined analysis**: Use both 18a and 18b results together
