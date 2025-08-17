# Material Selection System for PyKO Test 18a

## Overview

This directory contains a comprehensive material selection system for PyKO shock physics simulations. The system allows users to easily configure different material combinations for impactor and target without manually editing complex YAML files.

## 🎯 Quick Start

### 1. List Available Materials
```bash
python material_database.py --list-materials
```

### 2. Generate Custom Configuration
```bash
python material_database.py --impactor Aluminum --target Copper --min-velocity 100 --max-velocity 500 --velocity-steps 5
```

### 3. Run Example Configurations
```bash
python example_material_usage.py
```

## 📋 Available Materials

The system includes 8 common materials used in shock physics:

| Material | Density (kg/m³) | Sound Speed (m/s) | Spall Strength (MPa) | Description |
|----------|----------------|-------------------|---------------------|-------------|
| **Aluminum** | 2700 | 5200 | 276 | Lightweight, ductile metal |
| **Copper** | 8930 | 3900 | 1500 | High conductivity metal |
| **Steel** | 7900 | 4570 | 2000 | AISI 304 stainless steel |
| **Iron** | 7870 | 4600 | 1800 | Pure iron |
| **Titanium** | 4500 | 4780 | 1200 | Aerospace material |
| **Tungsten** | 19300 | 4020 | 3000 | Very dense refractory metal |
| **Gold** | 19320 | 3240 | 800 | Noble metal |
| **Lead** | 11340 | 2050 | 300 | Very soft, dense metal |

## 🔧 Command Line Usage

### Basic Syntax
```bash
python material_database.py [OPTIONS]
```

### Options
- `--list-materials`: Show all available materials
- `--impactor MATERIAL`: Set impactor material (default: Aluminum)
- `--target MATERIAL`: Set target material (default: Aluminum)
- `--impactor-thickness METERS`: Impactor thickness (default: 0.0001 m = 100 μm)
- `--target-thickness METERS`: Target thickness (default: 0.0002 m = 200 μm)
- `--impactor-cells N`: Number of cells for impactor (default: 100)
- `--target-cells N`: Number of cells for target (default: 200)
- `--min-velocity M/S`: Minimum impact velocity (default: 100 m/s)
- `--max-velocity M/S`: Maximum impact velocity (default: 300 m/s)
- `--velocity-steps N`: Number of velocity steps (default: 5)


### Examples

#### 1. Classic Aluminum on Copper
```bash
python material_database.py --impactor Aluminum --target Copper --min-velocity 100 --max-velocity 500 --velocity-steps 5
```

#### 2. High-Speed Steel on Steel
```bash
python material_database.py --impactor Steel --target Steel --min-velocity 500 --max-velocity 2000 --velocity-steps 8
```

#### 3. Dense Tungsten on Light Aluminum
```bash
python material_database.py --impactor Tungsten --target Aluminum --impactor-thickness 0.00005 --target-thickness 0.0003 --min-velocity 300 --max-velocity 1000
```



## 📁 Generated Files

For each configuration, the system generates two YAML files:

1. **`test18-with-interface-separation.yml`**: With spall physics enabled
2. **`test18-without-interface-separation.yml`**: With spall physics disabled

### File Structure
```
test18_a_multivelocity/
├── material_database.py          # Material database and generation tools
├── example_material_usage.py     # Example usage demonstrations
├── pko-test18_a_multivelocity.py # Main simulation script
├── test18-with-interface-separation.yml    # Generated config (with spall)
├── test18-without-interface-separation.yml # Generated config (no spall)
└── README_MATERIAL_SELECTION.md  # This file
```

## 🎯 Material Properties

Each material includes complete Mie-Gruneisen EOS parameters:

- **Density (ρ₀)**: Initial material density
- **Sound Speed (c₀)**: Bulk sound speed
- **Gruneisen Parameter (γ₀)**: Thermodynamic parameter
- **Hugoniot Slope (s₁)**: Linear Hugoniot parameter
- **Shear Modulus (G)**: Elastic shear modulus
- **Yield Strength (σy)**: Plastic yield strength
- **Spall Strength (σspall)**: Dynamic fracture threshold
- **Specific Heat (cv)**: Heat capacity at constant volume

## 🔬 Physics Considerations

### Material Selection Guidelines

1. **Density Contrast**: Higher density contrast leads to stronger shock waves
2. **Impedance Matching**: Similar acoustic impedances (ρ₀ × c₀) reduce reflection
3. **Spall Strength**: Higher spall strength materials resist fracture better
4. **Wave Speed**: Faster materials transmit shocks more quickly

### Recommended Combinations

| Application | Impactor | Target | Reason |
|-------------|----------|--------|--------|
| **Classic Shock** | Aluminum | Copper | Good impedance contrast |
| **Same Material** | Steel | Steel | No interface effects |
| **High Energy** | Tungsten | Lead | Extreme density contrast |
| **Aerospace** | Titanium | Aluminum | Lightweight combination |
| **Noble Metal** | Gold | Copper | High conductivity |

### Velocity Ranges

- **Low (100-500 m/s)**: Quasi-static compression, minimal spall
- **Medium (500-2000 m/s)**: Shock formation, moderate spall
- **High (2000+ m/s)**: Strong shocks, significant spall likely
- **Hypervelocity (5000+ m/s)**: Extreme conditions, extensive damage

## 🛠️ Advanced Usage

### Custom Material Properties

To add new materials, edit `material_database.py` and add entries to `MATERIAL_DATABASE`:

```python
'NewMaterial': {
    'density': 5000.0,        # kg/m³
    'sound_speed': 4000.0,    # m/s
    'gruneisen': 1.8,         # dimensionless
    'shear_modulus': 50.0E9,  # Pa
    'yield_strength': 300.0E6, # Pa
    'spall_strength': 1.0E9,  # Pa
    'specific_heat': 400.0,   # J/kg·K
    'hugoniot_slope': 1.4,    # dimensionless
    'description': 'Description of the material'
}
```

### Batch Processing

Create multiple configurations programmatically:

```python
from material_database import create_material_config_files

materials = ['Aluminum', 'Copper', 'Steel', 'Iron']
for impactor in materials:
    for target in materials:
        create_material_config_files(
            impactor_material=impactor,
            target_material=target,
            velocity_range=(100, 500),
            velocity_steps=5
        )
```

## 📊 Analysis and Results

After running simulations, analyze results using the main script:

```bash
python pko-test18_a_multivelocity.py
```

This generates:
- **FSV vs Time plots**: Free surface velocity evolution
- **Peak FSV vs Impact Velocity**: Velocity sweep analysis
- **Material property tables**: Configuration summary

## 🚨 Troubleshooting

### Common Issues

1. **Material Not Found**: Check spelling and use `--list-materials`
2. **Invalid Parameters**: Ensure positive values for thickness, velocity, cells
3. **File Generation Errors**: Check write permissions in directory
4. **Simulation Failures**: Reduce velocity range or increase mesh resolution

### Performance Tips

- **Coarse Mesh**: Start with fewer cells for quick testing
- **Velocity Range**: Use smaller ranges for initial exploration
- **Spall Physics**: Disable for faster simulations without fracture
- **Resolution**: Increase cells for better accuracy (slower computation)

## 📚 References

- **Mie-Gruneisen EOS**: Standard equation of state for shock physics
- **Spall Physics**: Dynamic fracture under tensile loading
- **Hugoniot Relations**: Shock wave equations of state
- **Material Properties**: Standard values from shock physics literature

## 🤝 Contributing

To add new materials or improve the system:

1. Add material properties to `MATERIAL_DATABASE`
2. Update documentation in this README
3. Test with example configurations
4. Verify physical consistency of parameters

---

**Note**: All material properties are based on standard shock physics literature and may need adjustment for specific applications or material conditions.
