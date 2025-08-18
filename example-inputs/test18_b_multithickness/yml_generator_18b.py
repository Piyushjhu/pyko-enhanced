#!/usr/bin/env python3
"""
YAML Generator for pyKO Test 18b: Multi-Thickness Analysis
==========================================================

This script generates YAML configuration files for multi-thickness analysis
where impact velocity is constant but target thickness is varied.

Key Features:
- Material selection (impactor and target)
- Strength model selection (HYDRO, VM, SG)
- Thickness sweep configuration
- Constant impact velocity
- Interactive parameter input

Author: pyKO Development Team
Date: 2024
"""

import yaml
import os
from typing import Dict, List, Any

# Import the material database for Test 18b
from material_database_18b import MATERIAL_DATABASE, STRENGTH_MODELS, get_available_materials

def generate_material_yaml(material_name: str, material_id: int, 
                          thickness: float, cells: int, xstart: float,
                          initial_velocity: float = 0.0, 
                          enable_spall: bool = True,
                          strength_model: str = 'VM') -> Dict[str, Any]:
    """
    Generate YAML configuration for a material based on database properties.
    
    Args:
        material_name (str): Name of the material from database
        material_id (int): Material ID (1 for flyer, 2 for target)
        thickness (float): Material thickness in meters
        cells (int): Number of mesh cells
        xstart (float): Starting x-position in meters
        initial_velocity (float): Initial velocity in m/s
        enable_spall (bool): Whether to enable spall physics
        strength_model (str): Strength model type ('HYDRO', 'VM', 'SG')
        
    Returns:
        dict: YAML configuration dictionary for the material
    """
    if material_name not in MATERIAL_DATABASE:
        raise ValueError(f"Material '{material_name}' not found in database")
    
    material_props = MATERIAL_DATABASE[material_name]
    
    # Generate strength model configuration
    str_config = generate_strength_model_config(material_props, strength_model)
    
    # Generate fracture configuration
    if enable_spall:
        frac_config = {
            'pfrac': material_props['spall_strength'],
            'nrhomin': 1.0
        }
    else:
        frac_config = {
            'pfrac': 1e20,  # Very high value to disable spall
            'nrhomin': 0.8
        }
    
    # Material configuration
    mat_config = {
        f'mat{material_id}': {
            'mesh': {
                'cells': cells,
                'xstart': xstart,
                'length': thickness
            },
            'init': {
                'up0': initial_velocity,
                'rho0': material_props['density'],
                'p0': 0.0,
                'e0': 0.0,
                't0': 298.0
            },
            'eos': {
                'name': f'{material_name} {"flyer" if material_id == 1 else "target"}',
                'type': 'MGR',
                'rhoref': material_props['density'],
                'c0': material_props['sound_speed'],
                's1': material_props['hugoniot_slope'],
                'gamma0': material_props['gruneisen'],
                'cv': material_props['specific_heat']
            },
            'str': str_config,
            'frac': frac_config
        }
    }
    
    return mat_config

def generate_strength_model_config(material_props: Dict[str, Any], strength_model: str) -> Dict[str, Any]:
    """
    Generate strength model configuration based on material properties and model type.
    
    Args:
        material_props (dict): Material properties from database
        strength_model (str): Strength model type ('HYDRO', 'VM', 'SG')
        
    Returns:
        dict: Strength model configuration
    """
    if strength_model == 'HYDRO':
        return {
            'type': 'HYDRO'
        }
    elif strength_model == 'VM':
        return {
            'type': 'VM',
            'gmod': material_props['shear_modulus'],
            'ys': material_props['yield_strength']
        }
    elif strength_model == 'SG':
        return {
            'type': 'SG',
            'Y0': material_props['sg_Y0'],
            'Ymax': material_props['sg_Ymax'],
            'beta': material_props['sg_beta'],
            'n': material_props['sg_n'],
            'b': material_props['sg_b'],
            'h': material_props['sg_h'],
            'Tm0': material_props['sg_Tm0'],
            'mu0': material_props['sg_mu0']
        }
    else:
        raise ValueError(f"Unknown strength model: {strength_model}. Available: {list(STRENGTH_MODELS.keys())}")

def generate_complete_yaml_config(impactor_material: str, target_material: str,
                                 impactor_thickness: float, target_thickness: float,
                                 impactor_cells: int, target_cells: int,
                                 thickness_sweep: Dict[str, Any],
                                 timing_params: Dict[str, float],
                                 enable_spall: bool = True,
                                 impactor_strength_model: str = 'VM',
                                 target_strength_model: str = 'VM',
                                 test_name: str = "Test 18b Multithickness") -> Dict[str, Any]:
    """
    Generate a complete YAML configuration file for Test 18b.
    
    Args:
        impactor_material (str): Name of the impactor material
        target_material (str): Name of the target material
        impactor_thickness (float): Impactor thickness in meters
        target_thickness (float): Target thickness in meters (will be varied)
        impactor_cells (int): Number of cells for impactor
        target_cells (int): Number of cells for target (will be varied)
        thickness_sweep (dict): Thickness sweep parameters
        timing_params (dict): Timing parameters (tstop, dtstart, dtoutput)
        enable_spall (bool): Whether to enable spall physics
        impactor_strength_model (str): Strength model for impactor
        target_strength_model (str): Strength model for target
        test_name (str): Name of the test
        
    Returns:
        dict: Complete YAML configuration
    """
    # Generate material configurations
    mat1_config = generate_material_yaml(
        impactor_material, 1, impactor_thickness, impactor_cells, 
        -impactor_thickness, thickness_sweep['constant_velocity'], enable_spall,
        impactor_strength_model
    )
    
    mat2_config = generate_material_yaml(
        target_material, 2, target_thickness, target_cells, 
        0.0, 0.0, enable_spall,
        target_strength_model
    )
    
    # Complete configuration
    config = {
        'name': test_name,
        'outputfilename': f'./pyko-{test_name.lower().replace(" ", "-")}-bin.dat',
        'outputformat': 'BIN',
        'thickness_sweep': thickness_sweep,
        'tstop': timing_params['tstop'],
        'dtstart': timing_params['dtstart'],
        'dtoutput': timing_params['dtoutput'],
        **mat1_config,
        **mat2_config,
        'boundaries': {
            'ibc': 'FREE',
            'ip0': 0.0,
            'obc': 'FREE',
            'op0': 0.0
        },
        'geometry': 'PLA',
        'gravity': 0.0,
        'pvoid': 0.0,
        
        # Analysis configuration parameters
        'min_tstop': 0.1e-06,  # Minimum simulation time (seconds)
        
        # Units definition
        'units': {
            'time': 'second',
            'length': 'meter',
            'velocity': 'meter/second',
            'density': 'kg/m^3',
            'mass': 'kg',
            'pressure': 'Pa',
            'temperature': 'K',
            'energy': 'J',
            'sp_energy': 'J/kg',
            'sp_entropy': 'J/kg/K',
            'sp_heat_cap': 'J/kg/K',
            'gravity': 'm/s^2',
            's2': 'second/meter'
        },
        
        # EOS table units
        'tableunits': {
            'density': 'g/cm^3',
            'temperature': 'K',
            'pressure': 'GPa',
            'sp_energy': 'MJ/kg',
            'hfree_energy': 'MJ/kg',
            'sp_entropy': 'MJ/K/kg',
            'sound_speed': 'cm/s',
            'sp_heat_cap': 'MJ/kg/K'
        },
        
        # Code units - internal units used by PyKO (Wilkins book units)
        'codeunits': {
            'time': 'microseconds',
            'length': 'cm',
            'mass': 'g',
            'density': 'g/cm^3',
            'relative_volume': 'dimensionless',
            'velocity': 'cm/microsecond',
            'pressure': 'megabar',
            'temperature': 'K',
            'energy': 'eu',
            'sp_energy': 'eu/g',
            'sp_heat_cap': 'eu/K/g',
            'sp_entropy': 'eu/K/g',
            'ie_perv0': 'eu/cm^3',
            'cv_perv0': 'eu/cm^3/K',
            'gravity': 'cm/microseconds^2',
            's2': 'microseconds/cm'
        },
        
        # Simulation control
        'nmat': 2
    }
    
    return config

def save_yaml_config(config: Dict[str, Any], filename: str):
    """
    Save YAML configuration to file with comprehensive comments.
    
    Args:
        config (dict): Configuration dictionary
        filename (str): Output filename
    """
    # Create commented YAML content
    yaml_content = []
    
    # Header comments
    yaml_content.extend([
        "# ==============================================================================",
        "# pyKO Test 18b: Multi-Thickness Analysis Configuration",
        "# ==============================================================================",
        "#",
        "# This file defines a thickness sweep analysis with the following setup:",
        f"# - Test Name: {config.get('name', 'Unknown')}",
        f"# - Impactor: {config.get('mat1', {}).get('eos', {}).get('name', 'Unknown')}",
        f"# - Target: {config.get('mat2', {}).get('eos', {}).get('name', 'Unknown')}",
        "# - Physics: Thickness sweep analysis with constant impact velocity",
        "#",
        "# ==============================================================================",
        "# SIMULATION PARAMETERS",
        "# ==============================================================================",
        ""
    ])
    
    # Basic simulation parameters
    yaml_content.extend([
        f"name: {config.get('name', 'Unknown')}  # Test name for identification",
        f"outputfilename: {config.get('outputfilename', 'unknown.dat')}  # Binary output file path",
        f"outputformat: {config.get('outputformat', 'BIN')}  # Output format (BIN = binary)",
        ""
    ])
    
    # Thickness sweep parameters
    if 'thickness_sweep' in config:
        ts = config['thickness_sweep']
        yaml_content.extend([
            "# ==============================================================================",
            "# THICKNESS SWEEP PARAMETERS",
            "# ==============================================================================",
            "# These parameters define the range of target thicknesses to simulate",
            "# The script will automatically run simulations for each thickness in the range",
            "thickness_sweep:",
            f"  constant_velocity: {float(ts.get('constant_velocity', 300.0))}  # Constant impact velocity (m/s)",
            f"  min_thickness: {float(ts.get('min_thickness', 50.0))}  # Minimum target thickness (μm)",
            f"  max_thickness: {float(ts.get('max_thickness', 500.0))}  # Maximum target thickness (μm)",
            f"  thickness_steps: {int(ts.get('thickness_steps', 5))}  # Number of thickness steps to simulate",
            ""
        ])
    
    # Timing parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# TIMING PARAMETERS",
        "# ==============================================================================",
        "# These control the simulation duration and output frequency",
        f"tstop: {float(config.get('tstop', 1e-6))}  # Total simulation time (seconds) - when simulation ends",
        f"dtstart: {float(config.get('dtstart', 1e-10))}  # Initial time step (seconds) - starting time step size",
        f"dtoutput: {float(config.get('dtoutput', 1e-9))}  # Output frequency (seconds) - how often to save data",
        "# Note: Timing should be sufficient for the thickest target",
        "",

        ""
    ])
    
    # Material 1 (Impactor)
    if 'mat1' in config:
        mat1 = config['mat1']
        yaml_content.extend([
            "# ==============================================================================",
            "# MATERIAL 1: IMPACTOR (FLYER)",
            "# ==============================================================================",
            "# This is the material that impacts the target",
            "mat1:",
            "  # Mesh definition - controls spatial discretization",
            "  mesh:",
            f"    cells: {int(mat1['mesh']['cells'])}  # Number of computational cells",
            f"    xstart: {float(mat1['mesh']['xstart'])}  # Starting x-position (meters)",
            f"    length: {float(mat1['mesh']['length'])}  # Material thickness (meters)",
            "",
            "  # Initial conditions - state at t=0",
            "  init:",
            f"    up0: {float(mat1['init']['up0'])}  # Initial velocity (m/s) - constant for all thicknesses",
            f"    rho0: {float(mat1['init']['rho0'])}  # Initial density (kg/m³)",
            f"    p0: {float(mat1['init']['p0'])}  # Initial pressure (Pa)",
            f"    e0: {float(mat1['init']['e0'])}  # Initial internal energy (J/kg)",
            f"    t0: {float(mat1['init']['t0'])}  # Initial temperature (K)",
            "",
            "  # Equation of State (EOS) - defines pressure-density-energy relationship",
            "  eos:",
            f"    name: {mat1['eos']['name']}  # Material name for identification",
            f"    type: {mat1['eos']['type']}  # EOS type (MGR = Mie-Gruneisen)",
            f"    rhoref: {float(mat1['eos']['rhoref'])}  # Reference density (kg/m³)",
            f"    c0: {float(mat1['eos']['c0'])}  # Bulk sound speed (m/s)",
            f"    s1: {float(mat1['eos']['s1'])}  # Hugoniot slope parameter (dimensionless)",
            f"    gamma0: {float(mat1['eos']['gamma0'])}  # Gruneisen parameter (dimensionless)",
            f"    cv: {float(mat1['eos']['cv'])}  # Specific heat capacity (J/kg/K)",
            "",
            "  # Strength model - defines material strength and plasticity",
            "  str:",
            f"    type: {mat1['str']['type']}  # Strength model (HYDRO/VM/SG)",
        ])
        
        # Add strength model specific parameters
        if mat1['str']['type'] == 'HYDRO':
            yaml_content.extend([
                "    # Hydrodynamic model - no strength parameters needed"
            ])
        elif mat1['str']['type'] == 'VM':
            yaml_content.extend([
                f"    gmod: {float(mat1['str']['gmod'])}  # Shear modulus (Pa)",
                f"    ys: {float(mat1['str']['ys'])}  # Yield strength (Pa)"
            ])
        elif mat1['str']['type'] == 'SG':
            yaml_content.extend([
                f"    Y0: {float(mat1['str']['Y0'])}  # Initial yield strength (Pa)",
                f"    Ymax: {float(mat1['str']['Ymax'])}  # Maximum yield strength (Pa)",
                f"    beta: {float(mat1['str']['beta'])}  # Strain rate sensitivity parameter",
                f"    n: {float(mat1['str']['n'])}  # Strain hardening exponent",
                f"    b: {float(mat1['str']['b'])}  # Strain rate parameter",
                f"    h: {float(mat1['str']['h'])}  # Strain hardening parameter",
                f"    Tm0: {float(mat1['str']['Tm0'])}  # Melting temperature (K)",
                f"    mu0: {float(mat1['str']['mu0'])}  # Initial shear modulus (Pa)"
            ])
        
        yaml_content.extend([
            "",
            "  # Fracture model - defines spall (tensile fracture) behavior",
            "  frac:",
            f"    pfrac: {float(mat1['frac']['pfrac'])}  # Fracture pressure (Pa) - tensile stress threshold for spall",
            f"    nrhomin: {float(mat1['frac']['nrhomin'])}  # Minimum density ratio - controls void formation during spall",
            ""
        ])
    
    # Material 2 (Target)
    if 'mat2' in config:
        mat2 = config['mat2']
        yaml_content.extend([
            "# ==============================================================================",
            "# MATERIAL 2: TARGET",
            "# ==============================================================================",
            "# This is the material that receives the impact (thickness will be varied)",
            "mat2:",
            "  # Mesh definition - controls spatial discretization",
            "  mesh:",
            f"    cells: {int(mat2['mesh']['cells'])}  # Number of computational cells (will be varied)",
            f"    xstart: {float(mat2['mesh']['xstart'])}  # Starting x-position (meters)",
            f"    length: {float(mat2['mesh']['length'])}  # Material thickness (meters) (will be varied)",
            "",
            "  # Initial conditions - state at t=0",
            "  init:",
            f"    up0: {float(mat2['init']['up0'])}  # Initial velocity (m/s) - target is stationary",
            f"    rho0: {float(mat2['init']['rho0'])}  # Initial density (kg/m³)",
            f"    p0: {float(mat2['init']['p0'])}  # Initial pressure (Pa)",
            f"    e0: {float(mat2['init']['e0'])}  # Initial internal energy (J/kg)",
            f"    t0: {float(mat2['init']['t0'])}  # Initial temperature (K)",
            "",
            "  # Equation of State (EOS) - defines pressure-density-energy relationship",
            "  eos:",
            f"    name: {mat2['eos']['name']}  # Material name for identification",
            f"    type: {mat2['eos']['type']}  # EOS type (MGR = Mie-Gruneisen)",
            f"    rhoref: {float(mat2['eos']['rhoref'])}  # Reference density (kg/m³)",
            f"    c0: {float(mat2['eos']['c0'])}  # Bulk sound speed (m/s)",
            f"    s1: {float(mat2['eos']['s1'])}  # Hugoniot slope parameter (dimensionless)",
            f"    gamma0: {float(mat2['eos']['gamma0'])}  # Gruneisen parameter (dimensionless)",
            f"    cv: {float(mat2['eos']['cv'])}  # Specific heat capacity (J/kg/K)",
            "",
            "  # Strength model - defines material strength and plasticity",
            "  str:",
            f"    type: {mat2['str']['type']}  # Strength model (HYDRO/VM/SG)",
        ])
        
        # Add strength model specific parameters
        if mat2['str']['type'] == 'HYDRO':
            yaml_content.extend([
                "    # Hydrodynamic model - no strength parameters needed"
            ])
        elif mat2['str']['type'] == 'VM':
            yaml_content.extend([
                f"    gmod: {float(mat2['str']['gmod'])}  # Shear modulus (Pa)",
                f"    ys: {float(mat2['str']['ys'])}  # Yield strength (Pa)"
            ])
        elif mat2['str']['type'] == 'SG':
            yaml_content.extend([
                f"    Y0: {float(mat2['str']['Y0'])}  # Initial yield strength (Pa)",
                f"    Ymax: {float(mat2['str']['Ymax'])}  # Maximum yield strength (Pa)",
                f"    beta: {float(mat2['str']['beta'])}  # Strain rate sensitivity parameter",
                f"    n: {float(mat2['str']['n'])}  # Strain hardening exponent",
                f"    b: {float(mat2['str']['b'])}  # Strain rate parameter",
                f"    h: {float(mat2['str']['h'])}  # Strain hardening parameter",
                f"    Tm0: {float(mat2['str']['Tm0'])}  # Melting temperature (K)",
                f"    mu0: {float(mat2['str']['mu0'])}  # Initial shear modulus (Pa)"
            ])
        
        yaml_content.extend([
            "",
            "  # Fracture model - defines spall (tensile fracture) behavior",
            "  frac:",
            f"    pfrac: {float(mat2['frac']['pfrac'])}  # Fracture pressure (Pa) - tensile stress threshold for spall",
            f"    nrhomin: {float(mat2['frac']['nrhomin'])}  # Minimum density ratio - controls void formation during spall",
            ""
        ])
    
    # Boundary conditions
    if 'boundaries' in config:
        bc = config['boundaries']
        yaml_content.extend([
            "# ==============================================================================",
            "# BOUNDARY CONDITIONS",
            "# ==============================================================================",
            "# These define how the simulation domain behaves at its edges",
            "boundaries:",
            f"  ibc: {bc.get('ibc', 'FREE')}  # Left boundary condition (FREE = free surface)",
            f"  ip0: {float(bc.get('ip0', 0.0))}  # Left boundary pressure (Pa)",
            f"  obc: {bc.get('obc', 'FREE')}  # Right boundary condition (FREE = free surface)",
            f"  op0: {float(bc.get('op0', 0.0))}  # Right boundary pressure (Pa)",
            ""
        ])
    
    # Geometry and other parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# GEOMETRY AND PHYSICS PARAMETERS",
        "# ==============================================================================",
        f"geometry: {config.get('geometry', 'PLA')}  # Geometry type (PLA = planar)",
        f"gravity: {float(config.get('gravity', 0.0))}  # Gravitational acceleration (m/s²)",
        f"pvoid: {float(config.get('pvoid', 0.0))}  # Void pressure (Pa)",
        ""
    ])
    
    # Units section
    yaml_content.extend([
        "# ==============================================================================",
        "# UNITS DEFINITION",
        "# ==============================================================================",
        "# User input parameters units - what units you specify values in",
        "# PyKO will automatically convert these to code units internally",
        "units:",
        "  time: second  # Time units for tstop, dtstart, dtoutput",
        "  length: meter  # Length units for mesh dimensions",
        "  velocity: meter/second  # Velocity units for up0",
        "  density: kg/m^3  # Density units for rho0, rhoref",
        "  mass: kg  # Mass units",
        "  pressure: Pa  # Pressure units for p0, pfrac, gmod, ys",
        "  temperature: K  # Temperature units for t0",
        "  energy: J  # Energy units",
        "  sp_energy: J/kg  # Specific energy units for e0",
        "  sp_entropy: J/kg/K  # Specific entropy units",
        "  sp_heat_cap: J/kg/K  # Specific heat capacity units for cv",
        "  gravity: m/s^2  # Gravitational acceleration units",
        "  s2: second/meter  # Hugoniot parameter units",
        "",
        "# EOS table units - for ANEOS tables (not used in this MGR configuration)",
        "tableunits:",
        "  density: g/cm^3  # Density units in EOS tables",
        "  temperature: K  # Temperature units in EOS tables",
        "  pressure: GPa  # Pressure units in EOS tables",
        "  sp_energy: MJ/kg  # Specific energy units in EOS tables",
        "  hfree_energy: MJ/kg  # Helmholtz free energy units",
        "  sp_entropy: MJ/K/kg  # Specific entropy units in EOS tables",
        "  sound_speed: cm/s  # Sound speed units in EOS tables",
        "  sp_heat_cap: MJ/kg/K  # Specific heat capacity units in EOS tables",
        "",
        "# Code units - internal units used by PyKO (Wilkins book units)",
        "# These are automatically converted from user units",
        "codeunits:",
        "  time: microseconds  # Internal time units",
        "  length: cm  # Internal length units",
        "  mass: g  # Internal mass units",
        "  density: g/cm^3  # Internal density units",
        "  relative_volume: dimensionless  # Volume ratio units",
        "  velocity: cm/microsecond  # Internal velocity units (10 km/s)",
        "  pressure: megabar  # Internal pressure units (100 GPa)",
        "  temperature: K  # Internal temperature units",
        "  energy: eu  # Internal energy units (10^12 ergs)",
        "  sp_energy: eu/g  # Internal specific energy units",
        "  sp_heat_cap: eu/K/g  # Internal specific heat capacity units",
        "  sp_entropy: eu/K/g  # Internal specific entropy units",
        "  ie_perv0: eu/cm^3  # Internal energy per volume units",
        "  cv_perv0: eu/cm^3/K  # Internal heat capacity per volume units",
        "  gravity: cm/microseconds^2  # Internal gravitational acceleration units",
        "  s2: microseconds/cm  # Internal Hugoniot parameter units",
        "",
        "# ==============================================================================",
        "# SIMULATION CONTROL",
        "# ==============================================================================",
        "nmat: 2  # Number of materials in simulation"
    ])
    
    # Force remove existing file to ensure clean overwrite
    try:
        if os.path.exists(filename):
            os.remove(filename)
            print(f"🗑️  Removed existing {filename}")
    except Exception as e:
        print(f"⚠️  Warning: Could not remove existing {filename}: {e}")
    
    # Write to file
    try:
        with open(filename, 'w') as file:
            file.write('\n'.join(yaml_content))
        print(f"✅ Successfully wrote {filename}")
    except Exception as e:
        print(f"❌ Error writing {filename}: {e}")
        print(f"   Current working directory: {os.getcwd()}")
        print(f"   File path: {os.path.abspath(filename)}")
        raise

def create_material_config_files(impactor_material: str, target_material: str,
                                impactor_thickness: float, target_thickness: float,
                                impactor_cells: int, target_cells: int,
                                constant_velocity: float,
                                min_thickness: float, max_thickness: float, thickness_steps: int,
                                enable_spall: bool = True,
                                impactor_strength_model: str = 'VM',
                                target_strength_model: str = 'VM'):
    """
    Create material configuration files for Test 18b.
    
    Args:
        impactor_material (str): Name of the impactor material
        target_material (str): Name of the target material
        impactor_thickness (float): Impactor thickness in meters
        target_thickness (float): Target thickness in meters (base value)
        impactor_cells (int): Number of cells for impactor
        target_cells (int): Number of cells for target (base value)
        constant_velocity (float): Constant impact velocity in m/s
        min_thickness (float): Minimum target thickness in μm
        max_thickness (float): Maximum target thickness in μm
        thickness_steps (int): Number of thickness steps
        enable_spall (bool): Whether to enable spall physics
        impactor_strength_model (str): Strength model for impactor
        target_strength_model (str): Strength model for target
    """
    # Thickness sweep parameters
    thickness_sweep = {
        'constant_velocity': constant_velocity,
        'min_thickness': min_thickness,
        'max_thickness': max_thickness,
        'thickness_steps': thickness_steps
    }
    
    # Timing parameters - ensure sufficient time for thickest target
    max_thickness_m = max_thickness / 1e6  # Convert to meters
    sound_speed = MATERIAL_DATABASE[target_material]['sound_speed']
    min_simulation_time = max_thickness_m / sound_speed * 3  # 3x shock transit time
    
    timing_params = {
        'tstop': max(0.5e-06, min_simulation_time),  # At least 0.5 μs or 3x shock transit
        'dtstart': 1.0e-10,    # Initial time step (0.1 ns)
        'dtoutput': 1.0e-09    # Output frequency (1 ns)
    }
    
    test_name = f"Test 18b {impactor_material}-{target_material} Multithickness"
    
    # Generate configuration
    config = generate_complete_yaml_config(
        impactor_material, target_material,
        impactor_thickness, target_thickness,
        impactor_cells, target_cells,
        thickness_sweep, timing_params,
        enable_spall,
        impactor_strength_model,
        target_strength_model,
        test_name
    )
    
    # Determine filename - save in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(script_dir, "test18b.yml")
    
    # Save configuration (this function will handle file overwriting)
    save_yaml_config(config, filename)
    print(f"✅ Generated {filename}")
    print(f"   Materials: {impactor_material} → {target_material}")
    print(f"   Strength models: {STRENGTH_MODELS[impactor_strength_model]['name']} → {STRENGTH_MODELS[target_strength_model]['name']}")
    print(f"   Impactor thickness: {impactor_thickness*1e6:.0f} μm")
    print(f"   Target thickness range: {min_thickness:.0f} - {max_thickness:.0f} μm ({thickness_steps} steps)")
    print(f"   Constant impact velocity: {constant_velocity:.0f} m/s")
    print(f"   Simulation time: {timing_params['tstop']*1e6:.1f} μs")

def main():
    """
    Interactive Test 18b Configuration Generator
    """
    print("🎯 pyKO Test 18b: Multi-Thickness Analysis Configuration Generator")
    print("=" * 70)
    
    # Show available materials first
    print("\n📋 Available Materials:")
    print("-" * 30)
    materials = get_available_materials()
    for i, material in enumerate(materials, 1):
        print(f"{i:2d}. {material}")
    
    print("\n" + "=" * 70)
    print("🔧 Configuration Setup")
    print("=" * 70)
    
    # Get impactor material
    while True:
        print(f"\nSelect impactor material (1-{len(materials)}):")
        try:
            choice = int(input("Enter number: ").strip())
            if 1 <= choice <= len(materials):
                impactor_material = materials[choice - 1]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(materials)}")
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Get target material
    while True:
        print(f"\nSelect target material (1-{len(materials)}):")
        try:
            choice = int(input("Enter number: ").strip())
            if 1 <= choice <= len(materials):
                target_material = materials[choice - 1]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(materials)}")
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Show available strength models
    print("\n🔧 Available Strength Models:")
    print("-" * 30)
    strength_models = list(STRENGTH_MODELS.keys())
    for i, model in enumerate(strength_models, 1):
        model_info = STRENGTH_MODELS[model]
        print(f"{i:2d}. {model} - {model_info['name']}")
        print(f"    {model_info['description']}")
    
    # Get impactor strength model
    while True:
        print(f"\nSelect impactor strength model (1-{len(strength_models)}):")
        try:
            choice = int(input("Enter number: ").strip())
            if 1 <= choice <= len(strength_models):
                impactor_strength_model = strength_models[choice - 1]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(strength_models)}")
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Get target strength model
    while True:
        print(f"\nSelect target strength model (1-{len(strength_models)}):")
        try:
            choice = int(input("Enter number: ").strip())
            if 1 <= choice <= len(strength_models):
                target_strength_model = strength_models[choice - 1]
                break
            else:
                print(f"❌ Please enter a number between 1 and {len(strength_models)}")
        except ValueError:
            print("❌ Please enter a valid number")
    
    # Get thicknesses
    print("\n📏 Thickness Configuration:")
    print("(Enter values in micrometers, e.g., 100 for 100 μm)")
    
    impactor_thickness = float(input(f"Impactor thickness (μm) [default: 100]: ").strip() or "100") / 1e6
    target_thickness = float(input(f"Base target thickness (μm) [default: 200]: ").strip() or "200") / 1e6
    
    # Get cell counts
    print("\n🔢 Mesh Resolution:")
    impactor_cells = int(input(f"Impactor cells [default: 100]: ").strip() or "100")
    target_cells = int(input(f"Base target cells [default: 200]: ").strip() or "200")
    
    # Get constant impact velocity
    print("\n⚡ Impact Velocity Configuration:")
    constant_velocity = float(input(f"Constant impact velocity (m/s) [default: 300]: ").strip() or "300")
    
    # Get thickness sweep parameters
    print("\n📊 Thickness Sweep Configuration:")
    min_thickness = float(input(f"Minimum target thickness (μm) [default: 50]: ").strip() or "50")
    max_thickness = float(input(f"Maximum target thickness (μm) [default: 500]: ").strip() or "500")
    thickness_steps = int(input(f"Number of thickness steps [default: 5]: ").strip() or "5")
    
    print("\n" + "=" * 70)
    print("📋 Configuration Summary")
    print("=" * 70)
    print(f"Impactor: {impactor_material} ({impactor_thickness*1e6:.0f} μm, {impactor_cells} cells)")
    print(f"  Strength model: {STRENGTH_MODELS[impactor_strength_model]['name']}")
    print(f"Target: {target_material} (thickness sweep: {min_thickness:.0f}-{max_thickness:.0f} μm)")
    print(f"  Strength model: {STRENGTH_MODELS[target_strength_model]['name']}")
    print(f"Constant impact velocity: {constant_velocity:.0f} m/s")
    print(f"Thickness steps: {thickness_steps}")
    
    # Confirm generation
    response = input("\nGenerate configuration? (y/n) [default: y]: ").strip().lower()
    if response in ['', 'y', 'yes']:
        create_material_config_files(
            impactor_material, target_material,
            impactor_thickness, target_thickness,
            impactor_cells, target_cells,
            constant_velocity,
            min_thickness, max_thickness, thickness_steps,
            True,  # Enable spall
            impactor_strength_model,
            target_strength_model
        )
        print("\n✅ Configuration generated successfully!")
        print(f"📁 File: test18b.yml")
        print(f"🚀 Ready to run: python pko-test18_b_multithickness.py")
    else:
        print("\n❌ Configuration generation cancelled.")

if __name__ == "__main__":
    main()
