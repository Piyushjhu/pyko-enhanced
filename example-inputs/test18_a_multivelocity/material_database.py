#!/usr/bin/env python3
"""
Material Database Module for PyKO Configuration
===============================================

This module provides a comprehensive database of common materials used in shock physics
simulations and functions to automatically generate YAML configurations based on user
material selections.

Author: PyKO Enhanced
Date: 2024
"""

import yaml
import os
import shutil
from typing import Dict, Any, Optional

# ===============================================================================
# MATERIAL DATABASE
# ===============================================================================

MATERIAL_DATABASE = {
    'Aluminum': {
        'density': 2700.0,        # kg/m³
        'sound_speed': 5200.0,    # m/s
        'gruneisen': 2.0,         # dimensionless
        'shear_modulus': 26.0E9,  # Pa
        'yield_strength': 207.0E6, # Pa
        'spall_strength': 276.0E6, # Pa
        'specific_heat': 896.0,   # J/kg·K
        'hugoniot_slope': 1.5,    # dimensionless
        'description': 'Lightweight, ductile metal with good shock properties'
    },
    
    'Copper': {
        'density': 8930.0,        # kg/m³
        'sound_speed': 3900.0,    # m/s
        'gruneisen': 1.99,        # dimensionless
        'shear_modulus': 46.0E9,  # Pa
        'yield_strength': 95.6E6, # Pa
        'spall_strength': 1.5E9,  # Pa
        'specific_heat': 385.0,   # J/kg·K
        'hugoniot_slope': 1.49,   # dimensionless
        'description': 'High conductivity metal with excellent shock properties'
    },
    
    'Steel': {
        'density': 7900.0,        # kg/m³
        'sound_speed': 4570.0,    # m/s
        'gruneisen': 2.17,        # dimensionless
        'shear_modulus': 77.0E9,  # Pa
        'yield_strength': 205.0E6, # Pa
        'spall_strength': 2.0E9,  # Pa
        'specific_heat': 502.0,   # J/kg·K
        'hugoniot_slope': 1.49,   # dimensionless
        'description': 'AISI 304 stainless steel - common structural material'
    },
    
    'Iron': {
        'density': 7870.0,        # kg/m³
        'sound_speed': 4600.0,    # m/s
        'gruneisen': 1.69,        # dimensionless
        'shear_modulus': 82.0E9,  # Pa
        'yield_strength': 180.0E6, # Pa
        'spall_strength': 1.8E9,  # Pa
        'specific_heat': 460.0,   # J/kg·K
        'hugoniot_slope': 1.92,   # dimensionless
        'description': 'Pure iron - fundamental material for shock studies'
    },
    
    'Titanium': {
        'density': 4500.0,        # kg/m³
        'sound_speed': 4780.0,    # m/s
        'gruneisen': 1.1,         # dimensionless
        'shear_modulus': 44.0E9,  # Pa
        'yield_strength': 830.0E6, # Pa
        'spall_strength': 1.2E9,  # Pa
        'specific_heat': 523.0,   # J/kg·K
        'hugoniot_slope': 1.14,   # dimensionless
        'description': 'High strength-to-weight ratio aerospace material'
    },
    
    'Tungsten': {
        'density': 19300.0,       # kg/m³
        'sound_speed': 4020.0,    # m/s
        'gruneisen': 1.62,        # dimensionless
        'shear_modulus': 161.0E9, # Pa
        'yield_strength': 550.0E6, # Pa
        'spall_strength': 3.0E9,  # Pa
        'specific_heat': 134.0,   # J/kg·K
        'hugoniot_slope': 1.24,   # dimensionless
        'description': 'Very dense refractory metal with high melting point'
    },
    
    'Gold': {
        'density': 19320.0,       # kg/m³
        'sound_speed': 3240.0,    # m/s
        'gruneisen': 2.97,        # dimensionless
        'shear_modulus': 27.0E9,  # Pa
        'yield_strength': 120.0E6, # Pa
        'spall_strength': 0.8E9,  # Pa
        'specific_heat': 129.0,   # J/kg·K
        'hugoniot_slope': 1.57,   # dimensionless
        'description': 'Noble metal with high density and ductility'
    },
    
    'Lead': {
        'density': 11340.0,       # kg/m³
        'sound_speed': 2050.0,    # m/s
        'gruneisen': 2.46,        # dimensionless
        'shear_modulus': 5.6E9,   # Pa
        'yield_strength': 18.0E6, # Pa
        'spall_strength': 0.3E9,  # Pa
        'specific_heat': 129.0,   # J/kg·K
        'hugoniot_slope': 1.46,   # dimensionless
        'description': 'Very soft, dense metal with low melting point'
    }
}

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

def get_available_materials() -> list:
    """
    Get list of available material names.
    
    Returns:
        list: List of available material names
    """
    return list(MATERIAL_DATABASE.keys())

def get_material_properties(material_name: str) -> Dict[str, Any]:
    """
    Get properties for a specific material.
    
    Args:
        material_name (str): Name of the material
        
    Returns:
        dict: Material properties dictionary
        
    Raises:
        ValueError: If material name is not found
    """
    if material_name not in MATERIAL_DATABASE:
        available = ', '.join(get_available_materials())
        raise ValueError(f"Material '{material_name}' not found. Available materials: {available}")
    
    return MATERIAL_DATABASE[material_name].copy()

def print_material_database():
    """
    Print the complete material database in a formatted table.
    """
    print("=" * 120)
    print("MATERIAL DATABASE FOR SHOCK PHYSICS SIMULATIONS")
    print("=" * 120)
    print(f"{'Material':<12} {'Density':<10} {'Sound':<8} {'Gruneisen':<10} {'Shear':<10} {'Yield':<10} {'Spall':<10} {'Description'}")
    print(f"{'Name':<12} {'(kg/m³)':<10} {'Speed':<8} {'Param':<10} {'Modulus':<10} {'Strength':<10} {'Strength':<10} {'(Brief)'}")
    print(f"{'':<12} {'':<10} {'(m/s)':<8} {'':<10} {'(GPa)':<10} {'(MPa)':<10} {'(MPa)':<10} {''}")
    print("-" * 120)
    
    for name, props in MATERIAL_DATABASE.items():
        print(f"{name:<12} {props['density']:<10.0f} {props['sound_speed']:<8.0f} "
              f"{props['gruneisen']:<10.2f} {props['shear_modulus']/1e9:<10.1f} "
              f"{props['yield_strength']/1e6:<10.0f} {props['spall_strength']/1e6:<10.0f} "
              f"{props['description']}")
    
    print("=" * 120)

def generate_material_yaml(material_name: str, material_id: int, 
                          thickness: float, cells: int, xstart: float,
                          initial_velocity: float = 0.0, 
                          enable_spall: bool = True) -> Dict[str, Any]:
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
        
    Returns:
        dict: YAML configuration dictionary for the material
    """
    props = get_material_properties(material_name)
    
    # Set spall parameters based on enable_spall flag
    if enable_spall:
        pfrac = props['spall_strength']
        nrhomin = 1.0  # No density reduction allowed
    else:
        pfrac = 1.0E20  # Effectively infinite (no spall)
        nrhomin = 0.8   # Allow some density reduction
    
    material_config = {
        f'mat{material_id}': {
            'mesh': {
                'cells': cells,
                'xstart': xstart,
                'length': thickness
            },
            'init': {
                'up0': initial_velocity,
                'rho0': props['density'],
                'p0': 0.0,
                'e0': 0.0,
                't0': 298.0
            },
            'eos': {
                'name': f'{material_name} {"flyer" if material_id == 1 else "target"}',
                'type': 'MGR',
                'rhoref': props['density'],
                'c0': props['sound_speed'],
                's1': props['hugoniot_slope'],
                'gamma0': props['gruneisen'],
                'cv': props['specific_heat']
            },
            'str': {
                'type': 'VM',
                'gmod': props['shear_modulus'],
                'ys': props['yield_strength']
            },
            'frac': {
                'pfrac': pfrac,
                'nrhomin': nrhomin
            }
        }
    }
    
    return material_config

def generate_complete_yaml_config(impactor_material: str, target_material: str,
                                 impactor_thickness: float, target_thickness: float,
                                 impactor_cells: int, target_cells: int,
                                 velocity_sweep: Dict[str, Any],
                                 timing_params: Dict[str, float],
                                 enable_spall: bool = True,
                                 test_name: str = "Custom Material Test") -> Dict[str, Any]:
    """
    Generate a complete YAML configuration file.
    
    Args:
        impactor_material (str): Name of impactor material
        target_material (str): Name of target material
        impactor_thickness (float): Impactor thickness in meters
        target_thickness (float): Target thickness in meters
        impactor_cells (int): Number of cells for impactor
        target_cells (int): Number of cells for target
        velocity_sweep (dict): Velocity sweep parameters
        timing_params (dict): Timing parameters (tstop, dtstart, dtoutput)
        enable_spall (bool): Whether to enable spall physics
        test_name (str): Name of the test
        
    Returns:
        dict: Complete YAML configuration
    """
    # Generate material configurations
    mat1_config = generate_material_yaml(
        impactor_material, 1, impactor_thickness, impactor_cells, 
        -impactor_thickness, velocity_sweep['min_velocity'], enable_spall
    )
    
    mat2_config = generate_material_yaml(
        target_material, 2, target_thickness, target_cells, 
        0.0, 0.0, enable_spall
    )
    
    # Complete configuration
    config = {
        'name': test_name,
        'outputfilename': f'./pyko-{test_name.lower().replace(" ", "-")}-bin.dat',
        'outputformat': 'BIN',
        'velocity_sweep': velocity_sweep,
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
        'density_threshold_factor': 0.1,  # Density threshold factor for FSV detection
        'right_half_factor': 0.5,  # Factor for right half domain analysis
        'debug_time_min': 0.04e-06,  # Debug time range start (seconds)
        'debug_time_max': 0.05e-06,  # Debug time range end (seconds)
        'debug_velocity_threshold': 200.0,  # Velocity threshold for debug output (m/s)
        'debug_fsv_count_threshold': 5,  # FSV count threshold for debug output
        'debug_initial_fsv_count': 3,  # Initial FSV count for debug output
        'debug_velocity_limit': 300.0,  # Velocity limit for debug output (m/s)
        'debug_rightmost_node_count': 5,  # Number of rightmost nodes to show in debug
        'interface_threshold': 0.01e-3,  # Interface detection threshold (meters)
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
        "# PyKO Shock Physics Simulation Configuration",
        "# ==============================================================================",
        "#",
        "# This file defines a shock physics simulation with the following setup:",
        f"# - Test Name: {config.get('name', 'Unknown')}",
        f"# - Impactor: {config.get('mat1', {}).get('eos', {}).get('name', 'Unknown')}",
        f"# - Target: {config.get('mat2', {}).get('eos', {}).get('name', 'Unknown')}",
        "# - Physics: Interface separation (spall) enabled",
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
    
    # Velocity sweep parameters
    if 'velocity_sweep' in config:
        vs = config['velocity_sweep']
        yaml_content.extend([
            "# ==============================================================================",
            "# VELOCITY SWEEP PARAMETERS",
            "# ==============================================================================",
            "# These parameters define the range of impact velocities to simulate",
            "# The script will automatically run simulations for each velocity in the range",
            "velocity_sweep:",
            f"  min_velocity: {vs.get('min_velocity', 100.0)}  # Minimum impact velocity (m/s)",
            f"  max_velocity: {vs.get('max_velocity', 1000.0)}  # Maximum impact velocity (m/s)",
            f"  velocity_steps: {vs.get('velocity_steps', 10)}  # Number of velocity steps to simulate",
            ""
        ])
    
    # Timing parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# TIMING PARAMETERS",
        "# ==============================================================================",
        "# These control the simulation duration and output frequency",
        f"tstop: {config.get('tstop', 1e-6)}  # Total simulation time (seconds) - when simulation ends",
        f"dtstart: {config.get('dtstart', 1e-10)}  # Initial time step (seconds) - starting time step size",
        f"dtoutput: {config.get('dtoutput', 1e-9)}  # Output frequency (seconds) - how often to save data",
        "# Note: Fixed timing used for all velocities - minimum 0.1 μs recommended",
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
            f"    cells: {mat1['mesh']['cells']}  # Number of computational cells",
            f"    xstart: {mat1['mesh']['xstart']}  # Starting x-position (meters)",
            f"    length: {mat1['mesh']['length']}  # Material thickness (meters)",
            "",
            "  # Initial conditions - state at t=0",
            "  init:",
            f"    up0: {mat1['init']['up0']}  # Initial velocity (m/s) - will be overridden by velocity sweep",
            f"    rho0: {mat1['init']['rho0']}  # Initial density (kg/m³)",
            f"    p0: {mat1['init']['p0']}  # Initial pressure (Pa)",
            f"    e0: {mat1['init']['e0']}  # Initial internal energy (J/kg)",
            f"    t0: {mat1['init']['t0']}  # Initial temperature (K)",
            "",
            "  # Equation of State (EOS) - defines pressure-density-energy relationship",
            "  eos:",
            f"    name: {mat1['eos']['name']}  # Material name for identification",
            f"    type: {mat1['eos']['type']}  # EOS type (MGR = Mie-Gruneisen)",
            f"    rhoref: {mat1['eos']['rhoref']}  # Reference density (kg/m³)",
            f"    c0: {mat1['eos']['c0']}  # Bulk sound speed (m/s)",
            f"    s1: {mat1['eos']['s1']}  # Hugoniot slope parameter (dimensionless)",
            f"    gamma0: {mat1['eos']['gamma0']}  # Gruneisen parameter (dimensionless)",
            f"    cv: {mat1['eos']['cv']}  # Specific heat capacity (J/kg/K)",
            "",
            "  # Strength model - defines material strength and plasticity",
            "  str:",
            f"    type: {mat1['str']['type']}  # Strength model (VM = Von Mises)",
            f"    gmod: {mat1['str']['gmod']}  # Shear modulus (Pa)",
            f"    ys: {mat1['str']['ys']}  # Yield strength (Pa)",
            "",
            "  # Fracture model - defines spall (tensile fracture) behavior",
            "  frac:",
            f"    pfrac: {mat1['frac']['pfrac']}  # Fracture pressure (Pa) - tensile stress threshold for spall",
            f"    nrhomin: {mat1['frac']['nrhomin']}  # Minimum density ratio - controls void formation during spall",
            ""
        ])
    
    # Material 2 (Target)
    if 'mat2' in config:
        mat2 = config['mat2']
        yaml_content.extend([
            "# ==============================================================================",
            "# MATERIAL 2: TARGET",
            "# ==============================================================================",
            "# This is the material that receives the impact",
            "mat2:",
            "  # Mesh definition - controls spatial discretization",
            "  mesh:",
            f"    cells: {mat2['mesh']['cells']}  # Number of computational cells",
            f"    xstart: {mat2['mesh']['xstart']}  # Starting x-position (meters)",
            f"    length: {mat2['mesh']['length']}  # Material thickness (meters)",
            "",
            "  # Initial conditions - state at t=0",
            "  init:",
            f"    up0: {mat2['init']['up0']}  # Initial velocity (m/s) - target is stationary",
            f"    rho0: {mat2['init']['rho0']}  # Initial density (kg/m³)",
            f"    p0: {mat2['init']['p0']}  # Initial pressure (Pa)",
            f"    e0: {mat2['init']['e0']}  # Initial internal energy (J/kg)",
            f"    t0: {mat2['init']['t0']}  # Initial temperature (K)",
            "",
            "  # Equation of State (EOS) - defines pressure-density-energy relationship",
            "  eos:",
            f"    name: {mat2['eos']['name']}  # Material name for identification",
            f"    type: {mat2['eos']['type']}  # EOS type (MGR = Mie-Gruneisen)",
            f"    rhoref: {mat2['eos']['rhoref']}  # Reference density (kg/m³)",
            f"    c0: {mat2['eos']['c0']}  # Bulk sound speed (m/s)",
            f"    s1: {mat2['eos']['s1']}  # Hugoniot slope parameter (dimensionless)",
            f"    gamma0: {mat2['eos']['gamma0']}  # Gruneisen parameter (dimensionless)",
            f"    cv: {mat2['eos']['cv']}  # Specific heat capacity (J/kg/K)",
            "",
            "  # Strength model - defines material strength and plasticity",
            "  str:",
            f"    type: {mat2['str']['type']}  # Strength model (VM = Von Mises)",
            f"    gmod: {mat2['str']['gmod']}  # Shear modulus (Pa)",
            f"    ys: {mat2['str']['ys']}  # Yield strength (Pa)",
            "",
            "  # Fracture model - defines spall (tensile fracture) behavior",
            "  frac:",
            f"    pfrac: {mat2['frac']['pfrac']}  # Fracture pressure (Pa) - tensile stress threshold for spall",
            f"    nrhomin: {mat2['frac']['nrhomin']}  # Minimum density ratio - controls void formation during spall",
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
            f"  ip0: {bc.get('ip0', 0.0)}  # Left boundary pressure (Pa)",
            f"  obc: {bc.get('obc', 'FREE')}  # Right boundary condition (FREE = free surface)",
            f"  op0: {bc.get('op0', 0.0)}  # Right boundary pressure (Pa)",
            ""
        ])
    
    # Geometry and other parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# GEOMETRY AND PHYSICS PARAMETERS",
        "# ==============================================================================",
        f"geometry: {config.get('geometry', 'PLA')}  # Geometry type (PLA = planar)",
        f"gravity: {config.get('gravity', 0.0)}  # Gravitational acceleration (m/s²)",
        f"pvoid: {config.get('pvoid', 0.0)}  # Void pressure (Pa)",
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
        ""
    ])
    
    # Final parameter
    yaml_content.extend([
        "# ==============================================================================",
        "# SIMULATION CONTROL",
        "# ==============================================================================",
        f"nmat: {config.get('nmat', 2)}  # Number of materials in simulation",
        "",
        "# ==============================================================================",
        "# END OF CONFIGURATION FILE",
        "# =============================================================================="
    ])
    
    # Write the commented YAML to file
    with open(filename, 'w') as f:
        f.write('\n'.join(yaml_content))

def create_material_config_files(impactor_material: str, target_material: str,
                                impactor_thickness: float = 0.0001,
                                target_thickness: float = 0.0002,
                                impactor_cells: int = 100,
                                target_cells: int = 200,
                                velocity_range: tuple = (100.0, 300.0),
                                velocity_steps: int = 5):
    """
    Create YAML configuration file with interface separation.
    
    Args:
        impactor_material (str): Name of impactor material
        target_material (str): Name of target material
        impactor_thickness (float): Impactor thickness in meters
        target_thickness (float): Target thickness in meters
        impactor_cells (int): Number of cells for impactor
        target_cells (int): Number of cells for target
        velocity_range (tuple): (min_velocity, max_velocity) in m/s
        velocity_steps (int): Number of velocity steps
    """
    # Velocity sweep parameters
    velocity_sweep = {
        'min_velocity': velocity_range[0],
        'max_velocity': velocity_range[1],
        'velocity_steps': velocity_steps
    }
    
    # Timing parameters (different for with/without spall)
    # Timing parameters (fixed timing for all velocities)
    timing_params = {
        'tstop': 0.150E-6,    # 0.15 μs (sufficient for most velocities)
        'dtstart': 0.0001E-6, # 0.0001 μs
        'dtoutput': 0.001E-6  # 0.001 μs (less frequent output)
    }
    test_name = f"Test 18a {impactor_material}-{target_material} Multivelocity"
    
    # Generate configuration
    config = generate_complete_yaml_config(
        impactor_material, target_material,
        impactor_thickness, target_thickness,
        impactor_cells, target_cells,
        velocity_sweep, timing_params,
        True, test_name  # Always enable spall/interface separation
    )
    
    # Determine filename (always with interface separation)
    filename = "test18.yml"
    
    # Save configuration
    save_yaml_config(config, filename)
    print(f"✅ Generated {filename}")
    print(f"   Materials: {impactor_material} → {target_material}")
    print(f"   Thicknesses: {impactor_thickness*1e6:.0f} μm → {target_thickness*1e6:.0f} μm")
    print(f"   Resolution: {impactor_cells} → {target_cells} cells")
    print(f"   Velocity range: {velocity_range[0]}-{velocity_range[1]} m/s ({velocity_steps} steps)")
    print(f"   Interface separation: Enabled")

# ===============================================================================
# MAIN FUNCTION FOR COMMAND LINE USE
# ===============================================================================

def main():
    """
    Interactive PyKO Material Configuration Generator
    """
    print("🎯 PyKO Material Configuration Generator")
    print("=" * 50)
    
    # Show available materials first
    print("\n📋 Available Materials:")
    print("-" * 30)
    materials = get_available_materials()
    for i, material in enumerate(materials, 1):
        print(f"{i:2d}. {material}")
    
    print("\n" + "=" * 50)
    print("🔧 Configuration Setup")
    print("=" * 50)
    
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
    
    # Get thicknesses
    print("\n📏 Thickness Configuration:")
    print("(Enter values in micrometers, e.g., 100 for 100 μm)")
    
    impactor_thickness = float(input(f"Impactor thickness (μm) [default: 100]: ").strip() or "100") / 1e6
    target_thickness = float(input(f"Target thickness (μm) [default: 200]: ").strip() or "200") / 1e6
    
    # Get cell counts
    print("\n🔢 Mesh Resolution:")
    impactor_cells = int(input(f"Impactor cells [default: 100]: ").strip() or "100")
    target_cells = int(input(f"Target cells [default: 200]: ").strip() or "200")
    
    # Get velocity sweep
    print("\n⚡ Velocity Sweep Configuration:")
    min_velocity = float(input(f"Minimum velocity (m/s) [default: 100]: ").strip() or "100")
    max_velocity = float(input(f"Maximum velocity (m/s) [default: 500]: ").strip() or "500")
    velocity_steps = int(input(f"Number of velocity steps [default: 5]: ").strip() or "5")
    
    # Confirm configuration
    print("\n" + "=" * 50)
    print("📋 Configuration Summary")
    print("=" * 50)
    print(f"Impactor: {impactor_material} ({impactor_thickness*1e6:.0f} μm, {impactor_cells} cells)")
    print(f"Target: {target_material} ({target_thickness*1e6:.0f} μm, {target_cells} cells)")
    print(f"Velocity range: {min_velocity}-{max_velocity} m/s ({velocity_steps} steps)")
    print(f"Interface separation: Enabled")
    
    confirm = input("\nGenerate configuration? (y/n) [default: y]: ").strip().lower()
    if confirm in ['', 'y', 'yes']:
        try:
            create_material_config_files(
                impactor_material, target_material,
                impactor_thickness, target_thickness,
                impactor_cells, target_cells,
                (min_velocity, max_velocity),
                velocity_steps
            )
            
            print("\n✅ Configuration generated successfully!")
            print("📁 File: test18.yml")
            print("🚀 Ready to run: python pko-test18_a_multivelocity.py")
            
        except ValueError as e:
            print(f"❌ Error: {e}")
            return 1
    else:
        print("❌ Configuration generation cancelled.")
    
    return 0

if __name__ == "__main__":
    exit(main())
