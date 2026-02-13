#!/usr/bin/env python3
"""
Material Database Module for PyKO Configuration - Test17 Version
===============================================================

This module provides a comprehensive database of common materials used in shock physics
simulations and functions to automatically generate YAML configurations based on user
material selections. Modified for Test17 with single velocity impact (no velocity sweeps).

Author: PyKO Enhanced
Date: 2024
"""

import yaml
import os
import shutil
from typing import Dict, Any, Optional

# ===============================================================================
# STRENGTH MODEL DATABASE
# ===============================================================================

STRENGTH_MODELS = {
    'HYDRO': {
        'name': 'Hydrodynamic',
        'description': 'No strength - material behaves as a fluid',
        'parameters': [],
        'type': 'HYDRO'
    },
    'VM': {
        'name': 'Von Mises',
        'description': 'Classic Von Mises yield criterion with constant yield strength',
        'parameters': ['gmod', 'ys'],
        'type': 'VM'
    },
    'SG': {
        'name': 'Steinberg-Guinan',
        'description': 'Advanced rate-dependent strength model with temperature and strain rate effects',
        'parameters': ['Y0', 'Ymax', 'beta', 'n', 'b', 'h', 'Tm0', 'mu0'],
        'type': 'SG'
    },
}

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
        'description': 'Lightweight, ductile metal with good shock properties',
        # Steinberg-Guinan parameters
        'sg_Y0': 0.29E9,          # Initial yield strength (Pa)
        'sg_Ymax': 0.68E9,        # Maximum yield strength (Pa)
        'sg_beta': 125.0,         # Strain rate sensitivity parameter
        'sg_n': 0.1,              # Strain hardening exponent
        'sg_b': 8.0,              # Strain rate parameter
        'sg_h': 6.2E-4,           # Strain hardening parameter
        'sg_Tm0': 1220.0,         # Melting temperature (K)
        'sg_mu0': 27.6E9,         # Initial shear modulus (Pa)
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
        'description': 'High conductivity metal with excellent shock properties',
        # Steinberg-Guinan parameters
        'sg_Y0': 0.12E9,          # Initial yield strength (Pa)
        'sg_Ymax': 0.64E9,        # Maximum yield strength (Pa)
        'sg_beta': 36.0,          # Strain rate sensitivity parameter
        'sg_n': 0.45,             # Strain hardening exponent
        'sg_b': 12.0,             # Strain rate parameter
        'sg_h': 0.025,            # Strain hardening parameter
        'sg_Tm0': 1358.0,         # Melting temperature (K)
        'sg_mu0': 47.7E9,         # Initial shear modulus (Pa)
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
        'description': 'AISI 304 stainless steel - common structural material',
        # Steinberg-Guinan parameters (estimated for AISI 304)
        'sg_Y0': 0.205E9,         # Initial yield strength (Pa)
        'sg_Ymax': 0.8E9,         # Maximum yield strength (Pa)
        'sg_beta': 80.0,          # Strain rate sensitivity parameter
        'sg_n': 0.2,              # Strain hardening exponent
        'sg_b': 10.0,             # Strain rate parameter
        'sg_h': 0.01,             # Strain hardening parameter
        'sg_Tm0': 1673.0,         # Melting temperature (K)
        'sg_mu0': 77.0E9,         # Initial shear modulus (Pa)
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
        'description': 'Pure iron - fundamental material for shock studies',
        # Steinberg-Guinan parameters (estimated for pure iron)
        'sg_Y0': 0.18E9,          # Initial yield strength (Pa)
        'sg_Ymax': 0.7E9,         # Maximum yield strength (Pa)
        'sg_beta': 70.0,          # Strain rate sensitivity parameter
        'sg_n': 0.25,             # Strain hardening exponent
        'sg_b': 9.0,              # Strain rate parameter
        'sg_h': 0.015,            # Strain hardening parameter
        'sg_Tm0': 1811.0,         # Melting temperature (K)
        'sg_mu0': 82.0E9,         # Initial shear modulus (Pa)
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
        'description': 'High strength-to-weight ratio aerospace material',
        # Steinberg-Guinan parameters (estimated for Ti-6Al-4V)
        'sg_Y0': 0.83E9,          # Initial yield strength (Pa)
        'sg_Ymax': 1.2E9,         # Maximum yield strength (Pa)
        'sg_beta': 60.0,          # Strain rate sensitivity parameter
        'sg_n': 0.15,             # Strain hardening exponent
        'sg_b': 7.0,              # Strain rate parameter
        'sg_h': 0.008,            # Strain hardening parameter
        'sg_Tm0': 1943.0,         # Melting temperature (K)
        'sg_mu0': 44.0E9,         # Initial shear modulus (Pa)
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
        'description': 'Very dense refractory metal with high melting point',
        # Steinberg-Guinan parameters (estimated for pure tungsten)
        'sg_Y0': 0.55E9,          # Initial yield strength (Pa)
        'sg_Ymax': 1.5E9,         # Maximum yield strength (Pa)
        'sg_beta': 100.0,         # Strain rate sensitivity parameter
        'sg_n': 0.1,              # Strain hardening exponent
        'sg_b': 15.0,             # Strain rate parameter
        'sg_h': 0.005,            # Strain hardening parameter
        'sg_Tm0': 3695.0,         # Melting temperature (K)
        'sg_mu0': 161.0E9,        # Initial shear modulus (Pa)
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
        'description': 'Noble metal with high density and ductility',
        # Steinberg-Guinan parameters (estimated for pure gold)
        'sg_Y0': 0.12E9,          # Initial yield strength (Pa)
        'sg_Ymax': 0.4E9,         # Maximum yield strength (Pa)
        'sg_beta': 50.0,          # Strain rate sensitivity parameter
        'sg_n': 0.3,              # Strain hardening exponent
        'sg_b': 6.0,              # Strain rate parameter
        'sg_h': 0.02,             # Strain hardening parameter
        'sg_Tm0': 1337.0,         # Melting temperature (K)
        'sg_mu0': 27.0E9,         # Initial shear modulus (Pa)
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
        'description': 'Very soft, dense metal with low melting point',
        # Steinberg-Guinan parameters (estimated for pure lead)
        'sg_Y0': 0.018E9,         # Initial yield strength (Pa)
        'sg_Ymax': 0.1E9,         # Maximum yield strength (Pa)
        'sg_beta': 30.0,          # Strain rate sensitivity parameter
        'sg_n': 0.5,              # Strain hardening exponent
        'sg_b': 4.0,              # Strain rate parameter
        'sg_h': 0.05,             # Strain hardening parameter
        'sg_Tm0': 600.0,          # Melting temperature (K)
        'sg_mu0': 5.6E9,          # Initial shear modulus (Pa)
    }
}

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

def get_available_materials() -> list:
    """Get list of available material names."""
    return list(MATERIAL_DATABASE.keys())

def get_material_properties(material_name: str) -> Dict[str, Any]:
    """Get properties for a specific material."""
    if material_name not in MATERIAL_DATABASE:
        available = ', '.join(get_available_materials())
        raise ValueError(f"Material '{material_name}' not found. Available materials: {available}")
    
    return MATERIAL_DATABASE[material_name].copy()

def print_material_database():
    """Print the complete material database in a formatted table."""
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
                          enable_spall: bool = True,
                          strength_model: str = 'VM') -> Dict[str, Any]:
    """Generate YAML configuration for a material based on database properties."""
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
            'str': generate_strength_model_config(props, strength_model),
            'frac': {
                'pfrac': pfrac,
                'nrhomin': nrhomin
            }
        }
    }
    
    return material_config

def generate_strength_model_config(material_props: Dict[str, Any], strength_model: str) -> Dict[str, Any]:
    """Generate strength model configuration based on material properties and model type."""
    if strength_model == 'HYDRO':
        return {'type': 'HYDRO'}
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
                                 impact_velocity: float,
                                 timing_params: Dict[str, float],
                                 enable_spall: bool = True,
                                 impactor_strength_model: str = 'VM',
                                 target_strength_model: str = 'VM',
                                 test_name: str = "Test17 Material Test") -> Dict[str, Any]:
    """Generate a complete YAML configuration file for single velocity impact."""
    # Generate material configurations
    mat1_config = generate_material_yaml(
        impactor_material, 1, impactor_thickness, impactor_cells, 
        -impactor_thickness, impact_velocity, enable_spall,
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
        'impact_velocity': impact_velocity,  # Single velocity for Test17
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
        'nmat': 2,
        # Required units section for PyKO
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
        }
    }
    
    return config

def save_yaml_config(config: Dict[str, Any], filename: str):
    """Save YAML configuration to file with comprehensive comments."""
    # Create commented YAML content
    yaml_content = []
    
    # Header comments
    yaml_content.extend([
        "# ==============================================================================",
        "# PyKO Test17 Shock Physics Simulation Configuration",
        "# ==============================================================================",
        "#",
        "# This file defines a shock physics simulation with the following setup:",
        f"# - Test Name: {config.get('name', 'Unknown')}",
        f"# - Impactor: {config.get('mat1', {}).get('eos', {}).get('name', 'Unknown')}",
        f"# - Target: {config.get('mat2', {}).get('eos', {}).get('name', 'Unknown')}",
        f"# - Impact Velocity: {config.get('impact_velocity', 0)} m/s",
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
        f"impact_velocity: {config.get('impact_velocity', 0)}  # Impact velocity (m/s)",
        ""
    ])
    
    # Timing parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# TIMING PARAMETERS",
        "# ==============================================================================",
        "# These control the simulation duration and output frequency",
        f"tstop: {float(config.get('tstop', 1e-6))}  # Total simulation time (seconds)",
        f"dtstart: {float(config.get('dtstart', 1e-10))}  # Initial time step (seconds)",
        f"dtoutput: {float(config.get('dtoutput', 1e-9))}  # Output frequency (seconds)",
        ""
    ])
    
    # Material configurations (simplified for brevity)
    for mat_id in [1, 2]:
        mat_key = f'mat{mat_id}'
        if mat_key in config:
            mat = config[mat_key]
            material_type = "IMPACTOR" if mat_id == 1 else "TARGET"
            yaml_content.extend([
                f"# ==============================================================================",
                f"# MATERIAL {mat_id}: {material_type}",
                f"# ==============================================================================",
                f"{mat_key}:",
                f"  mesh:",
                f"    cells: {int(mat['mesh']['cells'])}",
                f"    xstart: {float(mat['mesh']['xstart'])}",
                f"    length: {float(mat['mesh']['length'])}",
                f"  init:",
                f"    up0: {float(mat['init']['up0'])}",
                f"    rho0: {float(mat['init']['rho0'])}",
                f"    p0: {float(mat['init']['p0'])}",
                f"    e0: {float(mat['init']['e0'])}",
                f"    t0: {float(mat['init']['t0'])}",
                f"  eos:",
                f"    name: {mat['eos']['name']}",
                f"    type: {mat['eos']['type']}",
                f"    rhoref: {float(mat['eos']['rhoref'])}",
                f"    c0: {float(mat['eos']['c0'])}",
                f"    s1: {float(mat['eos']['s1'])}",
                f"    gamma0: {float(mat['eos']['gamma0'])}",
                f"    cv: {float(mat['eos']['cv'])}",
                f"  str:",
                f"    type: {mat['str']['type']}",
            ])
            
            # Add strength model specific parameters
            if mat['str']['type'] == 'VM':
                yaml_content.extend([
                    f"    gmod: {float(mat['str']['gmod'])}",
                    f"    ys: {float(mat['str']['ys'])}"
                ])
            elif mat['str']['type'] == 'SG':
                yaml_content.extend([
                    f"    Y0: {float(mat['str']['Y0'])}",
                    f"    Ymax: {float(mat['str']['Ymax'])}",
                    f"    beta: {float(mat['str']['beta'])}",
                    f"    n: {float(mat['str']['n'])}",
                    f"    b: {float(mat['str']['b'])}",
                    f"    h: {float(mat['str']['h'])}",
                    f"    Tm0: {float(mat['str']['Tm0'])}",
                    f"    mu0: {float(mat['str']['mu0'])}"
                ])
            
            yaml_content.extend([
                f"  frac:",
                f"    pfrac: {float(mat['frac']['pfrac'])}",
                f"    nrhomin: {float(mat['frac']['nrhomin'])}",
                ""
            ])
    
    # Boundary conditions and other parameters
    yaml_content.extend([
        "# ==============================================================================",
        "# BOUNDARY CONDITIONS",
        "# ==============================================================================",
        "boundaries:",
        f"  ibc: {config.get('boundaries', {}).get('ibc', 'FREE')}",
        f"  ip0: {float(config.get('boundaries', {}).get('ip0', 0.0))}",
        f"  obc: {config.get('boundaries', {}).get('obc', 'FREE')}",
        f"  op0: {float(config.get('boundaries', {}).get('op0', 0.0))}",
        "",
        "# ==============================================================================",
        "# GEOMETRY AND PHYSICS PARAMETERS",
        "# ==============================================================================",
        f"geometry: {config.get('geometry', 'PLA')}",
        f"gravity: {float(config.get('gravity', 0.0))}",
        f"pvoid: {float(config.get('pvoid', 0.0))}",
        f"nmat: {config.get('nmat', 2)}",
        "",
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
    
    # Write the commented YAML to file
    with open(filename, 'w') as f:
        f.write('\n'.join(yaml_content))

def create_material_config_files(impactor_material: str, target_material: str,
                                impactor_thickness: float = 0.0001,
                                target_thickness: float = 0.0002,
                                impactor_cells: int = 100,
                                target_cells: int = 200,
                                impact_velocity: float = 300.0,
                                impactor_strength_model: str = 'VM',
                                target_strength_model: str = 'VM'):
    """Create YAML configuration file for Test17 with single velocity impact."""
    # Timing parameters optimized for Test17
    timing_params = {
        'tstop': 0.50E-6,    # 0.5 μs (sufficient for spall studies)
        'dtstart': 0.0001E-6, # 0.0001 μs
        'dtoutput': 0.001E-6  # 0.001 μs
    }
    
    test_name = f"Test17 {impactor_material}-{target_material} Single Velocity"
    
    # Generate configuration
    config = generate_complete_yaml_config(
        impactor_material, target_material,
        impactor_thickness, target_thickness,
        impactor_cells, target_cells,
        impact_velocity, timing_params,
        True,  # Always enable spall/interface separation
        impactor_strength_model,
        target_strength_model,
        test_name
    )
    
    # Determine filename - save in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(script_dir, "test17-material-config.yml")
    
    # Save configuration
    save_yaml_config(config, filename)
    print(f"✅ Generated {filename}")
    print(f"   Materials: {impactor_material} → {target_material}")
    print(f"   Strength models: {STRENGTH_MODELS[impactor_strength_model]['name']} → {STRENGTH_MODELS[target_strength_model]['name']}")
    print(f"   Thicknesses: {impactor_thickness*1e6:.0f} μm → {target_thickness*1e6:.0f} μm")
    print(f"   Resolution: {impactor_cells} → {target_cells} cells")
    print(f"   Impact velocity: {impact_velocity} m/s")
    print(f"   Interface separation: Enabled")

# ===============================================================================
# MAIN FUNCTION FOR COMMAND LINE USE
# ===============================================================================

def main():
    """Interactive PyKO Material Configuration Generator for Test17"""
    print("🎯 PyKO Test17 Material Configuration Generator")
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
    
    # Show available strength models
    print("\n🔧 Available Strength Models:")
    print("-" * 30)
    strength_models = list(STRENGTH_MODELS.keys())
    for i, model in enumerate(strength_models, 1):
        model_info = STRENGTH_MODELS[model]
        print(f"{i:2d}. {model} - {model_info['name']}")
        print(f"    {model_info['description']}")
    
    # Get strength models
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
    
    # Get configuration parameters
    print("\n📏 Thickness Configuration:")
    print("(Enter values in micrometers, e.g., 100 for 100 μm)")
    
    impactor_thickness = float(input(f"Impactor thickness (μm) [default: 100]: ").strip() or "100") / 1e6
    target_thickness = float(input(f"Target thickness (μm) [default: 200]: ").strip() or "200") / 1e6
    
    print("\n🔢 Mesh Resolution:")
    impactor_cells = int(input(f"Impactor cells [default: 100]: ").strip() or "100")
    target_cells = int(input(f"Target cells [default: 200]: ").strip() or "200")
    
    print("\n⚡ Impact Velocity:")
    impact_velocity = float(input(f"Impact velocity (m/s) [default: 300]: ").strip() or "300")
    
    # Confirm configuration
    print("\n" + "=" * 50)
    print("📋 Configuration Summary")
    print("=" * 50)
    print(f"Impactor: {impactor_material} ({impactor_thickness*1e6:.0f} μm, {impactor_cells} cells)")
    print(f"  Strength model: {STRENGTH_MODELS[impactor_strength_model]['name']}")
    print(f"Target: {target_material} ({target_thickness*1e6:.0f} μm, {target_cells} cells)")
    print(f"  Strength model: {STRENGTH_MODELS[target_strength_model]['name']}")
    print(f"Impact velocity: {impact_velocity} m/s")
    print(f"Interface separation: Enabled")
    
    confirm = input("\nGenerate configuration? (y/n) [default: y]: ").strip().lower()
    if confirm in ['', 'y', 'yes']:
        try:
            create_material_config_files(
                impactor_material, target_material,
                impactor_thickness, target_thickness,
                impactor_cells, target_cells,
                impact_velocity,
                impactor_strength_model,
                target_strength_model
            )
            
            print("\n✅ Configuration generated successfully!")
            print("📁 File: test17-material-config.yml")
            print("🚀 Ready to run with PyKO")
            
        except ValueError as e:
            print(f"❌ Error: {e}")
            return 1
    else:
        print("❌ Configuration generation cancelled.")
    
    return 0

if __name__ == "__main__":
    exit(main())
