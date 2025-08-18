#!/usr/bin/env python3
"""
Material Database for pyKO Test 18b: Multi-Thickness Analysis
============================================================

This module contains material properties and strength model definitions
specifically for Test 18b multi-thickness analysis.

Author: pyKO Development Team
Date: 2024
"""

from typing import Dict, List, Any

# Available strength models for Test 18b
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
    }
}

# Material database with properties for Test 18b
MATERIAL_DATABASE = {
    'Aluminum': {
        # Basic properties
        'density': 2700.0,           # kg/m³
        'sound_speed': 5200.0,       # m/s
        'hugoniot_slope': 1.5,       # dimensionless
        'gruneisen': 2.0,            # dimensionless
        'specific_heat': 896.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 26.0e9,     # Pa
        'yield_strength': 207.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 276.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.29e9, 'sg_Ymax': 0.68e9, 'sg_beta': 125.0, 'sg_n': 0.1,
        'sg_b': 8.0, 'sg_h': 6.2e-4, 'sg_Tm0': 1220.0, 'sg_mu0': 27.6e9
    },
    
    'Copper': {
        # Basic properties
        'density': 8930.0,           # kg/m³
        'sound_speed': 3940.0,       # m/s
        'hugoniot_slope': 1.489,     # dimensionless
        'gruneisen': 1.96,           # dimensionless
        'specific_heat': 385.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 48.0e9,     # Pa
        'yield_strength': 120.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 200.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.12e9, 'sg_Ymax': 0.64e9, 'sg_beta': 36.0, 'sg_n': 0.45,
        'sg_b': 3.0, 'sg_h': 0.0, 'sg_Tm0': 1356.0, 'sg_mu0': 47.7e9
    },
    
    'Steel': {
        # Basic properties
        'density': 7850.0,           # kg/m³
        'sound_speed': 4570.0,       # m/s
        'hugoniot_slope': 1.49,      # dimensionless
        'gruneisen': 1.69,           # dimensionless
        'specific_heat': 460.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 80.0e9,     # Pa
        'yield_strength': 250.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 400.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.25e9, 'sg_Ymax': 0.65e9, 'sg_beta': 40.0, 'sg_n': 0.26,
        'sg_b': 3.0, 'sg_h': 0.0, 'sg_Tm0': 1811.0, 'sg_mu0': 81.8e9
    },
    
    'Iron': {
        # Basic properties
        'density': 7870.0,           # kg/m³
        'sound_speed': 4570.0,       # m/s
        'hugoniot_slope': 1.92,      # dimensionless
        'gruneisen': 1.69,           # dimensionless
        'specific_heat': 460.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 82.0e9,     # Pa
        'yield_strength': 200.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 300.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.20e9, 'sg_Ymax': 0.60e9, 'sg_beta': 35.0, 'sg_n': 0.25,
        'sg_b': 3.0, 'sg_h': 0.0, 'sg_Tm0': 1811.0, 'sg_mu0': 82.0e9
    },
    
    'Titanium': {
        # Basic properties
        'density': 4500.0,           # kg/m³
        'sound_speed': 5130.0,       # m/s
        'hugoniot_slope': 1.14,      # dimensionless
        'gruneisen': 1.1,            # dimensionless
        'specific_heat': 523.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 44.0e9,     # Pa
        'yield_strength': 830.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 500.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.83e9, 'sg_Ymax': 1.2e9, 'sg_beta': 80.0, 'sg_n': 0.1,
        'sg_b': 5.0, 'sg_h': 0.0, 'sg_Tm0': 1941.0, 'sg_mu0': 44.0e9
    },
    
    'Tungsten': {
        # Basic properties
        'density': 19300.0,          # kg/m³
        'sound_speed': 4030.0,       # m/s
        'hugoniot_slope': 1.237,     # dimensionless
        'gruneisen': 1.62,           # dimensionless
        'specific_heat': 134.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 160.0e9,    # Pa
        'yield_strength': 550.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 800.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.55e9, 'sg_Ymax': 1.5e9, 'sg_beta': 100.0, 'sg_n': 0.1,
        'sg_b': 5.0, 'sg_h': 0.0, 'sg_Tm0': 3695.0, 'sg_mu0': 160.0e9
    },
    
    'Gold': {
        # Basic properties
        'density': 19320.0,          # kg/m³
        'sound_speed': 3240.0,       # m/s
        'hugoniot_slope': 1.572,     # dimensionless
        'gruneisen': 2.97,           # dimensionless
        'specific_heat': 129.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 27.0e9,     # Pa
        'yield_strength': 100.0e6,   # Pa
        
        # Fracture properties
        'spall_strength': 150.0e6,   # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.10e9, 'sg_Ymax': 0.30e9, 'sg_beta': 20.0, 'sg_n': 0.1,
        'sg_b': 2.0, 'sg_h': 0.0, 'sg_Tm0': 1337.0, 'sg_mu0': 27.0e9
    },
    
    'Lead': {
        # Basic properties
        'density': 11340.0,          # kg/m³
        'sound_speed': 1960.0,       # m/s
        'hugoniot_slope': 1.46,      # dimensionless
        'gruneisen': 2.73,           # dimensionless
        'specific_heat': 129.0,      # J/kg/K
        
        # Strength properties (Von Mises)
        'shear_modulus': 5.6e9,      # Pa
        'yield_strength': 10.0e6,    # Pa
        
        # Fracture properties
        'spall_strength': 20.0e6,    # Pa
        
        # Steinberg-Guinan parameters
        'sg_Y0': 0.01e9, 'sg_Ymax': 0.05e9, 'sg_beta': 10.0, 'sg_n': 0.1,
        'sg_b': 1.0, 'sg_h': 0.0, 'sg_Tm0': 600.0, 'sg_mu0': 5.6e9
    }
}

def get_available_materials() -> List[str]:
    """
    Get list of available material names.
    
    Returns:
        List[str]: List of material names
    """
    return list(MATERIAL_DATABASE.keys())

def get_material_properties(material_name: str) -> Dict[str, Any]:
    """
    Get properties for a specific material.
    
    Args:
        material_name (str): Name of the material
        
    Returns:
        Dict[str, Any]: Material properties dictionary
        
    Raises:
        ValueError: If material not found
    """
    if material_name not in MATERIAL_DATABASE:
        raise ValueError(f"Material '{material_name}' not found. Available: {list(MATERIAL_DATABASE.keys())}")
    
    return MATERIAL_DATABASE[material_name]

def get_available_strength_models() -> List[str]:
    """
    Get list of available strength model types.
    
    Returns:
        List[str]: List of strength model types
    """
    return list(STRENGTH_MODELS.keys())

def get_strength_model_info(strength_model: str) -> Dict[str, Any]:
    """
    Get information about a specific strength model.
    
    Args:
        strength_model (str): Strength model type
        
    Returns:
        Dict[str, Any]: Strength model information
        
    Raises:
        ValueError: If strength model not found
    """
    if strength_model not in STRENGTH_MODELS:
        raise ValueError(f"Strength model '{strength_model}' not found. Available: {list(STRENGTH_MODELS.keys())}")
    
    return STRENGTH_MODELS[strength_model]
