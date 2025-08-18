#!/usr/bin/env python3
"""
pyKO Test 18b: Multi-Thickness Analysis
========================================

This script performs a thickness sweep analysis to study the effect of target thickness
on free surface velocity (FSV) response at constant impact velocity.

Key Features:
- Constant impact velocity (user-defined)
- Variable target thickness sweep
- FSV analysis for each thickness
- Comparative plotting of FSV vs time
- Peak FSV vs thickness correlation

Author: pyKO Development Team
Date: 2024
"""
# %%
import os
import sys
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import runpy
from typing import Dict, List, Any

# Add pyKO to path and import
sys.path.append('..')
from pyko import *
import pyko

########################################################################################################################
# LOAD MODULES AND SETUP
########################################################################################################################

# Load the modules setup
runpy.run_path(path_name='import-modules-simple.py')

########################################################################################################################
# LOAD THICKNESS SWEEP PARAMETERS
########################################################################################################################

# Load the configuration file to get thickness sweep parameters
filein = './test18_b_multithickness/test18b.yml'

print(f"📁 Input file: {filein}\n")

# Load configuration to get thickness sweep parameters
with open(filein, 'r') as file:
    config = yaml.safe_load(file)

# Extract thickness sweep parameters
thickness_sweep = config.get('thickness_sweep', {})
min_thickness = float(thickness_sweep.get('min_thickness', 50.0))  # μm
max_thickness = float(thickness_sweep.get('max_thickness', 500.0))  # μm
thickness_steps = int(thickness_sweep.get('thickness_steps', 5))

# Calculate thickness values
thicknesses = np.linspace(min_thickness, max_thickness, thickness_steps) / 1e6  # Convert to meters

# Get constant impact velocity
constant_velocity = float(config['mat1']['init']['up0'])

print("=== THICKNESS SWEEP PARAMETERS ===")
print(f"Constant impact velocity: {constant_velocity:.0f} m/s")
print(f"Thickness range: {min_thickness:.0f} to {max_thickness:.0f} μm")
print(f"Number of steps: {thickness_steps}")
print(f"Thicknesses to simulate: {thicknesses*1e6} μm")
print("=" * 40)

########################################################################################################################
# THICKNESS SWEEP ANALYSIS
########################################################################################################################

print("\n" + "="*80)
print("                    THICKNESS SWEEP ANALYSIS")
print("="*80)

# Store FSV data for all thicknesses
all_fsv_data = {}
all_peak_fsv = []
all_peak_times = []

print(f"Starting thickness sweep analysis...")
print(f"Will simulate {len(thicknesses)} different target thicknesses")
print("="*80)

# Loop through each thickness
for i, thickness in enumerate(thicknesses):
    # Create a new configuration for this thickness
    # Read the original YAML file fresh each time to avoid numpy object issues
    with open(filein, 'r') as file:
        current_config = yaml.safe_load(file)
    
    # Update the target thickness
    current_config['mat2']['mesh']['length'] = float(thickness)
    current_config['mat2']['mesh']['cells'] = int(thickness * 1e6)  # 1 cell per μm
    
    # Ensure timing parameters are floats
    current_config['tstop'] = float(current_config['tstop'])
    current_config['dtstart'] = float(current_config['dtstart'])
    current_config['dtoutput'] = float(current_config['dtoutput'])
    
    # Check if timing is sufficient for this thickness
    tstop = current_config['tstop']
    
    # Get minimum timing from YAML or use default, ensure it's a float
    min_tstop = float(current_config.get('min_tstop', 0.1e-06))  # Default 0.1 μs if not specified
    
    if tstop < min_tstop:
        print(f"\n❌ ERROR: Simulation time too short for thickness {thickness*1e6:.0f} μm!")
        print(f"   Current tstop: {tstop*1e6:.3f} μs")
        print(f"   Minimum required: {min_tstop*1e6:.3f} μs")
        print(f"   Please increase tstop in YAML file to at least {min_tstop*1e6:.3f} μs")
        continue
    
    print(f"\n🔄 Simulating thickness {i+1}/{len(thicknesses)}: {thickness*1e6:.0f} μm")
    print(f"    Timing: tstop={tstop*1e6:.3f} μs, dtoutput={current_config['dtoutput']*1e9:.1f} ns")
    
    # Update output filename to be unique for each thickness
    current_config['outputfilename'] = f'./test18_b_multithickness/pyko-test18b-thickness-{thickness*1e6:.0f}-bin.dat'
    
    # Write updated config to temporary file
    temp_yaml_file = f'./test18_b_multithickness/temp_thickness_{int(thickness*1e6)}.yml'
    
    # Validate configuration before writing
    print(f"    🔍 Validating configuration for thickness {thickness*1e6:.0f} μm...")
    print(f"    🔍 tstop: {current_config['tstop']*1e6:.3f} μs")
    print(f"    🔍 dtoutput: {current_config['dtoutput']*1e9:.1f} ns")
    print(f"    🔍 Target thickness: {current_config['mat2']['mesh']['length']*1e6:.1f} μm")
    print(f"    🔍 Target cells: {current_config['mat2']['mesh']['cells']}")
    print(f"    🔍 Impact velocity: {current_config['mat1']['init']['up0']:.1f} m/s")
    
    with open(temp_yaml_file, 'w') as file:
        yaml.dump(current_config, file, default_flow_style=False)
    
    print(f"    📝 Temporary YAML file created: {temp_yaml_file}")
    
    try:
        # Initialize and validate configuration
        print(f"    🔍 Initializing and validating PyKO configuration...")
        run = RunClass(fin=temp_yaml_file)
        run.checkinput()
        
        # Run pyKO with current thickness
        print(f"    🚀 Starting PyKO simulation for thickness {thickness*1e6:.0f} μm...")
        pyko.run(fin=temp_yaml_file, userdtstart=run.dtstart, verbose=False)
        print(f"    ✅ PyKO simulation completed successfully")
        
        # Load the output data - load ALL snapshots
        output_file = current_config['outputfilename']
        
        # Check if output file was created
        if not os.path.exists(output_file):
            print(f"    ❌ Output file not found: {output_file}")
            print(f"    🔍 Checking if PyKO simulation actually ran...")
            continue
            
        print(f"    📁 Output file found: {output_file}")
        pko_list = []
        
        with open(output_file, "rb") as f:
            while True:
                try:
                    pkodata = pickle.load(f)
                    # Convert to normal pandas DataFrame
                    df = pd.DataFrame({
                        "j": pkodata.j.magnitude,
                        "stepn": pkodata.stepn.magnitude,
                        "time": pkodata.time.magnitude,
                        "mat": pkodata.mat.magnitude,
                        "pos": pkodata.pos.magnitude,
                        "rho0": pkodata.rho0.magnitude,
                        "rho": pkodata.rho.magnitude,
                        "up": pkodata.up.magnitude,
                        "ie": pkodata.ie.magnitude,
                        "pres": pkodata.pres.magnitude,
                        "mass": pkodata.mass.magnitude,
                        "temp": pkodata.temp.magnitude,
                        "cs": pkodata.alocal.magnitude,
                        "sigmar": pkodata.sigmar.magnitude,
                        "sigmao": pkodata.sigmao.magnitude,
                        "etot": pkodata.etot.magnitude,
                        "dtminj": pkodata.dtminj.magnitude,
                    })
                    pko_list.append(df)
                except EOFError:
                    break
        
        # Combine all snapshots
        if pko_list:
            pko = pd.concat(pko_list, ignore_index=True)
            print(f"    Debug: Loaded {len(pko_list)} snapshots")
        else:
            print(f"    ❌ No data loaded from {output_file}")
            continue
        
        # Convert units based on pyKO documentation:
        # - Input units: User-specified (m/s in YAML)
        # - Code units: Internal (cm/μs for velocity)
        # - Output units: Same as input units (m/s)
        # So the output data should already be in m/s!
        
        pko['pos'] = pko['pos'] * 1e4  # m to μm
        pko['pres'] = pko['pres'] / 1e9  # Pa to GPa
        pko['time'] = pko['time'] * 1e6  # s to μs (microseconds)
        
        # Velocity should already be in m/s (input units)
        # Note: pko DataFrame already contains the velocity data from all snapshots
        
        print(f"    Debug: Raw velocity range: {pko['up'].min():.3f} to {pko['up'].max():.3f} m/s")
        
        # Get target material properties for analysis
        target_thickness = current_config['mat2']['mesh']['length']
        target_sound_speed = current_config['mat2']['eos']['c0']
        target_spall_strength = current_config['mat2']['frac']['pfrac'] / 1e6  # Convert to MPa
        target_spall_density_ratio = current_config['mat2']['frac']['nrhomin']
        
        print(f"    Debug: Position range: {pko['pos'].min():.3e} to {pko['pos'].max():.3e} μm")
        print(f"    Debug: Time range: {pko['time'].min():.3e} to {pko['time'].max():.3e} μs")
        print(f"    Debug: Found {len(pko_list)} unique time steps")
        
        # Convert position to μm for easier analysis
        pko['pos'] = pko['pos'] * 1e6  # Convert to μm
        print(f"    Debug: Position range: {pko['pos'].min():.3f} to {pko['pos'].max():.3f} μm")
        print(f"    Debug: Velocity range: {pko['up'].min():.3f} to {pko['up'].max():.3f} m/s")
        print(f"    Debug: Target thickness: {target_thickness*1e6:.0f} μm")
        print(f"    Debug: Sound speed: {target_sound_speed:.0f} m/s")
        print(f"    Debug: Target spall strength (pfrac): {target_spall_strength:.1f} MPa")
        print(f"    Debug: Target spall density ratio (nrhomin): {target_spall_density_ratio:.3f}")
        
        # Calculate expected shock arrival time
        expected_shock_arrival = target_thickness / target_sound_speed * 1e6  # Convert to μs
        print(f"    Debug: Expected shock arrival time: {expected_shock_arrival:.3f} μs")
        print(f"    Debug: Simulation time: {current_config['tstop']*1e6:.3f} μs")
        
        # Find the rightmost nodes (target free surface)
        target_material_id = 2  # Target material ID
        target_nodes = pko[pko['mat'] == target_material_id]
        
        if len(target_nodes) == 0:
            print(f"    ❌ No target material nodes found!")
            continue
        
        # Find the rightmost node (target free surface)
        rightmost_pos = target_nodes['pos'].max()
        target_free_surface_nodes = target_nodes[target_nodes['pos'] == rightmost_pos]
        
        print(f"    Debug: Rightmost 5 nodes:")
        rightmost_nodes = target_nodes.nlargest(5, 'pos')
        for idx, row in rightmost_nodes.iterrows():
            print(f"      Node {row['j']}: pos={row['pos']:.3f} μm, vel={row['up']:.3f} m/s, rho={row['rho']:.3f} kg/m³, mat={row['mat']:.1f}")
        
        print(f"    Debug t=0.000: Target free surface pos={rightmost_pos:.3f} μm")
        print(f"    Debug: Target nodes: {len(target_nodes)} out of {len(pko)} total")
        
        # Calculate FSV for each time step
        fsv_data = []
        unique_times = sorted(pko['time'].unique())
        
        for t in unique_times:
            time_data = pko[pko['time'] == t]
            target_time_data = time_data[time_data['mat'] == target_material_id]
            
            if len(target_time_data) == 0:
                continue
            
            # Find the rightmost node (target free surface)
            rightmost_node = target_time_data.loc[target_time_data['pos'].idxmax()]
            target_fsv = rightmost_node['up']
            
            fsv_data.append({
                'time': t,
                'fsv': target_fsv,
                'position': rightmost_node['pos']
            })
            
            # Debug output for critical times
            if t < 0.001 or t > 0.040:
                print(f"    Debug t={t:.3f}: Target Free Surface Velocity={target_fsv:.3f} m/s")
                print(f"    🔍 Material distribution: {dict(time_data['mat'].value_counts())}")
                print(f"    🔍 Target nodes density range: {target_time_data['rho'].min():.3f} to {target_time_data['rho'].max():.3f}")
                print(f"    🔍 Target rightmost node: pos={rightmost_node['pos']:.3f} μm")
        
        # Convert to DataFrame
        fsv_df = pd.DataFrame(fsv_data)
        
        if len(fsv_df) == 0:
            print(f"    ❌ No FSV data calculated!")
            continue
        
        # Find peak FSV
        peak_fsv = fsv_df['fsv'].max()
        peak_time = fsv_df.loc[fsv_df['fsv'].idxmax(), 'time']
        
        print(f"  ✅ Peak FSV: {peak_fsv:.1f} m/s at {peak_time:.3f} μs")
        
        # Store data for plotting
        thickness_key = f"{thickness*1e6:.0f}μm"
        all_fsv_data[thickness_key] = fsv_df
        all_peak_fsv.append(peak_fsv)
        all_peak_times.append(peak_time)
        
        # Clean up temporary file
        os.remove(temp_yaml_file)
        
    except Exception as e:
        print(f"  ❌ Error simulating thickness {thickness*1e6:.0f} μm: {str(e)}")
        print(f"  🔍 Error type: {type(e).__name__}")
        print(f"  🔍 Full error details: {str(e)}")
        continue

print(f"\n✅ Thickness sweep completed!")
print(f"Successfully simulated {len(all_fsv_data)} thicknesses")
print("="*80)

########################################################################################################################
# PLOTTING AND ANALYSIS
########################################################################################################################

if len(all_fsv_data) > 0:
    print("\n📊 Creating thickness sweep FSV plots...")
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Color map for different thicknesses
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_fsv_data)))
    
    # Plot 1: FSV vs Time for all thicknesses
    ax1.set_title('Free Surface Velocity vs Time for Different Target Thicknesses', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Time (μs)', fontsize=12)
    ax1.set_ylabel('Free Surface Velocity (m/s)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    
    for i, (thickness_key, fsv_df) in enumerate(all_fsv_data.items()):
        thickness_value = float(thickness_key.replace('μm', ''))
        ax1.plot(fsv_df['time'], fsv_df['fsv'], 
                color=colors[i], linewidth=2, 
                label=f'{thickness_value:.0f} μm Target')
    
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.set_xlim(0, max([df['time'].max() for df in all_fsv_data.values()]))
    
    # Plot 2: Peak FSV vs Target Thickness
    thickness_values = [float(key.replace('μm', '')) for key in all_fsv_data.keys()]
    
    ax2.set_title('Peak Free Surface Velocity vs Target Thickness', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Target Thickness (μm)', fontsize=12)
    ax2.set_ylabel('Peak FSV (m/s)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    # Plot data points
    ax2.scatter(thickness_values, all_peak_fsv, color='red', s=100, zorder=5)
    
    # Add trend line if we have enough points
    if len(thickness_values) > 2:
        z = np.polyfit(thickness_values, all_peak_fsv, 1)
        p = np.poly1d(z)
        ax2.plot(thickness_values, p(thickness_values), "r--", alpha=0.8, linewidth=2)
        
        # Calculate correlation
        correlation = np.corrcoef(thickness_values, all_peak_fsv)[0, 1]
        ax2.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=ax2.transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add data labels
    for i, (thickness, peak_fsv) in enumerate(zip(thickness_values, all_peak_fsv)):
        ax2.annotate(f'{peak_fsv:.1f}', 
                    (thickness, peak_fsv), 
                    textcoords="offset points", 
                    xytext=(0,10), 
                    ha='center', fontsize=9)
    
    plt.tight_layout()
    
    # Save plots
    plot_filename = './test18_b_multithickness/thickness_sweep_fsv.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"✅ Saved combined plot: {plot_filename}")
    
    # Create individual plots
    # Plot 1: FSV vs Time
    plt.figure(figsize=(10, 6))
    plt.title('Free Surface Velocity vs Time for Different Target Thicknesses', fontsize=14, fontweight='bold')
    plt.xlabel('Time (μs)', fontsize=12)
    plt.ylabel('Free Surface Velocity (m/s)', fontsize=12)
    plt.grid(True, alpha=0.3)
    
    for i, (thickness_key, fsv_df) in enumerate(all_fsv_data.items()):
        thickness_value = float(thickness_key.replace('μm', ''))
        plt.plot(fsv_df['time'], fsv_df['fsv'], 
                color=colors[i], linewidth=2, 
                label=f'{thickness_value:.0f} μm Target')
    
    plt.legend()
    plt.xlim(0, max([df['time'].max() for df in all_fsv_data.values()]))
    plt.tight_layout()
    plt.savefig('./test18_b_multithickness/fsv_vs_time.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved FSV vs Time plot: fsv_vs_time.png")
    
    # Plot 2: Peak FSV vs Thickness
    plt.figure(figsize=(10, 6))
    plt.title('Peak Free Surface Velocity vs Target Thickness', fontsize=14, fontweight='bold')
    plt.xlabel('Target Thickness (μm)', fontsize=12)
    plt.ylabel('Peak FSV (m/s)', fontsize=12)
    plt.grid(True, alpha=0.3)
    
    plt.scatter(thickness_values, all_peak_fsv, color='red', s=100, zorder=5)
    
    if len(thickness_values) > 2:
        z = np.polyfit(thickness_values, all_peak_fsv, 1)
        p = np.poly1d(z)
        plt.plot(thickness_values, p(thickness_values), "r--", alpha=0.8, linewidth=2)
        
        correlation = np.corrcoef(thickness_values, all_peak_fsv)[0, 1]
        plt.text(0.05, 0.95, f'Correlation: {correlation:.3f}', 
                transform=plt.gca().transAxes, fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    for i, (thickness, peak_fsv) in enumerate(zip(thickness_values, all_peak_fsv)):
        plt.annotate(f'{peak_fsv:.1f}', 
                    (thickness, peak_fsv), 
                    textcoords="offset points", 
                    xytext=(0,10), 
                    ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('./test18_b_multithickness/peak_fsv_vs_thickness.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved Peak FSV vs Thickness plot: peak_fsv_vs_thickness.png")
    
    # Print summary statistics
    print(f"\n📈 THICKNESS SWEEP SUMMARY:")
    print(f"Thickness range: {min_thickness:.0f} - {max_thickness:.0f} μm")
    print(f"Peak FSV range: {min(all_peak_fsv):.1f} - {max(all_peak_fsv):.1f} m/s")
    if len(thickness_values) > 2:
        correlation = np.corrcoef(thickness_values, all_peak_fsv)[0, 1]
        print(f"FSV-Thickness correlation: {correlation:.3f}")
    
    print(f"\n" + "="*80)
    print("MULTITHICKNESS ANALYSIS TEST 18b SUMMARY")
    print("="*80)
    print(f"Thickness sweep completed: {len(all_fsv_data)}/{len(thicknesses)} simulations successful")
    print(f"Constant impact velocity: {constant_velocity:.0f} m/s")
    print(f"Thickness range: {min_thickness:.0f} - {max_thickness:.0f} μm")
    print("="*80)

else:
    print("\n❌ No successful simulations to plot!")
    print("Check the configuration and try again.")

print("\n🎯 Test 18b Multi-Thickness Analysis Complete!")

# %%
