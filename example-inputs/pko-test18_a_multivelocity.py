# %%
# pyKO Test 18a: Multivelocity Analysis
# Combines spall functionality from Test 9 with interface separation analysis from Test 16
# 
# This test simulates:
# - Material 1 flyer impacting Material 2 target across a RANGE OF VELOCITIES (100-1000 m/s)
# - Both materials can have spall capability (strength models and fracture parameters from YAML)
# - Analysis includes spall detection and free surface velocity tracking
# - VELOCITY SWEEP: Automatically runs multiple simulations and plots FSV for all velocities
# - Fully configurable via YAML - no hardcoded material properties or analysis parameters

########################################################################################################################
# USER CONFIGURATION SWITCHES
########################################################################################################################

# Toggle switches for analysis modules (True = ON, False = OFF)
ENABLE_INTERFACE_ANALYSIS = True  # Interface analysis disabled - not needed for velocity sweep
ENABLE_FSV_ANALYSIS = True         # Set to True for velocity sweep FSV analysis
ENABLE_STRESS_ANALYSIS = False     # Set to False to skip stress analysis (not needed for velocity sweep)
ENABLE_SPALL_ANALYSIS = False      # Set to False to skip spall analysis (not needed for velocity sweep)

print("=== VELOCITY SWEEP ANALYSIS CONFIGURATION ===")
print(f"Interface Analysis: {'ENABLED' if ENABLE_INTERFACE_ANALYSIS else 'DISABLED'}")
print(f"Free Surface Velocity Analysis: {'ENABLED' if ENABLE_FSV_ANALYSIS else 'DISABLED'}")
print(f"Stress Analysis: {'ENABLED' if ENABLE_STRESS_ANALYSIS else 'DISABLED'}")
print(f"Spall Analysis: {'ENABLED' if ENABLE_SPALL_ANALYSIS else 'DISABLED'}")
print("============================================\n")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import runpy
import sys
import pickle
import yaml
import os

sys.path.insert(1, '/Users/piyushwanchoo/Documents/Post_Doc/DATA_ANALYSIS/Pyko_pw/pyko')
from pyko import *
import pyko

runpy.run_path(path_name='import-modules-simple.py')

########################################################################################################################
# LOAD VELOCITY SWEEP PARAMETERS
########################################################################################################################

# Load the configuration file to get velocity sweep parameters
filein = './test18_a_multivelocity/test18.yml'
print(f"📁 Input file: {filein}\n")

# Load configuration to get velocity sweep parameters
with open(filein, 'r') as file:
    config = yaml.safe_load(file)

# Extract velocity sweep parameters
velocity_sweep = config.get('velocity_sweep', {})
min_velocity = velocity_sweep.get('min_velocity', 100.0)
max_velocity = velocity_sweep.get('max_velocity', 1000.0)
velocity_steps = velocity_sweep.get('velocity_steps', 10)

# Calculate velocity array
if velocity_steps == 1:
    velocities = np.array([min_velocity])
    step_size = 0
else:
    velocities = np.linspace(min_velocity, max_velocity, velocity_steps)
    step_size = (max_velocity - min_velocity) / (velocity_steps - 1)

print("=== VELOCITY SWEEP PARAMETERS ===")
print(f"Velocity range: {min_velocity} to {max_velocity} m/s")
print(f"Number of steps: {velocity_steps}")
print(f"Step size: {step_size:.1f} m/s")
print(f"Velocities to simulate: {velocities}")
print("="*40)

########################################################################################################################
# VELOCITY SWEEP ANALYSIS
########################################################################################################################

print("\n" + "="*80)
print("                    VELOCITY SWEEP ANALYSIS")
print("="*80)

# Store FSV data for all velocities
all_fsv_data = {}
all_peak_fsv = []
all_peak_times = []

# Interface velocity tracking removed as requested

print(f"Starting velocity sweep analysis...")
print(f"Will simulate {len(velocities)} different impact velocities")
print("="*80)

# Loop through each velocity
for i, velocity in enumerate(velocities):
    # Create a new configuration for this velocity
    # Read the original YAML file fresh each time to avoid numpy object issues
    with open(filein, 'r') as file:
        current_config = yaml.safe_load(file)
    
    # Update the velocity
    current_config['mat1']['init']['up0'] = float(velocity)  # Ensure it's a Python float
    
    # Ensure timing parameters are floats
    current_config['tstop'] = float(current_config['tstop'])
    current_config['dtstart'] = float(current_config['dtstart'])
    current_config['dtoutput'] = float(current_config['dtoutput'])
    
    # Check if timing is sufficient for this velocity
    tstop = current_config['tstop']
    
    # Get minimum timing from YAML or use default, ensure it's a float
    min_tstop = float(current_config.get('min_tstop', 0.1e-06))  # Default 0.1 μs if not specified
    
    if tstop < min_tstop:
        print(f"\n❌ ERROR: Simulation time too short for velocity {velocity:.0f} m/s!")
        print(f"   Current tstop: {tstop*1e6:.3f} μs")
        print(f"   Minimum required: {min_tstop*1e6:.3f} μs")
        print(f"   Please increase tstop in YAML file to at least {min_tstop*1e6:.3f} μs")
        continue
    
    print(f"\n🔄 Simulating velocity {i+1}/{len(velocities)}: {velocity:.0f} m/s")
    print(f"    Timing: tstop={tstop*1e6:.3f} μs, dtoutput={current_config['dtoutput']*1e9:.1f} ns")
    
    # Update output filename to be unique for each velocity
    current_config['outputfilename'] = f'./test18_a_multivelocity/pyko-test18-velocity-{velocity:.0f}-bin.dat'
    
    # Write updated config to temporary file
    temp_yaml_file = f'./test18_a_multivelocity/temp_velocity_{int(velocity)}.yml'
    
    # Validate configuration before writing
    print(f"    🔍 Validating configuration for velocity {velocity:.0f} m/s...")
    print(f"    🔍 tstop: {current_config['tstop']*1e6:.3f} μs")
    print(f"    🔍 dtoutput: {current_config['dtoutput']*1e9:.1f} ns")
    print(f"    🔍 Impact velocity: {current_config['mat1']['init']['up0']:.1f} m/s")
    
    with open(temp_yaml_file, 'w') as file:
        yaml.dump(current_config, file, default_flow_style=False)
    
    print(f"    📝 Temporary YAML file created: {temp_yaml_file}")
    
    try:
        # Initialize and validate configuration
        print(f"    🔍 Initializing and validating PyKO configuration...")
        run = RunClass(fin=temp_yaml_file)
        run.checkinput()
        
        # Run pyKO with current velocity
        print(f"    🚀 Starting PyKO simulation for velocity {velocity:.0f} m/s...")
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
        print(f"    Debug: Position range: {pko['pos'].min():.3e} to {pko['pos'].max():.3e} μm")
        print(f"    Debug: Time range: {pko['time'].min():.3e} to {pko['time'].max():.3e} μs")
        
        # Calculate FSV for this velocity
        unique_times = np.unique(pko['time'])
        fsv_array = []
        time_array = []
        
        print(f"    Debug: Found {len(unique_times)} unique time steps")
        print(f"    Debug: Position range: {pko['pos'].min():.3f} to {pko['pos'].max():.3f} μm")
        print(f"    Debug: Velocity range: {pko['up'].min():.3f} to {pko['up'].max():.3f} m/s")
        print(f"    Debug: Target thickness: {current_config['mat2']['mesh']['length']*1e6:.0f} μm")
        print(f"    Debug: Sound speed: {current_config['mat2']['eos']['c0']:.0f} m/s")
        
        # Show spall parameters from YAML
        target_pfrac = current_config['mat2']['frac']['pfrac']
        target_nrhomin = current_config['mat2']['frac']['nrhomin']
        print(f"    Debug: Target spall strength (pfrac): {target_pfrac/1e6:.1f} MPa")
        print(f"    Debug: Target spall density ratio (nrhomin): {target_nrhomin:.3f}")
        
        # Calculate shock arrival time at free surface
        target_thickness = current_config['mat2']['mesh']['length']
        sound_speed = current_config['mat2']['eos']['c0']
        shock_arrival_time = target_thickness / sound_speed
        print(f"    Debug: Expected shock arrival time: {shock_arrival_time*1e6:.3f} μs")
        print(f"    Debug: Simulation time: {current_config['tstop']*1e6:.3f} μs")
        
        for t in unique_times:
            snapshot = pko[pko['time'] == t]
            if len(snapshot) > 0:
                # Find the target's free surface velocity (simple and clear)
                # Target material is Material 2, and its free surface is the rightmost node
                target_nodes = snapshot[snapshot['mat'] == 2]
                if len(target_nodes) > 0:
                    # Find the rightmost node of the target (free surface furthest from impact)
                    rightmost_target_idx = np.argmax(target_nodes['pos'])
                    fsv = target_nodes.iloc[rightmost_target_idx]['up']
                else:
                    # No target nodes found - skip this time step
                    print(f"    ⚠️  No target nodes found at t={t:.3f} μs")
                    continue
                
                # Additional debugging: Check the rightmost nodes to see what's happening (configurable)
                debug_rightmost_node_count = current_config.get('debug_rightmost_node_count', 5)  # Default 5
                debug_initial_fsv_count = current_config.get('debug_initial_fsv_count', 3)  # Default 3
                debug_velocity_limit = current_config.get('debug_velocity_limit', 300)  # Default 300 m/s
                if len(fsv_array) <= debug_initial_fsv_count or velocity <= debug_velocity_limit:
                    print(f"    Debug: Rightmost {debug_rightmost_node_count} nodes:")
                    rightmost_5 = snapshot.nlargest(debug_rightmost_node_count, 'pos')
                    for i, (idx, row) in enumerate(rightmost_5.iterrows()):
                        print(f"      Node {i+1}: pos={row['pos']:.3f} μm, vel={row['up']:.3f} m/s, rho={row['rho']:.3f} kg/m³, mat={row['mat']}")
                
                fsv_array.append(fsv)
                time_array.append(t)
                
                # Debug first few time steps and for low velocities (configurable)
                if len(fsv_array) <= debug_initial_fsv_count or velocity <= debug_velocity_limit:
                    print(f"    Debug t={t:.3f}: Target free surface pos={target_nodes.iloc[rightmost_target_idx]['pos']:.3f} μm")
                    print(f"    Debug: Target nodes: {len(target_nodes)} out of {len(snapshot)} total")
                    print(f"    Debug: Target Free Surface Velocity: {fsv:.3f} m/s")
                    print(f"    Debug: Impact velocity: {velocity:.1f} m/s")
                    print(f"    Debug: Target velocity range: {target_nodes['up'].min():.1f} to {target_nodes['up'].max():.1f} m/s")
                    print(f"    Debug: Target position range: {target_nodes['pos'].min():.3f} to {target_nodes['pos'].max():.3f} μm")
                
                # Additional debugging for the problematic time range (configurable)
                debug_time_min = current_config.get('debug_time_min', 0.04)  # Default 0.04 μs
                debug_time_max = current_config.get('debug_time_max', 0.05)  # Default 0.05 μs
                if debug_time_min <= t <= debug_time_max:
                    print(f"    🔍 CRITICAL t={t:.3f}: Target Free Surface Velocity={fsv:.3f} m/s")
                    print(f"    🔍 Material distribution: {snapshot['mat'].value_counts().to_dict()}")
                    print(f"    🔍 Target nodes density range: {target_nodes['rho'].min():.3f} to {target_nodes['rho'].max():.3f}")
                    print(f"    🔍 Target rightmost node: pos={target_nodes.iloc[rightmost_target_idx]['pos']:.3f} μm")
                    
                # Debug for low velocities to see what's happening (configurable)
                debug_velocity_threshold = current_config.get('debug_velocity_threshold', 200)  # Default 200 m/s
                debug_fsv_count_threshold = current_config.get('debug_fsv_count_threshold', 5)  # Default 5
                if velocity <= debug_velocity_threshold and len(fsv_array) <= debug_fsv_count_threshold:
                    print(f"    🔍 DEBUG t={t:.3f}: Target Free Surface Velocity={fsv:.3f} m/s")
                    print(f"    🔍 Target velocity range: {target_nodes['up'].min():.3f} to {target_nodes['up'].max():.3f} m/s")
                    print(f"    🔍 Target position range: {target_nodes['pos'].min():.3f} to {target_nodes['pos'].max():.3f} μm")
                    
                    # Check target free surface node details
                    rightmost_target = target_nodes.iloc[rightmost_target_idx]
                    print(f"    🔍 Target free surface: pos={rightmost_target['pos']:.3f} μm, vel={rightmost_target['up']:.3f} m/s, rho={rightmost_target['rho']:.3f} kg/m³")
                    print(f"    🔍 Density ratio: rho/rho0 = {rightmost_target['rho']/rightmost_target['rho0']:.3f}")
                    
                    # Check for spall damage using YAML parameters
                    target_nrhomin = current_config['mat2']['frac']['nrhomin']  # Get from YAML
                    if rightmost_target['rho']/rightmost_target['rho0'] < target_nrhomin:
                        print(f"    ⚠️  SPALL DETECTED at free surface! Density ratio = {rightmost_target['rho']/rightmost_target['rho0']:.3f} (threshold: {target_nrhomin:.3f})")
        
        if len(fsv_array) > 0:
            # Find peak FSV
            peak_idx = np.argmax(fsv_array)
            peak_fsv = fsv_array[peak_idx]
            peak_time = time_array[peak_idx]
            
            all_fsv_data[velocity] = {
                'times': time_array,
                'fsv': fsv_array,
                'peak_fsv': peak_fsv,
                'peak_time': peak_time
            }
            
            all_peak_fsv.append(peak_fsv)
            all_peak_times.append(peak_time)
            
            print(f"  ✅ Peak FSV: {peak_fsv:.1f} m/s at {peak_time:.3f} μs")
        else:
            print(f"  ⚠️  No FSV data found for velocity {velocity} m/s")
        
        # Interface velocity calculation removed as requested
            
    except Exception as e:
        print(f"  ❌ Error simulating velocity {velocity} m/s: {str(e)}")
        print(f"  🔍 Error type: {type(e).__name__}")
        print(f"  🔍 Full error details: {e}")
        continue
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_yaml_file):
            os.remove(temp_yaml_file)

print(f"\n✅ Velocity sweep completed!")
print(f"Successfully simulated {len(all_fsv_data)} velocities")
print("="*80)

########################################################################################################################
# PLOT VELOCITY SWEEP RESULTS
########################################################################################################################

if len(all_fsv_data) > 0:
    print("\n📊 Creating velocity sweep FSV plot...")
    
    # Create figure for FSV vs Time plot
    plt.figure(figsize=(12, 8))
    
    # Color map for different velocities
    colors_list = plt.cm.viridis(np.linspace(0, 1, len(all_fsv_data)))
    
    # Plot FSV curves for each velocity
    for i, (velocity, data) in enumerate(all_fsv_data.items()):
        plt.plot(data['times'], data['fsv'], 
                color=colors_list[i], 
                linewidth=2, 
                label=f'{velocity:.0f} m/s',
                alpha=0.8)
    
    # Add proper labels and title
    plt.xlabel('Time (μs)', fontsize=12)
    plt.ylabel('Free Surface Velocity (m/s)', fontsize=12)
    plt.title('Free Surface Velocity vs Time for Different Impact Velocities', fontsize=14, fontweight='bold')
    plt.legend(title='Impact Velocity', fontsize=10)
    plt.grid(True, alpha=0.3)  # Grid transparency
    plt.tight_layout()
    
    # Save the first plot
    plt.savefig('./test18_a_multivelocity/fsv_vs_time.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create second figure for Peak FSV vs Impact Velocity
    plt.figure(figsize=(10, 6))
    plt.plot(velocities[:len(all_peak_fsv)], all_peak_fsv, 'ro-', linewidth=2, markersize=8)
    plt.xlabel('Impact Velocity (m/s)', fontsize=12)
    plt.ylabel('Peak Free Surface Velocity (m/s)', fontsize=12)
    plt.title('Peak FSV vs Impact Velocity', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save the second plot
    plt.savefig('./test18_a_multivelocity/peak_fsv_vs_impact_velocity.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Interface velocity plotting removed as requested

    print("✅ Two plots generated:")
    print("   - fsv_vs_time.png: FSV vs Time for all impact velocities")
    print("   - peak_fsv_vs_impact_velocity.png: Peak FSV vs Impact Velocity correlation")
    
    # Summary statistics
    print(f"\n📈 VELOCITY SWEEP SUMMARY:")
    print(f"Velocity range: {min_velocity:.0f} - {max_velocity:.0f} m/s")
    print(f"Peak FSV range: {min(all_peak_fsv):.1f} - {max(all_peak_fsv):.1f} m/s")
    print(f"FSV amplification factor: {max(all_peak_fsv)/max_velocity:.2f}x")

else:
    print("❌ No FSV data available for plotting")

print("\n" + "="*60)
print("MULTIVELOCITY ANALYSIS TEST 18a SUMMARY")
print("="*60)
print(f"Velocity sweep completed: {len(all_fsv_data)}/{len(velocities)} simulations successful")
print(f"Velocity range: {min_velocity:.0f} - {max_velocity:.0f} m/s")
if len(all_peak_fsv) > 0:
    print(f"Peak FSV range: {min(all_peak_fsv):.1f} - {max(all_peak_fsv):.1f} m/s")
print("="*60)

# %%
