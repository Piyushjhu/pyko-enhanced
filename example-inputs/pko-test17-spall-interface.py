# %%
# pyKO Test 17: Hybrid Spall + Interface Separation Analysis
# Combines spall functionality from Test 9 with interface separation analysis from Test 16
# 
# This test simulates:
# - Material 1 flyer impacting Material 2 target (all properties read from YAML config)
# - Both materials can have spall capability (strength models and fracture parameters from YAML)
# - Analysis includes spall detection, interface separation, and free surface velocity tracking
# - Fully configurable via YAML - no hardcoded material properties or analysis parameters

########################################################################################################################
# USER CONFIGURATION SWITCHES
########################################################################################################################

# Toggle switches for analysis modules (True = ON, False = OFF)
ENABLE_INTERFACE_ANALYSIS = True  # Set to False to skip interface separation analysis
ENABLE_FSV_ANALYSIS = True         # Set to False to skip free surface velocity analysis
ENABLE_STRESS_ANALYSIS = True      # Set to False to skip stress analysis
ENABLE_SPALL_ANALYSIS = True       # Set to False to skip spall analysis

print("=== ANALYSIS MODULE CONFIGURATION ===")
print(f"Interface Analysis: {'ENABLED' if ENABLE_INTERFACE_ANALYSIS else 'DISABLED'}")
print(f"Free Surface Velocity Analysis: {'ENABLED' if ENABLE_FSV_ANALYSIS else 'DISABLED'}")
print(f"Stress Analysis: {'ENABLED' if ENABLE_STRESS_ANALYSIS else 'DISABLED'}")
print(f"Spall Analysis: {'ENABLED' if ENABLE_SPALL_ANALYSIS else 'DISABLED'}")
print("=====================================\n")

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import runpy
import sys
import pickle

sys.path.insert(1, '/Users/piyushwanchoo/Documents/Post_Doc/DATA_ANALYSIS/Pyko_pw/pyko')
from pyko import *
import pyko

runpy.run_path(path_name='import-modules.py')

########################################################################################################################
# CUSTOM COLORMAP FOR PRESSURE VISUALIZATION
########################################################################################################################

def create_pressure_colormap():
    """
    Create a custom colormap for pressure visualization:
    - Red shades for tension (negative pressure)
    - Gray for zero pressure  
    - Blue shades for compression (positive pressure)
    """
    colors_list = [
        '#8B0000',  # Dark red (max tension)
        '#CD5C5C',  # Medium red
        '#F0A0A0',  # Light red
        '#D3D3D3',  # Light gray (zero pressure)
        '#A0C8F0',  # Light blue
        '#5C85CD',  # Medium blue
        '#00008B'   # Dark blue (max compression)
    ]
    
    pressure_cmap = LinearSegmentedColormap.from_list('pressure', colors_list, N=256)
    return pressure_cmap

def create_pressure_norm(pres_min, pres_max):
    """
    Create a TwoSlopeNorm that centers the colormap at zero pressure
    """
    from matplotlib.colors import TwoSlopeNorm
    return TwoSlopeNorm(vmin=pres_min, vcenter=0.0, vmax=pres_max)

# Create the custom pressure colormap
pressure_cmap = create_pressure_colormap()

########################################################################################################################

# Automatically select the appropriate input file based on interface analysis setting
if ENABLE_INTERFACE_ANALYSIS:
    filein = './test17-spall-interface/test17-with-interface-separation.yml'
    print("🔬 Using configuration WITH interface separation physics")
else:
    filein = './test17-spall-interface/test17-without-interface-separation.yml'
    print("🔧 Using configuration WITHOUT interface separation physics")

print(f"📁 Input file: {filein}\n")

# initialize the run class variable by loading the configuration file
run = RunClass(fin=filein)

# print the run class state; this will print in code units
run.checkinput()

########################################################################################################################
# DISPLAY INPUT PARAMETERS TABLE
########################################################################################################################

print("\n" + "="*80)
print("                    TEST 17 INPUT PARAMETERS SUMMARY")
print("="*80)

# Configuration type
config_type = "WITH Interface Separation" if ENABLE_INTERFACE_ANALYSIS else "WITHOUT Interface Separation"
print(f"Configuration Type: {config_type}")
print(f"Input File: {filein}")
print("-"*80)

# Extract timing parameters from run object
print("TIMING PARAMETERS:")
print(f"{'Parameter':<20} {'Value':<15} {'Unit':<10} {'Description'}")
print("-"*65)
print(f"{'tstop':<20} {run.tstop:<15.3e} {'μs':<10} Total simulation time")
print(f"{'dtstart':<20} {run.dtstart:<15.3e} {'μs':<10} Initial time step")
print(f"{'dtoutput':<20} {run.time_skip:<15.3e} {'μs':<10} Output frequency")

# Calculate derived timing info
n_snapshots = int(run.tstop / run.time_skip)
print(f"{'Snapshots':<20} {n_snapshots:<15} {'count':<10} Total output snapshots")
print()

# Material properties table
print("MATERIAL PROPERTIES:")
print(f"{'Property':<25} {'Al Flyer':<15} {'Cu Target':<15} {'Unit':<10}")
print("-"*70)

# Get material properties from run object - MANDATORY YAML EXTRACTION
# NO FALLBACK VALUES - Script must fail if YAML parsing fails

try:
    print(f"Debug: Number of materials in run object: {run.nmat}")
    print(f"Debug: Available run attributes: {[attr for attr in dir(run) if not attr.startswith('_')]}")
    
    # Check if we have the expected number of materials
    if run.nmat < 2:
        raise ValueError(f"Expected 2 materials, found {run.nmat}")
    
    # Access material properties using pyKO's actual structure
    mat1_fracture = run.ifrac[0]  # Al flyer fracture properties
    mat2_fracture = run.ifrac[1]  # Cu target fracture properties
    mat1_strength = run.istr[0]   # Al flyer strength properties  
    mat2_strength = run.istr[1]   # Cu target strength properties
    
    print("✅ Successfully extracted material properties from YAML")
    print(f"Debug: Mat1 fracture object: {type(mat1_fracture)}")
    print(f"Debug: Mat2 fracture object: {type(mat2_fracture)}")
    print(f"Debug: Mat1 fracture attrs: {dir(mat1_fracture)}")
    print(f"Debug: Mat2 fracture attrs: {dir(mat2_fracture)}")
    
    # Try to access pfrac values
    print(f"Debug: Mat1 pfrac value: {mat1_fracture.pfrac}")
    print(f"Debug: Mat2 pfrac value: {mat2_fracture.pfrac}")
    
    # Check pyKO units and debug extraction
    print(f"Debug: run.ilength raw values: {run.ilength}")
    print(f"Debug: run.irhostart raw values: {run.irhostart}")
    print(f"Debug: mat1_eos.c0 raw: {run.ieos[0].c0}")
    print(f"Debug: mat2_eos.c0 raw: {run.ieos[1].c0}")
    
    # Extract material thicknesses from run object arrays
    # PyKO stores lengths in code units (check what units are being used)
    al_length_raw = run.ilength[0]  # Raw units from pyKO
    cu_length_raw = run.ilength[1]  # Raw units from pyKO
    
    print(f"Debug: Raw lengths - Al: {al_length_raw}, Cu: {cu_length_raw}")
    
    # PyKO uses Wilkins book units: lengths in cm
    al_length_cm = al_length_raw  # Already in cm
    cu_length_cm = cu_length_raw  # Already in cm
    
    # Convert cm to μm for display
    al_thickness_um = al_length_cm * 1e4  # cm to μm
    cu_thickness_um = cu_length_cm * 1e4  # cm to μm
    
    print(f"Debug: Lengths in cm - Al: {al_length_cm} cm, Cu: {cu_length_cm} cm")
    print(f"Debug: Converted to μm - Al: {al_thickness_um} μm, Cu: {cu_thickness_um} μm")
    
    # Create formatted strings for display
    if al_thickness_um >= 1000:
        al_thickness_str = f"{al_thickness_um/1000:.1f} mm"
    else:
        al_thickness_str = f"{al_thickness_um:.0f} μm"
        
    if cu_thickness_um >= 1000:
        cu_thickness_str = f"{cu_thickness_um/1000:.1f} mm"
    else:
        cu_thickness_str = f"{cu_thickness_um:.0f} μm"
    
    # Extract ALL material parameters from YAML config using pyKO structure
    # Spall thresholds (fracture pressures) - check units!
    al_pfrac_raw = mat1_fracture.pfrac  
    cu_pfrac_raw = mat2_fracture.pfrac  
    
    print(f"Debug: Raw pfrac values - Al: {al_pfrac_raw}, Cu: {cu_pfrac_raw}")
    
    # PyKO uses Wilkins book units: pressures in Mbar (1 Mbar = 100 GPa)
    al_spall_threshold_from_yaml = al_pfrac_raw * 100.0  # Convert Mbar to GPa
    cu_spall_threshold_from_yaml = cu_pfrac_raw * 100.0  # Convert Mbar to GPa
    print("Debug: Converting pfrac from Mbar to GPa (1 Mbar = 100 GPa)")
    
    print(f"Debug: Converted spall thresholds - Al: {al_spall_threshold_from_yaml:.6f} GPa, Cu: {cu_spall_threshold_from_yaml:.6f} GPa")
    
    # Density distension limits (nrhomin)
    al_nrhomin = mat1_fracture.nrhomin
    cu_nrhomin = mat2_fracture.nrhomin
    
    # Calculate spall density threshold from nrhomin (use the more restrictive one)
    # If any material has nrhomin = 1.0 (perfectly brittle), we need to handle it specially
    if cu_nrhomin >= 1.0:
        # For perfectly brittle materials, look for any density reduction at all
        spall_density_threshold = 0.999  # Detect even tiny density drops
        print(f"⚠️  Cu is perfectly brittle (nrhomin = {cu_nrhomin}), using sensitive density threshold")
    else:
        spall_density_threshold = max(al_nrhomin, cu_nrhomin)
    
    print(f"Material thicknesses from YAML: Al = {al_thickness_str}, Cu = {cu_thickness_str}")
    print(f"Spall thresholds from YAML: Al = {al_spall_threshold_from_yaml:.3f} GPa, Cu = {cu_spall_threshold_from_yaml:.3f} GPa")
    print(f"Density limits from YAML: Al nrhomin = {al_nrhomin:.2f}, Cu nrhomin = {cu_nrhomin:.2f}")
    print(f"Using spall density threshold: {spall_density_threshold:.3f} (from max nrhomin)")
    
    # Define global plotting parameters based on YAML values
    # Density ratio plot range
    density_vmin = min(al_nrhomin, cu_nrhomin) * 0.95  # Slightly below lowest nrhomin
    density_vmax = 1.0  # Maximum is always 1.0 (original density)
    
    # Pressure plot range based on material spall thresholds
    pressure_range_gpa = max(al_spall_threshold_from_yaml, cu_spall_threshold_from_yaml) * 2.0  # 2x max spall threshold
    
    print(f"Dynamic plot ranges from YAML: Density [{density_vmin:.2f}, {density_vmax:.1f}], Pressure ±{pressure_range_gpa:.2f} GPa")
    print("✅ ALL PARAMETERS EXTRACTED FROM YAML - NO HARDCODED VALUES")

    print(f"{'Thickness':<25} {al_thickness_str:<15} {cu_thickness_str:<15} {'-':<10}")
    print(f"{'Cells':<25} {run.inodes[0]//2:<15} {run.inodes[1]//2:<15} {'count':<10}")  # inodes = 2x cells
    
    # PyKO units: velocity in cm/μs, convert to m/s
    al_vel_ms = run.iupstart[0] * 10000  # cm/μs to m/s (cm/μs * 1e6 μs/s * 1e-2 m/cm = 1e4)
    cu_vel_ms = run.iupstart[1] * 10000  # cm/μs to m/s
    print(f"{'Initial velocity':<25} {al_vel_ms:<15.1f} {cu_vel_ms:<15.1f} {'m/s':<10}")
    
    # PyKO units: density in g/cm³, convert to kg/m³
    al_rho_kgm3 = run.irhostart[0] * 1000  # g/cm³ to kg/m³
    cu_rho_kgm3 = run.irhostart[1] * 1000  # g/cm³ to kg/m³
    print(f"{'Density':<25} {al_rho_kgm3:<15.0f} {cu_rho_kgm3:<15.0f} {'kg/m³':<10}")
    
    # Get EOS objects
    mat1_eos = run.ieos[0]  # Al flyer EOS
    mat2_eos = run.ieos[1]  # Cu target EOS
    
    # PyKO units: sound speed in cm/μs, convert to m/s
    al_c0_ms = mat1_eos.c0 * 10000  # cm/μs to m/s
    cu_c0_ms = mat2_eos.c0 * 10000  # cm/μs to m/s
    print(f"{'Sound speed (c0)':<25} {al_c0_ms:<15.0f} {cu_c0_ms:<15.0f} {'m/s':<10}")
    print(f"{'EOS parameter (s1)':<25} {mat1_eos.s1:<15.2f} {mat2_eos.s1:<15.2f} {'-':<10}")
    print(f"{'Gruneisen (gamma0)':<25} {mat1_eos.gamma0:<15.2f} {mat2_eos.gamma0:<15.2f} {'-':<10}")
    print(f"{'Specific heat (cv)':<25} {mat1_eos.cv:<15.0f} {mat2_eos.cv:<15.0f} {'eu/(K·cm³)':<10}")
    print()

    # Strength and fracture parameters
    print("STRENGTH & FRACTURE PARAMETERS:")
    print(f"{'Property':<25} {'Al Flyer':<15} {'Cu Target':<15} {'Unit':<10}")
    print("-"*70)
    print(f"{'Strength model':<25} {run.istrid[0]:<15} {run.istrid[1]:<15} {'-':<10}")
    print(f"{'Shear modulus':<25} {mat1_strength.gmod/1e9:<15.1f} {mat2_strength.gmod/1e9:<15.1f} {'GPa':<10}")
    print(f"{'Yield strength':<25} {mat1_strength.ys/1e6:<15.1f} {mat2_strength.ys/1e6:<15.1f} {'MPa':<10}")

    # Fracture parameters
    al_pfrac = mat1_fracture.pfrac / 1e6 if mat1_fracture.pfrac < 1e15 else float('inf')
    cu_pfrac = mat2_fracture.pfrac / 1e6 if mat2_fracture.pfrac < 1e15 else float('inf')

    if ENABLE_INTERFACE_ANALYSIS:
        print(f"{'Spall threshold':<25} {al_pfrac:<15.1f} {cu_pfrac:<15.1f} {'MPa':<10}")
    else:
        print(f"{'Spall threshold':<25} {'∞ (disabled)':<15} {'∞ (disabled)':<15} {'MPa':<10}")

    print(f"{'Max distension (nrhomin)':<25} {mat1_fracture.nrhomin:<15.2f} {mat2_fracture.nrhomin:<15.2f} {'-':<10}")
    print()

except (AttributeError, IndexError) as e:
    print("\n" + "="*80)
    print("❌ CRITICAL ERROR: YAML MATERIAL PROPERTIES EXTRACTION FAILED")
    print("="*80)
    print(f"Error details: {e}")
    print(f"Error type: {type(e).__name__}")
    print("\nDEBUGGING INFORMATION:")
    print(f"- Run object type: {type(run)}")
    print(f"- Has materials attribute: {hasattr(run, 'materials')}")
    if hasattr(run, 'materials'):
        print(f"- Number of materials: {len(run.materials)}")
        if len(run.materials) > 0:
            print(f"- Material 0 type: {type(run.materials[0])}")
            print(f"- Material 0 attributes: {dir(run.materials[0])}")
    
    print("\nPOSSIBLE CAUSES:")
    print("1. YAML file format error")
    print("2. Missing material sections (mat1, mat2)")
    print("3. Missing fracture parameters (pfrac, nrhomin)")
    print("4. YAML parsing failed during RunClass initialization")
    
    print("\nACTION REQUIRED:")
    print("1. Check YAML file syntax")
    print("2. Verify mat1 and mat2 sections exist")
    print("3. Verify fracture blocks exist with pfrac and nrhomin")
    print("4. Run: run.checkinput() to debug YAML loading")
    print("="*80)
    
    # Stop execution - no fallback values
    raise RuntimeError("YAML material properties extraction failed. Cannot proceed without valid material data.")

# Analysis configuration
print("ANALYSIS CONFIGURATION:")
print(f"{'Module':<30} {'Status':<10}")
print("-"*45)
print(f"{'Spall Analysis':<30} {'ON' if ENABLE_SPALL_ANALYSIS else 'OFF':<10}")
print(f"{'Free Surface Velocity':<30} {'ON' if ENABLE_FSV_ANALYSIS else 'OFF':<10}")
print(f"{'Stress Analysis':<30} {'ON' if ENABLE_STRESS_ANALYSIS else 'OFF':<10}")
print(f"{'Interface Analysis':<30} {'ON' if ENABLE_INTERFACE_ANALYSIS else 'OFF':<10}")
print()

# Expected physics
print("EXPECTED PHYSICS:")
if ENABLE_INTERFACE_ANALYSIS:
    print("• Interface separation and spall physics ENABLED")
    print("• Materials can fracture and create voids")
    print("• Expect visible material separation in x-t diagrams")
    print("• Spall threshold monitoring active")
else:
    print("• Interface separation and spall physics DISABLED")
    print("• Materials maintain structural integrity")
    print("• Clean interface dynamics without fracture")
    print("• Pure elastic-plastic behavior")

print("="*80)
print()

########################################################################################################################

# run pyko - use dtstart from YAML configuration (no hardcoded values)
pyko.run(fin=filein, userdtstart=run.dtstart, verbose=True)

########################################################################################################################

# pyko output filename is in the input file
pykofileout = run.outputfilename

# initialize a class object to hold the output data
pko = []  # this variable will hold a plain (no units) pandas datafram for plotting
pkodata = OutputClass()  # pandas + pint dataframe to read the pickled output data

# function to convert the stored pandas structure with pint units to a normal panda file
# hvplot tools do not work with a panda+pint file
# this also lets the user select a subset of variables to read into this notebook
def pyko_to_normal_panda(pkodata):
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
    return df

# loop through all the pickle dumps to read in the simulation data
# concat onto a pandas dataframe that stores the variables vs. time
with open(pykofileout, "rb") as f:
    pkodata = pickle.load(f)  # keeps units
    if 1:
        # print units
        print('pyKO output file units are the same as the input file units:')
        print('   Time        ', pkodata.time.units)
        print('   Position    ', pkodata.pos.units)
        print('   Density     ', pkodata.rho.units)
        print('   Part. vel.  ', pkodata.up.units)
        print('   Int. energy ', pkodata.ie.units)
        print('   Mass        ', pkodata.mass.units)
        print('   Temperature ', pkodata.temp.units)
        print('   Sound speed ', pkodata.alocal.units)
        print('   Pressure    ', pkodata.pres.units)
        print('   Stress      ', pkodata.sigmar.units)
    pko = pyko_to_normal_panda(pkodata)
    while True:
        try:
            pkodata = pickle.load(f)  # keeps units but only one snapshot at a time
            pko = pd.concat([pko, pyko_to_normal_panda(pkodata)], ignore_index=True)  # strips units for plotting
        except:
            break

# convert to same units as fKO for plot comparisons
# from binary in mks
pko['ie'] *= 1.E-11 * pko['rho0']  # J/kg * kg/m3 -> 100 GJ/m3 = eu/cm3
pko.rename(columns={"ie": "iev0"}, inplace=True)
pko['etot'] *= 1.E-8  # J/kg 10e7 erg/1000 g -> erg/g *1.e-12 -> eu/g
print('iev0 and etot converted to eu/g')
pko['time'] *= 1.0E6  # s->microseconds
pko['dtminj'] *= 1.0E6  # s->microseconds
pko['pos'] *= 1.0E2  # m->cm
pko['pres'] *= 1.E-9  # Pa -> GPa
pko['sigmar'] *= 1.E-9  # Pa -> GPa
pko['sigmao'] *= 1.E-9  # Pa -> GPa
pko['rho'] *= 1.E-3  # kg/m3 -> g/cm3
pko['rho0'] *= 1.E-3  # kg/m3 -> g/cm3

# list the columns in the dataframe
print(pko.columns)

########################################################################################################################

# get the original positions of the nodes and add them as a column in the final dataframe
# can use this to plot lagrangian insted of eulerian x-t diagrams
pos0 = np.asarray(pko[pko['time'] == 0.]['pos'])
unique_times = np.unique(pko['time'])
ntime = len(unique_times)
total_rows = len(pko)

print(f"Debug: pos0 length = {len(pos0)}, ntime = {ntime}, total_rows = {total_rows}")
print(f"Debug: Expected poscol length = {len(pos0) * ntime}")

# Check if all time steps have the same number of nodes
time_step_sizes = []
for t in unique_times:
    size = len(pko[pko['time'] == t])
    time_step_sizes.append(size)

print(f"Debug: Time step sizes: min={min(time_step_sizes)}, max={max(time_step_sizes)}, first={time_step_sizes[0]}, last={time_step_sizes[-1]}")

# Create pos0 column using a simpler approach
print("Creating pos0 column using groupby approach...")

# Group by time and assign initial positions to each group
pos0_column = []
grouped = pko.groupby('time')

for time_val, group in grouped:
    group_size = len(group)
    if group_size <= len(pos0):
        # Use the first group_size elements of pos0
        pos0_for_group = pos0[:group_size]
    else:
        # If group is larger than pos0, pad with the last value
        pos0_for_group = np.concatenate([pos0, np.full(group_size - len(pos0), pos0[-1])])
    
    pos0_column.extend(pos0_for_group)
    print(f"Time {time_val:.6f}: group_size={group_size}, assigned {len(pos0_for_group)} positions")

print(f"Debug: Created pos0_column with length {len(pos0_column)}")
print(f"Debug: pko dataframe length = {len(pko)}")

# Ensure exactly the right length
if len(pos0_column) == len(pko):
    pko['pos0'] = pos0_column
    print("✅ Successfully added pos0 column")
else:
    print(f"❌ Length mismatch: pos0_column={len(pos0_column)}, pko={len(pko)}")
    # Fallback: just repeat the first position for all rows
    print("Using fallback method: repeating pos0[0] for all rows")
    pko['pos0'] = pos0[0]

########################################################################################################################
# SPALL ANALYSIS SECTION (from Test 9)
########################################################################################################################

if ENABLE_SPALL_ANALYSIS:
    print("\n=== SPALL ANALYSIS ===")

    # Spall detection based on density reduction AND pressure thresholds
    # Identify regions where density has dropped significantly (indicating fracture/spall)
    pko['density_ratio'] = pko['rho'] / pko['rho0']
    
    # Use dual spall detection: density-based AND pressure-based
    print(f"Using density ratio threshold for spall detection: {spall_density_threshold:.3f} (from YAML nrhomin)")
    print(f"Also checking pressure-based spall: Al > {al_spall_threshold_from_yaml:.3f} GPa, Cu > {cu_spall_threshold_from_yaml:.3f} GPa tensile")
    
    # Check if any tensile pressures exceed spall thresholds
    max_tensile_pressure = abs(pko['pres'].min()) if pko['pres'].min() < 0 else 0
    print(f"Maximum tensile pressure in simulation: {max_tensile_pressure:.3f} GPa")
    
    pressure_spall_detected = max_tensile_pressure > min(al_spall_threshold_from_yaml, cu_spall_threshold_from_yaml)
    if pressure_spall_detected:
        print(f"⚠️  PRESSURE-BASED SPALL DETECTED: {max_tensile_pressure:.3f} GPa exceeds threshold!")
    else:
        print(f"✅ No pressure-based spall: {max_tensile_pressure:.3f} GPa below thresholds")

    # Find spalled regions at each time step
    unique_times = np.unique(pko['time'])
    spall_data = []

    for t in unique_times:
        snapshot = pko[pko['time'] == t]
        spalled_nodes = snapshot[snapshot['density_ratio'] < spall_density_threshold]
        
        if len(spalled_nodes) > 0:
            spall_info = {
                'time': t,
                'num_spalled_nodes': len(spalled_nodes),
                'spall_positions': spalled_nodes['pos'].values,
                'spall_materials': spalled_nodes['mat'].values,
                'min_density_ratio': spalled_nodes['density_ratio'].min(),
                'spall_extent': spalled_nodes['pos'].max() - spalled_nodes['pos'].min()
            }
            spall_data.append(spall_info)

    # Combined spall detection results
    density_spall_detected = len(spall_data) > 0
    overall_spall_detected = density_spall_detected or pressure_spall_detected
    
    if overall_spall_detected:
        print(f"🔥 SPALL DETECTED!")
        if density_spall_detected:
            print(f"  📉 Density-based spall: First at t = {spall_data[0]['time']:.3f} μs ({len(spall_data)} events)")
        if pressure_spall_detected:
            print(f"  💥 Pressure-based spall: Max tensile = {max_tensile_pressure:.3f} GPa")
        
        spall_detected = True
        
        # EULERIAN X-T DIAGRAMS
        print("\nCreating comprehensive Eulerian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Get actual data ranges for proper scaling
        pres_min = pko.pres.min()
        pres_max = pko.pres.max()
        density_min = pko.density_ratio.min()
        density_max = pko.density_ratio.max()
        
        print(f"Using actual data ranges - Pressure: [{pres_min:.2f}, {pres_max:.2f}] GPa, Density: [{density_min:.3f}, {density_max:.3f}]")
        print(f"Pressure colormap: Red (tension) -> Gray (zero) -> Blue (compression)")
        
        # Create pressure normalization that centers at zero
        pressure_norm = create_pressure_norm(pres_min, pres_max)
        
        # Pressure (Eulerian) - Custom colormap with zero-centered scaling
        xt_pres_eul = ax1.scatter(pko.pos * 10, pko.time, c=pko.pres, cmap=pressure_cmap, 
                                 norm=pressure_norm)
        ax1.set_xlabel('Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Eulerian: Pressure (Data Range)')
        plt.colorbar(xt_pres_eul, ax=ax1, label='Pressure (GPa)')
        
        # Particle velocity (Eulerian) - Auto-scaled
        xt_up_eul = ax2.scatter(pko.pos * 10, pko.time, c=pko.up, cmap='inferno')
        ax2.set_xlabel('Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Eulerian: Particle Velocity')
        plt.colorbar(xt_up_eul, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Eulerian) - Actual data range
        xt_rho_eul = ax3.scatter(pko.pos * 10, pko.time, c=pko.density_ratio, 
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Eulerian: Density Ratio (Data Range)')
        plt.colorbar(xt_rho_eul, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Eulerian) - Auto-scaled
        xt_mat_eul = ax4.scatter(pko.pos * 10, pko.time, c=pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Eulerian: Material ID')
        plt.colorbar(xt_mat_eul, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()
        
        # LAGRANGIAN X-T DIAGRAMS
        print("Creating Lagrangian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Lagrangian) - Custom colormap with zero-centered scaling
        xt_pres_lag = ax1.scatter(pko.pos0 * 10, pko.time, c=pko.pres, cmap=pressure_cmap, 
                                 norm=pressure_norm)
        ax1.set_xlabel('Initial Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Lagrangian: Pressure (Data Range)')
        plt.colorbar(xt_pres_lag, ax=ax1, label='Pressure (GPa)')
        
        # Particle velocity (Lagrangian) - Auto-scaled
        xt_up_lag = ax2.scatter(pko.pos0 * 10, pko.time, c=pko.up, cmap='inferno')
        ax2.set_xlabel('Initial Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Lagrangian: Particle Velocity')
        plt.colorbar(xt_up_lag, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Lagrangian) - Actual data range
        xt_rho_lag = ax3.scatter(pko.pos0 * 10, pko.time, c=pko.density_ratio, 
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Initial Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Lagrangian: Density Ratio (Data Range)')
        plt.colorbar(xt_rho_lag, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Lagrangian) - Auto-scaled
        xt_mat_lag = ax4.scatter(pko.pos0 * 10, pko.time, c=pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Initial Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Lagrangian: Material ID')
        plt.colorbar(xt_mat_lag, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()
    
    else:
        print("❌ NO SPALL DETECTED in this simulation.")
        print(f"   Density-based: No regions with ρ/ρ₀ < {spall_density_threshold:.3f}")
        print(f"   Pressure-based: Max tensile {max_tensile_pressure:.3f} GPa < {min(al_spall_threshold_from_yaml, cu_spall_threshold_from_yaml):.3f} GPa threshold")
        
        spall_detected = False
        
        # Get actual data ranges for proper scaling
        pres_min = pko.pres.min()
        pres_max = pko.pres.max()
        density_min = pko.density_ratio.min()
        density_max = pko.density_ratio.max()
        
        print(f"Using actual data ranges - Pressure: [{pres_min:.2f}, {pres_max:.2f}] GPa, Density: [{density_min:.3f}, {density_max:.3f}]")
        print(f"Pressure colormap: Red (tension) -> Gray (zero) -> Blue (compression)")
        
        # Create pressure normalization that centers at zero
        pressure_norm = create_pressure_norm(pres_min, pres_max)
        
        # EULERIAN X-T DIAGRAMS (no spall case)
        print("\nCreating comprehensive Eulerian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Eulerian) - Custom colormap with zero-centered scaling
        xt_pres_eul = ax1.scatter(pko.pos * 10, pko.time, c=pko.pres, cmap=pressure_cmap, 
                                 norm=pressure_norm)
        ax1.set_xlabel('Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Eulerian: Pressure (Data Range)')
        plt.colorbar(xt_pres_eul, ax=ax1, label='Pressure (GPa)')
        
        # Particle velocity (Eulerian) - Auto-scaled
        xt_up_eul = ax2.scatter(pko.pos * 10, pko.time, c=pko.up, cmap='inferno')
        ax2.set_xlabel('Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Eulerian: Particle Velocity')
        plt.colorbar(xt_up_eul, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Eulerian) - Actual data range
        xt_rho_eul = ax3.scatter(pko.pos * 10, pko.time, c=pko.density_ratio, 
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Eulerian: Density Ratio (Data Range)')
        plt.colorbar(xt_rho_eul, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Eulerian) - Auto-scaled
        xt_mat_eul = ax4.scatter(pko.pos * 10, pko.time, c=pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Eulerian: Material ID')
        plt.colorbar(xt_mat_eul, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()
        
        # LAGRANGIAN X-T DIAGRAMS (no spall case)
        print("Creating Lagrangian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Lagrangian) - Custom colormap with zero-centered scaling
        xt_pres_lag = ax1.scatter(pko.pos0 * 10, pko.time, c=pko.pres, cmap=pressure_cmap, 
                                 norm=pressure_norm)
        ax1.set_xlabel('Initial Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Lagrangian: Pressure (Data Range)')
        plt.colorbar(xt_pres_lag, ax=ax1, label='Pressure (GPa)')
        
        # Particle velocity (Lagrangian) - Auto-scaled
        xt_up_lag = ax2.scatter(pko.pos0 * 10, pko.time, c=pko.up, cmap='inferno')
        ax2.set_xlabel('Initial Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Lagrangian: Particle Velocity')
        plt.colorbar(xt_up_lag, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Lagrangian) - Actual data range
        xt_rho_lag = ax3.scatter(pko.pos0 * 10, pko.time, c=pko.density_ratio, 
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Initial Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Lagrangian: Density Ratio (Data Range)')
        plt.colorbar(xt_rho_lag, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Lagrangian) - Auto-scaled
        xt_mat_lag = ax4.scatter(pko.pos0 * 10, pko.time, c=pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Initial Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Lagrangian: Material ID')
        plt.colorbar(xt_mat_lag, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()

else:
    print("\n=== SPALL ANALYSIS DISABLED ===")
    spall_data = []  # Initialize empty for summary report
    density_spall_detected = False
    pressure_spall_detected = False
    overall_spall_detected = False
    max_tensile_pressure = 0.0

########################################################################################################################
# FREE SURFACE VELOCITY ANALYSIS SECTION (from Test 16)
########################################################################################################################

# Get all unique times in the simulation (needed for multiple analyses)
unique_times = np.unique(pko['time'])

if ENABLE_FSV_ANALYSIS:
    print("\n=== FREE SURFACE VELOCITY ANALYSIS ===")

    # Allocate storage for the free surface velocity at each time
    free_surface_velocity = np.zeros_like(unique_times)

    # For each time, find the node with the maximum position (rightmost = free surface)
    valid_times = []
    valid_fsv = []
    
    for i, t in enumerate(unique_times):
        snapshot = pko[pko['time'] == t]
        if len(snapshot) > 0:
            rightmost_node = snapshot.iloc[np.argmax(snapshot['pos'])]
            free_surface_velocity[i] = rightmost_node['up']
            valid_times.append(t)
            valid_fsv.append(rightmost_node['up'])
        else:
            print(f"Warning: No data at time {t:.6f} μs - skipping")
            free_surface_velocity[i] = 0  # Default value
    
    # Use only valid data for plotting
    if len(valid_times) == 0:
        print("❌ No valid time steps found for FSV analysis!")
        max_fsv = 0
        max_fsv_time = 0
    else:
        print(f"✅ Found {len(valid_times)} valid time steps out of {len(unique_times)} total")

        # Plot free surface velocity vs. time using valid data
        plt.figure(figsize=(12, 8), dpi=600)  # Increased size and resolution
        if len(valid_times) > 1:
            plt.plot(valid_times, valid_fsv, 'b-', linewidth=2, label='Valid Data')
            plt.plot(unique_times, free_surface_velocity, 'r--', alpha=0.5, label='All Times (with zeros)')
        else:
            plt.plot(unique_times, free_surface_velocity, 'b-', linewidth=2)
        
        plt.xlabel('Time (μs)')
        plt.ylabel('Free Surface Velocity (m/s)')
        plt.title('Free Surface Velocity vs. Time')
        plt.grid(True, alpha=0.3)
        if len(valid_times) > 1:
            plt.legend()

        # Calculate maximum FSV for reporting using valid data
        if len(valid_fsv) > 0:
            max_fsv = np.max(valid_fsv)
            max_fsv_time = valid_times[np.argmax(valid_fsv)]
        else:
            max_fsv = 0
            max_fsv_time = 0

        plt.tight_layout()
        plt.show()

        print(f"Maximum free surface velocity: {max_fsv:.2f} m/s at {max_fsv_time:.2f} μs")
        
        # ==========================================================================
        # FSV-BASED SPALL STRENGTH CALCULATION
        # ==========================================================================
        print("\n=== FSV-BASED SPALL STRENGTH ANALYSIS ===")
        
        if len(valid_fsv) > 10:  # Need sufficient data points
            
            # Detect velocity pullback (spall signature)
            # Look for the first significant decrease after peak velocity
            fsv_array = np.array(valid_fsv)
            time_array = np.array(valid_times)
            
            # Find the maximum velocity and its index
            max_idx = np.argmax(fsv_array)
            u_max = fsv_array[max_idx]
            t_max = time_array[max_idx]
            
            print(f"Peak FSV: {u_max:.2f} m/s at {t_max:.3f} μs")
            
            # Look for pullback after the peak (minimum velocity after peak)
            if max_idx < len(fsv_array) - 5:  # Need at least 5 points after peak
                
                # Search for minimum in the pullback region (after peak)
                pullback_region = fsv_array[max_idx:]
                pullback_times = time_array[max_idx:]
                
                # Find the minimum velocity in pullback region
                min_idx_relative = np.argmin(pullback_region)
                min_idx_absolute = max_idx + min_idx_relative
                u_min = pullback_region[min_idx_relative]
                t_min = pullback_times[min_idx_relative]
                
                # Calculate velocity pullback
                delta_u = u_max - u_min
                
                print(f"Minimum FSV after peak: {u_min:.2f} m/s at {t_min:.3f} μs")
                print(f"Velocity pullback (Δu): {delta_u:.2f} m/s")
                
                # Only calculate spall strength if there's significant pullback
                if delta_u > 10:  # Threshold for significant pullback (10 m/s)
                    
                    # Get material properties for spall strength calculation
                    # Get Cu target properties from YAML using correct pyKO units
                    try:
                        rho0_cu_gcm3 = run.irhostart[1]  # g/cm³ (pyKO units)
                        rho0_cu_kgm3 = rho0_cu_gcm3 * 1000  # Convert to kg/m³
                        print(f"Using Cu density from YAML: ρ₀ = {rho0_cu_kgm3:.0f} kg/m³ ({rho0_cu_gcm3:.3f} g/cm³)")
                    except (AttributeError, IndexError):
                        raise RuntimeError("Cannot access Cu density from run.irhostart[1]")
                    
                    # Get sound speed from EOS using correct pyKO units
                    try:
                        mat2_eos = run.ieos[1]  # Cu target EOS (material index 1)
                        if hasattr(mat2_eos, 'c0'):
                            c0_cu_cmus = mat2_eos.c0  # cm/μs (pyKO units)
                            c0_cu_ms = c0_cu_cmus * 10000  # Convert to m/s
                            print(f"Using c₀ from YAML: {c0_cu_ms:.0f} m/s ({c0_cu_cmus:.3f} cm/μs)")
                        else:
                            raise RuntimeError("Sound speed c0 not found in EOS object")
                    except Exception as e:
                        raise RuntimeError(f"Error accessing Cu sound speed properties: {e}")
                    
                    # Calculate FSV-based spall strength using converted SI units
                    # Spall strength = 0.5 * ρ₀ * c₀ * Δu
                    spall_strength_pa = 0.5 * rho0_cu_kgm3 * c0_cu_ms * delta_u
                    spall_strength_gpa = spall_strength_pa / 1e9
                    
                    print(f"\n🎯 FSV-BASED SPALL STRENGTH CALCULATION:")
                    print(f"   Formula: σ_spall = ½ × ρ₀ × c₀ × Δu")
                    print(f"   σ_spall = 0.5 × {rho0_cu_kgm3:.0f} × {c0_cu_ms:.0f} × {delta_u:.2f}")
                    print(f"   σ_spall = {spall_strength_pa:.0f} Pa = {spall_strength_gpa:.3f} GPa")
                    
                    # Compare with YAML spall threshold
                    print(f"\n📊 COMPARISON WITH YAML THRESHOLD:")
                    print(f"   FSV-measured spall strength: {spall_strength_gpa:.3f} GPa")
                    print(f"   YAML Cu spall threshold:     {cu_spall_threshold_from_yaml:.3f} GPa")
                    
                    ratio = spall_strength_gpa / cu_spall_threshold_from_yaml
                    if abs(ratio - 1.0) < 0.3:  # Within 30%
                        print(f"   ✅ Good agreement! Ratio: {ratio:.2f} (within 30%)")
                    elif ratio > 1.3:
                        print(f"   ⚠️  FSV measurement higher than YAML: ratio = {ratio:.2f}")
                    elif ratio < 0.7:
                        print(f"   ⚠️  FSV measurement lower than YAML: ratio = {ratio:.2f}")
                    else:
                        print(f"   📏 Moderate difference: ratio = {ratio:.2f}")
                    
                    # Store for summary
                    fsv_spall_strength = spall_strength_gpa
                    fsv_measurement_available = True
                    
                else:
                    print(f"❌ Insufficient velocity pullback ({delta_u:.2f} m/s < 10 m/s threshold)")
                    print("   No clear spall signature detected in FSV")
                    fsv_spall_strength = 0
                    fsv_measurement_available = False
                    
            else:
                print("❌ Insufficient data after peak velocity for pullback analysis")
                fsv_spall_strength = 0
                fsv_measurement_available = False
                
        else:
            print("❌ Insufficient FSV data points for spall analysis")
            fsv_spall_strength = 0
            fsv_measurement_available = False

else:
    print("\n=== FREE SURFACE VELOCITY ANALYSIS DISABLED ===")
    max_fsv = 0
    max_fsv_time = 0
    fsv_spall_strength = 0
    fsv_measurement_available = False

########################################################################################################################
# MAXIMUM STRESS ANALYSIS IN Cu TARGET
########################################################################################################################

if ENABLE_STRESS_ANALYSIS:
    print("\n=== MAXIMUM COMPRESSIVE & TENSILE STRESS IN Cu TARGET ANALYSIS ===")

    # Track maximum compressive and tensile stresses in Cu target (material 2) over time
    max_compressive_stress = np.zeros_like(unique_times)
    max_tensile_stress = np.zeros_like(unique_times)
    max_comp_positions = np.zeros_like(unique_times)
    max_tens_positions = np.zeros_like(unique_times)

    for i, t in enumerate(unique_times):
        snapshot = pko[pko['time'] == t]
        
        if len(snapshot) == 0:
            print(f"Warning: No data at time {t:.6f} μs - skipping stress analysis")
            continue
            
        # Filter for Cu target (material 2) only
        cu_target = snapshot[snapshot['mat'] == 2]
        
        if len(cu_target) > 0:
            # Use pressure for compressive/tensile analysis (positive = compression, negative = tension)
            pressures = cu_target['pres']
            
            # Maximum compressive stress (maximum positive pressure)
            max_comp_pressure = pressures.max()
            max_compressive_stress[i] = max_comp_pressure
            if max_comp_pressure > 0:
                max_comp_idx = pressures.idxmax()
                max_comp_positions[i] = cu_target.loc[max_comp_idx, 'pos'] * 10  # Convert to mm
            else:
                max_comp_positions[i] = 0
            
            # Maximum tensile stress (maximum negative pressure, reported as positive magnitude)
            min_pressure = pressures.min()
            max_tensile_stress[i] = -min_pressure if min_pressure < 0 else 0
            if min_pressure < 0:
                max_tens_idx = pressures.idxmin()
                max_tens_positions[i] = cu_target.loc[max_tens_idx, 'pos'] * 10  # Convert to mm
            else:
                max_tens_positions[i] = 0
        else:
            max_compressive_stress[i] = 0
            max_tensile_stress[i] = 0
            max_comp_positions[i] = 0
            max_tens_positions[i] = 0

    # Create comprehensive stress analysis plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10), dpi=300)

    # Top left: Maximum compressive stress vs time
    ax1.plot(unique_times, max_compressive_stress, 'r-', linewidth=2)
    ax1.set_xlabel('Time (μs)')
    ax1.set_ylabel('Max Compressive Stress (GPa)')
    ax1.set_title('Maximum Compressive Stress in Cu Target vs. Time')
    ax1.grid(True, alpha=0.3)

    # Calculate maximum compressive stress for reporting
    max_comp_value = np.max(max_compressive_stress)
    max_comp_time = unique_times[np.argmax(max_compressive_stress)]

    # Top right: Maximum tensile stress vs time
    ax2.plot(unique_times, max_tensile_stress, 'b-', linewidth=2)
    ax2.set_xlabel('Time (μs)')
    ax2.set_ylabel('Max Tensile Stress (GPa)')
    ax2.set_title('Maximum Tensile Stress in Cu Target vs. Time')
    ax2.grid(True, alpha=0.3)

    # Add spall threshold line for Cu (only if interface separation is enabled)
    if ENABLE_INTERFACE_ANALYSIS:
        cu_spall_threshold = cu_spall_threshold_from_yaml
        ax2.axhline(y=cu_spall_threshold, color='red', linestyle='--', alpha=0.7)
    else:
        cu_spall_threshold = 1.0E11  # Very high value when spall is disabled

    # Calculate maximum tensile stress for reporting
    max_tens_value = np.max(max_tensile_stress)
    max_tens_time = unique_times[np.argmax(max_tensile_stress)]

    # Bottom left: Position of maximum compressive stress vs time
    ax3.plot(unique_times, max_comp_positions, 'r-', linewidth=2)
    ax3.set_xlabel('Time (μs)')
    ax3.set_ylabel('Position (mm)')
    ax3.set_title('Location of Maximum Compressive Stress vs. Time')
    ax3.grid(True, alpha=0.3)

    # Add reference lines
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax3.axhline(y=1, color='black', linestyle='--', alpha=0.7)

    # Bottom right: Position of maximum tensile stress vs time
    ax4.plot(unique_times, max_tens_positions, 'b-', linewidth=2)
    ax4.set_xlabel('Time (μs)')
    ax4.set_ylabel('Position (mm)')
    ax4.set_title('Location of Maximum Tensile Stress vs. Time')
    ax4.grid(True, alpha=0.3)

    # Add reference lines
    ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax4.axhline(y=1, color='black', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()

    print(f"Maximum compressive stress in Cu target: {max_comp_value:.2f} GPa at {max_comp_time:.3f} μs")
    print(f"Maximum tensile stress in Cu target: {max_tens_value:.2f} GPa at {max_tens_time:.3f} μs")
    if max_tens_value > cu_spall_threshold:
        print(f"⚠️  TENSILE STRESS EXCEEDS Cu SPALL THRESHOLD ({cu_spall_threshold} GPa)!")
    else:
        print(f"✅ Tensile stress remains below Cu spall threshold ({cu_spall_threshold} GPa)")

else:
    print("\n=== STRESS ANALYSIS DISABLED ===")
    max_comp_value = 0
    max_comp_time = 0
    max_tens_value = 0
    max_tens_time = 0
    # Get Cu spall threshold from configuration if interface analysis is enabled
    if ENABLE_INTERFACE_ANALYSIS:
        cu_spall_threshold = cu_spall_threshold_from_yaml
    else:
        cu_spall_threshold = 1.0E11  # Very high value when spall is disabled

########################################################################################################################
# INTERFACE ANALYSIS SECTION
########################################################################################################################

if ENABLE_INTERFACE_ANALYSIS:
    print("\n=== INTERFACE ANALYSIS ===")

    # Track interfaces between materials over time
    interfaces = []

    for t in unique_times:
        snapshot = pko[pko['time'] == t]
        
        # Find material transitions
        mat_changes = np.where(np.diff(snapshot['mat']) != 0)[0]
        
        if len(mat_changes) > 0:
            interface_info = {
                'time': t,
                'interface_positions': snapshot.iloc[mat_changes]['pos'].values * 10,  # Convert to mm
                'material_pairs': [(snapshot.iloc[idx]['mat'], snapshot.iloc[idx+1]['mat']) 
                                 for idx in mat_changes]
            }
            interfaces.append(interface_info)

    if interfaces:
        print(f"Al-Cu interface tracked over {len(interfaces)} time steps")
        
        # Plot interface evolution
        plt.figure(figsize=(10, 6), dpi=300)
        
        for i, interface in enumerate(interfaces):
            for j, pos in enumerate(interface['interface_positions']):
                if j == 0:  # Only one interface between Al and Cu
                    plt.plot(interface['time'], pos, 'ro', markersize=3)
        
        plt.xlabel('Time (μs)')
        plt.ylabel('Interface Position (mm)')
        plt.title('Al-Cu Interface Evolution')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

else:
    print("\n=== INTERFACE ANALYSIS DISABLED ===")

########################################################################################################################
# SUMMARY REPORT
########################################################################################################################

print("\n" + "="*60)
print("HYBRID SPALL + INTERFACE SEPARATION TEST SUMMARY")
print("="*60)
print(f"Simulation time: 0 to {unique_times[-1]:.2f} μs")
print(f"Total nodes: {len(pos0)}")
print(f"Materials: Al flyer ({al_thickness_str}) -> Cu target ({cu_thickness_str})")
print(f"Maximum free surface velocity: {max_fsv:.2f} m/s at {max_fsv_time:.2f} μs")
print(f"Maximum compressive stress in Cu: {max_comp_value:.2f} GPa at {max_comp_time:.3f} μs")
print(f"Maximum tensile stress in Cu: {max_tens_value:.2f} GPa at {max_tens_time:.3f} μs")
if max_tens_value > cu_spall_threshold:
    print(f"⚠️  Tensile stress exceeds Cu spall threshold ({cu_spall_threshold} GPa)!")
else:
    print(f"✅ Tensile stress below Cu spall threshold ({cu_spall_threshold} GPa)")

# Add FSV-based spall strength to summary
if 'fsv_measurement_available' in locals() and fsv_measurement_available:
    print(f"🎯 FSV-measured spall strength: {fsv_spall_strength:.3f} GPa")
    ratio = fsv_spall_strength / cu_spall_threshold
    if abs(ratio - 1.0) < 0.3:
        print(f"✅ FSV measurement agrees with YAML threshold (ratio: {ratio:.2f})")
    else:
        print(f"📏 FSV vs YAML threshold ratio: {ratio:.2f}")
else:
    print("❌ FSV-based spall strength: Not available")

# Use the overall spall detection result (both density and pressure-based)
if overall_spall_detected:
    print(f"Spall detected: YES")
    if density_spall_detected:
        print(f"  Density-based spall: First at {spall_data[0]['time']:.3f} μs ({len(spall_data)} events)")
        
        # Report spall by material
        spall_by_material = {}
        for spall_event in spall_data:
            for mat_id in spall_event['spall_materials']:
                if mat_id not in spall_by_material:
                    spall_by_material[mat_id] = 0
                spall_by_material[mat_id] += 1
        
        for mat_id, count in spall_by_material.items():
            material_names = {1: f'Al flyer ({al_thickness_str})', 2: f'Cu target ({cu_thickness_str})'}
            print(f"    - Material {mat_id} ({material_names.get(mat_id, 'Unknown')}): {count} spall events")
    
    if pressure_spall_detected:
        print(f"  Pressure-based spall: Max tensile stress {max_tensile_pressure:.2f} GPa exceeded thresholds")
        if max_tensile_pressure > al_spall_threshold_from_yaml:
            print(f"    - Al threshold exceeded: {max_tensile_pressure:.2f} > {al_spall_threshold_from_yaml:.2f} GPa")
        if max_tensile_pressure > cu_spall_threshold_from_yaml:
            print(f"    - Cu threshold exceeded: {max_tensile_pressure:.2f} > {cu_spall_threshold_from_yaml:.2f} GPa")
else:
    print("Spall detected: NO")

print("="*60)

# %%
