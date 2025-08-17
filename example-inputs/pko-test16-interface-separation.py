# %%
# import modules
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import runpy
import sys

sys.path.insert(1, '/Users/piyushwanchoo/Documents/Post_Doc/DATA_ANALYSIS/Pyko_pw/pyko')
from pyko import *
import pyko

runpy.run_path(path_name='import-modules.py')

########################################################################################################################

# path to the input file
filein = './test16/test16-interface-separation.yml'

# initialize the run class variable by loading the configuration file
run = RunClass(fin=filein)

# print the run class state; this will print in code units
run.checkinput()

########################################################################################################################

# run pyko
pyko.run(fin=filein, userdtstart=0.001, verbose=True)

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
ntime = len(np.unique(pko['time']))
poscol = np.tile(pos0, ntime)
pko['pos0'] = poscol

# make eulerian x-t diagrams for the pressure and particle velocity
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,10))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5), dpi=300)
xt_pres = ax1.scatter(pko.pos * 10, pko.time, c=pko.pres, cmap='coolwarm_r', norm=colors.CenteredNorm())
# xt_pres = ax1.scatter(pko.pos*10, pko.time, c=pko.mat, cmap=cm.viridis)
ax1.set_xlabel('Position (mm)')
ax1.set_ylabel('Time (us)')
plt.colorbar(xt_pres, label='Pressure (GPa)')

xt_up = ax2.scatter(pko.pos * 10, pko.time, c=pko.up, cmap='inferno')
ax2.set_xlabel('Position (mm)')
ax2.set_ylabel('Time (us)')
plt.colorbar(xt_up, label='Particle Velocity (m/s)')

# xt_rho = ax3.scatter(pko.pos * 10, pko.time, c=pko.rho / pko.rho0, cmap=cm.coolwarm_r,
#                      norm=colors.CenteredNorm(vcenter=1.0))
# ax3.set_xlabel('Position (mm)')
# ax3.set_ylabel('Time (us)')
# plt.colorbar(xt_rho, label=r'$\rho/\rho_0$')

# xt_temp = ax4.scatter(pko.pos * 10, pko.time, c=pko.temp, cmap=cm.inferno)
# ax4.set_xlabel('Position (mm)')
# ax4.set_ylabel('Time (us)')
# plt.colorbar(xt_temp, label='Temperature (K)')

plt.tight_layout()
plt.show()


# --- Free Surface Velocity (FSV) Calculation and Plotting ---

# Get all unique times in the simulation
unique_times = np.unique(pko['time'])

# Allocate storage for the free surface velocity at each time
free_surface_velocity = np.zeros_like(unique_times)

# For each time, find the node with the maximum position (rightmost = free surface)
for i, t in enumerate(unique_times):
    snapshot = pko[pko['time'] == t]
    rightmost_node = snapshot.iloc[np.argmax(snapshot['pos'])]
    free_surface_velocity[i] = rightmost_node['up']

# Plot free surface velocity vs. time
plt.figure(dpi=300)
plt.plot(unique_times, free_surface_velocity, label='Free Surface Velocity')
plt.xlabel('Time (μs)')
plt.ylabel('Free Surface Velocity (m/s)')
plt.title('Free Surface Velocity vs. Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Optionally, print the maximum FSV and when it occurs
max_fsv = np.max(free_surface_velocity)
max_fsv_time = unique_times[np.argmax(free_surface_velocity)]
print(f"Maximum free surface velocity: {max_fsv:.2f} m/s at {max_fsv_time:.2f} μs")

# %%
