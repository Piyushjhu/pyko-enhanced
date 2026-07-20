# %%
# Animated Lagrangian x-t characteristic diagram with synchronized FSV plot.
# Standalone — runs pyKO from test17-material-config.yml and needs no other script.
#
# Left panel : Lagrangian x-t diagram.
#              Blue lines  = compressive wave fronts  (pressure > 0 GPa)
#              Red  lines  = rarefaction / tensile     (pressure < 0 GPa)
#              No filled colour — contour lines only.
# Right panel: Free-surface velocity (FSV) vs time.  Points appear as each
#              time step is revealed by the animation sweep.

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT   = os.path.abspath(os.path.join(_SCRIPT_DIR, '..'))
sys.path.insert(1, _REPO_ROOT)

from pyko import *
import pyko

# ── User-tuneable settings ─────────────────────────────────────────────────────
CONFIG = os.path.join(_SCRIPT_DIR, 'test17-spall-interface', 'test17-material-config.yml')

N_COMP_LEVELS = 14    # contour lines on the compressive side (blue)
N_TENS_LEVELS = 12    # contour lines on the tensile/rarefaction side (red)
STRESS_CLIP   = 0.01  # GPa  — levels below this magnitude are suppressed (noise guard)
FRAME_SKIP    = 3     # reveal every N-th time step (raise to reduce frame count / speed up)
INTERVAL_MS   = 80    # milliseconds between frames (only affects interactive backends)

# ── Run simulation ─────────────────────────────────────────────────────────────
print("Loading configuration …")
run = RunClass(fin=CONFIG)
run.checkinput()

print("Running pyKO …")
pyko.run(fin=CONFIG, userdtstart=run.dtstart, verbose=False)

# ── Load pickle output ─────────────────────────────────────────────────────────
def _snap_to_df(d):
    return pd.DataFrame({
        "j":      d.j.magnitude,
        "time":   d.time.magnitude,
        "mat":    d.mat.magnitude,
        "pos":    d.pos.magnitude,
        "rho0":   d.rho0.magnitude,
        "up":     d.up.magnitude,
        "pres":   d.pres.magnitude,
        "sigmar": d.sigmar.magnitude,
    })

snaps = []
with open(run.outputfilename, "rb") as fh:
    while True:
        try:
            snaps.append(_snap_to_df(pickle.load(fh)))
        except EOFError:
            break

pko = pd.concat(snaps, ignore_index=True)

# Unit conversions (consistent with main test17 script)
pko['time']   *= 1e6    # s  → μs
pko['pos']    *= 1e2    # m  → cm
pko['pres']   *= 1e-9   # Pa → GPa  (positive = compression)
pko['sigmar'] *= 1e-9   # Pa → GPa
pko['rho0']   *= 1e-3   # kg/m³ → g/cm³
# 'up' stays in m/s

# ── Lagrangian (initial) positions ─────────────────────────────────────────────
t0   = pko['time'].min()
pos0 = pko.loc[pko['time'] == t0, 'pos'].values.copy()   # cm, length = n_nodes
n_nodes = len(pos0)

pos0_col = []
for _, grp in pko.groupby('time', sort=True):
    n = len(grp)
    row = pos0[:n] if n <= n_nodes else np.append(pos0, np.full(n - n_nodes, pos0[-1]))
    pos0_col.extend(row)
pko['pos0'] = pos0_col   # cm

# ── Identify free surface of the target (mat == 2) ────────────────────────────
mat2_mask = pko['mat'] == 2
if mat2_mask.any():
    fs_pos0_cm = pko.loc[mat2_mask, 'pos0'].max()
else:
    fs_pos0_cm = pko['pos0'].max()   # fallback: global rightmost node

# Also find the flyer–target interface for the boundary marker
if (pko['mat'] == 1).any() and (pko['mat'] == 2).any():
    flyer_right_cm  = pko.loc[pko['mat'] == 1, 'pos0'].max()
    target_left_cm  = pko.loc[pko['mat'] == 2, 'pos0'].min()
    interface_cm    = 0.5 * (flyer_right_cm + target_left_cm)
else:
    interface_cm = None

# ── Build sorted time array and FSV series ────────────────────────────────────
unique_times = np.sort(pko['time'].unique())
n_times      = len(unique_times)

# FSV: particle velocity at free surface node at each time
fsv_t, fsv_v = [], []
for t, grp in pko.groupby('time', sort=True):
    row = grp[np.isclose(grp['pos0'], fs_pos0_cm, atol=1e-6)]
    if not row.empty:
        fsv_t.append(t)
        fsv_v.append(row['up'].iat[0])
fsv_t = np.array(fsv_t)   # μs
fsv_v = np.array(fsv_v)   # m/s

# ── Build 2-D stress grid: shape (n_times, n_nodes) ──────────────────────────
stress_grid = np.full((n_times, n_nodes), np.nan)
for i, (_, grp) in enumerate(pko.groupby('time', sort=True)):
    v = grp['pres'].values          # GPa, positive = compression
    n = min(len(v), n_nodes)
    stress_grid[i, :n] = v[:n]

# Axes coordinates for the contour plots
X_mm = pos0 * 10       # cm → mm
T_us = unique_times    # μs

# Fixed contour levels (computed once so colours are consistent across frames)
abs_max     = max(np.nanmax(np.abs(stress_grid)), STRESS_CLIP * 2)
comp_levels = np.linspace(STRESS_CLIP, abs_max, N_COMP_LEVELS)   # blue
tens_levels = np.linspace(-abs_max, -STRESS_CLIP, N_TENS_LEVELS) # red

# Frames to render
frame_indices = list(range(2, n_times + 1, max(1, FRAME_SKIP)))

# ── Figure layout ─────────────────────────────────────────────────────────────
plt.rcParams.update({'font.size': 13})

fig = plt.figure(figsize=(15, 6.5), dpi=120)
fig.suptitle('Lagrangian x-t Wave Characteristics  |  Free Surface Velocity',
             fontsize=14, fontweight='bold')

ax_xt  = fig.add_subplot(1, 2, 1)
ax_fsv = fig.add_subplot(1, 2, 2)

# x-t panel
ax_xt.set_xlabel('Lagrangian Position $X_0$ (mm)')
ax_xt.set_ylabel('Time (μs)')
ax_xt.set_title('blue = compression    red = rarefaction / tension', fontsize=11)
ax_xt.set_xlim(X_mm.min(), X_mm.max())
ax_xt.set_ylim(T_us.min(), T_us.max())
ax_xt.grid(True, alpha=0.25)

# Material boundary markers (drawn once, static)
if interface_cm is not None:
    ax_xt.axvline(interface_cm * 10, color='dimgray', lw=1.2, ls=':', alpha=0.8,
                  label='flyer | target')
ax_xt.axvline(fs_pos0_cm * 10, color='green', lw=1.4, ls='--', alpha=0.8,
              label='free surface')
ax_xt.legend(fontsize=10, loc='upper left')

# FSV panel
ax_fsv.set_xlabel('Time (μs)')
ax_fsv.set_ylabel('Free Surface Velocity (m/s)')
ax_fsv.set_title(f'Target free surface  ($X_0$ = {fs_pos0_cm * 10:.2f} mm)', fontsize=11)
ax_fsv.set_xlim(T_us.min(), T_us.max())
v_pad = max(50, (fsv_v.max() - fsv_v.min()) * 0.12) if len(fsv_v) else 50
ax_fsv.set_ylim(min(fsv_v.min() - v_pad, -20) if len(fsv_v) else -50,
                (fsv_v.max() + v_pad)           if len(fsv_v) else 500)
ax_fsv.axhline(0, color='k', lw=0.8, ls='--', alpha=0.5)
ax_fsv.grid(True, alpha=0.25)

# Animated line elements
time_bar, = ax_xt.plot([], [], color='dimgray', lw=1.0, ls='--', alpha=0.7)
fsv_line,  = ax_fsv.plot([], [], color='black', lw=2.0)
fsv_dot,   = ax_fsv.plot([], [], 'o', color='black', ms=7)

# Container for contour collections so we can remove them each frame
_contours = []

fig.tight_layout(rect=[0, 0, 1, 0.94])

# ── Animation callbacks ────────────────────────────────────────────────────────
def init():
    time_bar.set_data([], [])
    fsv_line.set_data([], [])
    fsv_dot.set_data([], [])
    return time_bar, fsv_line, fsv_dot


def update(i_end):
    """Draw characteristics up to time-step i_end and update FSV."""

    # Remove contours from the previous frame (handles both old and new matplotlib API)
    for cs in _contours:
        try:
            cs.remove()                        # matplotlib >= 3.8
        except (AttributeError, ValueError):
            for col in cs.collections:         # matplotlib < 3.8 fallback
                col.remove()
    _contours.clear()

    sub_stress = stress_grid[:i_end]
    sub_T      = T_us[:i_end]

    if i_end >= 2 and not np.all(np.isnan(sub_stress)):
        XX, TT = np.meshgrid(X_mm, sub_T)
        s_max = np.nanmax(sub_stress)
        s_min = np.nanmin(sub_stress)

        # Compressive characteristics — blue
        if s_max > STRESS_CLIP:
            lvls = comp_levels[comp_levels <= s_max]
            if len(lvls):
                cs = ax_xt.contour(XX, TT, sub_stress, levels=lvls,
                                   colors='steelblue', linewidths=0.9, alpha=0.9)
                _contours.append(cs)

        # Tensile / rarefaction characteristics — red
        if s_min < -STRESS_CLIP:
            lvls = tens_levels[tens_levels >= s_min]
            if len(lvls):
                cs = ax_xt.contour(XX, TT, sub_stress, levels=lvls,
                                   colors='firebrick', linewidths=0.9, alpha=0.9)
                _contours.append(cs)

    # Sweep line at current time
    t_now = T_us[i_end - 1]
    time_bar.set_data([X_mm.min(), X_mm.max()], [t_now, t_now])

    # FSV: reveal all points up to t_now
    mask = fsv_t <= t_now
    if mask.any():
        fsv_line.set_data(fsv_t[mask], fsv_v[mask])
        fsv_dot.set_data([fsv_t[mask][-1]], [fsv_v[mask][-1]])

    return time_bar, fsv_line, fsv_dot


ani = animation.FuncAnimation(
    fig, update,
    frames=frame_indices,
    init_func=init,
    interval=INTERVAL_MS,
    blit=False,
    repeat=True,
)

print(f"Rendering {len(frame_indices)} animation frames …")

# In Jupyter / VS Code interactive Python, FuncAnimation + plt.show() produces
# only a static blank frame.  Use to_jshtml() so the animation plays inline.
try:
    from IPython.display import HTML, display as ipy_display
    plt.rcParams['animation.html'] = 'jshtml'
    plt.close(fig)                  # suppress the duplicate static figure
    ipy_display(HTML(ani.to_jshtml(fps=max(1, 1000 // INTERVAL_MS))))
except ImportError:
    # Plain script / terminal: use the event loop directly
    plt.show()

# %%
