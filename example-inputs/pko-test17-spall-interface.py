# %%
# pyKO Test 17: Hybrid Spall + Interface Separation Analysis
# Combines spall functionality from Test 9 with interface separation analysis from Test 16
# 
# This test simulates:
# - Material 1 flyer impacting Material 2 target (all properties read from YAML config)
# - Both materials can have spall capability (strength models and fracture parameters from YAML)
# - Analysis includes spall detection, interface separation, and free surface velocity tracking
# - Fully configurable via YionAML - no hardcoded material properties or analysis parameters

########################################################################################################################
# USER CONFIGURATION SWITCHES
########################################################################################################################

# Toggle switches for analysis modules (True = ON, False = OFF)
ENABLE_INTERFACE_ANALYSIS = True  # Set to False to skip interface separation analysis
ENABLE_FSV_ANALYSIS = True         # Set to False to skip free surface velocity analysis
ENABLE_STRESS_ANALYSIS = True      # Set to False to skip stress analysis
ENABLE_SPALL_ANALYSIS = True      # Set to False to skip spall analysis

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
import os
import sys
import pickle

# Paths relative to this file so the script runs from any cwd (e.g. repo root or example-inputs)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..'))

# Set global font size for all plots
plt.rcParams.update({'font.size': 20})

sys.path.insert(1, _REPO_ROOT)
from pyko import *
import pyko
# Not using example-inputs/import-modules.py: that pulls hvplot/holoviews for notebooks; this script only needs pyko + mpl.

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
        '#FFFFFF',  # White (zero pressure)
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
# X-T DIAGRAM RENDERING
########################################################################################################################

def xt_pcolormesh(ax, xseries, cseries, **kwargs):
    """Render an x-t diagram as a per-cell pcolormesh instead of scatter markers.

    Scatter draws one discrete circular marker per (cell, output-time), which
    aliases the sharp shock fronts into a staircase of 'teeth' and leaves gaps
    between output rows. pcolormesh tiles the cells edge-to-edge, giving clean
    wave fronts and revealing the low-amplitude reflected/release waves. Unlike
    tricontourf it does NOT interpolate: every tile is one true cell x one
    output dump, so the plot resolution is exactly the data resolution
    (n_cells wide by n_output_times tall). To make the tiles finer, increase
    the mesh 'cells' (x-resolution) and/or decrease 'dtoutput' (t-resolution)
    in the YAML config. Figure dpi only changes how many screen pixels each tile
    spans, not how much information is shown.

    Assumes every output time has the same number of cells in the same
    Lagrangian order. This holds for pyKO binary output: spall fractures insert
    void/boundary cells that are NOT written to the file, so the written zones
    are a fixed, mass-conserving, index-aligned set (verified: per-cell mass is
    identical across all frames). Falls back to scatter if the data is ragged.

    x is scaled by 10 (cm -> mm) to match the original scatter calls.
    """
    x = np.asarray(xseries, dtype=float); c = np.asarray(cseries, dtype=float)
    t = np.asarray(pko['time'], dtype=float)
    times = np.unique(t); ntime = len(times); ncell = len(t) // ntime if ntime else 0
    if ntime == 0 or ntime * ncell != len(t):
        # ragged grid (unexpected) -> keep the old scatter behaviour
        return ax.scatter(x * 10, t, c=c, **kwargs)
    order = np.argsort(t, kind='stable')  # group rows by time, preserve cell order
    X = (x[order] * 10).reshape(ntime, ncell)
    C = c[order].reshape(ntime, ncell)
    if np.allclose(X, X[0]):
        # Lagrangian frame: x is the same at every time -> use 1-D separable
        # (monotonic) coordinates so pcolormesh gets exact cell edges, no warning.
        return ax.pcolormesh(X[0], times, C, shading='nearest', **kwargs)
    # Eulerian frame: the mesh moves, so x is a full 2-D field (a cell's
    # position is not monotonic in time -> pcolormesh emits a benign
    # cell-edge warning; the deforming quad mesh is still the correct picture).
    Y = t[order].reshape(ntime, ncell)
    return ax.pcolormesh(X, Y, C, shading='nearest', **kwargs)


def plot_spall_zoom(pko, pressure_cmap, pressure_norm, density_min, density_max, dpi=300):
    """Second figure: a high-resolution Lagrangian x-t zoom on the spall region.

    The spall plane is found automatically as the point of peak tension (the
    most negative Lagrangian stress) in the whole space-time field. A window is
    framed around the strongly-tensile zone (cells reaching >= 50% of the peak
    tension), so the zoom follows the physics if the setup changes rather than
    being hard-coded. Two panels are shown at the same location -- stress (where
    the two release waves cross and tension peaks) and density ratio (where the
    fracture opens and rho/rho0 drops) -- reusing xt_pcolormesh so it is the
    exact per-cell view at full data resolution.

    Returns without plotting if the data never goes into tension.
    """
    xmm = np.asarray(pko['pos0'], dtype=float) * 10.0   # mm, matches xt_pcolormesh x-axis
    t   = np.asarray(pko['time'], dtype=float)           # microseconds
    p   = np.asarray(pko['pres'], dtype=float)           # GPa
    ip = int(np.argmin(p)); peak = p[ip]
    if peak >= 0:
        print("Spall zoom: no tension in the data; skipping zoom figure.")
        return
    sel = p <= 0.5 * peak                                # strongly-tensile cells
    xs = xmm[sel]; ts = t[sel]
    xpad = 0.15 * (xs.max() - xs.min()) + 0.003
    tpad = 0.15 * (ts.max() - ts.min()) + 0.002
    xlim = (xs.min() - xpad, xs.max() + xpad)
    ylim = (ts.min() - tpad, ts.max() + tpad)

    print(f"Creating high-resolution spall-region zoom "
          f"(peak tension {peak:.2f} GPa at x0={xmm[ip]:.3f} mm, t={t[ip]*1e3:.1f} ns)...")
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(14, 8), dpi=dpi)
    fig.suptitle(f'Spall region zoom - peak tension {peak:.2f} GPa '
                 f'@ x0={xmm[ip]:.3f} mm, t={t[ip]*1e3:.1f} ns', fontsize=18)

    mA = xt_pcolormesh(axa, pko.pos0, pko.pres, cmap=pressure_cmap, norm=pressure_norm)
    axa.set_title('Lagrangian: Stress (GPa)')
    plt.colorbar(mA, ax=axa, label='Stress (GPa)')
    mB = xt_pcolormesh(axb, pko.pos0, pko.density_ratio,
                       cmap='coolwarm', vmin=density_min, vmax=density_max)
    axb.set_title('Lagrangian: Density Ratio')
    plt.colorbar(mB, ax=axb, label=r'$\rho/\rho_0$')
    for ax in (axa, axb):
        ax.set_xlim(*xlim); ax.set_ylim(*ylim)
        ax.set_xlabel('Initial Position (mm)'); ax.set_ylabel('Time (μs)')
        ax.grid(True, alpha=0.2, which='both')
        # crosshair marks the peak-tension (spall) point
        ax.axvline(xmm[ip], color='k', ls=':', lw=0.9, alpha=0.6)
        ax.axhline(t[ip], color='k', ls=':', lw=0.9, alpha=0.6)
    plt.tight_layout()
    plt.show()


def _compute_xt_fans(pko, run=None):
    """Compute the x-t rarefaction-fan data: (a) fan edges traced from pyKO's own
    sound-speed field and (b) analytic Hugoniot tracking of the shock and the
    plastic/elastic release fronts. Returns a dict for _draw_xt_fans, or None.

    (a) SIMULATION fans -- a release fan is bounded by two characteristics. In
    initial-position coordinates h = x0 a characteristic obeys
        dh/dt = +/- c * (rho / rho0)                                (Lagrangian)
    with c = pyKO's local sound speed. A head (shock arrival) and a tail (bottom
    of the release dip) are integrated into the material for each free surface.

    (b) ANALYTIC tracking -- impedance-match the flyer/target for the impact
    particle velocity up, then track: shock Us = c0 + s*up; plastic (bulk)
    release c0 + 2s*up; and the elastic release head at the shocked-state
    longitudinal Lagrangian speed sqrt(c_bulk_E^2 + 4G/3rho)*(rho/rho0). The
    elastic head is fast but limited to ~HEL amplitude.
    """
    from scipy.interpolate import RegularGridInterpolator
    from scipy.ndimage import uniform_filter1d
    # ---- gridded fields (ntime, ncell) on the fixed Lagrangian mesh ----------
    t_all = np.asarray(pko['time'], dtype=float)
    times = np.unique(t_all); nt = len(times)
    ncell = len(t_all) // nt if nt else 0
    if nt < 3 or nt * ncell != len(t_all):
        print("x-t fans: unexpected/ragged data; skipping."); return None
    order = np.argsort(t_all, kind='stable')
    h = np.asarray(pko['pos0'], dtype=float)[order].reshape(nt, ncell)[0] * 10.0   # mm, fixed
    P  = np.asarray(pko['pres'], dtype=float)[order].reshape(nt, ncell)            # GPa
    cs = np.asarray(pko['cs'],   dtype=float)[order].reshape(nt, ncell) * 1e-3     # m/s -> mm/us
    rr = (np.asarray(pko['rho'], dtype=float) /
          np.asarray(pko['rho0'], dtype=float))[order].reshape(nt, ncell)
    Vlag = cs * rr                                       # mm/us, Lagrangian char. speed
    Vlag[0] = Vlag[1]                                    # cs = 0 at t = 0 -> use next frame
    si = np.argsort(h); h = h[si]; P = P[:, si]; Vlag = Vlag[:, si]   # h strictly ascending
    speed = RegularGridInterpolator((times, h), Vlag, bounds_error=False, fill_value=None)

    def launch_times(cell):
        """(head, tail) launch times at a free-surface cell, read from its own
        compression history: head = shock arrival, tail = bottom of release dip."""
        p = P[:, cell]; pc = float(np.nanmax(p))
        if pc <= 0:
            return times[0], times[0]
        ps = uniform_filter1d(p, 5)                      # light smoothing vs cell noise
        above = np.where(ps > 0.1 * pc)[0]
        t_head = times[above[0]] if above.size else times[0]
        ipk = int(np.argmax(ps))
        rise = np.where(np.diff(ps[ipk:]) > 0)[0]        # first up-tick after the peak
        k_tail = ipk + (rise[0] + 1 if rise.size else nt - 1 - ipk)
        return t_head, times[min(k_tail, nt - 1)]

    def trace(h0, t0, sign):
        """Integrate dh/dt = sign * c*rho/rho0 from (h0, t0) into the material."""
        hs = [h0]; ts = [t0]; hc = h0
        for k in range(int(np.searchsorted(times, t0)), nt - 1):
            v = float(speed((times[k], np.clip(hc, h[0], h[-1]))))
            hc = hc + sign * v * (times[k + 1] - times[k])
            if hc < h[0] or hc > h[-1]:
                break
            hs.append(hc); ts.append(times[k + 1])
        return np.asarray(hs), np.asarray(ts)

    cf = int(np.argmin(h)); ct = int(np.argmax(h))       # flyer FS (left), target FS (right)
    h_ff, h_tf = h[cf], h[ct]
    fh, ft = launch_times(cf); th, tt = launch_times(ct)
    fan_fh = trace(h_ff, fh, +1); fan_ft = trace(h_ff, ft, +1)   # flyer release -> +h
    fan_th = trace(h_tf, th, -1); fan_tt = trace(h_tf, tt, -1)   # target release -> -h
    print(f"Rarefaction fans (sim): flyer surface head@{fh:.4f}/tail@{ft:.4f} us, "
          f"target surface head@{th:.4f}/tail@{tt:.4f} us")

    # ---- (b) analytic Hugoniot tracking (calculation-based) ------------------
    analytic = None
    if run is not None and getattr(run, 'nmat', 0) >= 2:
        try:
            from scipy.optimize import brentq
            jf, jt = 0, run.nmat - 1                      # flyer (left), target (right)
            # pyKO code units -> SI: c0,up in cm/us (x1e4 = m/s); rho0 in g/cm3; length cm (x10 = mm)
            cF = run.ieos[jf].c0 * 1e4; sF = run.ieos[jf].s1; rF = run.ieos[jf].rho0
            cT = run.ieos[jt].c0 * 1e4; sT = run.ieos[jt].s1; rT = run.ieos[jt].rho0
            # shear modulus / yield (code pressure unit = megabar -> x1e11 = Pa); may be absent
            GF = getattr(run.istr[jf], 'gmod', 0.0) * 1e11; YF = getattr(run.istr[jf], 'ys', 0.0) * 1e11
            GT = getattr(run.istr[jt], 'gmod', 0.0) * 1e11; YT = getattr(run.istr[jt], 'ys', 0.0) * 1e11
            Vimp = float(run.iupstart[jf]) * 1e4          # flyer impact velocity (m/s)
            Lf = abs(run.ilength[jf]) * 10.0              # flyer thickness (mm)
            hff = run.ixstart[jf] * 10.0                  # flyer free surface (mm)

            def elastic_lag(c0, s, G, rho0_si, u):
                """Shocked-state longitudinal-elastic Lagrangian speed (m/s):
                a_L = sqrt(c_bulk_Euler^2 + 4G/3rho) * (rho/rho0), with the bulk
                (Hugoniot) speed and the shear term from the strength model."""
                Us_ = c0 + s*u; comp = Us_/(Us_-u); cbE = (c0 + 2*s*u)/comp
                rho = rho0_si*comp
                return np.sqrt(cbE**2 + (4.0/3.0)*G/rho)*comp
            if Vimp > 0:
                # impedance match: rho0_T*(cT+sT*u)*u = rho0_F*(cF+sF*(V-u))*(V-u)
                up = brentq(lambda u: rT*(cT+sT*u)*u - rF*(cF+sF*(Vimp-u))*(Vimp-u), 1.0, Vimp-1.0)
                upF = Vimp - up                           # flyer shock up-jump
                UsT = (cT + sT*up)   * 1e-3               # mm/us
                relT = (cT + 2*sT*up) * 1e-3
                UsF = (cF + sF*upF)  * 1e-3
                relF = (cF + 2*sF*upF) * 1e-3
                t1 = Lf / UsF                             # shock reaches flyer FS
                tI = t1 + Lf / relF                       # release back at interface
                t_sfs = h_tf / UsT                        # shock reaches target FS
                # would the PLASTIC (bulk) release catch the shock? h_rel = h_shock
                t_catch = relT * tI / (relT - UsT) if relT > UsT else np.inf
                h_catch = UsT * t_catch
                # elastic (longitudinal) release head: fast precursor, HEL-limited amplitude
                aeF = elastic_lag(cF, sF, GF, rF*1e3, upF) * 1e-3   # mm/us
                aeT = elastic_lag(cT, sT, GT, rT*1e3, up)  * 1e-3
                tIe = t1 + Lf/aeF                          # elastic release back at interface
                te_catch = aeT * tIe / (aeT - UsT) if aeT > UsT else np.inf
                he_catch = UsT * te_catch
                e_before = he_catch < h_tf                 # catches shock before breakout?
                he_end = min(he_catch, h_tf); te_end = min(te_catch, tIe + h_tf/aeT)
                # HEL ~ Y*(K0 + 4G/3)/(2G); K0 = rho0*c0^2
                K0T = (rT*1e3)*cT**2
                HELT = (YT*(K0T + 4/3*GT)/(2*GT))/1e9 if GT > 0 else 0.0
                Pshock = rT*(cT+sT*up)*up * 1e-3          # GPa (rho0 g/cm3 * m/s^2 -> /1e3? see note)
                analytic = dict(up=up, UsT=UsT*1e3, relT=relT*1e3, aeT=aeT*1e3, aeF=aeF*1e3,
                                t1=t1, tI=tI, t_sfs=t_sfs, h_catch=h_catch, t_catch=t_catch,
                                he_catch=he_catch, te_catch=te_catch,
                                shock=(np.array([0.0, h_tf]), np.array([0.0, t_sfs])),
                                rel_f=(np.array([hff, 0.0]), np.array([t1, tI])),
                                # plastic release runs from interface up to the free surface
                                rel_t=(np.array([0.0, relT*(min(tI + h_tf/relT, times.max()) - tI)]),
                                       np.array([tI, min(tI + h_tf/relT, times.max())])),
                                # elastic head: flyer FS -> interface -> catch point (or free surface)
                                elastic=(np.array([hff, 0.0, he_end]), np.array([t1, tIe, te_end])),
                                e_catch_pt=(he_catch, te_catch) if e_before else None,
                                HELT=HELT)
                catches = h_catch < h_tf and t_catch < t_sfs
                print(f"Analytic: up={up:.0f} m/s, Us(target)={UsT*1e3:.0f} m/s.")
                print(f"  Plastic release (C+2S*up={relT*1e3:.0f}) catches shock at "
                      f"h={h_catch:.2f} mm -> {'CATCHES in target' if catches else 'NEVER in target'} "
                      f"(target={h_tf:.2f} mm).")
                print(f"  Elastic head   (long. {aeT*1e3:.0f}) would overtake shock at "
                      f"h={he_catch:.2f} mm -> "
                      f"{'OVERTAKES before breakout' if e_before else f'target only {h_tf:.2f} mm, shock breaks out first'}; "
                      f"amplitude ~HEL ({HELT:.2f} GPa) << shock, so peak stress ~unattenuated either way.")
        except Exception as e:
            print(f"Analytic tracking skipped: {e}")

    amax = float(np.nanmax(np.abs(P)))
    return dict(h=h, times=times, P=P, nt=nt, ncell=ncell,
                fan_fh=fan_fh, fan_ft=fan_ft, fan_th=fan_th, fan_tt=fan_tt,
                h_ff=h_ff, h_tf=h_tf, analytic=analytic, amax=amax)


def _draw_xt_fans(ax, d, pressure_cmap, transpose=False, legend=True, fs=7, legend_loc='upper left', norm=None):
    """Draw the x-t rarefaction-fan diagram (data from _compute_xt_fans) onto ax.
    If transpose, the axes are swapped so time is horizontal and initial position
    is vertical (used for stacking above a free-surface-velocity trace)."""
    from matplotlib.colors import AsinhNorm
    h, times, P = d['h'], d['times'], d['P']; analytic = d['analytic']
    if norm is None:
        norm = AsinhNorm(linear_width=0.1 * d['amax'], vmin=-d['amax'], vmax=d['amax'])
    if transpose:
        qm = ax.pcolormesh(times, h, P.T, cmap=pressure_cmap, norm=norm, shading='nearest')
    else:
        qm = ax.pcolormesh(h, times, P, cmap=pressure_cmap, norm=norm, shading='nearest')
    def L(pos, tim):                 # (x, y) with the chosen orientation
        return (tim, pos) if transpose else (pos, tim)
    def posline(val, **kw):          # a constant-position reference line
        (ax.axhline if transpose else ax.axvline)(val, **kw)
    ax.plot(*L(*d['fan_fh']), color='#0044ff', lw=2.0, label='flyer rarefaction: head (sim)')
    ax.plot(*L(*d['fan_ft']), color='#0044ff', lw=2.0, ls='--', label='flyer rarefaction: tail (sim)')
    ax.plot(*L(*d['fan_th']), color='#00a000', lw=2.0, label='target rarefaction: head (sim)')
    ax.plot(*L(*d['fan_tt']), color='#00a000', lw=2.0, ls='--', label='target rarefaction: tail (sim)')
    if analytic is not None:
        hr = np.r_[analytic['rel_f'][0], analytic['rel_t'][0]]
        tr_ = np.r_[analytic['rel_f'][1], analytic['rel_t'][1]]
        ax.plot(*L(analytic['shock'][0], analytic['shock'][1]), 'k--', lw=1.8, label='shock (C+S·up)')
        ax.plot(*L(hr, tr_), 'k-', lw=1.8, label='plastic release (C+2S·up)')
        ax.plot(*L(analytic['elastic'][0], analytic['elastic'][1]), color='#ff8800', lw=2.0,
                label=f"elastic release head (~HEL {analytic['HELT']:.2f} GPa)")
        if analytic['e_catch_pt'] is not None:
            xy = L(np.array([analytic['e_catch_pt'][0]]), np.array([analytic['e_catch_pt'][1]]))
            ax.plot(xy[0], xy[1], 'o', color='#ff8800', ms=8, mec='k', zorder=5)
    posline(0.0,        color='0.35', ls='-.', lw=1.2, label='flyer / target interface')
    posline(d['h_ff'],  color='0.55', ls=':',  lw=1.2, label='flyer free surface')
    posline(d['h_tf'],  color='k',    ls='--', lw=1.2, label='target free surface')
    prange = (d['h_ff'] - 0.05 * (d['h_tf'] - d['h_ff']), d['h_tf'] * 1.03)
    if transpose:
        ax.set_xlim(0, times.max()); ax.set_ylim(*prange); ax.invert_yaxis()
        ax.set_ylabel('Initial Position (mm)')
    else:
        ax.set_xlim(*prange); ax.set_ylim(0, times.max())
        ax.set_xlabel('Initial Position (mm)'); ax.set_ylabel('Time (μs)')
    if legend:
        ax.legend(loc=legend_loc, fontsize=fs)
    return qm


def plot_rarefaction_fans(pko, pressure_cmap, run=None, dpi=200):
    """Lagrangian x-t stress diagram: sim-traced rarefaction fans plus analytic
    Hugoniot tracking of the shock and the plastic/elastic release fronts."""
    d = _compute_xt_fans(pko, run)
    if d is None:
        return
    fig, ax = plt.subplots(figsize=(10, 8), dpi=dpi)
    _draw_xt_fans(ax, d, pressure_cmap, transpose=False)
    ax.set_title('Rarefaction fans (sim) vs analytic shock/release tracking')
    plt.tight_layout()
    plt.show()


def plot_xt_fsv_stack(pko, pressure_cmap, run=None, dpi=200):
    """Stacked figure relating the x-t diagram to the free-surface velocity:
    top = x-t stress map + wave tracking (time horizontal, position vertical with
    the target free surface at the bottom); bottom = free-surface velocity vs
    time on the SAME time axis. Vertical guide lines mark the wave arrivals that
    produce the FSV features (shock breakout -> rise, spall -> pullback)."""
    d = _compute_xt_fans(pko, run)
    if d is None:
        return
    # free-surface velocity = up of the right-most (free-surface) node at each time
    times = np.sort(pko['time'].unique())
    fsv = np.array([pko.loc[g.index[np.argmax(g['pos'].values)], 'up']
                    for _, g in pko.groupby('time')])
    tg = np.array([t for t, _ in pko.groupby('time')])
    o = np.argsort(tg); tg = tg[o]; fsv = fsv[o]
    # guide-line times: shock breakout, FSV peak, spall pull-back minimum
    guides = []
    if d['analytic'] is not None and np.isfinite(d['analytic']['t_sfs']):
        guides.append(('breakout', d['analytic']['t_sfs']))
    if fsv.size > 5:
        ipk = int(np.argmax(fsv)); guides.append(('peak', tg[ipk]))
        seg = fsv[ipk:]
        if seg.size > 2:
            imin = ipk + int(np.argmin(seg))         # spall pull-back = velocity minimum after peak
            if imin > ipk:
                guides.append(('pull-back', tg[imin]))
    # constrained_layout auto-spaces the panels so labels never collide
    fig = plt.figure(figsize=(11, 13), dpi=dpi, constrained_layout=True)
    gs = fig.add_gridspec(4, 1, height_ratios=[0.10, 3, 1.1, 1.1])
    axC = fig.add_subplot(gs[0])                       # stress colorbar band (top)
    axT = fig.add_subplot(gs[1])
    axB = fig.add_subplot(gs[2], sharex=axT)
    axL = fig.add_subplot(gs[3]); axL.axis('off')     # dedicated legend panel
    qm = _draw_xt_fans(axT, d, pressure_cmap, transpose=True, legend=False)
    axT.tick_params(labelbottom=False)                # time labels only on the FSV panel below
    cb = fig.colorbar(qm, cax=axC, orientation='horizontal')
    _am = d['amax']                                   # clean ticks (asinh auto-ticks crowd at 0)
    _ticks = [t for t in (-_am, -1.0, 0.0, 1.0, _am) if -_am <= t <= _am]
    cb.set_ticks(_ticks); cb.ax.set_xticklabels([f'{t:.1f}' for t in _ticks])
    cb.set_label('Stress (GPa) — pyKO  (tension <0<  compression)', fontsize=11)
    axC.xaxis.set_ticks_position('top'); axC.xaxis.set_label_position('top')
    axB.plot(tg, fsv, 'b-', lw=2, label='free surface velocity')
    axB.set_xlabel('Time (μs)'); axB.set_ylabel('Free surface\nvelocity (m/s)')
    axB.grid(True, alpha=0.3)
    axB.set_xlim(0, times.max())
    for lbl, t in guides:                            # vertical guides across both panels
        axT.axvline(t, color='0.5', ls='--', lw=1.0, alpha=0.8)
        axB.axvline(t, color='0.5', ls='--', lw=1.0, alpha=0.8)
        axB.annotate(lbl, xy=(t, 1.0), xycoords=('data', 'axes fraction'),
                     xytext=(2, -2), textcoords='offset points', fontsize=9,
                     rotation=90, va='top', ha='left', color='0.35')
    # collect all handles/labels (x-t panel + FSV) into the legend-only panel
    hs, ls_ = axT.get_legend_handles_labels()
    hb, lb = axB.get_legend_handles_labels()
    axL.legend(hs + hb, ls_ + lb, loc='center', ncol=3, fontsize=12,
               frameon=False, handlelength=2.5, columnspacing=1.8)
    plt.show()

########################################################################################################################

# Automatically select the appropriate input file based on interface analysis setting
if ENABLE_INTERFACE_ANALYSIS:
    # Use the new material database configuration
    filein = os.path.join(_SCRIPT_DIR, 'test17-spall-interface', 'test17-material-config.yml')
    print("🔬 Using material database configuration WITH interface separation physics")
else:
    # Fallback to original configuration without interface separation
    filein = os.path.join(_SCRIPT_DIR, 'test17-spall-interface', 'test17-without-interface-separation.yml')
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

# Get material names from YAML configuration
try:
    # Debug: Check what's available in the run object
    print("🔍 DEBUGGING MATERIAL NAME EXTRACTION:")
    print(f"Debug: run.ieos[0] type: {type(run.ieos[0])}")
    print(f"Debug: run.ieos[0] attributes: {dir(run.ieos[0])}")
    print(f"Debug: run.ieos[1] type: {type(run.ieos[1])}")
    print(f"Debug: run.ieos[1] attributes: {dir(run.ieos[1])}")
    
    # Try multiple ways to get material names
    mat1_name = "Material 1"
    mat2_name = "Material 2"
    
    # Method 1: Try to get from EOS name attribute
    if hasattr(run.ieos[0], 'name'):
        mat1_name = run.ieos[0].name
        print(f"Debug: Found mat1 name in EOS: {mat1_name}")
    elif hasattr(run.ieos[0], 'eosname'):
        mat1_name = run.ieos[0].eosname
        print(f"Debug: Found mat1 name in eosname: {mat1_name}")
    
    if hasattr(run.ieos[1], 'name'):
        mat2_name = run.ieos[1].name
        print(f"Debug: Found mat2 name in EOS: {mat2_name}")
    elif hasattr(run.ieos[1], 'eosname'):
        mat2_name = run.ieos[1].eosname
        print(f"Debug: Found mat2 name in eosname: {mat2_name}")
    
    # Method 2: Try to get from run object attributes
    if hasattr(run, 'material_names'):
        if len(run.material_names) >= 2:
            mat1_name = run.material_names[0]
            mat2_name = run.material_names[1]
            print(f"Debug: Found names in run.material_names: {mat1_name}, {mat2_name}")
    
    # Method 3: Try to get from YAML config directly
    try:
        import yaml
        with open(filein, 'r') as f:
            yaml_config = yaml.safe_load(f)
        
        if 'mat1' in yaml_config and 'eos' in yaml_config['mat1'] and 'name' in yaml_config['mat1']['eos']:
            mat1_name = yaml_config['mat1']['eos']['name']
            print(f"Debug: Found mat1 name in YAML: {mat1_name}")
        
        if 'mat2' in yaml_config and 'eos' in yaml_config['mat2'] and 'name' in yaml_config['mat2']['eos']:
            mat2_name = yaml_config['mat2']['eos']['name']
            print(f"Debug: Found mat2 name in YAML: {mat2_name}")

        # Define thickness strings for summary/reporting
        try:
            al_len_m = float(yaml_config['mat1']['mesh']['length'])
            cu_len_m = float(yaml_config['mat2']['mesh']['length'])
            al_thickness_str = f"{al_len_m*1e6:.0f} μm"
            cu_thickness_str = f"{cu_len_m*1e6:.0f} μm"
        except Exception as _:
            # Fallback placeholders; will be replaced later if available
            al_thickness_str = "n/a"
            cu_thickness_str = "n/a"
    except Exception as yaml_error:
        print(f"Debug: Could not read YAML directly: {yaml_error}")
    
    # Clean up names (remove "flyer" or "target" suffixes)
    mat1_name = mat1_name.replace(" flyer", "").replace(" target", "")
    mat2_name = mat2_name.replace(" flyer", "").replace(" target", "")
    
    print(f"✅ Final material names: {mat1_name}, {mat2_name}")
    print(f"Material 1 (Impactor): {mat1_name}")
    print(f"Material 2 (Target): {mat2_name}")
    print()
    
    print(f"{'Property':<25} {mat1_name:<15} {mat2_name:<15} {'Unit':<10}")
    print("-"*70)
except Exception as e:
    print(f"Warning: Could not extract material names from YAML: {e}")
    print(f"{'Property':<25} {'Material 1':<15} {'Material 2':<15} {'Unit':<10}")
    print("-"*70)
    mat1_name = "Material 1"
    mat2_name = "Material 2"

# Get material properties from run object - MANDATORY YAML EXTRACTION
# NO FALLBACK VALUES - Script must fail if YAML parsing fails

try:
    print(f"Debug: Number of materials in run object: {run.nmat}")
    print(f"Debug: Available run attributes: {[attr for attr in dir(run) if not attr.startswith('_')]}")
    
    # Check if we have the expected number of materials
    if run.nmat < 2:
        raise ValueError(f"Expected 2 materials, found {run.nmat}")
    
    # Access material properties using pyKO's actual structure
    mat1_fracture = run.ifrac[0]  # Material 1 fracture properties
    mat2_fracture = run.ifrac[1]  # Material 2 fracture properties
    mat1_strength = run.istr[0]   # Material 1 strength properties  
    mat2_strength = run.istr[1]   # Material 2 strength properties
    
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
    mat1_length_raw = run.ilength[0]  # Raw units from pyKO
    mat2_length_raw = run.ilength[1]  # Raw units from pyKO
    
    print(f"Debug: Raw lengths - {mat1_name}: {mat1_length_raw}, {mat2_name}: {mat2_length_raw}")
    
    # PyKO uses Wilkins book units: lengths in cm
    mat1_length_cm = mat1_length_raw  # Already in cm
    mat2_length_cm = mat2_length_raw  # Already in cm
    
    # Convert cm to μm for display
    mat1_thickness_um = mat1_length_cm * 1e4  # cm to μm
    mat2_thickness_um = mat2_length_cm * 1e4  # cm to μm
    
    print(f"Debug: Lengths in cm - {mat1_name}: {mat1_length_cm} cm, {mat2_name}: {mat2_length_cm} cm")
    print(f"Debug: Converted to μm - {mat1_name}: {mat1_thickness_um} μm, {mat2_name}: {mat2_thickness_um} μm")
    
    # Create formatted strings for display
    if mat1_thickness_um >= 1000:
        mat1_thickness_str = f"{mat1_thickness_um/1000:.1f} mm"
    else:
        mat1_thickness_str = f"{mat1_thickness_um:.0f} μm"
        
    if mat2_thickness_um >= 1000:
        mat2_thickness_str = f"{mat2_thickness_um/1000:.1f} mm"
    else:
        mat2_thickness_str = f"{mat2_thickness_um:.0f} μm"
    
    # Extract ALL material parameters from YAML config using pyKO structure
    # Spall thresholds (fracture pressures) - check units!
    mat1_pfrac_raw = mat1_fracture.pfrac  
    mat2_pfrac_raw = mat2_fracture.pfrac  
    
    print(f"Debug: Raw pfrac values - {mat1_name}: {mat1_pfrac_raw}, {mat2_name}: {mat2_pfrac_raw}")
    
    # PyKO uses Wilkins book units: pressures in Mbar (1 Mbar = 100 GPa)
    mat1_spall_threshold_from_yaml = mat1_pfrac_raw * 100.0  # Convert Mbar to GPa
    mat2_spall_threshold_from_yaml = mat2_pfrac_raw * 100.0  # Convert Mbar to GPa
    print("Debug: Converting pfrac from Mbar to GPa (1 Mbar = 100 GPa)")
    
    print(f"Debug: Converted spall thresholds - {mat1_name}: {mat1_spall_threshold_from_yaml:.6f} GPa, {mat2_name}: {mat2_spall_threshold_from_yaml:.6f} GPa")
    
    # Density distension limits (nrhomin)
    mat1_nrhomin = mat1_fracture.nrhomin
    mat2_nrhomin = mat2_fracture.nrhomin
    
    # Calculate spall density threshold from nrhomin (use the more restrictive one)
    # If any material has nrhomin = 1.0 (perfectly brittle), we need to handle it specially
    if mat2_nrhomin >= 1.0:
        # For perfectly brittle materials, look for any density reduction at all
        spall_density_threshold = 0.999  # Detect even tiny density drops
        print(f"⚠️  {mat2_name} is perfectly brittle (nrhomin = {mat2_nrhomin}), using sensitive density threshold")
    else:
        spall_density_threshold = max(mat1_nrhomin, mat2_nrhomin)
    
    print(f"Material thicknesses from YAML: {mat1_name} = {mat1_thickness_str}, {mat2_name} = {mat2_thickness_str}")
    print(f"Spall thresholds from YAML: {mat1_name} = {mat1_spall_threshold_from_yaml:.3f} GPa, {mat2_name} = {mat2_spall_threshold_from_yaml:.3f} GPa")
    print(f"Density limits from YAML: {mat1_name} nrhomin = {mat1_nrhomin:.2f}, {mat2_name} nrhomin = {mat2_nrhomin:.2f}")
    print(f"Using spall density threshold: {spall_density_threshold:.3f} (from max nrhomin)")
    
    # Define global plotting parameters based on YAML values
    # Density ratio plot range
    density_vmin = min(mat1_nrhomin, mat2_nrhomin) * 0.95  # Slightly below lowest nrhomin
    density_vmax = 1.0  # Maximum is always 1.0 (original density)
    
    # Pressure plot range based on material spall thresholds
    pressure_range_gpa = max(mat1_spall_threshold_from_yaml, mat2_spall_threshold_from_yaml) * 2.0  # 2x max spall threshold
    
    print(f"Dynamic plot ranges from YAML: Density [{density_vmin:.2f}, {density_vmax:.1f}], Pressure ±{pressure_range_gpa:.2f} GPa")
    print("✅ ALL PARAMETERS EXTRACTED FROM YAML - NO HARDCODED VALUES")

    print(f"{'Thickness':<25} {mat1_thickness_str:<15} {mat2_thickness_str:<15} {'-':<10}")
    print(f"{'Cells':<25} {run.inodes[0]//2:<15} {run.inodes[1]//2:<15} {'count':<10}")  # inodes = 2x cells
    
    # PyKO units: velocity in cm/μs, convert to m/s
    mat1_vel_ms = run.iupstart[0] * 10000  # cm/μs to m/s (cm/μs * 1e6 μs/s * 1e-2 m/cm = 1e4)
    mat2_vel_ms = run.iupstart[1] * 10000  # cm/μs to m/s
    print(f"{'Initial velocity':<25} {mat1_vel_ms:<15.1f} {mat2_vel_ms:<15.1f} {'m/s':<10}")
    
    # PyKO units: density in g/cm³, convert to kg/m³
    mat1_rho_kgm3 = run.irhostart[0] * 1000  # g/cm³ to kg/m³
    mat2_rho_kgm3 = run.irhostart[1] * 1000  # g/cm³ to kg/m³
    print(f"{'Density':<25} {mat1_rho_kgm3:<15.0f} {mat2_rho_kgm3:<15.0f} {'kg/m³':<10}")
    
    # Get EOS objects
    mat1_eos = run.ieos[0]  # Material 1 EOS
    mat2_eos = run.ieos[1]  # Material 2 EOS
    
    # PyKO units: sound speed in cm/μs, convert to m/s
    mat1_c0_ms = mat1_eos.c0 * 10000  # cm/μs to m/s
    mat2_c0_ms = mat2_eos.c0 * 10000  # cm/μs to m/s
    print(f"{'Sound speed (c0)':<25} {mat1_c0_ms:<15.0f} {mat2_c0_ms:<15.0f} {'m/s':<10}")
    print(f"{'EOS parameter (s1)':<25} {mat1_eos.s1:<15.2f} {mat2_eos.s1:<15.2f} {'-':<10}")
    print(f"{'Gruneisen (gamma0)':<25} {mat1_eos.gamma0:<15.2f} {mat2_eos.gamma0:<15.2f} {'-':<10}")
    print(f"{'Specific heat (cv)':<25} {mat1_eos.cv:<15.0f} {mat2_eos.cv:<15.0f} {'eu/(K·cm³)':<10}")
    print()

    # Strength and fracture parameters
    print("STRENGTH & FRACTURE PARAMETERS:")
    print(f"{'Property':<25} {mat1_name:<15} {mat2_name:<15} {'Unit':<10}")
    print("-"*70)
    print(f"{'Strength model':<25} {run.istrid[0]:<15} {run.istrid[1]:<15} {'-':<10}")
    print(f"{'Shear modulus':<25} {mat1_strength.gmod/1e9:<15.1f} {mat2_strength.gmod/1e9:<15.1f} {'GPa':<10}")
    print(f"{'Yield strength':<25} {mat1_strength.ys/1e6:<15.1f} {mat2_strength.ys/1e6:<15.1f} {'MPa':<10}")

    # Fracture parameters
    mat1_pfrac = mat1_fracture.pfrac / 1e6 if mat1_fracture.pfrac < 1e15 else float('inf')
    mat2_pfrac = mat2_fracture.pfrac / 1e6 if mat2_fracture.pfrac < 1e15 else float('inf')

    if ENABLE_INTERFACE_ANALYSIS:
        print(f"{'Spall threshold':<25} {mat1_pfrac:<15.1f} {mat2_pfrac:<15.1f} {'MPa':<10}")
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

# run pyko - let pyKO compute a stable CFL-limited initial time step
# (min(dx/cs)/10). Do NOT force userdtstart=run.dtstart: the YAML dtstart
# (1e-9 s) is ~50x larger than the stable step for this ~1 um mesh, which
# makes the opening step blow up (impact face -> ~390 GPa, ~8600 m/s) and
# drives runaway spall. See pyko.py makegrid initial-timestep guess.
pyko.run(fin=filein, verbose=True)

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

    # Prepare spall thresholds from YAML/run objects
    # Prefer values from the run object's fracture models; fall back to YAML if available
    def _get_attr(obj, name):
        try:
            return getattr(obj, name)
        except Exception:
            return None

    # Defaults if nothing is found (kept conservative to avoid false positives)
    al_spall_threshold_from_yaml = None
    cu_spall_threshold_from_yaml = None
    spall_density_threshold = None

    # Prefer thresholds from mandatory extraction (pfrac in code units Mbar → GPa via *100)
    try:
        al_spall_threshold_from_yaml = mat1_spall_threshold_from_yaml
        cu_spall_threshold_from_yaml = mat2_spall_threshold_from_yaml
    except NameError:
        pass

    # Fall back: YAML pfrac is in Pa → GPa
    try:
        if 'yaml_config' in globals() or 'yaml_config' in locals():
            cfg = yaml_config
            if al_spall_threshold_from_yaml is None:
                al_spall_threshold_from_yaml = abs(float(cfg['mat1']['frac'].get('pfrac', 0.0))) / 1e9
            if cu_spall_threshold_from_yaml is None:
                cu_spall_threshold_from_yaml = abs(float(cfg['mat2']['frac'].get('pfrac', 0.0))) / 1e9
            if spall_density_threshold is None:
                spall_density_threshold = float(cfg['mat2']['frac'].get('nrhomin', 0.9))
    except Exception:
        pass

    # nrhomin from fracture object if still missing
    try:
        nrhomin_val = _get_attr(mat2_fracture, 'nrhomin')
        if spall_density_threshold is None and nrhomin_val is not None:
            spall_density_threshold = float(nrhomin_val)
    except Exception:
        pass

    # Final safe fallbacks if still missing
    if al_spall_threshold_from_yaml is None:
        # 0.276 GPa default (276 MPa) commonly used for Al demo configs
        al_spall_threshold_from_yaml = 0.276
    if cu_spall_threshold_from_yaml is None:
        cu_spall_threshold_from_yaml = 0.276
    if spall_density_threshold is None:
        spall_density_threshold = 0.9

    # Spall detection based on density reduction AND pressure thresholds
    # Identify regions where density has dropped significantly (indicating fracture/spall)
    pko['density_ratio'] = pko['rho'] / pko['rho0']
    
    # Use dual spall detection: density-based AND pressure-based
    print(f"Using density ratio threshold for spall detection: {spall_density_threshold:.3f} (from YAML nrhomin)")
    print(
        f"Also checking pressure-based spall: {mat1_name} > {al_spall_threshold_from_yaml:.3f} GPa, "
        f"{mat2_name} > {cu_spall_threshold_from_yaml:.3f} GPa tensile"
    )
    
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
        
        print(f"Using actual data ranges - Stress: [{pres_min:.2f}, {pres_max:.2f}] GPa, Density: [{density_min:.3f}, {density_max:.3f}]")
        print(f"Stress colormap: Red (tension) -> Gray (zero) -> Blue (compression)")

        # Symmetric stress normalization: equal colormap range on both sides of zero
        _stress_abs_max = max(abs(pres_min), abs(pres_max))
        pressure_norm = colors.Normalize(vmin=-_stress_abs_max, vmax=_stress_abs_max)

        # Stress (Eulerian) - Custom colormap with symmetric zero-centered scaling
        xt_pres_eul = xt_pcolormesh(ax1, pko.pos, pko.pres, cmap=pressure_cmap,
                                 norm=pressure_norm)
        ax1.set_xlabel('Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Eulerian: Stress (GPa, symmetric)')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.minorticks_on()
        ax1.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_pres_eul, ax=ax1, label='Stress (GPa)')
        
        # Particle velocity (Eulerian) - Auto-scaled
        xt_up_eul = xt_pcolormesh(ax2, pko.pos, pko.up, cmap='inferno')
        ax2.set_xlabel('Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Eulerian: Particle Velocity')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.minorticks_on()
        ax2.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_up_eul, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Eulerian) - Actual data range
        xt_rho_eul = xt_pcolormesh(ax3, pko.pos, pko.density_ratio,
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Eulerian: Density Ratio (Data Range)')
        ax3.grid(True, alpha=0.3, which='both')
        ax3.minorticks_on()
        ax3.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_rho_eul, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Eulerian) - Auto-scaled
        xt_mat_eul = xt_pcolormesh(ax4, pko.pos, pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Eulerian: Material ID')
        ax4.grid(True, alpha=0.3, which='both')
        ax4.minorticks_on()
        ax4.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_mat_eul, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()
        
        # LAGRANGIAN X-T DIAGRAMS
        print("Creating Lagrangian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Lagrangian) - Custom colormap with zero-centered scaling
        xt_pres_lag = xt_pcolormesh(ax1, pko.pos0, pko.pres, cmap=pressure_cmap,
                                 norm=pressure_norm)
        ax1.set_xlabel('Initial Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Lagrangian: Stress (GPa, symmetric)')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.minorticks_on()
        ax1.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_pres_lag, ax=ax1, label='Stress (GPa)')
        
        # Particle velocity (Lagrangian) - Auto-scaled
        xt_up_lag = xt_pcolormesh(ax2, pko.pos0, pko.up, cmap='inferno')
        ax2.set_xlabel('Initial Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Lagrangian: Particle Velocity')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.minorticks_on()
        ax2.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_up_lag, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Lagrangian) - Actual data range
        xt_rho_lag = xt_pcolormesh(ax3, pko.pos0, pko.density_ratio,
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Initial Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Lagrangian: Density Ratio (Data Range)')
        ax3.grid(True, alpha=0.3, which='both')
        ax3.minorticks_on()
        ax3.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_rho_lag, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Lagrangian) - Auto-scaled
        xt_mat_lag = xt_pcolormesh(ax4, pko.pos0, pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Initial Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Lagrangian: Material ID')
        ax4.grid(True, alpha=0.3, which='both')
        ax4.minorticks_on()
        ax4.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_mat_lag, ax=ax4, label='Material ID')

        plt.tight_layout()
        plt.show()

        # Second figure: high-resolution zoom onto the spall region
        plot_spall_zoom(pko, pressure_cmap, pressure_norm, density_min, density_max)

        # Rarefaction fans from each free surface, traced through pyKO's sound-speed
        # field (method of characteristics); where they overlap -> spall.
        plot_rarefaction_fans(pko, pressure_cmap, run)

        # Stacked x-t diagram + free-surface velocity (related on a shared time axis)
        plot_xt_fsv_stack(pko, pressure_cmap, run)

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
        
        print(f"Using actual data ranges - Stress: [{pres_min:.2f}, {pres_max:.2f}] GPa, Density: [{density_min:.3f}, {density_max:.3f}]")
        print(f"Stress colormap: Red (tension) -> Gray (zero) -> Blue (compression)")

        # Symmetric stress normalization: equal colormap range on both sides of zero
        _stress_abs_max = max(abs(pres_min), abs(pres_max))
        pressure_norm = colors.Normalize(vmin=-_stress_abs_max, vmax=_stress_abs_max)

        # EULERIAN X-T DIAGRAMS (no spall case)
        print("\nCreating comprehensive Eulerian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Eulerian) - Custom colormap with zero-centered scaling
        xt_pres_eul = xt_pcolormesh(ax1, pko.pos, pko.pres, cmap=pressure_cmap,
                                 norm=pressure_norm)
        ax1.set_xlabel('Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Eulerian: Stress (GPa, symmetric)')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.minorticks_on()
        ax1.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_pres_eul, ax=ax1, label='Stress (GPa)')
        
        # Particle velocity (Eulerian) - Auto-scaled
        xt_up_eul = xt_pcolormesh(ax2, pko.pos, pko.up, cmap='inferno')
        ax2.set_xlabel('Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Eulerian: Particle Velocity')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.minorticks_on()
        ax2.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_up_eul, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Eulerian) - Actual data range
        xt_rho_eul = xt_pcolormesh(ax3, pko.pos, pko.density_ratio,
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Eulerian: Density Ratio (Data Range)')
        ax3.grid(True, alpha=0.3, which='both')
        ax3.minorticks_on()
        ax3.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_rho_eul, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Eulerian) - Auto-scaled
        xt_mat_eul = xt_pcolormesh(ax4, pko.pos, pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Eulerian: Material ID')
        ax4.grid(True, alpha=0.3, which='both')
        ax4.minorticks_on()
        ax4.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_mat_eul, ax=ax4, label='Material ID')
        
        plt.tight_layout()
        plt.show()
        
        # LAGRANGIAN X-T DIAGRAMS (no spall case)
        print("Creating Lagrangian x-t diagrams...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
        
        # Pressure (Lagrangian) - Custom colormap with zero-centered scaling
        xt_pres_lag = xt_pcolormesh(ax1, pko.pos0, pko.pres, cmap=pressure_cmap,
                                 norm=pressure_norm)
        ax1.set_xlabel('Initial Position (mm)')
        ax1.set_ylabel('Time (μs)')
        ax1.set_title('Lagrangian: Stress (GPa, symmetric)')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.minorticks_on()
        ax1.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_pres_lag, ax=ax1, label='Stress (GPa)')
        
        # Particle velocity (Lagrangian) - Auto-scaled
        xt_up_lag = xt_pcolormesh(ax2, pko.pos0, pko.up, cmap='inferno')
        ax2.set_xlabel('Initial Position (mm)')
        ax2.set_ylabel('Time (μs)')
        ax2.set_title('Lagrangian: Particle Velocity')
        ax2.grid(True, alpha=0.3, which='both')
        ax2.minorticks_on()
        ax2.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_up_lag, ax=ax2, label='Particle Velocity (m/s)')
        
        # Density ratio (Lagrangian) - Actual data range
        xt_rho_lag = xt_pcolormesh(ax3, pko.pos0, pko.density_ratio,
                                cmap='coolwarm', vmin=density_min, vmax=density_max)
        ax3.set_xlabel('Initial Position (mm)')
        ax3.set_ylabel('Time (μs)')
        ax3.set_title('Lagrangian: Density Ratio (Data Range)')
        ax3.grid(True, alpha=0.3, which='both')
        ax3.minorticks_on()
        ax3.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_rho_lag, ax=ax3, label=r'$\rho/\rho_0$')
        
        # Material ID (Lagrangian) - Auto-scaled
        xt_mat_lag = xt_pcolormesh(ax4, pko.pos0, pko.mat, cmap='viridis', alpha=0.7)
        ax4.set_xlabel('Initial Position (mm)')
        ax4.set_ylabel('Time (μs)')
        ax4.set_title('Lagrangian: Material ID')
        ax4.grid(True, alpha=0.3, which='both')
        ax4.minorticks_on()
        ax4.grid(True, alpha=0.1, which='minor')
        plt.colorbar(xt_mat_lag, ax=ax4, label='Material ID')

        plt.tight_layout()
        plt.show()

        # Wave-front race: even without spall, show whether the flyer release
        # overtook the shock before it reached the target free surface.
        plot_rarefaction_fans(pko, pressure_cmap, run)

        # Stacked x-t diagram + free-surface velocity (related on a shared time axis)
        plot_xt_fsv_stack(pko, pressure_cmap, run)

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
        plt.grid(True, alpha=0.3, which='both')
        plt.minorticks_on()
        plt.grid(True, alpha=0.1, which='minor')
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
            
            # Enhanced FSV analysis diagnostics
            print(f"Debug: Total FSV data points: {len(fsv_array)}")
            print(f"Debug: Peak occurs at index {max_idx} of {len(fsv_array)} points")
            print(f"Debug: Data points after peak: {len(fsv_array) - max_idx}")
            
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
                    
                    # Get material properties for spall strength calculation (target = mat 2)
                    try:
                        rho0_cu_gcm3 = run.irhostart[1]  # g/cm³ (pyKO units)
                        rho0_cu_kgm3 = rho0_cu_gcm3 * 1000  # Convert to kg/m³
                        print(f"Using {mat2_name} density from run: ρ₀ = {rho0_cu_kgm3:.0f} kg/m³ ({rho0_cu_gcm3:.3f} g/cm³)")
                    except (AttributeError, IndexError):
                        raise RuntimeError(f"Cannot access {mat2_name} density from run.irhostart[1]")
                    
                    # Get sound speed from EOS using correct pyKO units
                    try:
                        mat2_eos = run.ieos[1]  # Target EOS (material index 1)
                        if hasattr(mat2_eos, 'c0'):
                            c0_cu_cmus = mat2_eos.c0  # cm/μs (pyKO units)
                            c0_cu_ms = c0_cu_cmus * 10000  # Convert to m/s
                            print(f"Using {mat2_name} c₀ from run: {c0_cu_ms:.0f} m/s ({c0_cu_cmus:.3f} cm/μs)")
                        else:
                            raise RuntimeError("Sound speed c0 not found in EOS object")
                    except Exception as e:
                        raise RuntimeError(f"Error accessing {mat2_name} sound speed properties: {e}")
                    
                    # Calculate FSV-based spall strength using converted SI units
                    # Spall strength = 0.5 * ρ₀ * c₀ * Δu
                    spall_strength_pa = 0.5 * rho0_cu_kgm3 * c0_cu_ms * delta_u
                    spall_strength_gpa = spall_strength_pa / 1e9
                    
                    print(f"\n🎯 FSV-BASED SPALL STRENGTH CALCULATION:")
                    print(f"   Formula: σ_spall = ½ × ρ₀ × c₀ × Δu")
                    print(f"   σ_spall = 0.5 × {rho0_cu_kgm3:.0f} × {c0_cu_ms:.0f} × {delta_u:.2f}")
                    print(f"   σ_spall = {spall_strength_pa:.0f} Pa = {spall_strength_gpa:.3f} GPa")
                    
                    # Compare with YAML spall threshold
                    # cu_spall_threshold_from_yaml is set by the spall analysis block;
                    # fall back to mat2_spall_threshold_from_yaml when spall analysis is disabled.
                    if 'cu_spall_threshold_from_yaml' not in globals():
                        try:
                            cu_spall_threshold_from_yaml = mat2_spall_threshold_from_yaml
                        except NameError:
                            cu_spall_threshold_from_yaml = 0.276  # GPa conservative default
                    print(f"\n📊 COMPARISON WITH YAML THRESHOLD:")
                    print(f"   FSV-measured spall strength: {spall_strength_gpa:.3f} GPa")
                    print(f"   YAML {mat2_name} spall threshold:     {cu_spall_threshold_from_yaml:.3f} GPa")
                    
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
# MAXIMUM STRESS ANALYSIS IN TARGET (MAT 2)
########################################################################################################################

if ENABLE_STRESS_ANALYSIS:
    print(f"\n=== MAXIMUM COMPRESSIVE & TENSILE STRESS IN {mat2_name} TARGET (MAT 2) ===")

    # Track maximum compressive and tensile stresses in target (material 2) over time
    max_compressive_stress = np.zeros_like(unique_times)
    max_tensile_stress = np.zeros_like(unique_times)
    max_comp_positions = np.zeros_like(unique_times)
    max_tens_positions = np.zeros_like(unique_times)

    for i, t in enumerate(unique_times):
        snapshot = pko[pko['time'] == t]
        
        if len(snapshot) == 0:
            print(f"Warning: No data at time {t:.6f} μs - skipping stress analysis")
            continue
            
        # Filter for target (material 2) only
        mat2_target = snapshot[snapshot['mat'] == 2]
        
        if len(mat2_target) > 0:
            # Use pressure for compressive/tensile analysis (positive = compression, negative = tension)
            pressures = mat2_target['pres']
            
            # Maximum compressive stress (maximum positive pressure)
            max_comp_pressure = pressures.max()
            max_compressive_stress[i] = max_comp_pressure
            if max_comp_pressure > 0:
                max_comp_idx = pressures.idxmax()
                max_comp_positions[i] = mat2_target.loc[max_comp_idx, 'pos'] * 10  # Convert to mm
            else:
                max_comp_positions[i] = 0
            
            # Maximum tensile stress (maximum negative pressure, reported as positive magnitude)
            min_pressure = pressures.min()
            max_tensile_stress[i] = -min_pressure if min_pressure < 0 else 0
            if min_pressure < 0:
                max_tens_idx = pressures.idxmin()
                max_tens_positions[i] = mat2_target.loc[max_tens_idx, 'pos'] * 10  # Convert to mm
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
    ax1.set_title(f'Maximum Compressive Stress in {mat2_name} Target vs. Time')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.minorticks_on()
    ax1.grid(True, alpha=0.1, which='minor')

    # Calculate maximum compressive stress for reporting
    max_comp_value = np.max(max_compressive_stress)
    max_comp_time = unique_times[np.argmax(max_compressive_stress)]

    # Top right: Maximum tensile stress vs time
    ax2.plot(unique_times, max_tensile_stress, 'b-', linewidth=2)
    ax2.set_xlabel('Time (μs)')
    ax2.set_ylabel('Max Tensile Stress (GPa)')
    ax2.set_title(f'Maximum Tensile Stress in {mat2_name} Target vs. Time')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.minorticks_on()
    ax2.grid(True, alpha=0.1, which='minor')

    # Add spall threshold line for mat2 target (only if interface separation is enabled)
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
    ax3.grid(True, alpha=0.3, which='both')
    ax3.minorticks_on()
    ax3.grid(True, alpha=0.1, which='minor')

    # Add reference lines
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax3.axhline(y=1, color='black', linestyle='--', alpha=0.7)

    # Bottom right: Position of maximum tensile stress vs time
    ax4.plot(unique_times, max_tens_positions, 'b-', linewidth=2)
    ax4.set_xlabel('Time (μs)')
    ax4.set_ylabel('Position (mm)')
    ax4.set_title('Location of Maximum Tensile Stress vs. Time')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.minorticks_on()
    ax4.grid(True, alpha=0.1, which='minor')

    # Add reference lines
    ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.7)
    ax4.axhline(y=1, color='black', linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()

    print(f"Maximum compressive stress in {mat2_name} target: {max_comp_value:.2f} GPa at {max_comp_time:.3f} μs")
    print(f"Maximum tensile stress in {mat2_name} target: {max_tens_value:.2f} GPa at {max_tens_time:.3f} μs")
    if max_tens_value > cu_spall_threshold:
        print(f"⚠️  TENSILE STRESS EXCEEDS {mat2_name} SPALL THRESHOLD ({cu_spall_threshold:.3f} GPa)!")
    else:
        print(f"✅ Tensile stress remains below {mat2_name} spall threshold ({cu_spall_threshold:.3f} GPa)")

else:
    print("\n=== STRESS ANALYSIS DISABLED ===")
    max_comp_value = 0
    max_comp_time = 0
    max_tens_value = 0
    max_tens_time = 0
    # Get mat2 spall threshold from configuration if interface analysis is enabled
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
        print(f"{mat1_name}-{mat2_name} interface tracked over {len(interfaces)} time steps")
        
        # Plot interface evolution
        plt.figure(figsize=(10, 6), dpi=300)
        
        for i, interface in enumerate(interfaces):
            for j, pos in enumerate(interface['interface_positions']):
                if j == 0:  # Primary material interface
                    plt.plot(interface['time'], pos, 'ro', markersize=3)
        
        plt.xlabel('Time (μs)')
        plt.ylabel('Interface Position (mm)')
        plt.title(f'{mat1_name}-{mat2_name} Interface Evolution')
        plt.grid(True, alpha=0.3, which='both')
        plt.minorticks_on()
        plt.grid(True, alpha=0.1, which='minor')
        plt.tight_layout()
        plt.show()

else:
    print("\n=== INTERFACE ANALYSIS DISABLED ===")

########################################################################################################################
# SUMMARY REPORT
########################################################################################################################

# Ensure thickness strings are defined for summary
try:
    _ = al_thickness_str
except NameError:
    al_thickness_str = "n/a"
try:
    _ = cu_thickness_str
except NameError:
    cu_thickness_str = "n/a"

print("\n" + "="*60)
print("HYBRID SPALL + INTERFACE SEPARATION TEST SUMMARY")
print("="*60)
print(f"Simulation time: 0 to {unique_times[-1]:.2f} μs")
print(f"Total nodes: {len(pos0)}")
try:
    _sum_mat = f"{mat1_name} ({mat1_thickness_str}) -> {mat2_name} ({mat2_thickness_str})"
except NameError:
    _sum_mat = f"Impactor ({al_thickness_str}) -> Target ({cu_thickness_str})"
print(f"Materials: {_sum_mat}")
print(f"Maximum free surface velocity: {max_fsv:.2f} m/s at {max_fsv_time:.2f} μs")
print(f"Maximum compressive stress in {mat2_name}: {max_comp_value:.2f} GPa at {max_comp_time:.3f} μs")
print(f"Maximum tensile stress in {mat2_name}: {max_tens_value:.2f} GPa at {max_tens_time:.3f} μs")
if max_tens_value > cu_spall_threshold:
    print(f"⚠️  Tensile stress exceeds {mat2_name} spall threshold ({cu_spall_threshold:.3f} GPa)!")
else:
    print(f"✅ Tensile stress below {mat2_name} spall threshold ({cu_spall_threshold:.3f} GPa)")

# Add FSV-based spall strength to summary
if 'fsv_measurement_available' in locals() and fsv_measurement_available:
    print(f"🎯 FSV-measured spall strength: {fsv_spall_strength:.3f} GPa")
    ratio = (
        fsv_spall_strength / cu_spall_threshold
        if cu_spall_threshold is not None and cu_spall_threshold > 1e-30
        else float('nan')
    )
    if ratio != ratio:  # NaN
        print("📏 FSV vs YAML threshold ratio: N/A (mat2 spall threshold unset or ~0)")
    elif abs(ratio - 1.0) < 0.3:
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
            try:
                material_names = {
                    1: f'{mat1_name} ({mat1_thickness_str})',
                    2: f'{mat2_name} ({mat2_thickness_str})',
                }
            except NameError:
                material_names = {
                    1: f'Impactor ({al_thickness_str})',
                    2: f'Target ({cu_thickness_str})',
                }
            print(f"    - Material {mat_id} ({material_names.get(mat_id, 'Unknown')}): {count} spall events")
    
    if pressure_spall_detected:
        print(f"  Pressure-based spall: Max tensile stress {max_tensile_pressure:.2f} GPa exceeded thresholds")
        if max_tensile_pressure > al_spall_threshold_from_yaml:
            print(f"    - {mat1_name} threshold exceeded: {max_tensile_pressure:.2f} > {al_spall_threshold_from_yaml:.2f} GPa")
        if max_tensile_pressure > cu_spall_threshold_from_yaml:
            print(f"    - {mat2_name} threshold exceeded: {max_tensile_pressure:.2f} > {cu_spall_threshold_from_yaml:.2f} GPa")
else:
    print("Spall detected: NO")

print("="*60)

# %%
