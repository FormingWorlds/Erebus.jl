#!/usr/bin/env python3
"""
Generate 6-panel summary benchmark figure for Core Formation and Planetesimal
Differentiation in Erebus.jl.
"""

import os
import sys
import json
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap, BoundaryNorm

# Interra visual palette
STRATA = {
    'gold': '#F0BA5E',
    'amber': '#DE7037',
    'magma': '#A23E3C',
    'plum': '#6A2A4D',
    'cobalt': '#3E5B86',
    'ink': '#15101C'
}
NEUTRALS = {
    'paper': '#FAF5E6',
    'mist': '#C8BCA9',
    'graphite': '#3A3140',
    'bone': '#E6DAB6'
}

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'mathtext.fontset': 'stixsans',
    'axes.edgecolor': NEUTRALS['mist'],
    'axes.linewidth': 1.0,
    'axes.labelcolor': NEUTRALS['graphite'],
    'axes.facecolor': '#FFFFFF',
    'xtick.color': NEUTRALS['graphite'],
    'ytick.color': NEUTRALS['graphite'],
    'grid.color': NEUTRALS['mist'],
    'grid.linestyle': ':',
    'grid.linewidth': 0.8,
    'legend.frameon': True,
    'legend.facecolor': '#FFFFFF',
    'legend.edgecolor': NEUTRALS['mist'],
    'figure.facecolor': '#FFFFFF',
})

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "output_files", "core_formation_benchmark_data.json")
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")

R_PLANET_KM = 50.0
BOX_SIZE_KM = 65.0
N_GRID = 129


def classify_differentiation_regime(T_2d, phi_fe_2d, Fm_2d, phi_ice_2d=None):
    """
    Classify internal state into physical regimes:
    0: Primordial cold icy mix (T < 273.15 K or phi_ice > 0.05)
    1: Dehydrated crystalline rock (273.15 <= T < 1213 K)
    2: Molten Fe-FeS percolation mush (1213 <= T < 1416 K, phi_fe > 0.02)
    3: Silicate magma ocean suspension (T >= 1416 K or Fm >= 0.40)
    4: Segregated metallic iron core (phi_fe >= 0.50)
    """
    regime = np.zeros_like(T_2d)
    if phi_ice_2d is not None:
        regime[(T_2d < 273.15) | (phi_ice_2d > 0.05)] = 0
    else:
        regime[T_2d < 273.15] = 0
    regime[(T_2d >= 273.15) & (T_2d < 1213.0)] = 1
    regime[(T_2d >= 1213.0) & (T_2d < 1416.0)] = 2
    regime[(T_2d >= 1416.0) | (Fm_2d >= 0.40)] = 3
    regime[phi_fe_2d >= 0.50] = 4
    return regime


def generate_benchmark_figure(data_path=DATA_PATH, assets_dir=ASSETS_DIR, output_dir=OUTPUT_FILES_DIR):
    """Generate 6-panel core formation benchmark figure."""
    os.makedirs(assets_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data file {data_path} not found.")

    with open(data_path, "r") as f:
        data = json.load(f)

    ref_data = data["reference"]

    x_coords = np.linspace(-BOX_SIZE_KM, BOX_SIZE_KM, N_GRID)
    y_coords = np.linspace(-BOX_SIZE_KM, BOX_SIZE_KM, N_GRID)
    X, Y = np.meshgrid(x_coords, y_coords)
    R = np.sqrt(X**2 + Y**2)

    r_1d = np.array(ref_data["r"]) / 1000.0
    snapshots_T = ref_data["snapshot_T"]
    snapshots_phi = ref_data["snapshot_phi_fe"]
    snapshots_phi_ice = ref_data.get("snapshot_phi_ice", {})
    snapshots_Fm = ref_data["snapshot_Fm"]

    snap_keys = sorted(snapshots_T.keys(), key=lambda k: float(k))

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.6), dpi=200, constrained_layout=True)

    t_ref_myr = np.array(ref_data["times_yr"]) / 1.0e6
    t_noheat_myr = np.array(data["no_heating"]["times_yr"]) / 1.0e6
    r_km = r_1d

    # Panel (a): 2D Snapshot of Fully Differentiated Body at t = 3.0 Ma
    ax_a = axes[0, 0]
    ax_a.text(0.04, 0.93, '(a)', transform=ax_a.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    snap_3ma_key = min(snap_keys, key=lambda k: abs(float(k) - 3.0e6))
    T_3ma_1d = np.array(snapshots_T[snap_3ma_key])
    phi_3ma_1d = np.array(snapshots_phi[snap_3ma_key])
    Fm_3ma_1d = np.array(snapshots_Fm[snap_3ma_key])
    phi_ice_3ma_1d = np.array(snapshots_phi_ice.get(snap_3ma_key, np.zeros_like(T_3ma_1d)))

    T_3ma_2d = np.interp(R, r_1d, T_3ma_1d, right=150.0)
    phi_3ma_2d = np.interp(R, r_1d, phi_3ma_1d, right=0.0)
    Fm_3ma_2d = np.interp(R, r_1d, Fm_3ma_1d, right=0.0)
    phi_ice_3ma_2d = np.interp(R, r_1d, phi_ice_3ma_1d, right=0.0)

    regime_3ma = classify_differentiation_regime(T_3ma_2d, phi_3ma_2d, Fm_3ma_2d, phi_ice_3ma_2d)
    regime_3ma[R > R_PLANET_KM] = np.nan

    regime_colors = ['#88B7D5', '#9E8B7D', STRATA['amber'], STRATA['magma'], STRATA['gold']]
    regime_cmap = ListedColormap(regime_colors)
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], 5)
    cf_a = ax_a.imshow(regime_3ma, extent=[-BOX_SIZE_KM, BOX_SIZE_KM, -BOX_SIZE_KM, BOX_SIZE_KM],
                       origin='lower', cmap=regime_cmap, norm=norm, interpolation='nearest')
    ax_a.add_patch(Circle((0.0, 0.0), R_PLANET_KM, fill=False, edgecolor=NEUTRALS['graphite'], lw=1.2))
    ax_a.set_aspect('equal')
    ax_a.set_title(r'Differentiated Planetesimal ($t = 3.0$ Ma)', fontsize=11, fontweight='bold')
    ax_a.set_xlabel('x [km]', fontsize=9.5)
    ax_a.set_ylabel('y [km]', fontsize=9.5)
    cb_a = fig.colorbar(cf_a, ax=ax_a, ticks=[0, 1, 2, 3, 4], fraction=0.046, pad=0.04)
    cb_a.ax.set_yticklabels(['Primordial Ice', 'Rock', 'Percolation', 'Magma Ocean', 'Metallic Core'], fontsize=7.5)

    # Panel (b): Central Temperature & Dissipation Heating Impact
    ax_b = axes[0, 1]
    ax_b.text(0.04, 0.93, '(b)', transform=ax_b.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax_b.plot(t_ref_myr, ref_data["T_core_hist"], '-', color=STRATA['magma'], linewidth=2.2, label=r'Coupled ($Q_\mathrm{seg}$ ON)')
    ax_b.plot(t_noheat_myr, data["no_heating"]["T_core_hist"], '--', color=STRATA['cobalt'], linewidth=2.0, label=r'No seg heating ($Q_\mathrm{seg}$ OFF)')
    ax_b.axhline(1800.0, color=STRATA['gold'], linestyle=':', linewidth=1.4, label=r'Silicate Liquidus ($1800$ K)')
    ax_b.axhline(1416.0, color=STRATA['ink'], linestyle='--', linewidth=1.4, label=r'Silicate Solidus ($1416$ K)')
    ax_b.axhline(1213.0, color=STRATA['amber'], linestyle='-.', linewidth=1.4, label=r'Fe-FeS Eutectic ($1213$ K)')
    ax_b.axhline(273.15, color='cyan', linestyle='--', linewidth=1.2, label=r'Ice Melting ($273.15$ K)')
    ax_b.set_xlim(0, 3.5)
    ax_b.set_ylim(100, 2500)
    ax_b.set_xlabel('Time [Ma]', fontsize=9.5)
    ax_b.set_ylabel(r'Central Temperature $T_\mathrm{core}$ [K]', fontsize=9.5)
    ax_b.set_title('Core Thermal Runaway and Dissipation', fontsize=11, fontweight='bold')
    ax_b.grid(True)
    ax_b.legend(loc='lower right', fontsize=8.0)

    # Panel (c): Core Formation & Magma Ocean Front Growth
    ax_c = axes[0, 2]
    ax_c.text(0.04, 0.93, '(c)', transform=ax_c.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax_c.plot(t_ref_myr, ref_data["R_core_hist"], '-', color=STRATA['gold'], linewidth=2.4, label=r'Metallic Core Radius ($\phi_\mathrm{fe} \geq 0.50$)')
    ax_c.plot(t_ref_myr, ref_data["R_magma_hist"], '-', color=STRATA['magma'], linewidth=2.0, label=r'Magma Ocean Boundary ($F_m \geq 0.40$)')
    ax_c.axhline(R_PLANET_KM, color=NEUTRALS['graphite'], linestyle=':', linewidth=1.2, label='Planet Surface (50 km)')
    ax_c.set_xlim(0, 3.5)
    ax_c.set_ylim(0, 55)
    ax_c.set_xlabel('Time [Ma]', fontsize=9.5)
    ax_c.set_ylabel('Radius [km]', fontsize=9.5)
    ax_c.set_title('Differentiation Fronts Timeline', fontsize=11, fontweight='bold')
    ax_c.grid(True)
    ax_c.legend(loc='lower right', fontsize=8.0)

    # Panel (d): Radial Metal Concentration Evolution phi_fe(r)
    ax_d = axes[1, 0]
    ax_d.text(0.04, 0.93, '(d)', transform=ax_d.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    plot_times = [0.0, 1.0e6, 1.5e6, 3.0e6]
    time_colors = [STRATA['cobalt'], STRATA['amber'], STRATA['magma'], STRATA['gold']]
    for pt, col in zip(plot_times, time_colors):
        k_near = min(snap_keys, key=lambda k: abs(float(k) - pt))
        phi_prof = np.array(snapshots_phi[k_near])
        ax_d.plot(r_km, phi_prof, '-', color=col, linewidth=2.0, label=f'$t = {float(k_near)/1e6:.1f}$ Ma')
    ax_d.axhline(0.12, color=NEUTRALS['graphite'], linestyle=':', linewidth=1.2, label=r'Initial uniform $\phi_0 = 0.12$')
    ax_d.axhline(0.65, color=STRATA['gold'], linestyle='--', linewidth=1.2, label=r'Packing ceiling $\phi_\mathrm{pack} = 0.65$')
    ax_d.set_xlim(0, 50.0)
    ax_d.set_ylim(0, 0.72)
    ax_d.set_xlabel('Radius $r$ [km]', fontsize=9.5)
    ax_d.set_ylabel(r'Metal Volume Fraction $\phi_\mathrm{fe}$ [-]', fontsize=9.5)
    ax_d.set_title('Mantle Depletion & Core Ponding Profiles', fontsize=11, fontweight='bold')
    ax_d.grid(True)
    ax_d.legend(loc='center right', fontsize=8.0)

    # Panel (e): Segregation Transport Regimes
    ax_e = axes[1, 1]
    ax_e.text(0.04, 0.93, '(e)', transform=ax_e.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    t_perc_myr = np.array(data["perc_only"]["times_yr"]) / 1.0e6
    t_settle_myr = np.array(data["settle_only"]["times_yr"]) / 1.0e6
    ax_e.plot(t_ref_myr, np.maximum(ref_data["v_seg_peak_hist"], 1.0e-12), '-', color=STRATA['magma'], linewidth=2.2, label='Coupled Hermite Transition')
    ax_e.plot(t_perc_myr, np.maximum(data["perc_only"]["v_seg_peak_hist"], 1.0e-12), '--', color=STRATA['amber'], linewidth=2.0, label='Porous Percolation only')
    ax_e.plot(t_settle_myr, np.maximum(data["settle_only"]["v_seg_peak_hist"], 1.0e-12), ':', color=STRATA['cobalt'], linewidth=2.0, label='Stokes Settling only')
    ax_e.set_yscale('log')
    ax_e.set_xlim(0.5, 3.5)
    ax_e.set_ylim(1.0e-10, 1.0e-1)
    ax_e.set_xlabel('Time [Ma]', fontsize=9.5)
    ax_e.set_ylabel(r'Peak Segregation Velocity $v_\mathrm{seg}$ [m/s]', fontsize=9.5)
    ax_e.set_title('Velocity Dynamics by Transport Regime', fontsize=11, fontweight='bold')
    ax_e.grid(True)
    ax_e.legend(loc='center right', fontsize=8.0)

    # Panel (f): Droplet Size Physics Sensitivity
    ax_f = axes[1, 2]
    ax_f.text(0.04, 0.93, '(f)', transform=ax_f.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    t_fixed_myr = np.array(data["drop_fixed"]["times_yr"]) / 1.0e6
    t_turb_myr = np.array(data["drop_turb"]["times_yr"]) / 1.0e6
    ax_f.plot(t_ref_myr, ref_data["R_core_hist"], '-', color=STRATA['gold'], linewidth=2.2, label=r'Weber Equilibrium Mean ($d \propto \sqrt{\sigma/g}$)')
    ax_f.plot(t_fixed_myr, data["drop_fixed"]["R_core_hist"], '--', color=STRATA['cobalt'], linewidth=2.0, label=r'Fixed Droplet ($d = 1.0$ cm)')
    ax_f.plot(t_turb_myr, data["drop_turb"]["R_core_hist"], ':', color=STRATA['plum'], linewidth=2.0, label=r'Turbulent Dynamic Breakup ($d \propto v^{-2}$)')
    ax_f.set_xlim(0.8, 3.5)
    ax_f.set_ylim(0, 30)
    ax_f.set_xlabel('Time [Ma]', fontsize=9.5)
    ax_f.set_ylabel('Core Radius [km]', fontsize=9.5)
    ax_f.set_title('Core Growth Sensitivity to Droplet Physics', fontsize=11, fontweight='bold')
    ax_f.grid(True)
    ax_f.legend(loc='lower right', fontsize=8.0)

    fig.savefig(os.path.join(assets_dir, "core_formation_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(assets_dir, "core_formation_benchmark.pdf"))
    fig.savefig(os.path.join(output_dir, "core_formation_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(output_dir, "core_formation_benchmark.pdf"))
    plt.close(fig)
    print(f"Saved benchmark summary figure to {assets_dir}/core_formation_benchmark.png")


def main():
    """Main entry point for generating core formation benchmark figures."""
    generate_benchmark_figure()


if __name__ == "__main__":
    main()
