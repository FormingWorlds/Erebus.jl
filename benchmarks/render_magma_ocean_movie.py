#!/usr/bin/env python3
"""
Render 2D movie frames and generate benchmark assets for the Planetesimal Magma Ocean
cooling model in Erebus.jl.
"""

import os
import sys
import json
import subprocess
from multiprocessing import Pool
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import Normalize, LogNorm

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

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "output_files", "soft_turbulence_benchmark_data.json")
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")
FRAME_DIR = "/tmp/erebus_magma_frames"

os.makedirs(ASSETS_DIR, exist_ok=True)
os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)
os.makedirs(FRAME_DIR, exist_ok=True)

with open(DATA_PATH, "r") as f:
    data = json.load(f)

# 2D Grid Setup
R_PLANET_KM = 50.0
BOX_SIZE_KM = 70.0
N_GRID = 129
x_coords = np.linspace(-BOX_SIZE_KM, BOX_SIZE_KM, N_GRID)
y_coords = np.linspace(-BOX_SIZE_KM, BOX_SIZE_KM, N_GRID)
X, Y = np.meshgrid(x_coords, y_coords)
R = np.sqrt(X**2 + Y**2)

r_1d = np.array(data["res_on"]["r"]) / 1000.0  # km
snapshots_T = data["res_on"]["snapshot_T"]
snapshots_Fm = data["res_on"]["snapshot_Fm"]
snapshots_k = data["res_on"]["snapshot_k"]

snap_keys = sorted(snapshots_T.keys(), key=lambda k: float(k))

def render_single_frame(task):
    idx, key = task
    t_yr = float(key)
    t_kyr = t_yr / 1000.0
    out_png = os.path.join(FRAME_DIR, f"frame_{idx:04d}.png")

    T_1d = np.array(snapshots_T[key])
    Fm_1d = np.array(snapshots_Fm[key])
    k_1d = np.array(snapshots_k[key])

    # Interpolate 1D spherical profile onto 2D Cartesian grid
    T_2d = np.interp(R, r_1d, T_1d, right=300.0)
    Fm_2d = np.interp(R, r_1d, Fm_1d, right=0.0)
    k_2d = np.interp(R, r_1d, k_1d, right=3.0)

    # Outside planetesimal: mask with NaN
    mask_outside = R > R_PLANET_KM
    T_2d[mask_outside] = np.nan
    Fm_2d[mask_outside] = np.nan
    k_2d[mask_outside] = np.nan

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0), dpi=120, constrained_layout=True)
    fig.suptitle(f"Erebus 2D Planetesimal Magma Ocean Benchmark (128x128)  |  Time = {t_kyr:5.1f} kyr",
                 fontsize=14, fontweight='bold', color=NEUTRALS['graphite'])

    def format_ax(ax, title, has_ylabel=True):
        ax.set_facecolor('white')
        ax.add_patch(Circle((0.0, 0.0), R_PLANET_KM, fill=False, edgecolor=NEUTRALS['graphite'], lw=1.2, ls='-'))
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=12, fontweight='bold', pad=8)
        ax.set_xlabel("x [km]", fontsize=10)
        if has_ylabel:
            ax.set_ylabel("y [km]", fontsize=10)
        ax.set_xlim(-BOX_SIZE_KM, BOX_SIZE_KM)
        ax.set_ylim(-BOX_SIZE_KM, BOX_SIZE_KM)

    # Panel 1: Temperature T
    ax = axes[0]
    format_ax(ax, "Temperature $T$", has_ylabel=True)
    t_levels = np.linspace(200.0, 1900.0, 35)
    cf1 = ax.contourf(X, Y, T_2d, levels=t_levels, cmap='viridis', extend='both')
    ax.contour(X, Y, T_2d, levels=[1400.0, 1800.0], colors=['white', STRATA['gold']], linewidths=1.2, linestyles=['--', ':'])
    fig.colorbar(cf1, ax=ax, ticks=[300, 600, 900, 1200, 1500, 1800], label="Temperature [K]", fraction=0.046, pad=0.04)

    # Panel 2: Melt Fraction Fm
    ax = axes[1]
    format_ax(ax, "Silicate Melt Fraction $F_m$", has_ylabel=False)
    fm_levels = np.linspace(0.0, 1.0, 21)
    cf2 = ax.contourf(X, Y, Fm_2d, levels=fm_levels, cmap='magma', vmin=0.0, vmax=1.0)
    ax.contour(X, Y, Fm_2d, levels=[0.30, 0.40, 0.50], colors=['white', STRATA['gold'], 'cyan'], linewidths=1.2, linestyles=[':', '-', ':'])
    fig.colorbar(cf2, ax=ax, ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], label="Melt Fraction [-]", fraction=0.046, pad=0.04)

    # Panel 3: Effective Thermal Conductivity keff
    ax = axes[2]
    format_ax(ax, "Effective Conductivity $k_\mathrm{eff}$", has_ylabel=False)
    log_k_levels = np.logspace(0, 5, 26)
    cf3 = ax.contourf(X, Y, k_2d, levels=log_k_levels, norm=LogNorm(vmin=1.0, vmax=1.0e5), cmap='inferno')
    fig.colorbar(cf3, ax=ax, ticks=[1e0, 1e1, 1e2, 1e3, 1e4, 1e5], label=r"$k_\mathrm{eff}$ [$\mathrm{W/(m\cdot K)}$]", fraction=0.046, pad=0.04)

    fig.savefig(out_png)
    plt.close(fig)
    return idx

def generate_video():
    print(f"Rendering {len(snap_keys)} frames to {FRAME_DIR}...")
    tasks = list(enumerate(snap_keys))
    with Pool() as pool:
        for idx in pool.imap_unordered(render_single_frame, tasks):
            if idx % 20 == 0:
                print(f"  Rendered frame {idx}/{len(snap_keys)}")

    mp4_out = os.path.join(ASSETS_DIR, "magma_ocean_cooling_128.mp4")
    gif_out = os.path.join(ASSETS_DIR, "magma_ocean_cooling_128.gif")

    env = dict(os.environ)
    env["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/Cellar/x265/4.2/lib"

    print("Encoding MP4 video...")
    cmd_mp4 = [
        "ffmpeg", "-y", "-framerate", "20",
        "-i", os.path.join(FRAME_DIR, "frame_%04d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
        mp4_out
    ]
    subprocess.run(cmd_mp4, check=True, env=env)
    print(f"Saved {mp4_out}")

    print("Encoding GIF animation...")
    palette_png = "/tmp/erebus_magma_palette.png"
    cmd_pal = [
        "ffmpeg", "-y", "-i", os.path.join(FRAME_DIR, "frame_%04d.png"),
        "-vf", "fps=15,scale=960:-1:flags=lanczos,palettegen=stats_mode=diff",
        palette_png
    ]
    subprocess.run(cmd_pal, check=True, env=env)

    cmd_gif = [
        "ffmpeg", "-y", "-framerate", "15",
        "-i", os.path.join(FRAME_DIR, "frame_%04d.png"),
        "-i", palette_png,
        "-lavfi", "fps=15,scale=960:-1:flags=lanczos [x]; [x][1:v] paletteuse=dither=bayer:bayer_scale=3",
        gif_out
    ]
    subprocess.run(cmd_gif, check=True, env=env)
    print(f"Saved {gif_out}")

def generate_summary_figure():
    print("Generating 6-panel summary figure...")
    fig, axes = plt.subplots(2, 3, figsize=(16, 10.0), constrained_layout=True)

    res_off = data["res_off"]
    res_on = data["res_on"]

    r_km = np.array(res_off["r"]) / 1000.0
    t_off_kyr = np.array(res_off["times_yr"]) / 1000.0
    t_on_kyr = np.array(res_on["times_yr"]) / 1000.0

    # Panel (a): 2D Snapshot of Magma Ocean at t = 15 kyr
    ax_a = axes[0, 0]
    ax_a.text(0.04, 0.93, '(a)', transform=ax_a.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    T_15_1d = np.array(res_on["snapshot_T"]["15000.0"])
    T_15_2d = np.interp(R, r_1d, T_15_1d, right=300.0)
    T_15_2d[R > R_PLANET_KM] = np.nan
    cf_a = ax_a.contourf(X, Y, T_15_2d, levels=np.linspace(200, 1900, 35), cmap='viridis')
    ax_a.contour(X, Y, T_15_2d, levels=[1400.0, 1800.0], colors=['white', STRATA['gold']], linewidths=1.2, linestyles=['--', ':'])
    ax_a.add_patch(Circle((0.0, 0.0), R_PLANET_KM, fill=False, edgecolor=NEUTRALS['graphite'], lw=1.2))
    ax_a.set_aspect('equal')
    ax_a.set_title(r'2D Thermal Field at $t = 15$ kyr', fontsize=12, fontweight='bold')
    ax_a.set_xlabel('x [km]')
    ax_a.set_ylabel('y [km]')
    fig.colorbar(cf_a, ax=ax_a, fraction=0.046, pad=0.04, label='Temperature [K]')

    # Panel (b): Central Temperature Evolution T_core(t)
    ax_b = axes[0, 1]
    ax_b.text(0.04, 0.93, '(b)', transform=ax_b.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax_b.plot(t_off_kyr, res_off["T_core_hist"], '-', color=STRATA['cobalt'], linewidth=2.2, label='Conduction only (Soft Turb. OFF)')
    ax_b.plot(t_on_kyr, res_on["T_core_hist"], '-', color=STRATA['magma'], linewidth=2.2, label='Regularized Soft Turb. ON')
    ax_b.axhline(1800.0, color=STRATA['gold'], linestyle=':', linewidth=1.5, label='Liquidus $T_l = 1800$ K')
    ax_b.axhline(1400.0, color=STRATA['ink'], linestyle='--', linewidth=1.5, label='Solidus $T_s = 1400$ K')
    ax_b.set_xlim(0, 50.0)
    ax_b.set_ylim(1300, 1900)
    ax_b.set_xlabel('Time [kyr]')
    ax_b.set_ylabel(r'Central Temperature $T_\mathrm{core}$ [K]')
    ax_b.set_title('Core Thermal Quenching', fontsize=12, fontweight='bold')
    ax_b.grid(True)
    ax_b.legend(loc='lower left', fontsize=8.5)

    # Panel (c): Solidification Front Retreat r_melt(t)
    ax_c = axes[0, 2]
    ax_c.text(0.04, 0.93, '(c)', transform=ax_c.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax_c.plot(t_off_kyr, np.array(res_off["melt_radius_hist"]) / 1000.0, '-', color=STRATA['cobalt'], linewidth=2.2, label='Conduction only (Soft Turb. OFF)')
    ax_c.plot(t_on_kyr, np.array(res_on["melt_radius_hist"]) / 1000.0, '-', color=STRATA['magma'], linewidth=2.2, label='Regularized Soft Turb. ON')
    ax_c.set_xlim(0, 50.0)
    ax_c.set_ylim(0, 35.0)
    ax_c.set_xlabel('Time [kyr]')
    ax_c.set_ylabel(r'Magma Ocean Radius ($F_m \geq 0.40$) [km]')
    ax_c.set_title('Magma Ocean Crystallization Front', fontsize=12, fontweight='bold')
    ax_c.grid(True)
    ax_c.legend(loc='upper right', fontsize=8.5)

    # Panel (d): Radial Temperature Profiles T(r)
    ax_d = axes[1, 0]
    ax_d.text(0.04, 0.93, '(d)', transform=ax_d.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    times_to_plot = ["0.0", "5000.0", "15000.0", "50000.0"]
    time_colors = [STRATA['gold'], STRATA['amber'], STRATA['magma'], STRATA['plum']]
    for t_key, col in zip(times_to_plot, time_colors):
        t_val = float(t_key) / 1000.0
        ax_d.plot(r_km, res_off["snapshot_T"][t_key], ':', color=col, linewidth=1.6)
        ax_d.plot(r_km, res_on["snapshot_T"][t_key], '-', color=col, linewidth=2.2, label=f'$t = {t_val:.0f}$ kyr')
    ax_d.plot([], [], '-', color=NEUTRALS['graphite'], linewidth=2.0, label='Soft Turb. ON')
    ax_d.plot([], [], ':', color=NEUTRALS['graphite'], linewidth=1.6, label='Soft Turb. OFF')
    ax_d.set_xlim(0, 50)
    ax_d.set_ylim(200, 1900)
    ax_d.set_xlabel('Radius $r$ [km]')
    ax_d.set_ylabel('Temperature $T(r)$ [K]')
    ax_d.set_title('Radial Temperature Profiles', fontsize=12, fontweight='bold')
    ax_d.grid(True)
    ax_d.legend(loc='lower left', fontsize=8.0)

    # Panel (e): Effective Thermal Conductivity keff(r)
    ax_e = axes[1, 1]
    ax_e.text(0.04, 0.93, '(e)', transform=ax_e.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    for t_key, col in zip(["0.0", "5000.0", "15000.0"], [STRATA['gold'], STRATA['amber'], STRATA['magma']]):
        t_val = float(t_key) / 1000.0
        ax_e.plot(r_km, res_on["snapshot_k"][t_key], '-', color=col, linewidth=2.2, label=f'$t = {t_val:.0f}$ kyr')
    ax_e.axhline(3.0, color=NEUTRALS['graphite'], linestyle=':', linewidth=1.5, label=r'$k_\mathrm{cond} = 3.0\ \mathrm{W/(m\cdot K)}$')
    ax_e.set_yscale('log')
    ax_e.set_xlim(0, 50)
    ax_e.set_ylim(1.0, 5e5)
    ax_e.set_xlabel('Radius $r$ [km]')
    ax_e.set_ylabel(r'$k_\mathrm{eff}$ [$\mathrm{W/(m\cdot K)}$]')
    ax_e.set_title('Convective Conductivity Profiles', fontsize=12, fontweight='bold')
    ax_e.grid(True)
    ax_e.legend(loc='upper right', fontsize=8.5)

    # Panel (f): Surface Heat Flux Evolution q_surf(t)
    ax_f = axes[1, 2]
    ax_f.text(0.04, 0.93, '(f)', transform=ax_f.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    q_off = np.array(res_off["q_surf_hist"])
    q_on = np.array(res_on["q_surf_hist"])
    ax_f.plot(t_off_kyr, q_off, '-', color=STRATA['cobalt'], linewidth=2.2, label='Conduction only (Soft Turb. OFF)')
    ax_f.plot(t_on_kyr, q_on, '-', color=STRATA['magma'], linewidth=2.2, label='Regularized Soft Turb. ON')
    ax_f.set_xlim(0, 50.0)
    ax_f.set_ylim(0, max(np.max(q_on), np.max(q_off)) * 1.1)
    ax_f.set_xlabel('Time [kyr]')
    ax_f.set_ylabel(r'Surface Heat Flux $q_\mathrm{surf}$ [$\mathrm{W/m^2}$]')
    ax_f.set_title('Planetary Heat Loss', fontsize=12, fontweight='bold')
    ax_f.grid(True)
    ax_f.legend(loc='upper right', fontsize=8.5)

    summary_png_assets = os.path.join(ASSETS_DIR, "magma_ocean_cooling_benchmark.png")
    summary_png_output = os.path.join(OUTPUT_FILES_DIR, "magma_ocean_cooling_benchmark.png")
    summary_pdf_output = os.path.join(OUTPUT_FILES_DIR, "magma_ocean_cooling_benchmark.pdf")
    fig.savefig(summary_png_assets, dpi=200)
    fig.savefig(summary_png_output, dpi=300)
    fig.savefig(summary_pdf_output)
    plt.close(fig)
    print(f"Saved {summary_png_assets}")

def generate_convergence_figure():
    print("Generating resolution convergence figure across 32, 64, 128, 256...")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    res_32 = data["res_grid32"]
    res_64 = data["res_grid64"]
    res_128 = data["res_grid128"]
    res_256 = data["res_grid256"]

    r_32_km = np.array(res_32["r"]) / 1000.0
    r_64_km = np.array(res_64["r"]) / 1000.0
    r_128_km = np.array(res_128["r"]) / 1000.0
    r_256_km = np.array(res_256["r"]) / 1000.0

    # Panel (a): Radial Temperature at t = 15 kyr across Resolutions
    ax1.text(0.04, 0.93, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax1.plot(r_32_km, res_32["snapshot_T"]["15000.0"], ':', color=STRATA['cobalt'], linewidth=2.0, label=r'$N_r = 32$ ($\Delta r = 1.56$ km)')
    ax1.plot(r_64_km, res_64["snapshot_T"]["15000.0"], '--', color=STRATA['gold'], linewidth=2.0, label=r'$N_r = 64$ ($\Delta r = 0.78$ km)')
    ax1.plot(r_128_km, res_128["snapshot_T"]["15000.0"], '-', color=STRATA['amber'], linewidth=2.0, label=r'$N_r = 128$ ($\Delta r = 0.39$ km)')
    ax1.plot(r_256_km, res_256["snapshot_T"]["15000.0"], '-', color=STRATA['magma'], linewidth=1.5, label=r'$N_r = 256$ ($\Delta r = 0.20$ km)')
    ax1.set_xlim(0, 50)
    ax1.set_ylim(250, 1900)
    ax1.set_xlabel('Planetesimal Radius $r$ [km]')
    ax1.set_ylabel('Temperature $T(r)$ [K] at $t = 15$ kyr')
    ax1.set_title('Spatial Convergence: $T(r)$ at 15 kyr', fontsize=12, fontweight='bold')
    ax1.grid(True)
    ax1.legend(loc='lower left', fontsize=8.5)

    # Panel (b): Melt Fraction Profile Fm(r) at t = 15 kyr
    ax2.text(0.04, 0.93, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax2.plot(r_32_km, res_32["snapshot_Fm"]["15000.0"], ':', color=STRATA['cobalt'], linewidth=2.0, label=r'$N_r = 32$')
    ax2.plot(r_64_km, res_64["snapshot_Fm"]["15000.0"], '--', color=STRATA['gold'], linewidth=2.0, label=r'$N_r = 64$')
    ax2.plot(r_128_km, res_128["snapshot_Fm"]["15000.0"], '-', color=STRATA['amber'], linewidth=2.0, label=r'$N_r = 128$')
    ax2.plot(r_256_km, res_256["snapshot_Fm"]["15000.0"], '-', color=STRATA['magma'], linewidth=1.5, label=r'$N_r = 256$')
    ax2.axhline(0.40, color=NEUTRALS['graphite'], linestyle=':', linewidth=1.2, label=r'$\phi_\mathrm{crit} = 0.40$')
    ax2.set_xlim(0, 40)
    ax2.set_ylim(-0.05, 1.05)
    ax2.set_xlabel('Planetesimal Radius $r$ [km]')
    ax2.set_ylabel('Melt Fraction $F_m(r)$ at $t = 15$ kyr')
    ax2.set_title('Spatial Convergence: Melt Profile $F_m(r)$', fontsize=12, fontweight='bold')
    ax2.grid(True)
    ax2.legend(loc='lower left', fontsize=8.5)

    # Panel (c): Convergence Error vs Cell Width Delta r
    ax3.text(0.04, 0.93, '(c)', transform=ax3.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    t_256_core = res_256["snapshot_T"]["15000.0"][0]
    dr_vals = [50.0 / 32, 50.0 / 64, 50.0 / 128]
    err_vals = [
        abs(res_32["snapshot_T"]["15000.0"][0] - t_256_core) / t_256_core,
        abs(res_64["snapshot_T"]["15000.0"][0] - t_256_core) / t_256_core,
        abs(res_128["snapshot_T"]["15000.0"][0] - t_256_core) / t_256_core,
    ]
    ax3.loglog(dr_vals, err_vals, 'o-', color=STRATA['magma'], linewidth=2.2, markersize=7, label=r'Core $T$ relative error vs $N_r = 256$')
    dr_ref = np.array([0.3, 1.6])
    ax3.loglog(dr_ref, err_vals[-1] * (dr_ref / dr_vals[-1])**2, 'k--', linewidth=1.5, label=r'Second-order $\mathcal{O}(\Delta r^2)$ slope')
    ax3.set_xlabel(r'Grid Cell Size $\Delta r$ [km]')
    ax3.set_ylabel('Relative Error [-]')
    ax3.set_title('Numerical Order of Convergence', fontsize=12, fontweight='bold')
    ax3.grid(True, which="both")
    ax3.legend(loc='lower right', fontsize=8.5)

    conv_png_assets = os.path.join(ASSETS_DIR, "magma_ocean_grid_convergence.png")
    fig.savefig(conv_png_assets, dpi=200)
    plt.close(fig)
    print(f"Saved {conv_png_assets}")

def generate_regularization_figure():
    print("Generating conductivity regularization comparison figure...")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    Fm = np.array(data["Fm_vals"])
    dT = np.array(data["dT_vals"])
    w_T = np.array(data["w_T_vals"])

    # Panel (a): keff vs Fm
    ax1.text(0.04, 0.93, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax1.plot(Fm, data["k_i2elvis_curves"]["100"], '--', color=NEUTRALS['mist'], linewidth=2.0,
             label=r'Raw i2elvis ($F_m=0.40$ step)')
    colors_eta = {'10': STRATA['gold'], '100': STRATA['amber'], '1000': STRATA['magma']}
    labels_eta = {
        '10': r'Erebus ($\eta_\mathrm{fluid} = 10\ \mathrm{Pa\cdot s}$)',
        '100': r'Erebus ($\eta_\mathrm{fluid} = 100\ \mathrm{Pa\cdot s}$, baseline)',
        '1000': r'Erebus ($\eta_\mathrm{fluid} = 1000\ \mathrm{Pa\cdot s}$)'
    }
    for key in ['10', '100', '1000']:
        ax1.plot(Fm, data["k_erebus_curves"][key], '-', color=colors_eta[key], linewidth=2.2, label=labels_eta[key])
    ax1.axvspan(0.30, 0.50, color=STRATA['gold'], alpha=0.12, label=r'Transition window $[0.30, 0.50]$')
    ax1.set_yscale('log')
    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(1.0, 5e6)
    ax1.set_xlabel(r'Silicate Melt Fraction $F_m$ [-]')
    ax1.set_ylabel(r'$k_\mathrm{eff}$ [$\mathrm{W/(m\cdot K)}$]')
    ax1.set_title('Thermal Conductivity Regularization', fontsize=12, fontweight='bold')
    ax1.grid(True)
    ax1.legend(loc='lower right', fontsize=8.5)

    # Panel (b): Logarithmic derivative d(log10 k)/dFm
    ax2.text(0.04, 0.93, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    for key in ['10', '100', '1000']:
        ax2.plot(Fm, data["dk_dFm_erebus"][key], '-', color=colors_eta[key], linewidth=2.2, label=labels_eta[key])
    ax2.axvline(0.40, color=STRATA['ink'], linestyle=':', linewidth=2.0, label=r'i2elvis Dirac $\delta$-spike at $0.40$')
    ax2.set_xlim(0.20, 0.60)
    ax2.set_ylim(-2, 45)
    ax2.set_xlabel(r'Silicate Melt Fraction $F_m$ [-]')
    ax2.set_ylabel(r'$\mathrm{d}(\log_{10} k_\mathrm{eff}) / \mathrm{d}F_m$ [-]')
    ax2.set_title('Logarithmic Derivative Smoothness ($C^1$)', fontsize=12, fontweight='bold')
    ax2.grid(True)
    ax2.legend(loc='upper right', fontsize=8.5)

    # Panel (c): Surface temperature difference gating
    ax3.text(0.04, 0.93, '(c)', transform=ax3.transAxes, fontsize=12, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    ax3.plot(dT, w_T, '-', color=STRATA['cobalt'], linewidth=2.4, label=r'Gating factor $w_T = \mathrm{clamp}(\Delta T / \Delta T_\mathrm{min}, 0, 1)$')
    ax3.axvline(10.0, color=STRATA['amber'], linestyle='--', linewidth=1.5, label=r'$\Delta T_\mathrm{min} = 10\ \mathrm{K}$ anchor')
    ax3.axhspan(0, 1, color=STRATA['cobalt'], alpha=0.08)
    ax3.set_xlim(0.0, 30.0)
    ax3.set_ylim(-0.05, 1.05)
    ax3.set_xlabel(r'Temperature Difference $\Delta T = T - T_\mathrm{surface}$ [K]')
    ax3.set_ylabel(r'Thermal Boundary Gate $w_T$ [-]')
    ax3.set_title('Surface Singularity Guard', fontsize=12, fontweight='bold')
    ax3.grid(True)
    ax3.legend(loc='lower right', fontsize=8.5)

    reg_png_assets = os.path.join(ASSETS_DIR, "magma_ocean_regularization.png")
    fig.savefig(reg_png_assets, dpi=200)
    plt.close(fig)
    print(f"Saved {reg_png_assets}")

if __name__ == "__main__":
    generate_video()
    generate_summary_figure()
    generate_convergence_figure()
    generate_regularization_figure()
    print("All benchmark assets generated successfully.")
