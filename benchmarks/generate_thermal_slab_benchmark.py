#!/usr/bin/env python3
"""
Generate benchmark figure for 2D analytical thermal slab diffusion in Erebus.jl.

Validates the transient thermal conduction solver against the exact 2D analytical
cosine decay solution in an insulated box with homogeneous Neumann boundaries.
"""

import os
import sys
import json
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

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

DATA_PATH = os.path.join(
    os.path.dirname(__file__), "..", "output_files", "thermal_slab_benchmark_data.json"
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def generate_benchmark_figure():
    """Generate 2-panel thermal slab benchmark figure."""
    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(
            f"Benchmark data file not found at {DATA_PATH}. "
            "Run benchmarks/export_thermal_slab_benchmark.jl first."
        )

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    times_Myr = data["times_Myr"]
    x_km = np.array(data["x_centerline_km"])
    profiles_ana = [np.array(p) for p in data["profiles_analytical"]]
    profiles_num = [np.array(p) for p in data["profiles_numerical"]]
    dx_km = np.array(data["dx_km"])
    l2_errors = np.array(data["l2_errors"])
    linf_errors = np.array(data["linf_errors"])

    colors = [STRATA['cobalt'], STRATA['plum'], STRATA['magma'], STRATA['amber']]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=200)

    # Panel (a): Centerline temperature profiles
    ax_a = axes[0]
    ax_a.text(
        0.04, 0.92, '(a)', transform=ax_a.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )

    for idx, (t_val, col) in enumerate(zip(times_Myr, colors)):
        t_label = f"$t = {t_val:.1f}$ Ma" if t_val > 0.0 else "$t = 0$"
        ax_a.plot(
            x_km, profiles_ana[idx], '-', color=col, linewidth=1.8,
            label=f"Analytical ({t_label})"
        )
        # Downsample marker points to avoid clutter
        step_pts = max(1, len(x_km) // 10)
        ax_a.plot(
            x_km[::step_pts], profiles_num[idx][::step_pts], 'o', color=col,
            markersize=5.0, markeredgewidth=1.0, markeredgecolor='white',
            label=f"Numerical ({t_label})"
        )

    ax_a.set_xlabel('Horizontal Coordinate $x$ [km]', fontsize=10.0)
    ax_a.set_ylabel('Temperature $T$ [K]', fontsize=10.0)
    ax_a.set_xlim(0.0, 100.0)
    ax_a.set_ylim(240.0, 360.0)
    ax_a.set_title('Centerline Temperature Profiles ($y = L_y/2$)', fontsize=11.0, fontweight='bold')
    ax_a.grid(True)
    ax_a.legend(loc='lower right', fontsize=7.5, ncol=2)

    # Panel (b): Grid convergence and error vs dx
    ax_b = axes[1]
    ax_b.text(
        0.04, 0.92, '(b)', transform=ax_b.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )

    ax_b.loglog(
        dx_km, linf_errors, 's-', color=STRATA['plum'], linewidth=1.8, markersize=6.0,
        label=r'$L_\infty$ relative error'
    )
    ax_b.loglog(
        dx_km, l2_errors, 'o-', color=STRATA['cobalt'], linewidth=1.8, markersize=6.0,
        label=r'$L_2$ relative error'
    )
    ax_b.axhline(
        1.0e-3, color=STRATA['amber'], linestyle='--', linewidth=1.5,
        label=r'Target tolerance ($10^{-3}$)'
    )

    ax_b.set_xlabel(r'Grid Spacing $\Delta x$ [km]', fontsize=10.0)
    ax_b.set_ylabel('Relative Error', fontsize=10.0)
    ax_b.set_xlim(1.0, 10.0)
    ax_b.set_ylim(1.0e-4, 5.0e-3)
    ax_b.set_title('Spatial Grid Convergence ($t = 0.1\\,\\tau$)', fontsize=11.0, fontweight='bold')
    ax_b.grid(True, which="both")
    ax_b.legend(loc='lower right', fontsize=8.5)

    fig.tight_layout()

    fig.savefig(os.path.join(ASSETS_DIR, "thermal_slab_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(ASSETS_DIR, "thermal_slab_benchmark.pdf"))
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "thermal_slab_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "thermal_slab_benchmark.pdf"))
    plt.close(fig)
    print(f"Saved benchmark figure to {ASSETS_DIR}/thermal_slab_benchmark.png")


def main():
    """Main entry point for generating thermal slab benchmark figures."""
    generate_benchmark_figure()


if __name__ == "__main__":
    main()
