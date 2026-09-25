#!/usr/bin/env python3
"""
Generate benchmark figure for Jeans kinetic atmospheric effusion in Erebus.jl.

Validates atomic hydrogen kinetic escape flux against the analytical
Maxwell-Boltzmann exobase effusion integral across exobase temperatures.
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
    os.path.dirname(__file__), "..", "output_files", "jeans_effusion_benchmark_data.json"
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def generate_benchmark_figure():
    """Generate 2-panel Jeans effusion benchmark figure."""
    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(
            f"Benchmark data file not found at {DATA_PATH}. "
            "Run benchmarks/export_jeans_effusion_benchmark.jl first."
        )

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    T_sweep = np.array(data["T_sweep"])
    lambda_sweep = np.array(data["lambda_sweep"])
    phi_ref = np.array(data["phi_ref"])
    phi_code = np.array(data["phi_code"])
    rel_errors = np.array(data["rel_errors"])
    T_hot = data["T_hot"]
    flux_eff_hot = data["flux_effusion_only"]
    flux_hydro_hot = data["flux_hydro_active"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=200)

    # Panel (a): Effusion particle flux vs Exobase Temperature
    ax_a = axes[0]
    ax_a.text(
        0.04, 0.92, '(a)', transform=ax_a.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )
    ax_a.plot(
        T_sweep, phi_ref, '-', color=STRATA['cobalt'], linewidth=2.2,
        label=r'Maxwell-Boltzmann: $\frac{n v_{\mathrm{th}}}{2\sqrt{\pi}}(1+\lambda)e^{-\lambda}$'
    )
    # Plot subsampled points for clarity
    step = max(1, len(T_sweep) // 15)
    ax_a.plot(
        T_sweep[::step], phi_code[::step], 'o', color=STRATA['magma'], markersize=6.0,
        markeredgewidth=1.2, markeredgecolor='white', label='Erebus.jl (kinetic mode)'
    )
    # Highlight branch switch point
    ax_a.plot(
        [T_hot], [flux_hydro_hot], 's', color=STRATA['amber'], markersize=7.0,
        markeredgewidth=1.2, markeredgecolor='white', label=r'Hydrodynamic branch ($\lambda < 2$)'
    )

    ax_a.set_xlabel(r'Exobase Temperature $T_{\mathrm{exo}}$ [K]', fontsize=10.0)
    ax_a.set_ylabel(r'Escape Number Flux $\Phi$ [$\mathrm{m}^{-2}\,\mathrm{s}^{-1}$]', fontsize=10.0)
    ax_a.set_yscale('log')
    ax_a.set_xlim(150.0, 3100.0)
    ax_a.grid(True)
    ax_a.legend(loc='lower right', fontsize=8.5)

    # Panel (b): Relative Error vs Exobase Temperature
    ax_b = axes[1]
    ax_b.text(
        0.04, 0.92, '(b)', transform=ax_b.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )
    # Guard against exact zeros in log scale
    err_plot = np.maximum(rel_errors, 1.0e-17)
    ax_b.semilogy(
        T_sweep, err_plot, '.-', color=STRATA['plum'], linewidth=1.5,
        markersize=4.0, label=r'Relative Error $|\Phi_{\mathrm{code}} - \Phi_{\mathrm{ref}}| / \Phi_{\mathrm{ref}}$'
    )
    ax_b.axhline(
        1.0e-6, color=STRATA['magma'], linestyle='--', linewidth=1.2,
        label=r'Quality Threshold ($10^{-6}$)'
    )

    ax_b.set_xlabel(r'Exobase Temperature $T_{\mathrm{exo}}$ [K]', fontsize=10.0)
    ax_b.set_ylabel(r'Relative Difference', fontsize=10.0)
    ax_b.set_xlim(150.0, 3100.0)
    ax_b.set_ylim(1.0e-18, 1.0e-4)
    ax_b.grid(True)
    ax_b.legend(loc='upper right', fontsize=8.5)

    fig.tight_layout()

    out_asset = os.path.join(ASSETS_DIR, "jeans_effusion_benchmark.png")
    out_files = os.path.join(OUTPUT_FILES_DIR, "jeans_effusion_benchmark.png")
    fig.savefig(out_asset, dpi=200, bbox_inches='tight')
    fig.savefig(out_files, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Figure saved to: {out_asset}")
    print(f"Figure saved to: {out_files}")


if __name__ == "__main__":
    generate_benchmark_figure()
