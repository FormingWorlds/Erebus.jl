#!/usr/bin/env python3
"""
Generate benchmark figure for magma ocean volatile degassing in Erebus.jl.

Validates retained volatile inventory in Lagrangian melt markers against the
analytical Burnham (1979) / Dixon et al. (1995) solubility law across silicate melt fractions.
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
    os.path.dirname(__file__), "..", "output_files", "degassing_benchmark_data.json"
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def generate_benchmark_figure():
    """Generate 2-panel degassing benchmark figure."""
    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(
            f"Benchmark data file not found at {DATA_PATH}. "
            "Run benchmarks/export_degassing_benchmark.jl first."
        )

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    F_ref = np.array(data["F_ref"])
    ret_ref = np.array(data["retained_ref_wtpct"])
    F_sim = np.array(data["F_sim"])
    ret_sim = np.array(data["retained_sim_wtpct"])
    rel_errors = np.array(data["rel_errors"])
    p_H2O_MPa = data["p_H2O_MPa"]
    water_As = data["water_As"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=200)

    # Panel (a): Retained water vs melt fraction
    ax_a = axes[0]
    ax_a.text(
        0.04, 0.92, '(a)', transform=ax_a.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )
    ax_a.plot(
        F_ref, ret_ref, '-', color=STRATA['cobalt'], linewidth=2.2,
        label=r'Burnham (1979): $w = F \cdot A_s \sqrt{p_{\mathrm{H}_2\mathrm{O}}}$'
    )
    ax_a.plot(
        F_sim, ret_sim, 'o', color=STRATA['magma'], markersize=6.5, markeredgewidth=1.2,
        markeredgecolor='white', label='Lagrangian markers'
    )
    ax_a.set_xlabel('Silicate Melt Fraction $F$', fontsize=10.0)
    ax_a.set_ylabel(r'Retained Water $w_{\mathrm{H}_2\mathrm{O}}$ [wt%]', fontsize=10.0)
    ax_a.set_xlim(0.05, 1.05)
    ax_a.set_ylim(0.0, max(ret_ref) * 1.15)
    ax_a.set_title('Retained Water vs Melt Fraction', fontsize=11.0, fontweight='bold')
    ax_a.grid(True)
    ax_a.legend(loc='lower right', fontsize=8.5)

    # Panel (b): Relative error vs tolerance
    ax_b = axes[1]
    ax_b.text(
        0.04, 0.92, '(b)', transform=ax_b.transAxes, fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9)
    )
    # Plot machine precision baseline or measured error
    clipped_err = np.maximum(rel_errors, 1.0e-16)
    ax_b.semilogy(
        F_sim, clipped_err, 's-', color=STRATA['plum'], linewidth=1.8, markersize=5.5,
        label='Simulation relative residual'
    )
    ax_b.axhline(
        1.0e-6, color=STRATA['amber'], linestyle='--', linewidth=1.5,
        label=r'Verification threshold ($10^{-6}$)'
    )
    ax_b.set_xlabel('Silicate Melt Fraction $F$', fontsize=10.0)
    ax_b.set_ylabel('Relative Error', fontsize=10.0)
    ax_b.set_xlim(0.05, 1.05)
    ax_b.set_ylim(1.0e-17, 1.0e-4)
    ax_b.set_title('Analytical Solution Agreement', fontsize=11.0, fontweight='bold')
    ax_b.grid(True)
    ax_b.legend(loc='upper right', fontsize=8.5)

    fig.tight_layout()

    fig.savefig(os.path.join(ASSETS_DIR, "degassing_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(ASSETS_DIR, "degassing_benchmark.pdf"))
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "degassing_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "degassing_benchmark.pdf"))
    plt.close(fig)
    print(f"Saved benchmark figure to {ASSETS_DIR}/degassing_benchmark.png")


def main():
    """Main entry point for generating degassing benchmark figures."""
    generate_benchmark_figure()


if __name__ == "__main__":
    main()
