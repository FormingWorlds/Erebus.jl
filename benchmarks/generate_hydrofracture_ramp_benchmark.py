#!/usr/bin/env python3
"""
Generate benchmark figure for hydrofracture regularisation ramp and Darcy under-relaxation in Erebus.jl.

Validates:
(a) C^1 continuous overpressure regularisation ramp s(x) vs hard kink law
(b) Regularised derivative ds/dx showing smooth transition over width delta
(c) Geometric convergence of Darcy resistance under-relaxation across plastic iterations
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
    os.path.dirname(__file__), "..", "output_files", "hydrofracture_ramp_benchmark_data.json"
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def generate_benchmark_figure():
    """Generate 3-panel hydrofracture stability and limiters benchmark figure."""
    if not os.path.exists(DATA_PATH):
        print(f"Data file not found at {DATA_PATH}. Running exporter...")
        import subprocess
        subprocess.run(
            ["julia", "--project=.", "benchmarks/export_hydrofracture_ramp_benchmark.jl"],
            check=True,
            cwd=os.path.join(os.path.dirname(__file__), "..")
        )

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    x_sweep = np.array(data["x_sweep"])
    ramp_curves = data["ramp_curves"]
    k_iters = np.array(data["k_iters"])
    relaxation_curves = data["relaxation_curves"]

    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(14.0, 4.2), dpi=200)

    delta_colors = {
        "0.0": STRATA['ink'],
        "0.02": STRATA['cobalt'],
        "0.05": STRATA['amber'],
        "0.1": STRATA['magma'],
    }
    delta_styles = {
        "0.0": (':', 2.0, "Kink law ($\\delta = 0$)"),
        "0.02": ('-', 1.8, "$C^1$ ramp ($\\delta = 0.02$)"),
        "0.05": ('-', 2.0, "$C^1$ ramp ($\\delta = 0.05$)"),
        "0.1": ('--', 1.8, "$C^1$ ramp ($\\delta = 0.10$)"),
    }

    # Panel (a): Regularised overpressure s(x)
    for key, (ls, lw, label) in delta_styles.items():
        if key in ramp_curves:
            s_vals = np.array(ramp_curves[key]["s"])
            ax_a.plot(x_sweep, s_vals, ls, color=delta_colors[key], linewidth=lw, label=label)

    ax_a.axvline(0.0, color=NEUTRALS['graphite'], linestyle='--', linewidth=0.8, alpha=0.6)
    ax_a.set_xlabel("Normalized Overpressure $x = (-P_\\mathrm{eff} - \\sigma_t) / \\sigma_t$", fontsize=9.5)
    ax_a.set_ylabel("Regularised Overpressure $s(x)$", fontsize=9.5)
    ax_a.set_title("(a) $C^1$ Overpressure Function", fontsize=10.5, fontweight='bold', loc='left')
    ax_a.set_xlim(-0.08, 0.25)
    ax_a.set_ylim(-0.01, 0.25)
    ax_a.grid(True)
    ax_a.legend(loc="upper left", fontsize=8.0)

    # Panel (b): Derivative ds/dx
    for key, (ls, lw, label) in delta_styles.items():
        if key in ramp_curves:
            ds_vals = np.array(ramp_curves[key]["ds"])
            ax_b.plot(x_sweep, ds_vals, ls, color=delta_colors[key], linewidth=lw, label=label)

    ax_b.axvline(0.0, color=NEUTRALS['graphite'], linestyle='--', linewidth=0.8, alpha=0.6)
    ax_b.set_xlabel("Normalized Overpressure $x = (-P_\\mathrm{eff} - \\sigma_t) / \\sigma_t$", fontsize=9.5)
    ax_b.set_ylabel("Derivative $s'(x) = \\mathrm{d}s / \\mathrm{d}x$", fontsize=9.5)
    ax_b.set_title("(b) Regularised Derivative Continuity", fontsize=10.5, fontweight='bold', loc='left')
    ax_b.set_xlim(-0.08, 0.25)
    ax_b.set_ylim(-0.05, 1.15)
    ax_b.grid(True)
    ax_b.legend(loc="center right", fontsize=8.0)

    # Panel (c): Darcy resistance under-relaxation convergence
    theta_colors = {
        "0.1": STRATA['plum'],
        "0.3": STRATA['amber'],
        "0.5": STRATA['cobalt'],
        "1.0": STRATA['ink'],
    }
    theta_styles = {
        "0.1": ('--', 1.8, "$\\theta = 0.1$"),
        "0.3": ('-', 2.2, "$\\theta = 0.3$ (recommended)"),
        "0.5": ('-.', 1.8, "$\\theta = 0.5$"),
        "1.0": (':', 2.0, "$\\theta = 1.0$ (unrelaxed)"),
    }

    for key, (ls, lw, label) in theta_styles.items():
        if key in relaxation_curves:
            res_vals = np.array(relaxation_curves[key])
            ax_c.semilogy(k_iters, res_vals, ls, color=theta_colors[key], linewidth=lw, label=label)

    ax_c.set_xlabel("Plastic Iteration $k$", fontsize=9.5)
    ax_c.set_ylabel("Normalized Resistance Error $|r^{(k)} - r_\\mathrm{new}| / |r^{(0)} - r_\\mathrm{new}|$", fontsize=9.5)
    ax_c.set_title("(c) Darcy Under-Relaxation Convergence", fontsize=10.5, fontweight='bold', loc='left')
    ax_c.set_xlim(0, 20)
    ax_c.set_ylim(1.0e-16, 2.0)
    ax_c.grid(True, which="both")
    ax_c.legend(loc="lower left", fontsize=8.0)

    fig.tight_layout()

    fig.savefig(os.path.join(ASSETS_DIR, "hydrofracture_ramp_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(ASSETS_DIR, "hydrofracture_ramp_benchmark.pdf"))
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "hydrofracture_ramp_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "hydrofracture_ramp_benchmark.pdf"))
    plt.close(fig)
    print(f"Saved benchmark figure to {ASSETS_DIR}/hydrofracture_ramp_benchmark.png")


def main():
    """Main entry point for generating hydrofracture ramp benchmark figures."""
    generate_benchmark_figure()


if __name__ == "__main__":
    main()
