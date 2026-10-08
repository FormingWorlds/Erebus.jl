#!/usr/bin/env python3
"""
Generate benchmark figure for differentiated planetesimal self-gravity in Erebus.jl.

Validates the 3D enclosed-mass gravity formulation against exact analytical solutions
for a two-layer differentiated body and demonstrates the limitation of the 2D Poisson solver.
"""

import os
import json
import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

# Neutral scientific palette
PALETTE = {
    "slate_blue": "#1F4E79",
    "crimson": "#C0392B",
    "green": "#27AE60",
    "amber": "#D97706",
    "gold": "#D4AC0D",
    "purple": "#6C3483",
    "teal": "#16A085",
    "charcoal": "#1A1A1A",
    "dark_gray": "#222222",
    "grid_gray": "#E0E0E0",
    "off_white": "#F5F5F5",
    "white": "#FFFFFF",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "stixsans",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "axes.labelcolor": "#222222",
        "axes.facecolor": "#FFFFFF",
        "xtick.color": "#222222",
        "ytick.color": "#222222",
        "grid.color": "#E0E0E0",
        "grid.linestyle": ":",
        "grid.linewidth": 0.6,
        "legend.frameon": True,
        "legend.facecolor": "#FFFFFF",
        "legend.edgecolor": "#CCCCCC",
        "figure.facecolor": "#FFFFFF",
    }
)

DATA_PATH = os.path.join(
    os.path.dirname(__file__),
    "..",
    "output_files",
    "gravity_two_layer_benchmark_data.json",
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def generate_benchmark_figure():
    """Generate 2-panel gravity benchmark figure."""
    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(
            f"Benchmark data file not found at {DATA_PATH}. "
            "Run benchmarks/export_gravity_two_layer_benchmark.jl first."
        )

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    r_eval_km = np.array(data["r_eval_km"])
    g_ana = np.array(data["g_analytical"])
    g_p2d = np.array(data["g_poisson2d"])
    r_bins_km = np.array(data["r_bins_km"])
    g_bins = np.array(data["g_bins"])
    R_p = data["R_planet_km"]
    R_c = data["R_core_km"]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=200)

    # Panel (a): Radial gravity profiles
    ax_a = axes[0]
    ax_a.plot(
        r_eval_km, g_ana, "-", color=PALETTE["charcoal"], linewidth=2.0, label="3D Analytical"
    )
    ax_a.plot(
        r_bins_km[::4],
        g_bins[::4],
        "o",
        color=PALETTE["slate_blue"],
        markersize=4.0,
        markeredgewidth=0.8,
        markeredgecolor="white",
        label="Enclosed Mass (:enclosed_mass)",
    )
    ax_a.plot(
        r_eval_km,
        g_p2d,
        "--",
        color=PALETTE["amber"],
        linewidth=1.8,
        label="2D Poisson (:poisson2d)",
    )

    # Structural boundaries
    ax_a.axvline(R_c, color=PALETTE["purple"], linestyle=":", linewidth=1.2)
    ax_a.axvline(R_p, color=PALETTE["crimson"], linestyle=":", linewidth=1.2)
    ax_a.text(
        R_c - 1.5,
        0.058,
        r"$r_c$",
        color=PALETTE["purple"],
        fontsize=10,
        ha="right",
        va="top",
    )
    ax_a.text(
        R_p - 1.5,
        0.058,
        r"$R$",
        color=PALETTE["crimson"],
        fontsize=10,
        ha="right",
        va="top",
    )

    ax_a.set_xlabel("Radial Distance $r$ [km]", fontsize=10.0)
    ax_a.set_ylabel("Gravitational Acceleration $g$ [m/s$^2$]", fontsize=10.0)
    ax_a.set_xlim(0.0, 70.0)
    ax_a.set_ylim(0.0, 0.065)
    ax_a.set_title(
        "(a) Radial Gravity in Differentiated Planetesimal",
        fontsize=11.0,
        fontweight="bold",
        loc="left",
    )
    ax_a.grid(True)
    ax_a.legend(loc="lower right", fontsize=8.5)

    # Panel (b): Core-excess gravity outside core (rc <= r <= R)
    ax_b = axes[1]
    const_G = 6.6743e-11
    rho_m = data["rho_mantle"]

    mask = (r_eval_km >= R_c) & (r_eval_km <= R_p)
    r_sub = r_eval_km[mask]
    g_excess_ana = g_ana[mask] - (4.0 / 3.0) * np.pi * const_G * rho_m * (
        r_sub * 1000.0
    )
    g_excess_p2d = g_p2d[mask] - (4.0 / 3.0) * np.pi * const_G * rho_m * (
        r_sub * 1000.0
    )

    mask_bins = (r_bins_km >= R_c) & (r_bins_km <= R_p)
    r_bins_sub = r_bins_km[mask_bins]
    g_excess_enc = g_bins[mask_bins] - (4.0 / 3.0) * np.pi * const_G * rho_m * (
        r_bins_sub * 1000.0
    )

    ax_b.plot(
        r_sub,
        g_excess_ana * 1000.0,
        "-",
        color=PALETTE["charcoal"],
        linewidth=2.0,
        label="3D Analytical",
    )
    ax_b.plot(
        r_bins_sub[::2],
        g_excess_enc[::2] * 1000.0,
        "o",
        color=PALETTE["slate_blue"],
        markersize=4.0,
        markeredgewidth=0.8,
        markeredgecolor="white",
        label="Enclosed Mass (:enclosed_mass)",
    )
    ax_b.plot(
        r_sub,
        g_excess_p2d * 1000.0,
        "--",
        color=PALETTE["amber"],
        linewidth=1.8,
        label=r"2D Poisson ($R/r_c \approx 2\times$)",
    )

    ax_b.axvline(R_c, color=PALETTE["purple"], linestyle=":", linewidth=1.2)
    ax_b.axvline(R_p, color=PALETTE["crimson"], linestyle=":", linewidth=1.2)

    ax_b.set_xlabel("Radial Distance $r$ [km]", fontsize=10.0)
    ax_b.set_ylabel(
        r"Core-Excess Gravity $g_\mathrm{excess}$ [mm/s$^2$]", fontsize=10.0
    )
    ax_b.set_xlim(24.0, 52.0)
    ax_b.set_ylim(0.0, 32.0)
    ax_b.set_xticks([25.0, 30.0, 35.0, 40.0, 45.0, 50.0])
    ax_b.set_xticklabels([r"$r_c = 25$", "30", "35", "40", "45", r"$R = 50$"])
    ax_b.set_title(
        r"(b) Core-Excess Gravity in Mantle ($r_c \leq r \leq R$)",
        fontsize=11.0,
        fontweight="bold",
        loc="left",
    )
    ax_b.grid(True)
    ax_b.legend(loc="upper right", fontsize=8.5)

    fig.tight_layout()

    fig.savefig(os.path.join(ASSETS_DIR, "gravity_two_layer_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(ASSETS_DIR, "gravity_two_layer_benchmark.pdf"))
    fig.savefig(
        os.path.join(OUTPUT_FILES_DIR, "gravity_two_layer_benchmark.png"), dpi=200
    )
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "gravity_two_layer_benchmark.pdf"))
    plt.close(fig)
    print(f"Saved benchmark figure to {ASSETS_DIR}/gravity_two_layer_benchmark.png")


def main():
    """Main entry point for generating gravity benchmark figures."""
    generate_benchmark_figure()


if __name__ == "__main__":
    main()
