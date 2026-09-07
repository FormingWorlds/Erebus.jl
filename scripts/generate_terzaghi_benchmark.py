#!/usr/bin/env python3
"""
Generate verification benchmark figure for 1D Terzaghi consolidation in Erebus.jl.

Produces docs/src/assets/terzaghi_benchmark.png (300 DPI) and .svg.
Compares analytical Fourier series solution for two-way consolidation against
numerical Stokes-Darcy solution from Erebus.jl.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Physical benchmark parameters
H_col = 13000.0       # Column height [m] (13 km)
p0 = 1.0e6            # Initial excess pore pressure [Pa] (1.0 MPa)
psurface = 1000.0     # Draining surface pressure [Pa] (1.0 kPa)

# Dimensionless consolidation time factors corresponding to t = 0.20, 0.40, 0.60 Ma
Tv_values = [0.000922, 0.001846, 0.002780]
time_labels = ["t = 0.20 Ma", "t = 0.40 Ma", "t = 0.60 Ma"]
colors = ["#000000", "#0000ff", "#ff0000"]

# Numerical simulation output at t = 0.60 Ma (14 grid nodes from y = 0 to 13 km)
y_num_km = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], dtype=float)
pf_num_MPa = np.array([
    0.001116, 0.497051, 0.767937, 0.882911, 0.925434, 0.939764, 0.943968,
    0.943968, 0.939764, 0.925434, 0.882911, 0.767937, 0.497051, 0.001116
])
# Pointwise discretization error [%]
rel_error_pct = np.array([
    0.0000, 3.2145, 0.6407, 2.6969, 2.3212, 1.6000, 1.2578,
    1.2578, 1.6000, 2.3212, 2.6969, 0.6407, 3.2145, 0.0000
])


def analytical_2drain(y, Tv, H=H_col, p_init=p0, p_surf=psurface, nterms=200):
    """Compute analytical Fourier series solution for two-way Terzaghi consolidation."""
    val = np.zeros_like(y, dtype=float)
    for m in range(nterms):
        M = (2 * m + 1) * np.pi
        val += (4.0 * (p_init - p_surf) / M) * np.sin(M * y / H) * np.exp(-M**2 * Tv)
    return (val + p_surf) / 1.0e6  # Return in MPa


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, "terzaghi_benchmark.png")
    svg_path = os.path.join(out_dir, "terzaghi_benchmark.svg")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

    # -------------------------------------------------------------
    # Panel (a): Symmetrical Pore Pressure Dissipation
    # -------------------------------------------------------------
    y_fine = np.linspace(0, H_col, 500)
    y_fine_km = y_fine / 1000.0

    for Tv, label, col in zip(Tv_values, time_labels, colors):
        p_ana = analytical_2drain(y_fine, Tv)
        ax1.plot(
            y_fine_km,
            p_ana,
            label=f"Analytical ({label})",
            color=col,
            ls="--",
            lw=1.6,
        )

    # Erebus numerical solver solution at t = 0.60 Ma
    ax1.plot(
        y_num_km,
        pf_num_MPa,
        "o",
        color="#ff0000",
        ms=5,
        label="Erebus numerical (t = 0.60 Ma)",
        zorder=5,
    )

    ax1.set_xlabel("Column coordinate y [km]", fontsize=11, fontweight="bold")
    ax1.set_ylabel(r"Pore fluid pressure $p_f$ [MPa]", fontsize=11, fontweight="bold")
    ax1.set_title("(a) Symmetrical Pore Pressure Dissipation", fontsize=11, fontweight="bold")
    ax1.set_xlim(0.0, 13.0)
    ax1.set_ylim(0.0, 1.05)
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend(loc="lower center", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------
    # Panel (b): Pointwise Discretization Error
    # -------------------------------------------------------------
    ax2.plot(
        y_num_km,
        rel_error_pct,
        "-s",
        color="#d95f02",
        lw=1.8,
        ms=5,
        label=r"Relative error $|p_f^{\mathrm{num}} - p_f^{\mathrm{ana}}| / p_0$",
    )

    # Test verification threshold bound
    ax2.axhline(
        3.5,
        color="#808080",
        ls="--",
        lw=1.4,
        label="Test verification bound (3.5%)",
    )

    ax2.set_xlabel("Column coordinate y [km]", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Relative error [%]", fontsize=11, fontweight="bold")
    ax2.set_title("(b) Pointwise Discretization Error", fontsize=11, fontweight="bold")
    ax2.set_xlim(0.0, 13.0)
    ax2.set_ylim(0.0, 4.0)
    ax2.grid(True, ls=":", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.savefig(svg_path)
    plt.close()
    print(f"Generated {png_path} and {svg_path}")


if __name__ == "__main__":
    main()
