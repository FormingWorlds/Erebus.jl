#!/usr/bin/env python3
"""
Generate benchmark figure for short-lived radiogenic isotope decay in Erebus.jl.

Validates 26Al and 60Fe specific radiogenic power decay curves and cumulative energy
release against analytical solutions and literature half-lives (Tang & Dauphas 2012).
"""

import os
import sys
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

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'mathtext.fontset': 'stixsans',
    'axes.edgecolor': '#333333',
    'axes.linewidth': 0.8,
    'axes.labelcolor': '#222222',
    'axes.facecolor': '#FFFFFF',
    'xtick.color': '#222222',
    'ytick.color': '#222222',
    'grid.color': '#E0E0E0',
    'grid.linestyle': ':',
    'grid.linewidth': 0.6,
    'legend.frameon': True,
    'legend.facecolor': '#FFFFFF',
    'legend.edgecolor': '#CCCCCC',
    'figure.facecolor': '#FFFFFF',
})

DATA_PATH = os.path.join(
    os.path.dirname(__file__), "..", "output_files", "radiogenic_decay_benchmark_data.json"
)
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")


def main():
    """Generate 2-panel radiogenic decay benchmark figure."""
    if not os.path.exists(DATA_PATH):
        print(f"Benchmark data file not found at: {DATA_PATH}")
        print("Run benchmarks/export_radiogenic_decay_benchmark.jl first.")
        return

    with open(DATA_PATH, "r") as f:
        data = json.load(f)

    t_Myr = np.array(data["t_Myr"])
    Q_al = np.array(data["Q_al_W_kg"])
    Q_fe = np.array(data["Q_fe_W_kg"])
    t_half_al = data["t_half_al_Myr"]
    t_half_fe = data["t_half_fe_Myr"]
    tau_al = t_half_al / np.log(2.0)
    tau_fe = t_half_fe / np.log(2.0)

    # Analytical curves
    Q_al_ana = data["Q0_al"] * np.exp(-t_Myr / tau_al)
    Q_fe_ana = data["Q0_fe"] * np.exp(-t_Myr / tau_fe)

    # Cumulative specific energy release [MJ / kg]
    # E(t) = integral_0^t Q(t') dt'
    sec_per_Myr = 1.0e6 * 31_540_000.0
    E_al_MJ = (data["Q0_al"] * (tau_al * sec_per_Myr) * (1.0 - np.exp(-t_Myr / tau_al))) * 1.0e-6
    E_fe_MJ = (data["Q0_fe"] * (tau_fe * sec_per_Myr) * (1.0 - np.exp(-t_Myr / tau_fe))) * 1.0e-6
    E_tot_MJ = E_al_MJ + E_fe_MJ

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # --- Panel (a): Specific Power Decay ---
    ax1.plot(t_Myr, Q_al, color=PALETTE["slate_blue"], lw=2.2, label=r'$^{26}\mathrm{Al}$ (Erebus)')
    ax1.plot(t_Myr, Q_al_ana, color=PALETTE["charcoal"], lw=1.2, ls='--', label=r'$^{26}\mathrm{Al}$ analytical')
    ax1.plot(t_Myr, Q_fe, color=PALETTE["amber"], lw=2.2, label=r'$^{60}\mathrm{Fe}$ (Erebus)')
    ax1.plot(t_Myr, Q_fe_ana, color=PALETTE["crimson"], lw=1.2, ls='--', label=r'$^{60}\mathrm{Fe}$ analytical')

    ax1.axvline(t_half_al, color=PALETTE["slate_blue"], ls=':', lw=0.9, alpha=0.7)
    ax1.axvline(t_half_fe, color=PALETTE["amber"], ls=':', lw=0.9, alpha=0.7)

    ax1.set_yscale('log')
    ax1.set_xlim(0.0, 10.0)
    ax1.set_ylim(1.0e-15, 1.0e-6)
    ax1.set_xlabel('Elapsed Time since CAIs [Myr]', fontsize=11)
    ax1.set_ylabel(r'Specific Radiogenic Power $Q(t)$ [$\mathrm{W\,kg^{-1}}$]', fontsize=11)
    ax1.set_title('(a) Radiogenic Isotope Power Decay', fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, which='both', alpha=0.4)
    ax1.legend(loc='upper right', fontsize=9.5)

    # --- Panel (b): Cumulative Energy Release ---
    ax2.plot(t_Myr, E_tot_MJ, color=PALETTE["purple"], lw=2.2, label=r'Total ($^{26}\mathrm{Al} + {}^{60}\mathrm{Fe}$)')
    ax2.plot(t_Myr, E_al_MJ, color=PALETTE["slate_blue"], lw=1.8, ls='-', label=r'$^{26}\mathrm{Al}$ cumulative')
    ax2.plot(t_Myr, E_fe_MJ, color=PALETTE["amber"], lw=1.8, ls='-', label=r'$^{60}\mathrm{Fe}$ cumulative')

    ax2.set_xlim(0.0, 10.0)
    ax2.set_ylim(0.0, max(E_tot_MJ) * 1.1)
    ax2.set_xlabel('Elapsed Time since CAIs [Myr]', fontsize=11)
    ax2.set_ylabel(r'Cumulative Specific Energy [$\mathrm{MJ\,kg^{-1}}$]', fontsize=11)
    ax2.set_title('(b) Cumulative Radiogenic Energy Release', fontsize=12, fontweight='bold', pad=10)
    ax2.grid(True, alpha=0.4)
    ax2.legend(loc='lower right', fontsize=9.5)

    fig.tight_layout()

    os.makedirs(ASSETS_DIR, exist_ok=True)
    os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

    fig.savefig(os.path.join(ASSETS_DIR, "radiogenic_decay_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(ASSETS_DIR, "radiogenic_decay_benchmark.pdf"))
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "radiogenic_decay_benchmark.png"), dpi=200)
    fig.savefig(os.path.join(OUTPUT_FILES_DIR, "radiogenic_decay_benchmark.pdf"))
    plt.close(fig)

    print("Saved radiogenic_decay_benchmark figures to assets and output_files directories.")


if __name__ == "__main__":
    main()
