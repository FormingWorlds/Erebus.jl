#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for multi-phase HCNSPO
volatile mixtures, ammonia-water freezing point depression, multi-snowline
disk condensation, and refractory organic pyrolysis in Erebus.jl.
"""

import json
import os
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
    'cream': '#F2EAD3',
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

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")
os.makedirs(ASSETS_DIR, exist_ok=True)
os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)


def compute_freezing_curve(x_nh3_arr, lambda_nh3=None, T_pure=273.15, T_eutectic=176.0, x_eutectic=0.33):
    """Compute equilibrium freezing point temperature array."""
    if lambda_nh3 is None:
        lambda_nh3 = (T_pure - T_eutectic) / x_eutectic
    t_freeze = np.zeros_like(x_nh3_arr)
    for i, x in enumerate(x_nh3_arr):
        t_dep = T_pure - lambda_nh3 * x
        t_freeze[i] = max(t_dep, T_eutectic)
    return t_freeze


def evaluate_condensation_profile(T_arr):
    """Evaluate multi-snowline ice condensation fractions across temperatures."""
    T_conds = {
        'H2O': 160.0,
        'NH3': 135.0,
        'CO2': 75.0,
        'H2S': 75.0,
        'CH4': 45.0,
        'CO': 25.0,
        'N2': 18.0
    }
    nominal_mass_fracs = {
        'H2O': 0.82,
        'CO2': 0.08,
        'NH3': 0.04,
        'CH4': 0.02,
        'CO': 0.02,
        'N2': 0.01,
        'H2S': 0.01
    }

    condensed = {sp: np.zeros_like(T_arr) for sp in nominal_mass_fracs}
    for sp, Tc in T_conds.items():
        if sp in condensed:
            condensed[sp] = np.where(T_arr <= Tc, nominal_mass_fracs[sp], 0.0)
    return condensed, T_conds


def evaluate_pyrolysis_profile(T_arr, T_pyro=600.0, f_graphite=0.60, DeltaT=250.0):
    """Evaluate refractory organic pyrolysis products across temperatures."""
    C_refr_rem = np.zeros_like(T_arr)
    C_graphite = np.zeros_like(T_arr)
    C_gas = np.zeros_like(T_arr)
    N_refr_rem = np.zeros_like(T_arr)
    N_gas = np.zeros_like(T_arr)

    for i, T in enumerate(T_arr):
        if T <= T_pyro:
            C_refr_rem[i] = 1.0
            C_graphite[i] = 0.0
            C_gas[i] = 0.0
            N_refr_rem[i] = 1.0
            N_gas[i] = 0.0
        else:
            xi = min(0.80, (T - T_pyro) / DeltaT)
            C_refr_rem[i] = 1.0 - xi
            C_graphite[i] = f_graphite * xi
            C_gas[i] = (1.0 - f_graphite) * xi
            N_refr_rem[i] = 1.0 - xi
            N_gas[i] = xi

    return C_refr_rem, C_graphite, C_gas, N_refr_rem, N_gas


def main():
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 10.5))
    plt.subplots_adjust(hspace=0.32, wspace=0.28, left=0.08, right=0.96, top=0.94, bottom=0.07)

    # -------------------------------------------------------------------------
    # Panel (a): Freezing Point Depression & Eutectic Floor
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    x_nh3 = np.linspace(0.0, 0.40, 300)
    t_freeze = compute_freezing_curve(x_nh3)

    ax_a.plot(x_nh3 * 100.0, t_freeze, color=STRATA['cobalt'], lw=2.5,
              label="Equilibrium liquidus $T_{\\mathrm{freeze}}$")
    ax_a.axhline(273.15, color=NEUTRALS['graphite'], ls="--", lw=1.5,
                 label="Pure $\\mathrm{H_2O}$ melting (273.15 K)")
    ax_a.axhline(176.0, color=STRATA['magma'], ls=":", lw=1.8,
                 label="Eutectic floor $T_{\\mathrm{eutectic}}$ (176 K)")

    ax_a.fill_between(x_nh3 * 100.0, t_freeze, 273.15, color=STRATA['cobalt'], alpha=0.15)

    ax_a.set_xlim(0.0, 40.0)
    ax_a.set_ylim(160.0, 290.0)
    ax_a.set_xlabel("Ammonia concentration $X_{\\mathrm{NH_3}}$ [wt%]")
    ax_a.set_ylabel("Equilibrium freezing point $T$ [K]")
    ax_a.set_title("(a) Ammonia-water freezing depression and eutectic floor", loc="left", pad=8)
    ax_a.grid(True)
    ax_a.legend(loc="upper right", framealpha=0.92, fontsize=9.5)

    # -------------------------------------------------------------------------
    # Panel (b): Multi-Snowline Volatile Condensation Sequence
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    t_disk = np.linspace(10.0, 200.0, 500)
    condensed, t_conds = evaluate_condensation_profile(t_disk)

    colors = {
        'H2O': STRATA['cobalt'],
        'CO2': STRATA['amber'],
        'NH3': STRATA['plum'],
        'CH4': STRATA['gold'],
        'CO': STRATA['magma'],
        'N2': STRATA['ink']
    }

    for sp in ['H2O', 'CO2', 'NH3', 'CH4', 'CO', 'N2']:
        ax_b.plot(t_disk, condensed[sp] * 100.0, lw=2.0, color=colors[sp],
                  label=f"{sp} ($T_{{\\mathrm{{cond}}}} = {int(t_conds[sp])}$ K)")

    ax_b.set_xlim(10.0, 300.0)
    ax_b.set_ylim(-2.0, 125.0)
    ax_b.set_xlabel("Midplane disk temperature $T_{\\mathrm{disk}}$ [K]")
    ax_b.set_ylabel("Condensed ice fraction [wt% of accreted ice]")
    ax_b.set_title("(b) Multi-snowline volatile ice condensation sequence", loc="left", pad=8)
    ax_b.grid(True)
    # Legend placed in upper right where T > 165 K and y > 85 wt% (completely empty)
    ax_b.legend(loc="upper right", framealpha=0.92, fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (c): Refractory Element Delivery across Snowlines
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    t_disk_wide = np.linspace(10.0, 300.0, 400)
    c_refr = np.full_like(t_disk_wide, 0.60)
    c_ice_cond = np.zeros_like(t_disk_wide)
    c_ice_cond[t_disk_wide <= 75.0] += 0.20   # CO2 ice
    c_ice_cond[t_disk_wide <= 45.0] += 0.10   # CH4 ice
    c_ice_cond[t_disk_wide <= 25.0] += 0.10   # CO ice
    s_refr = np.full_like(t_disk_wide, 0.89)
    s_ice_cond = np.where(t_disk_wide <= 75.0, 0.11, 0.0)

    ax_c.plot(t_disk_wide, c_refr * 100.0, color=STRATA['ink'], lw=2.2,
              label="Refractory Carbon (Bergin+ 2026)")
    ax_c.plot(t_disk_wide, (c_refr + c_ice_cond) * 100.0, color=STRATA['ink'], lw=1.8, ls="--",
              label="Total Carbon (Refractory + Ice)")
    ax_c.plot(t_disk_wide, s_refr * 100.0, color=STRATA['amber'], lw=2.2,
              label="Refractory FeS Sulfur (Kama+ 2019)")
    ax_c.plot(t_disk_wide, (s_refr + s_ice_cond) * 100.0, color=STRATA['amber'], lw=1.8, ls="--",
              label="Total Sulfur (FeS + H2S Ice)")

    ax_c.set_xlim(10.0, 300.0)
    ax_c.set_ylim(35.0, 115.0)
    ax_c.set_xlabel("Midplane disk temperature $T_{\\mathrm{disk}}$ [K]")
    ax_c.set_ylabel("Delivered element fraction [% of cosmic budget]")
    ax_c.set_title("(c) Refractory carbon and sulfur delivery across disk orbits", loc="left", pad=8)
    ax_c.grid(True)
    ax_c.legend(loc="lower left", framealpha=0.92, fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (d): Refractory IOM Pyrolysis & Thermal Breakdown
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    t_body = np.linspace(300.0, 1100.0, 400)
    c_rem, c_graph, c_gas, n_rem, n_gas = evaluate_pyrolysis_profile(t_body)

    ax_d.plot(t_body, c_rem * 100.0, color=STRATA['cobalt'], lw=2.2,
              label="IOM Carbon ($C_{\\mathrm{refr}}$)")
    ax_d.plot(t_body, c_graph * 100.0, color=NEUTRALS['graphite'], lw=2.2,
              label="Graphite residue ($C_{\\mathrm{graph}}$)")
    ax_d.plot(t_body, c_gas * 100.0, color=STRATA['amber'], lw=2.2,
              label="Devolatilized gas ($C_{\\mathrm{gas}}$)")
    ax_d.plot(t_body, n_gas * 100.0, color=STRATA['magma'], lw=1.8, ls=":",
              label="Devolatilized gas ($N_{\\mathrm{gas}}$)")

    # Mass conservation verification curve
    c_tot = (c_rem + c_graph + c_gas) * 100.0
    ax_d.plot(t_body, c_tot, color=NEUTRALS['graphite'], ls="--", lw=1.0, alpha=0.6,
              label="Total carbon balance (100%)")

    ax_d.set_xlim(300.0, 1100.0)
    ax_d.set_ylim(-5.0, 155.0)
    ax_d.set_xlabel("Planetesimal interior temperature $T$ [K]")
    ax_d.set_ylabel("Phase partition [% of initial inventory]")
    ax_d.set_title("(d) Refractory organic pyrolysis and devolatilization", loc="left", pad=8)
    ax_d.grid(True)
    # Legend placed in upper right where y in [105, 150] (strictly above all curves)
    ax_d.legend(loc="upper right", framealpha=0.92, fontsize=8.5)

    # Save benchmark figures
    png_path = os.path.join(ASSETS_DIR, "volatile_mixtures_benchmark.png")
    pdf_path = os.path.join(ASSETS_DIR, "volatile_mixtures_benchmark.pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Generated {png_path} and {pdf_path}")

    # Export benchmark summary dataset
    json_path = os.path.join(OUTPUT_FILES_DIR, "volatile_mixtures_benchmark_data.json")
    benchmark_data = {
        "ammonia_water_freezing": {
            "pure_melting_point_K": 273.15,
            "eutectic_temperature_K": 176.0,
            "depression_slope_K_per_fraction": float((273.15 - 176.0) / 0.33),
            "eutectic_nh3_fraction": 0.33,
            "tested_nh3_fractions": [0.0, 0.05, 0.15, 0.33, 0.40],
            "computed_freezing_points_K": compute_freezing_curve(np.array([0.0, 0.05, 0.15, 0.33, 0.40])).tolist()
        },
        "multi_snowlines": {
            "H2O_condensation_T_K": 160.0,
            "NH3_condensation_T_K": 135.0,
            "CO2_condensation_T_K": 75.0,
            "H2S_condensation_T_K": 75.0,
            "CH4_condensation_T_K": 45.0,
            "CO_condensation_T_K": 25.0,
            "N2_condensation_T_K": 18.0
        },
        "refractory_fractions": {
            "carbon_refractory_fraction": 0.60,
            "sulfur_refractory_fraction": 0.89,
            "nitrogen_refractory_fraction": 0.10,
            "phosphorus_refractory_fraction": 0.98,
            "pyrolysis_threshold_T_K": 600.0,
            "graphite_residue_yield": 0.60,
            "carbon_mass_conservation_residual": float(np.max(np.abs(c_tot - 100.0)))
        }
    }
    with open(json_path, "w") as f:
        json.dump(benchmark_data, f, indent=2)
    print(f"Generated benchmark summary: {json_path}")


if __name__ == "__main__":
    main()
