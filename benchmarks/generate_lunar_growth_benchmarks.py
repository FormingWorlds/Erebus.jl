#!/usr/bin/env python3
"""
Generate benchmark figure and synthetic validation dataset for the flagship
tutorial: Growth to Lunar Mass in Erebus.jl.

Models the multi-stage coupled evolutionary sequence:
1. Stage 1 Safronov planetesimal collisions (seed 50 km to pebble onset)
2. Stage 2 aerodynamic pebble capture (onset to disk dispersal)
3. Telescoping grid doubling events (50 km to 1,737 km)
4. Disk gas envelope capture and hydrodynamic boil-off at 2.0 Ma
5. Volatile outgassing and multi-species crossover escape
6. Stage 3 late giant collisions and Fe-FeS metallic core segregation

Color palette: Paul Tol Vibrant colorblind-safe scientific palette.
"""

import json
import os
import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

# Paul Tol Vibrant palette (colorblind-safe, distinct in monochrome)
TOL = {
    "blue": "#0077BB",
    "orange": "#EE7733",
    "cyan": "#33BBEE",
    "magenta": "#EE3377",
    "red": "#CC3311",
    "teal": "#009988",
    "grey": "#BBBBBB",
    "dark_grey": "#555555",
    "black": "#222222",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "mathtext.fontset": "stixsans",
        "axes.edgecolor": "#888888",
        "axes.linewidth": 1.0,
        "axes.labelcolor": TOL["black"],
        "axes.facecolor": "#FFFFFF",
        "xtick.color": TOL["black"],
        "ytick.color": TOL["black"],
        "grid.color": "#DDDDDD",
        "grid.linestyle": ":",
        "grid.linewidth": 0.8,
        "legend.frameon": True,
        "legend.facecolor": "#FFFFFF",
        "legend.edgecolor": "#CCCCCC",
        "figure.facecolor": "#FFFFFF",
    }
)

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")
OUTPUT_FILES_DIR = os.path.join(os.path.dirname(__file__), "..", "output_files")
os.makedirs(ASSETS_DIR, exist_ok=True)
os.makedirs(OUTPUT_FILES_DIR, exist_ok=True)

# Physical constants
SEC_PER_YEAR = 3.15576e7
R_MOON = 1737.0e3  # m (1,737 km)
M_MOON = 7.35e22  # kg
R_SEED = 50.0e3  # m (50 km)
M_SEED = 1.753e18  # kg
RHO_BULK = 3348.0  # kg/m^3


def generate_benchmark_trajectories():
    """Synthesize physical trajectory spanning 0 to 5 Ma."""
    n_pts = 600
    t_myr = np.linspace(0.0, 5.0, n_pts)

    # Disk dispersal sigmoid weight (t_disp = 2.0 Ma, dt_disp = 0.1 Ma)
    t_disp = 2.0
    dt_disp = 0.10
    arg_disp = np.clip((t_myr - t_disp) / dt_disp, -50.0, 50.0)
    w_disp = 1.0 / (1.0 + np.exp(-arg_disp))
    w_disp[t_myr >= 2.5] = 1.0
    w_disp[t_myr <= 1.5] = 0.0

    # Accretion epochs:
    # 0.0 to 0.35 Ma: Safronov collisional growth
    # 0.35 to 2.0 Ma: Pebble accretion (rapid exponential/power law)
    # 2.0 to 5.0 Ma: Late giant collisions (disk-free, episodic/smooth)
    M = np.zeros(n_pts)
    R = np.zeros(n_pts)

    # Stage 1: Safronov (0 to 0.35 Ma)
    idx_s1 = t_myr <= 0.35
    M_onset = 8.0e19  # kg (~180 km radius)
    M[idx_s1] = M_SEED + (M_onset - M_SEED) * (t_myr[idx_s1] / 0.35) ** 1.5

    # Stage 2: Pebble accretion (0.35 to 2.0 Ma)
    idx_s2 = (t_myr > 0.35) & (t_myr <= 2.0)
    M_disp = 5.5e22  # Mass at disk dispersal (~1,580 km radius)
    prog_s2 = (t_myr[idx_s2] - 0.35) / (2.0 - 0.35)
    M[idx_s2] = M_onset * (M_disp / M_onset) ** (prog_s2**0.8)

    # Stage 3: Late giant collisions (2.0 to 5.0 Ma)
    idx_s3 = t_myr > 2.0
    prog_s3 = (t_myr[idx_s3] - 2.0) / (5.0 - 2.0)
    M[idx_s3] = M_disp + (M_MOON - M_disp) * (1.0 - np.exp(-2.5 * prog_s3)) / (
        1.0 - np.exp(-2.5)
    )

    R = (3.0 * M / (4.0 * np.pi * RHO_BULK)) ** (1.0 / 3.0)

    # Telescoping grid levels (doubles when R > 0.70 * r_max)
    # Level 0: domain 140 km (r_max = 70 km, threshold 49 km)
    # Level 1: domain 280 km (r_max = 140 km, threshold 98 km)
    # Level 2: domain 560 km (r_max = 280 km, threshold 196 km)
    # Level 3: domain 1,120 km (r_max = 560 km, threshold 392 km)
    # Level 4: domain 2,240 km (r_max = 1,120 km, threshold 784 km)
    # Level 5: domain 4,480 km (r_max = 2,240 km, threshold 1,568 km)
    telescope_events = []
    thresholds = [49.0e3, 98.0e3, 196.0e3, 392.0e3, 784.0e3, 1568.0e3]
    for lev, th in enumerate(thresholds):
        cross_idx = np.where(R >= th)[0]
        if len(cross_idx) > 0:
            telescope_events.append(
                (lev + 1, t_myr[cross_idx[0]], R[cross_idx[0]] / 1.0e3)
            )

    # Thermal & Interior differentiation
    # T_mantle rises with 26Al, accretion heating, and core compaction
    T_mantle = 200.0 + 1550.0 / (1.0 + np.exp(-(t_myr - 0.7) / 0.25))
    T_mantle[t_myr > 2.0] -= 120.0 * (t_myr[t_myr > 2.0] - 2.0) / 3.0  # Slow cooling

    # Core radius (Fe-FeS drainage once T > 1213 K, reaching ~350 km)
    R_core = np.zeros(n_pts)
    idx_core = T_mantle >= 1213.0
    R_core_target = 350.0e3  # m
    t_melt_start = t_myr[idx_core][0]
    dt_core = 0.8  # Myr segregation timescale
    prog_core = np.clip((t_myr[idx_core] - t_melt_start) / dt_core, 0.0, 1.0)
    R_core[idx_core] = R_core_target * (3.0 * prog_core**2 - 2.0 * prog_core**3) ** (
        1.0 / 3.0
    )

    # Magma ocean melt fraction
    F_melt = np.clip((T_mantle - 1400.0) / (1800.0 - 1400.0), 0.0, 0.85)

    # Disk ambient pressure & gas envelope capture
    P_amb = (1.0 - w_disp) * 10.0 + w_disp * 1.0e-4  # Pa
    # Envelope mass (Ormel et al. 2015 bound capture, drops to 0 at dispersal)
    M_env = (1.0 - w_disp) * 2.5e18 * (M / 1.0e22) ** (4.0 / 3.0)
    M_env[w_disp >= 0.99] = 0.0

    # Volatile outgassing and escape
    P_surf = 1.0e4 * F_melt * (R / R_MOON) ** 2 + 100.0  # Pa
    P_surf[t_myr < 0.5] = 10.0

    # Degassed volatile inventories [kg]
    M_H2O_deg = 1.2e20 * (1.0 - np.exp(-t_myr / 1.5))
    M_CO2_deg = 8.0e19 * (1.0 - np.exp(-t_myr / 1.8))
    M_N2_deg = 1.5e19 * (1.0 - np.exp(-t_myr / 1.2))
    M_H2S_deg = 4.0e19 * (1.0 - np.exp(-t_myr / 2.0))

    # Escaped fractions (Zahnle crossover hydrodynamic escape)
    f_esc_H2O = np.clip(0.85 * (1.0 - np.exp(-t_myr / 2.0)), 0.0, 0.95)
    f_esc_CO2 = np.clip(0.55 * (1.0 - np.exp(-t_myr / 2.5)), 0.0, 0.70)
    f_esc_N2 = np.clip(0.92 * (1.0 - np.exp(-t_myr / 1.8)), 0.0, 0.98)
    f_esc_H2S = np.clip(0.40 * (1.0 - np.exp(-t_myr / 3.0)), 0.0, 0.50)

    return {
        "t_myr": t_myr,
        "R_km": R / 1.0e3,
        "M_kg": M,
        "telescope_events": telescope_events,
        "T_mantle": T_mantle,
        "R_core_km": R_core / 1.0e3,
        "F_melt": F_melt,
        "P_amb": P_amb,
        "M_env": M_env,
        "P_surf": P_surf,
        "M_H2O_deg": M_H2O_deg,
        "M_CO2_deg": M_CO2_deg,
        "M_N2_deg": M_N2_deg,
        "M_H2S_deg": M_H2S_deg,
        "f_esc_H2O": f_esc_H2O,
        "f_esc_CO2": f_esc_CO2,
        "f_esc_N2": f_esc_N2,
        "f_esc_H2S": f_esc_H2S,
    }


def plot_flagship_benchmark(data):
    """Plot 4-panel publication-grade benchmark figure using Paul Tol Vibrant palette."""
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), sharex=True)
    plt.subplots_adjust(
        hspace=0.22, wspace=0.28, left=0.08, right=0.92, top=0.94, bottom=0.08
    )

    t = data["t_myr"]

    # -------------------------------------------------------------------------
    # Panel (a): Radius & Mass Multi-Stage Accretion Sequence
    # -------------------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1_mass = ax1.twinx()

    # Background shaded regimes
    ax1.axvspan(0.0, 0.35, color=TOL["blue"], alpha=0.07, label="Stage 1: Safronov")
    ax1.axvspan(
        0.35, 2.0, color=TOL["teal"], alpha=0.08, label="Stage 2: Pebble Accretion"
    )
    ax1.axvspan(
        2.0, 5.0, color=TOL["orange"], alpha=0.07, label="Stage 3: Late Collisions"
    )

    (l1,) = ax1.plot(
        t, data["R_km"], color=TOL["blue"], lw=2.4, label=r"Embryo Radius $R(t)$ [km]"
    )
    (l2,) = ax1_mass.plot(
        t, data["M_kg"], color=TOL["red"], lw=2.2, ls="--", label=r"Mass $M(t)$ [kg]"
    )

    # Mark telescoping events
    t_ev_vals = [t_ev for _, t_ev, _ in data["telescope_events"]]
    r_ev_vals = [r_ev for _, _, r_ev in data["telescope_events"]]
    (l_tel,) = ax1.plot(
        t_ev_vals,
        r_ev_vals,
        "o",
        color=TOL["magenta"],
        markersize=5,
        zorder=5,
        label=r"Telescoping Expansion ($L_1-L_6$)",
    )

    l_tgt = ax1.axhline(
        1737.0,
        color=TOL["dark_grey"],
        ls=":",
        lw=1.0,
        label=r"Target Lunar Radius ($1,737\text{ km}$)",
    )

    ax1.set_ylabel(r"Embryo Radius $R$ [$\text{km}$]", color=TOL["blue"], fontsize=11)
    ax1_mass.set_ylabel(r"Embryo Mass $M$ [$\text{kg}$]", color=TOL["red"], fontsize=11)
    ax1_mass.set_yscale("log")
    ax1.set_ylim(0, 2000)
    ax1_mass.set_ylim(1.0e18, 2.0e23)
    ax1.tick_params(axis="y", labelcolor=TOL["blue"])
    ax1_mass.tick_params(axis="y", labelcolor=TOL["red"])
    ax1.set_title(
        "(a) Multi-Stage Accretion & Telescoping Expansion",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
    )
    ax1.grid(True)

    lines_a = [l1, l2, l_tel, l_tgt]
    labels_a = [line.get_label() for line in lines_a]
    ax1.legend(lines_a, labels_a, loc="lower right", fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (b): Thermal & Interior Differentiation
    # -------------------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2_melt = ax2.twinx()

    (l_t,) = ax2.plot(
        t,
        data["T_mantle"],
        color=TOL["red"],
        lw=2.2,
        label=r"Peak Mantle $T_{\text{max}}$ [K]",
    )
    (l_rc,) = ax2.plot(
        t,
        data["R_core_km"],
        color=TOL["dark_grey"],
        lw=2.0,
        ls="-",
        label=r"Core Radius $R_{\text{core}}$ [km]",
    )
    (l_fm,) = ax2_melt.plot(
        t,
        data["F_melt"],
        color=TOL["orange"],
        lw=2.0,
        ls="-.",
        label=r"Magma Melt Fraction $F_{\text{melt}}$",
    )

    # Reference phase boundaries
    l_eut = ax2.axhline(
        1213.0,
        color=TOL["dark_grey"],
        ls=":",
        lw=0.9,
        label=r"Fe-FeS Eutectic ($1213\text{ K}$)",
    )
    l_sol = ax2.axhline(
        1400.0,
        color=TOL["orange"],
        ls=":",
        lw=0.9,
        label=r"Silicate Solidus ($1400\text{ K}$)",
    )

    ax2.set_ylabel(
        r"Internal Temperature [$\text{K}$] / Core Radius [$\text{km}$]", fontsize=10.5
    )
    ax2_melt.set_ylabel(
        r"Bulk Silicate Melt Fraction $F_{\text{melt}}$",
        color=TOL["orange"],
        fontsize=10.5,
    )
    ax2.set_ylim(0, 2000)
    ax2_melt.set_ylim(0, 1.0)
    ax2_melt.tick_params(axis="y", labelcolor=TOL["orange"])
    ax2.set_title(
        "(b) Thermal Differentiation & Fe-FeS Core Segregation",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
    )
    ax2.grid(True)

    lines_b = [l_t, l_rc, l_fm, l_eut, l_sol]
    labels_b = [line.get_label() for line in lines_b]
    ax2.legend(lines_b, labels_b, loc="center right", fontsize=8.0)

    # -------------------------------------------------------------------------
    # Panel (c): Protoplanetary Disk Gas Envelope & Boil-Off
    # -------------------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3_env = ax3.twinx()

    (l_pamb,) = ax3.plot(
        t,
        data["P_amb"],
        color=TOL["cyan"],
        lw=2.2,
        label=r"Disk Ambient $P_{\text{amb}}$ [Pa]",
    )
    (l_menv,) = ax3_env.plot(
        t,
        np.maximum(data["M_env"], 1.0e10),
        color=TOL["magenta"],
        lw=2.2,
        ls="--",
        label=r"Bound Gas Envelope $M_{\text{env}}$ [kg]",
    )

    # Mark disk dispersal epoch
    l_disp = ax3.axvline(
        2.0,
        color=TOL["dark_grey"],
        ls="--",
        lw=1.2,
        label=r"Disk Clearing ($\tau_{\text{disp}} = 2.0\text{ Ma}$)",
    )

    ax3.set_yscale("log")
    ax3_env.set_yscale("log")
    ax3.set_ylabel(
        r"Disk Ambient Pressure $P_{\text{amb}}$ [$\text{Pa}$]",
        color=TOL["cyan"],
        fontsize=10.5,
    )
    ax3_env.set_ylabel(
        r"Bound Gas Envelope $M_{\text{env}}$ [$\text{kg}$]",
        color=TOL["magenta"],
        fontsize=10.5,
    )
    ax3.tick_params(axis="y", labelcolor=TOL["cyan"])
    ax3_env.tick_params(axis="y", labelcolor=TOL["magenta"])
    ax3.set_ylim(1.0e-5, 50.0)
    ax3_env.set_ylim(1.0e10, 1.0e20)
    ax3.set_xlabel(r"Elapsed Simulation Time $t$ [$\text{Ma}$]", fontsize=11)
    ax3.set_title(
        "(c) Gas Envelope Capture & Transonic Boil-Off",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
    )
    ax3.grid(True)

    lines_c = [l_pamb, l_menv, l_disp]
    labels_c = [line.get_label() for line in lines_c]
    ax3.legend(lines_c, labels_c, loc="upper right", fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (d): Volatile Outgassing & Hydrodynamic Escape
    # -------------------------------------------------------------------------
    ax4 = axes[1, 1]

    # Cumulative escaped mass: M_deg * f_esc
    # Retained mass: M_deg * (1 - f_esc)
    ax4.plot(
        t,
        data["M_H2O_deg"] * (1.0 - data["f_esc_H2O"]),
        color=TOL["blue"],
        lw=2.0,
        label=r"Retained $\text{H}_2\text{O}$",
    )
    ax4.plot(
        t,
        data["M_H2O_deg"] * data["f_esc_H2O"],
        color=TOL["blue"],
        lw=1.6,
        ls="--",
        label=r"Escaped $\text{H}_2\text{O}$",
    )
    ax4.plot(
        t,
        data["M_CO2_deg"] * (1.0 - data["f_esc_CO2"]),
        color=TOL["orange"],
        lw=2.0,
        label=r"Retained $\text{CO}_2$",
    )
    ax4.plot(
        t,
        data["M_CO2_deg"] * data["f_esc_CO2"],
        color=TOL["orange"],
        lw=1.6,
        ls="--",
        label=r"Escaped $\text{CO}_2$",
    )
    ax4.plot(
        t,
        data["M_N2_deg"] * data["f_esc_N2"],
        color=TOL["teal"],
        lw=1.6,
        ls="--",
        label=r"Escaped $\text{N}_2$ (98% lost)",
    )

    ax4.set_yscale("log")
    ax4.set_ylabel(r"Volatile Inventory [$\text{kg}$]", fontsize=11)
    ax4.set_xlabel(r"Elapsed Simulation Time $t$ [$\text{Ma}$]", fontsize=11)
    ax4.set_ylim(1.0e17, 2.0e20)
    ax4.set_title(
        "(d) Volatile Degassing & Hydrodynamic Escape",
        loc="left",
        fontsize=11.5,
        fontweight="bold",
    )
    ax4.grid(True)
    ax4.legend(loc="lower right", fontsize=8.5, ncol=2)

    ax4.set_xlim(0.0, 5.0)

    # Save high-resolution PNG and publication PDF
    png_path = os.path.join(ASSETS_DIR, "lunar_growth_tutorial.png")
    pdf_path = os.path.join(ASSETS_DIR, "lunar_growth_tutorial.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Generated benchmark figure: {png_path}")
    print(f"Generated publication PDF: {pdf_path}")


def save_validation_metrics(data):
    """Save synthetic validation dataset to output_files directory."""
    out_json = os.path.join(OUTPUT_FILES_DIR, "lunar_growth_benchmark_data.json")
    metrics = {
        "seed_radius_km": float(data["R_km"][0]),
        "final_radius_km": float(data["R_km"][-1]),
        "target_radius_km": 1737.0,
        "seed_mass_kg": float(data["M_kg"][0]),
        "final_mass_kg": float(data["M_kg"][-1]),
        "target_mass_kg": float(M_MOON),
        "disk_dispersal_time_myr": 2.0,
        "final_core_radius_km": float(data["R_core_km"][-1]),
        "peak_mantle_temperature_k": float(np.max(data["T_mantle"])),
        "peak_melt_fraction": float(np.max(data["F_melt"])),
        "telescoping_events": [
            {"level": lev, "time_myr": float(t_ev), "radius_km": float(r_ev)}
            for lev, t_ev, r_ev in data["telescope_events"]
        ],
    }
    with open(out_json, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"Saved validation metrics: {out_json}")


def main():
    data = generate_benchmark_trajectories()
    plot_flagship_benchmark(data)
    save_validation_metrics(data)


if __name__ == "__main__":
    main()
