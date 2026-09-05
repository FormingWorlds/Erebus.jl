#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for protoplanetary disk temperature evolution
and water snowline migration across orbital distance and host star mass.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Time array: 0 to 5 Myr
t_Myr = np.concatenate([
    np.linspace(0.0, 0.5, 600),
    np.linspace(0.505, 5.0, 500)
])
t_Myr = np.unique(t_Myr)

# Base physical constants & fiducial parameters (Sun-like star at 1 AU)
T_cloud = 30.0         # Molecular cloud temperature [K]
T_irr_1AU = 150.0      # Irradiation temperature at 1 AU for 1 M_sun [K]
T_peak_1AU = 520.0     # Peak midplane temperature at 1 AU for 1 M_sun [K]
t_peak_1AU = 0.12      # Peak time at 1 AU for 1 M_sun [Myr]
t_visc_0 = 0.25        # Viscous decay timescale for 1 M_sun [Myr]
gamma = 1.4            # Viscous power-law exponent
alpha = 2.0            # Infall rise power exponent

# Power-law scalings with orbital distance r [AU] and stellar mass M_star [M_sun]
q_irr = 3.0 / 7.0      # Flared disk irradiation: r^(-3/7)
q_visc = 0.75          # Viscous heating: r^(-3/4)
p_r_t = 0.25           # Radial peak time delay: r^(1/4)

p_M_irr = 0.25         # Stellar mass scaling for irradiation: M_star^0.25
p_M_visc = 0.30        # Viscous heating scaling with M_star: M_star^0.30
p_M_t = 0.40           # Infall duration scaling with cloud mass: M_cloud^0.40
p_M_visc_decay = 0.30  # Disk viscous lifetime scaling with mass: M_star^0.30

def disk_temperature_2A(t, r_au=2.5, M_star=1.0):
    """
    Evaluate Model 2 (Class 0 to II) disk temperature for t [Myr], r [AU], and M_star [M_sun].
    """
    T_irr = T_irr_1AU * (M_star**p_M_irr) * (r_au**(-q_irr))
    T_peak = T_peak_1AU * (M_star**p_M_visc) * (r_au**(-q_visc))
    t_peak = t_peak_1AU * (M_star**p_M_t) * (r_au**p_r_t)
    tau_star = 0.8 * t_peak

    x = np.maximum(1e-12, t / t_peak)
    norm = 1.0 + (alpha / gamma)
    f_acc = norm * (x**alpha) / (1.0 + (alpha / gamma) * (x**(alpha + gamma)))
    g_star = 1.0 - np.exp(-t / tau_star)
    T_eff_irr4 = T_cloud**4 + (T_irr**4 - T_cloud**4) * g_star

    T4 = T_eff_irr4 + (T_peak**4 - T_irr**4) * f_acc
    return np.maximum(T_cloud, T4**0.25)

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "disk_temperature_multidistance_multimass.png")

    fig, axs = plt.subplots(2, 2, figsize=(14, 10), dpi=200)

    # Panel 1: Different Orbital Distances around Sun-like Star (M_star = 1.0)
    radii = [0.5, 1.0, 1.5, 2.0, 2.5, 3.5, 5.0]
    cmap_r = plt.cm.plasma(np.linspace(0.1, 0.9, len(radii)))

    for r, col in zip(radii, cmap_r):
        T_curve = disk_temperature_2A(t_Myr, r_au=r, M_star=1.0)
        axs[0, 0].plot(t_Myr, T_curve, label=f"r = {r:.1f} AU", color=col, lw=2.0)

    axs[0, 0].axhline(273.15, color="blue", ls="--", lw=1.2, alpha=0.8, label="Ice Melting (273 K)")
    axs[0, 0].axhline(170.0, color="cyan", ls=":", lw=1.2, alpha=0.9, label="Snowline (170 K)")
    axs[0, 0].set_title(r"(A) Different Orbital Distances ($M_\star = 1.0\,M_\odot$)", fontsize=11, fontweight="bold")
    axs[0, 0].set_xlabel("Time after collapse onset [Myr]", fontsize=10)
    axs[0, 0].set_ylabel("Midplane Temperature [K]", fontsize=10)
    axs[0, 0].set_xlim(0, 3.0)
    axs[0, 0].set_ylim(0, 650)
    axs[0, 0].grid(True, alpha=0.3)
    axs[0, 0].legend(fontsize=8, loc="upper right")

    # Panel 2: Early Phase Zoom (0 - 0.6 Myr) around Sun-like Star
    for r, col in zip(radii, cmap_r):
        T_curve = disk_temperature_2A(t_Myr, r_au=r, M_star=1.0)
        axs[0, 1].plot(t_Myr, T_curve, label=f"r = {r:.1f} AU", color=col, lw=2.0)

    axs[0, 1].axhline(273.15, color="blue", ls="--", lw=1.2, alpha=0.8)
    axs[0, 1].axhline(170.0, color="cyan", ls=":", lw=1.2, alpha=0.9)
    axs[0, 1].set_title(r"(B) Early Accretion and Cooling (0 - 0.6 Myr)", fontsize=11, fontweight="bold")
    axs[0, 1].set_xlabel("Time after collapse onset [Myr]", fontsize=10)
    axs[0, 1].set_ylabel("Midplane Temperature [K]", fontsize=10)
    axs[0, 1].set_xlim(0, 0.6)
    axs[0, 1].set_ylim(0, 650)
    axs[0, 1].grid(True, alpha=0.3)

    # Panel 3: Different Stellar Masses at r = 2.5 AU
    masses = [0.1, 0.3, 1.0, 2.0]
    mass_labels = [
        r"$0.1\,M_\odot$ (Very low-mass M-dwarf)",
        r"$0.3\,M_\odot$ (Mid M-dwarf)",
        r"$1.0\,M_\odot$ (Solar analogue)",
        r"$2.0\,M_\odot$ (Intermediate mass)"
    ]
    cmap_m = ["#7570b3", "#e7298a", "#d95f02", "#1b9e77"]

    for m, lab, col in zip(masses, mass_labels, cmap_m):
        T_curve = disk_temperature_2A(t_Myr, r_au=2.5, M_star=m)
        axs[1, 0].plot(t_Myr, T_curve, label=lab, color=col, lw=2.2)

    axs[1, 0].axhline(273.15, color="blue", ls="--", lw=1.2, alpha=0.8, label="Ice Melting (273 K)")
    axs[1, 0].axhline(170.0, color="cyan", ls=":", lw=1.2, alpha=0.9, label="Snowline (170 K)")
    axs[1, 0].set_title(r"(C) Different Stellar Masses at Planetesimal Zone ($r = 2.5$ AU)", fontsize=11, fontweight="bold")
    axs[1, 0].set_xlabel("Time after collapse onset [Myr]", fontsize=10)
    axs[1, 0].set_ylabel("Midplane Temperature [K]", fontsize=10)
    axs[1, 0].set_xlim(0, 3.0)
    axs[1, 0].set_ylim(0, 350)
    axs[1, 0].grid(True, alpha=0.3)
    axs[1, 0].legend(fontsize=8, loc="upper right")

    # Panel 4: Water Snowline Migration over Time (r_snow where T = 170 K)
    r_grid = np.linspace(0.2, 8.0, 300)
    t_snow_eval = np.linspace(0.01, 3.0, 400)

    for m, lab, col in zip(masses, mass_labels, cmap_m):
        r_snow = []
        for t_val in t_snow_eval:
            T_profile = [disk_temperature_2A(t_val, r_au=r_val, M_star=m) for r_val in r_grid]
            idx = np.where(np.array(T_profile) <= 170.0)[0]
            if len(idx) > 0 and idx[0] > 0:
                i0 = idx[0] - 1
                i1 = idx[0]
                r_cross = r_grid[i0] + (170.0 - T_profile[i0]) * (r_grid[i1] - r_grid[i0]) / (T_profile[i1] - T_profile[i0])
                r_snow.append(r_cross)
            else:
                r_snow.append(np.nan)
        axs[1, 1].plot(t_snow_eval, r_snow, label=lab, color=col, lw=2.2)

    axs[1, 1].set_title(r"(D) Water Snowline ($T=170$ K) Outward Expansion & Inward Retreat", fontsize=11, fontweight="bold")
    axs[1, 1].set_xlabel("Time after collapse onset [Myr]", fontsize=10)
    axs[1, 1].set_ylabel("Water Snowline Location [AU]", fontsize=10)
    axs[1, 1].set_xlim(0, 3.0)
    axs[1, 1].set_ylim(0, 6.0)
    axs[1, 1].grid(True, alpha=0.3)
    axs[1, 1].legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    plt.savefig(out_fig)
    plt.close()
    print("Saved diagnostic figure to:", out_fig)

if __name__ == "__main__":
    main()
