#!/usr/bin/env python3
"""
Generate diagnostic figures for atmospheric Jeans kinetic escape parameterization.

Illustrates:
(a) Jeans parameter lambda vs planetesimal radius for H2O, N2, CO2.
(b) Normalized Jeans flux (1 + lambda) * exp(-lambda) vs lambda.
(c) Escape timescale vs planetesimal radius.
(d) Coupled atmospheric mass evolution and retention vs escape across planetesimal sizes.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Physical constants
G = 6.67430e-11  # m^3 / (kg s^2)
kB = 1.380649e-23  # J / K
NA = 6.02214076e23  # mol^-1

# Molecular masses [kg]
MASSES = {
    "H2O": 18.01528e-3 / NA,
    "N2": 28.01340e-3 / NA,
    "CO2": 44.00950e-3 / NA,
}


def compute_escape_velocity(M, r):
    return np.sqrt(2.0 * G * M / r)


def compute_thermal_velocity(T, m):
    return np.sqrt(2.0 * kB * T / m)


def compute_jeans_parameter(M, r, T, m):
    return (G * M * m) / (kB * T * r)


def compute_jeans_flux_ratio(lam):
    return (1.0 + lam) * np.exp(-lam)


def generate_benchmark_figure(output_path):
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), dpi=300)
    plt.subplots_adjust(hspace=0.28, wspace=0.26)

    rho_bulk = 2500.0  # kg / m^3
    T_exo = 200.0  # K

    # -------------------------------------------------------------------------
    # Panel (a): Jeans Parameter vs Radius
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    R_km = np.linspace(10.0, 2000.0, 500)
    R_m = R_km * 1.0e3
    M_planet = (4.0 / 3.0) * np.pi * (R_m ** 3) * rho_bulk

    colors = {"H2O": "#1f77b4", "N2": "#2ca02c", "CO2": "#d62728"}
    labels = {
        "H2O": r"$\mathrm{H}_2\mathrm{O}$ (18 g/mol)",
        "N2": r"$\mathrm{N}_2$ (28 g/mol)",
        "CO2": r"$\mathrm{CO}_2$ (44 g/mol)",
    }

    for sp, m_sp in MASSES.items():
        lam = compute_jeans_parameter(M_planet, R_m, T_exo, m_sp)
        ax_a.plot(R_km, lam, color=colors[sp], lw=2.2, label=labels[sp])

    ax_a.axhline(1.0, color="gray", ls="--", lw=1.2, alpha=0.8, label=r"$\lambda = 1$ (Effusion boundary)")
    ax_a.axhline(10.0, color="black", ls=":", lw=1.2, alpha=0.8, label=r"$\lambda = 10$ (Retention threshold)")
    ax_a.set_yscale("log")
    ax_a.set_xlabel("Planetesimal Radius $R$ [km]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Jeans Parameter $\\lambda$", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Jeans Parameter vs Body Size ($T = 200$ K)", fontsize=12, fontweight="bold")
    ax_a.set_xlim(10, 2000)
    ax_a.set_ylim(1e-3, 1e2)
    ax_a.grid(True, ls=":", alpha=0.6)
    ax_a.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Kinetic Flux Ratio vs Lambda
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    lam_range = np.linspace(0.0, 20.0, 500)
    flux_ratio = compute_jeans_flux_ratio(lam_range)

    ax_b.plot(lam_range, flux_ratio, color="#800080", lw=2.5, label=r"$\Phi_{\mathrm{Jeans}} / \Phi_{\mathrm{eff}} = (1 + \lambda) e^{-\lambda}$")
    ax_b.axvline(1.0, color="gray", ls="--", lw=1.2, alpha=0.8)
    ax_b.axvline(10.0, color="black", ls=":", lw=1.2, alpha=0.8)
    ax_b.set_yscale("log")
    ax_b.set_xlabel(r"Jeans Parameter $\lambda$", fontsize=11, fontweight="bold")
    ax_b.set_ylabel(r"Escape Flux Ratio $\Phi / \Phi_{\mathrm{eff}}$", fontsize=11, fontweight="bold")
    ax_b.set_title("(b) Gravitational Escape Flux Suppression", fontsize=12, fontweight="bold")
    ax_b.set_xlim(0, 20)
    ax_b.set_ylim(1e-8, 2.0)
    ax_b.grid(True, ls=":", alpha=0.6)
    ax_b.legend(loc="upper right", fontsize=9.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Characteristic Escape Timescale vs Radius
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    sec_per_year = 365.25 * 86400.0

    for sp in ["H2O", "CO2"]:
        m_sp = MASSES[sp]
        v_th = compute_thermal_velocity(T_exo, m_sp)
        lam = compute_jeans_parameter(M_planet, R_m, T_exo, m_sp)
        effusion_factor = (1.0 + lam) * np.exp(-lam)
        # k_escape [s^-1]
        k_esc = (v_th * lam / (2.0 * np.sqrt(np.pi) * R_m)) * effusion_factor
        tau_years = np.where(k_esc > 1e-30, 1.0 / (k_esc * sec_per_year), 1e12)
        ax_c.plot(R_km, tau_years, color=colors[sp], lw=2.2, label=f"{labels[sp]}")

    ax_c.set_yscale("log")
    ax_c.set_xlabel("Planetesimal Radius $R$ [km]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel(r"Characteristic Escape Timescale $\tau_{\mathrm{esc}}$ [yr]", fontsize=11, fontweight="bold")
    ax_c.set_title("(c) Atmospheric Loss Timescale vs Body Size", fontsize=12, fontweight="bold")
    ax_c.set_xlim(10, 2000)
    ax_c.set_ylim(1e-4, 1e12)
    ax_c.grid(True, ls=":", alpha=0.6)
    ax_c.legend(loc="upper left", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (d): Dynamic Atmospheric Inventory Evolution
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    time_yr = np.linspace(0.0, 5.0, 500)
    time_s = time_yr * sec_per_year
    dt_step = time_s[1] - time_s[0]

    # Three body sizes:
    # R = 50 km (rapid effusion, complete escape)
    # R = 1350 km (intermediate body, dynamic partitioning, tau ~ 1.7 yr)
    # R = 2000 km (planetary embryo, strong gravitational retention)
    vent_flux_rate = 100.0  # kg / s steady venting

    bodies = [
        {"R": 50.0, "col": "#1f77b4", "label_atm": None, "label_esc": r"Escaped to Space ($R=50$ km)"},
        {"R": 1350.0, "col": "#2ca02c", "label_atm": r"Atmosphere ($R=1350$ km)", "label_esc": r"Escaped to Space ($R=1350$ km)"},
        {"R": 2000.0, "col": "#d62728", "label_atm": r"Atmosphere ($R=2000$ km)", "label_esc": None},
    ]

    for b in bodies:
        R_val = b["R"] * 1.0e3
        M_val = (4.0 / 3.0) * np.pi * (R_val ** 3) * rho_bulk
        m_sp = MASSES["H2O"]
        v_th = compute_thermal_velocity(T_exo, m_sp)
        lam = compute_jeans_parameter(M_val, R_val, T_exo, m_sp)
        effusion_factor = (1.0 + lam) * np.exp(-lam)
        k_esc = (v_th * lam / (2.0 * np.sqrt(np.pi) * R_val)) * effusion_factor

        M_atm = np.zeros_like(time_s)
        M_escaped = np.zeros_like(time_s)
        cur_M = 0.0
        cur_esc = 0.0

        for i in range(1, len(time_s)):
            x = k_esc * dt_step
            if x < 1e-6:
                int_fac = dt_step * (1.0 - 0.5 * x + (x**2) / 6.0)
                next_M = cur_M * np.exp(-x) + vent_flux_rate * int_fac
            elif x > 40.0:
                next_M = vent_flux_rate / k_esc
            else:
                next_M = cur_M * np.exp(-x) + (vent_flux_rate / k_esc) * (1.0 - np.exp(-x))
            delta_esc = max(0.0, cur_M + vent_flux_rate * dt_step - next_M)
            cur_esc += delta_esc
            cur_M = next_M
            M_atm[i] = cur_M
            M_escaped[i] = cur_esc

        if b["label_atm"]:
            ax_d.plot(time_yr, M_atm * 1e-9, color=b["col"], lw=2.2, ls="-", label=b["label_atm"])
        if b["label_esc"]:
            ls_e = "--" if b["R"] == 50.0 else ":"
            lw_e = 2.0 if b["R"] == 50.0 else 1.8
            ax_d.plot(time_yr, M_escaped * 1e-9, color=b["col"], lw=lw_e, ls=ls_e, label=b["label_esc"])

    ax_d.set_xlabel("Time [yr]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel(r"Volatile Inventory [$10^9$ kg]", fontsize=11, fontweight="bold")
    ax_d.set_title(r"(d) Inventory Partitioning (Venting = $100$ kg/s)", fontsize=12, fontweight="bold")
    ax_d.set_xlim(0, 5)
    ax_d.grid(True, ls=":", alpha=0.6)
    ax_d.legend(loc="upper left", fontsize=8.5, framealpha=0.9)

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Generated benchmark figure: {output_path}")


if __name__ == "__main__":
    out_png = os.path.join(
        os.path.dirname(__file__), "..", "docs", "src", "assets", "jeans_escape_benchmark.png"
    )
    generate_benchmark_figure(out_png)
