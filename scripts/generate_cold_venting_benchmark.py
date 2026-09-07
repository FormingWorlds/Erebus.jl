#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for cold surface venting,
ice cold-trap Clausius-Clapeyron vapor pressure, disk dispersal transitions,
and coupled sublimation latent cooling in Erebus.jl.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Physical constants
L_sun = 3.828e26       # Solar luminosity [W]
sigma_sb = 5.670374419e-8 # Stefan-Boltzmann constant [W/(m^2 K^4)]
AU_m = 1.495978707e11   # 1 AU [m]
P0_triple = 611.66     # Triple-point vapor pressure [Pa]
T0_triple = 273.16     # Triple-point temperature [K]
L_sub = 2.83e6         # Latent heat of ice sublimation [J/kg]
Rv = 461.5             # Specific gas constant for water vapor [J/(kg K)]
rho_f = 1000.0         # Reference water fluid density [kg/m^3]

def compute_dispersal_weight(t_yr, t_disp_yr=3.0e6, dt_disp_yr=0.1e6):
    """Sigmoid dispersal weight w_disp in [0, 1]."""
    arg = -(t_yr - t_disp_yr) / dt_disp_yr
    arg = np.clip(arg, -100.0, 100.0)
    return 1.0 / (1.0 + np.exp(arg))

def compute_t_eq(r_au=2.7, albedo=0.06):
    """Solar radiative equilibrium temperature [K]."""
    d = r_au * AU_m
    return ((1.0 - albedo) * L_sun / (16.0 * np.pi * sigma_sb * d**2))**0.25

def compute_p_sat_ice(T_K):
    """Clausius-Clapeyron ice sublimation vapor pressure [Pa]."""
    arg = -(L_sub / Rv) * (1.0 / np.maximum(T_K, 1.0) - 1.0 / T0_triple)
    arg = np.clip(arg, -300.0, 50.0)
    return P0_triple * np.exp(arg)

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "cold_surface_venting_benchmark.png")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=300)
    plt.subplots_adjust(hspace=0.32, wspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Disk Dispersal & Ambient Conditions Evolution
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    t_myr = np.linspace(0.0, 5.0, 500)
    t_yr = t_myr * 1e6
    t_disp_myr = 3.0
    dt_disp_myr = 0.1
    w_disp = compute_dispersal_weight(t_yr, t_disp_myr * 1e6, dt_disp_myr * 1e6)

    # Disk temperature decay + dispersal to T_eq
    T_disk_init = 240.0
    T_disk_t = T_disk_init * (1.0 + t_myr / 0.5)**(-0.2)
    T_eq = compute_t_eq(r_au=2.7, albedo=0.06) # ~166.8 K
    T_amb = (1.0 - w_disp) * T_disk_t + w_disp * T_eq

    P_disk = 10.0 # Pa
    P_space = 1.0e-4 # Pa
    P_amb = (1.0 - w_disp) * P_disk + w_disp * P_space

    line1 = ax_a.plot(t_myr, T_amb, color="#1f77b4", lw=2.5, label="Ambient Temperature $T_{\\mathrm{amb}}(t)$")
    ax_a.axhline(T_eq, color="#1f77b4", ls="--", lw=1.2, alpha=0.7, label=f"Solar $T_{{\\mathrm{{eq}}}}$ ({T_eq:.1f} K)")
    ax_a.set_xlabel("Time [Myr]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Ambient Temperature [K]", fontsize=11, fontweight="bold", color="#1f77b4")
    ax_a.tick_params(axis="y", labelcolor="#1f77b4")
    ax_a.set_xlim(0, 5)
    ax_a.grid(True, ls=":", alpha=0.6)

    ax_a2 = ax_a.twinx()
    line2 = ax_a2.plot(t_myr, P_amb, color="#d62728", lw=2.0, ls="-.", label="Ambient Pressure $P_{\\mathrm{amb}}(t)$")
    ax_a2.set_ylabel("Ambient Pressure [Pa]", fontsize=11, fontweight="bold", color="#d62728")
    ax_a2.set_yscale("log")
    ax_a2.tick_params(axis="y", labelcolor="#d62728")

    # Combined legend
    lines = line1 + [ax_a.lines[1]] + line2
    labels = [l.get_label() for l in lines]
    ax_a.legend(lines, labels, loc="center right", fontsize=9, framealpha=0.9)
    ax_a.set_title("(a) Disk Dispersal & Ambient Evolution ($d = 2.7$ AU)", fontsize=12, fontweight="bold")

    # -------------------------------------------------------------------------
    # Panel (b): Water Ice Clausius-Clapeyron Cold Trap
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    T_eval = np.linspace(100.0, 273.16, 400)
    P_sat = compute_p_sat_ice(T_eval)

    ax_b.plot(T_eval, P_sat, color="#2ca02c", lw=2.5, label="Ice Sublimation $P_{\\mathrm{sat,ice}}(T)$")
    ax_b.axhline(P_disk, color="#d62728", ls=":", lw=1.5, label="Nebular Gas $P_{\\mathrm{disk}} = 10$ Pa")
    ax_b.axhline(1.0, color="#ff7f0e", ls=":", lw=1.5, label="Intermediate $P_{\\mathrm{amb}} = 1$ Pa")
    ax_b.axhline(P_space, color="#9467bd", ls=":", lw=1.5, label="Vacuum Floor $P_{\\mathrm{space}} = 10^{-4}$ Pa")

    # Mark crossover points
    T_cross_disk = T_eval[np.argmin(np.abs(P_sat - P_disk))]
    ax_b.scatter([T_cross_disk], [P_disk], color="#d62728", s=60, zorder=5)
    ax_b.annotate(f"Free Venting\n$T > {T_cross_disk:.0f}$ K",
                  xy=(T_cross_disk, P_disk), xytext=(T_cross_disk - 35, P_disk * 8),
                  arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2),
                  fontsize=9, fontweight="bold", color="#d62728")

    ax_b.set_yscale("log")
    ax_b.set_ylim(1e-12, 1e3)
    ax_b.set_xlim(100, 280)
    ax_b.set_xlabel("Surface Temperature $T_{\\mathrm{surf}}$ [K]", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Equilibrium Vapor Pressure [Pa]", fontsize=11, fontweight="bold")
    ax_b.grid(True, ls=":", alpha=0.6)
    ax_b.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax_b.set_title("(b) Ice Cold-Trap Equilibrium Vapor Pressure", fontsize=12, fontweight="bold")

    # -------------------------------------------------------------------------
    # Panel (c): Leaky Robin Darcy Venting Flux
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    delta_P = np.linspace(-0.5, 5.0, 300) # Overpressure in MPa
    delta_P_Pa = delta_P * 1e6
    delta_x = 4375.0 # Grid resolution [m]
    eta_f = 1e-3     # Fluid viscosity [Pa s]

    permeabilities = [1e-13, 1e-12, 1e-11, 1e-10]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    for k_v, col in zip(permeabilities, colors):
        # Flux q_vent = (k_vent / (eta_f * delta)) * max(0, Delta P)
        q_flux = (k_v / (eta_f * delta_x)) * np.maximum(0.0, delta_P_Pa) * 1e6 # um/s
        ax_c.plot(delta_P, q_flux, lw=2.0, color=col, label=f"$k_{{\\mathrm{{vent}}}} = 10^{{{int(np.log10(k_v))}}}$ m$^2$")

    ax_c.axvline(0.0, color="gray", ls="--", lw=1.0)
    ax_c.annotate("Zero Flux Gate\n($P_f \\leq P_{\\mathrm{vent}}$)",
                  xy=(0.0, 5.0), xytext=(-0.45, 12.0),
                  arrowprops=dict(arrowstyle="->", color="gray", lw=1.0),
                  fontsize=9, color="#555555")

    ax_c.set_xlim(-0.5, 5.0)
    ax_c.set_xlabel("Pore Fluid Overpressure $\\Delta P = P_f - P_{\\mathrm{vent}}$ [MPa]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Venting Discharge Flux $q_{\\mathrm{vent}}$ [$\\mu$m/s]", fontsize=11, fontweight="bold")
    ax_c.grid(True, ls=":", alpha=0.6)
    ax_c.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax_c.set_title("(c) Leaky Robin Boundary Drainage Flux", fontsize=12, fontweight="bold")

    # -------------------------------------------------------------------------
    # Panel (d): Marker Porosity Drainage & Latent Heat Sink
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    time_vent_kyr = np.linspace(0.0, 100.0, 400)
    time_vent_s = time_vent_kyr * 1000.0 * 3.15576e7
    phi_0 = 0.20
    phi_min = 1.0e-4

    # Exponential porosity drawdown to floor
    tau_drain = 25.0 * 1000.0 * 3.15576e7 # 25 kyr timescale
    phi_t = phi_min + (phi_0 - phi_min) * np.exp(-time_vent_s / tau_drain)
    S_vent = (phi_0 - phi_min) / tau_drain * np.exp(-time_vent_s / tau_drain)
    Q_lat_mag = L_sub * rho_f * S_vent # W/m^3

    ax_d.plot(time_vent_kyr, phi_t * 100.0, color="#1f77b4", lw=2.5, label="Boundary Porosity $\\phi(t)$")
    ax_d.axhline(phi_min * 100.0, color="#1f77b4", ls=":", lw=1.2, label=f"Porosity Floor ({phi_min*100:.2f}%)")
    ax_d.set_xlabel("Venting Elapsed Time [kyr]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel("Boundary Porosity [%]", fontsize=11, fontweight="bold", color="#1f77b4")
    ax_d.tick_params(axis="y", labelcolor="#1f77b4")
    ax_d.set_ylim(0, 22)
    ax_d.grid(True, ls=":", alpha=0.6)

    ax_d2 = ax_d.twinx()
    ax_d2.plot(time_vent_kyr, Q_lat_mag * 1e3, color="#9467bd", lw=2.0, ls="-.", label="Sublimation Cooling $|Q_{\\mathrm{lat}}|$")
    ax_d2.set_ylabel("Sublimation Heat Sink $|Q_{\\mathrm{lat}}|$ [mW/m$^3$]", fontsize=11, fontweight="bold", color="#9467bd")
    ax_d2.tick_params(axis="y", labelcolor="#9467bd")

    # Combined legend
    lines_d = ax_d.lines + ax_d2.lines
    labels_d = [l.get_label() for l in lines_d]
    ax_d.legend(lines_d, labels_d, loc="center right", fontsize=9, framealpha=0.9)
    ax_d.set_title("(d) Porosity Drainage & Latent Sublimation Cooling", fontsize=12, fontweight="bold")

    plt.savefig(out_fig, bbox_inches="tight")
    plt.close()
    print(f"Figure successfully generated: {out_fig}")

if __name__ == "__main__":
    main()
