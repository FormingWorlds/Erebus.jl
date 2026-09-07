#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for cold lid hydrofracture breaching,
cryogenic pore ice permeability sealing, and episodic venting dynamics in Erebus.jl.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def compute_ice_sealed_permeability(k0, T_K, T_freeze=273.15, dT_seal=10.0, k_min_ratio=1.0e-6):
    """Compute cryogenic ice-sealed permeability [m^2]."""
    k = np.full_like(T_K, float(k0))
    subfreeze = T_K < T_freeze
    if np.any(subfreeze):
        ratio = (1.0 - k_min_ratio) * np.exp(-(T_freeze - T_K[subfreeze]) / dT_seal) + k_min_ratio
        k[subfreeze] = k0 * ratio
    return k

def compute_hydrofracture_permeability(k0, Peff_Pa, sigma_t_Pa=1.0e7, kappa_frac=1.0e3, gamma=1.0, kmax=1.0e-9):
    """Compute hydrofracture-enhanced permeability [m^2]."""
    overpressure = np.maximum(0.0, -Peff_Pa - sigma_t_Pa)
    norm_op = overpressure / sigma_t_Pa
    factor = 1.0 + kappa_frac * (norm_op ** gamma)
    return np.clip(k0 * factor, k0, kmax)

def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "hydrofracture_venting_benchmark.png")

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), dpi=300)
    plt.subplots_adjust(hspace=0.34, wspace=0.38)

    # -------------------------------------------------------------------------
    # Panel (a): Cryogenic Pore Ice Permeability Sealing
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    T_range = np.linspace(50.0, 300.0, 500)
    k0_ref = 1.0e-11 # m^2

    dt_seal_list = [5.0, 10.0, 20.0]
    colors_a = ["#1f77b4", "#2ca02c", "#ff7f0e"]
    for dt_s, col in zip(dt_seal_list, colors_a):
        k_sealed = compute_ice_sealed_permeability(k0_ref, T_range, dT_seal=dt_s, k_min_ratio=1.0e-6)
        ax_a.plot(T_range, k_sealed / k0_ref, color=col, lw=2.2,
                  label=f"$\\Delta T_{{\\mathrm{{seal}}}} = {dt_s:.0f}$ K")

    ax_a.axvline(273.15, color="black", ls="--", lw=1.2, alpha=0.8,
                 label="Water Freezing ($T_{\\mathrm{freeze}} = 273.15$ K)")
    ax_a.axhline(1.0e-6, color="gray", ls=":", lw=1.2, alpha=0.8,
                 label="Cryogenic Floor ($k_{\\mathrm{min}}/k_0 = 10^{-6}$)")
    ax_a.set_yscale("log")
    ax_a.set_xlabel("Surface Temperature $T$ [K]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Permeability Ratio $k_{\\mathrm{eff}} / k_0$ [-]", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Cryogenic Pore Ice Sealing vs Temperature", fontsize=12, fontweight="bold")
    ax_a.set_xlim(50, 300)
    ax_a.set_ylim(5.0e-7, 2.0)
    ax_a.grid(True, which="both", ls=":", alpha=0.6)
    ax_a.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Surface Boundary Regime Diagram (T_surf vs P_eff)
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    T_grid = np.linspace(100.0, 320.0, 250)
    Peff_grid_MPa = np.linspace(-30.0, 10.0, 250)
    TT, PP = np.meshgrid(T_grid, Peff_grid_MPa)

    # Effective permeability across phase space
    sigma_t_val = 10.0 # MPa
    k_eff_map = np.zeros_like(TT)
    for i in range(TT.shape[0]):
        for j in range(TT.shape[1]):
            t_val = TT[i, j]
            peff_val = PP[i, j] * 1.0e6
            if peff_val <= -sigma_t_val * 1.0e6:
                # Breached
                k_eff_map[i, j] = compute_hydrofracture_permeability(
                    k0_ref, peff_val, sigma_t_Pa=sigma_t_val*1.0e6, kappa_frac=1.0e3, gamma=1.0, kmax=1.0e-9
                )
            else:
                # Unbreached with ice sealing
                k_eff_map[i, j] = compute_ice_sealed_permeability(
                    k0_ref, np.array([t_val]), T_freeze=273.15, dT_seal=10.0, k_min_ratio=1.0e-6
                )[0]

    cs = ax_b.contourf(TT, PP, np.log10(k_eff_map), levels=np.linspace(-17, -9, 17),
                       cmap="viridis", alpha=0.85)
    cbar = fig.colorbar(cs, ax=ax_b, pad=0.03)
    cbar.set_label(r"$\log_{10}(k_{\mathrm{eff}}\ [\mathrm{m}^2])$", fontsize=10, fontweight="bold")

    # Boundary lines
    ax_b.axhline(-sigma_t_val, color="red", ls="--", lw=2.0,
                 label=r"Tensile Rupture ($P_{\mathrm{eff}} = -\sigma_t = -10$ MPa)")
    ax_b.axvline(273.15, color="cyan", ls="--", lw=2.0,
                 label=r"Cryogenic Boundary ($T = 273.15$ K)")

    # Annotate quadrants with well-bounded positions
    ax_b.text(180, 4.0, "Cryogenically Sealed\nIntact Lid", color="white",
              fontsize=9.0, fontweight="bold", ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.6))
    ax_b.text(180, -20.0, "Breached Cold Lid\nHydrofracture Vent", color="white",
              fontsize=9.0, fontweight="bold", ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="darkred", alpha=0.6))
    ax_b.text(290, 4.0, "Warm Darcy\nPermeable Sink", color="white",
              fontsize=9.0, fontweight="bold", ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="navy", alpha=0.6))
    ax_b.text(290, -20.0, "Warm Breached\nHydrofracture", color="white",
              fontsize=9.0, fontweight="bold", ha="center", va="center",
              bbox=dict(boxstyle="round,pad=0.3", facecolor="darkgreen", alpha=0.6))

    ax_b.set_xlabel("Surface Temperature $T_{\\mathrm{surf}}$ [K]", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Effective Pressure $P_{\\mathrm{eff}} = P_t - P_f$ [MPa]", fontsize=11, fontweight="bold")
    ax_b.set_title("(b) Venting Regime & Effective Permeability Map", fontsize=12, fontweight="bold")
    ax_b.set_xlim(100, 320)
    ax_b.set_ylim(-30, 10)
    ax_b.legend(loc="lower left", fontsize=8.0, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Episodic Rupture, Venting Pulse & Resealing Cycle
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    time_hr = np.linspace(0.0, 120.0, 600)
    P_litho = 5.0 # MPa
    sigma_t_sc = 10.0 # MPa
    P_crit = P_litho + sigma_t_sc # 15 MPa

    # Generate saw-tooth pressure pulses
    period = 30.0 # hours
    t_mod = time_hr % period
    P_f = np.zeros_like(time_hr)
    for idx, (t, tm) in enumerate(zip(time_hr, t_mod)):
        if tm <= 24.0:
            P_f[idx] = 6.0 + 10.0 * (tm / 24.0) # ramps from 6 to 16 MPa
        else:
            P_f[idx] = 16.0 - 10.0 * ((tm - 24.0) / 6.0) # drops back to 6 MPa

    P_eff = P_litho - P_f # MPa
    is_breached = P_eff <= -sigma_t_sc

    ax_c.plot(time_hr, P_f, color="#d62728", lw=2.2, label=r"Fluid Pressure $P_f(t)$")
    ax_c.axhline(P_litho, color="black", ls=":", lw=1.5, label=r"Lithostatic Pressure $P_t = 5$ MPa")
    ax_c.axhline(P_crit, color="darkred", ls="--", lw=1.8, label=r"Tensile Rupture Limit $P_t + \sigma_t = 15$ MPa")

    ax_c.set_xlabel("Time [hours]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Pressure [MPa]", fontsize=11, fontweight="bold", color="#d62728")
    ax_c.tick_params(axis="y", labelcolor="#d62728")
    ax_c.set_xlim(0, 120)
    ax_c.set_ylim(0, 20)
    ax_c.grid(True, ls=":", alpha=0.6)

    ax_c2 = ax_c.twinx()
    k_cycle = np.zeros_like(P_eff)
    for idx, (peff_v, br) in enumerate(zip(P_eff, is_breached)):
        if br:
            k_cycle[idx] = compute_hydrofracture_permeability(
                k0_ref, peff_v * 1.0e6, sigma_t_Pa=sigma_t_sc * 1.0e6, kappa_frac=1.0e3, gamma=1.0, kmax=1.0e-9
            )
        else:
            k_cycle[idx] = compute_ice_sealed_permeability(
                k0_ref, np.array([150.0]), T_freeze=273.15, dT_seal=10.0, k_min_ratio=1.0e-6
            )[0]

    ax_c2.plot(time_hr, k_cycle, color="#2ca02c", lw=2.0, ls="-.", label=r"Effective Permeability $k_{\mathrm{eff}}(t)$")
    ax_c2.set_ylabel(r"Permeability $k_{\mathrm{eff}}$ [$\mathrm{m}^2$]", fontsize=11, fontweight="bold", color="#2ca02c")
    ax_c2.set_yscale("log")
    ax_c2.set_ylim(1.0e-18, 1.0e-8)
    ax_c2.tick_params(axis="y", labelcolor="#2ca02c")

    lines_c = [ax_c.lines[0], ax_c.lines[1], ax_c.lines[2], ax_c2.lines[0]]
    labels_c = [l.get_label() for l in lines_c]
    ax_c.legend(lines_c, labels_c, loc="upper right", fontsize=8.0, framealpha=0.9)
    ax_c.set_title("(c) Episodic Overpressure Rupture & Resealing", fontsize=12, fontweight="bold")

    # -------------------------------------------------------------------------
    # Panel (d): Episodic Venting Flux & Cumulative Released Mass
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    eta_f = 1.0e-3 # Pa s
    dx = 4375.0 # m
    P_vent = 10.0 # Pa
    q_vent = (k_cycle / eta_f) * np.maximum(0.0, (P_f * 1.0e6) - P_vent) / dx # m/s

    rho_water = 1000.0
    mass_flux = rho_water * q_vent # kg/(m^2 s)
    dt_sec = (time_hr[1] - time_hr[0]) * 3600.0
    cum_mass = np.cumsum(mass_flux * dt_sec) # kg/m^2

    k_unsealed = k0_ref
    q_unsealed = (k_unsealed / eta_f) * np.maximum(0.0, (P_f * 1.0e6) - P_vent) / dx
    cum_mass_unsealed = np.cumsum(rho_water * q_unsealed * dt_sec)

    ax_d.plot(time_hr, mass_flux, color="#1f77b4", lw=2.0, label="Pulsed Venting Flux (Hydrofracture Gated)")
    ax_d.set_xlabel("Time [hours]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel(r"Surface Mass Flux [$\mathrm{kg}/(\mathrm{m}^2\cdot\mathrm{s})$]",
                    fontsize=11, fontweight="bold", color="#1f77b4")
    ax_d.tick_params(axis="y", labelcolor="#1f77b4")
    ax_d.set_yscale("log")
    ax_d.set_xlim(0, 120)
    ax_d.set_ylim(1.0e-9, 1.0e-1)
    ax_d.grid(True, which="both", ls=":", alpha=0.6)

    ax_d2 = ax_d.twinx()
    ax_d2.plot(time_hr, cum_mass, color="#9467bd", lw=2.5, label="Cumulative Vented Mass (Gated)")
    ax_d2.plot(time_hr, cum_mass_unsealed, color="gray", lw=1.5, ls="--", label="Cumulative Mass (Ungated Leakage)")
    ax_d2.set_ylabel(r"Cumulative Mass [$\mathrm{kg}/\mathrm{m}^2$]",
                     fontsize=11, fontweight="bold", color="#9467bd")
    ax_d2.tick_params(axis="y", labelcolor="#9467bd")

    lines_d = [ax_d.lines[0], ax_d2.lines[0], ax_d2.lines[1]]
    labels_d = [l.get_label() for l in lines_d]
    ax_d.legend(lines_d, labels_d, loc="center left", fontsize=8.0, framealpha=0.9)
    ax_d.set_title("(d) Pulsed Cryovolcanic Venting vs Steady Leakage", fontsize=12, fontweight="bold")

    plt.savefig(out_fig, dpi=300)
    plt.close()
    print(f"Generated benchmark figure: {out_fig}")

if __name__ == "__main__":
    main()
