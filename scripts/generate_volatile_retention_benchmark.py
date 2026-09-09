#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for volatile retention floors in
nominally anhydrous minerals (NAMs) and low-temperature venting drainage
coupling in Erebus.jl.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compute_retention_floor(T_K, C_floor, T_sol=1400.0, dT=200.0, law="nams_exponential"):
    """
    Compute temperature-dependent volatile retention floor.

    Parameters
    ----------
    T_K : float or numpy.ndarray
        Temperature in Kelvin.
    C_floor : float
        Nominal subsolidus retention floor in ppmw.
    T_sol : float, optional
        Reference solidus temperature in Kelvin.
    dT : float, optional
        Melt transition width in Kelvin.
    law : str, optional
        Retention law: "constant_floor", "linear_melt_blend", or "nams_exponential".

    Returns
    -------
    float or numpy.ndarray
        Effective retention floor in ppmw.
    """
    T = np.asarray(T_K)
    if law == "constant_floor":
        return np.full_like(T, C_floor, dtype=float)
    elif law == "linear_melt_blend":
        phi_melt = np.clip((T - T_sol) / max(dT, 1.0), 0.0, 1.0)
        return C_floor * (1.0 - phi_melt)
    elif law == "nams_exponential":
        dT_sup = np.maximum(0.0, T - T_sol)
        return C_floor * np.exp(-dT_sup / max(dT, 1.0))
    else:
        raise ValueError(f"Unknown retention law: {law}")


def compute_decompression_solubility(P_Pa, C_bulk_ppm=500.0, C_floor_ppm=50.0, with_retention=True):
    """
    Compute dissolved volatile concentration during isothermal decompression.

    Parameters
    ----------
    P_Pa : numpy.ndarray
        Pressure range in Pascals.
    C_bulk_ppm : float, optional
        Bulk initial volatile concentration in ppmw.
    C_floor_ppm : float, optional
        Retention floor in ppmw.
    with_retention : bool, optional
        Whether retention floor clamping is active.

    Returns
    -------
    numpy.ndarray
        Dissolved volatile concentration in ppmw.
    """
    # Water Henry/Burnham style solubility: C_eq [wt%] = 0.40 * sqrt(P_MPa)
    # Convert to ppmw: C_eq [ppmw] = 0.40 * sqrt(P_MPa) * 10000
    P_MPa = np.maximum(0.0, P_Pa) * 1.0e-6
    S_eq_ppm = 0.40 * np.sqrt(P_MPa) * 10000.0

    if not with_retention:
        return np.minimum(C_bulk_ppm, S_eq_ppm)

    C_mob = max(0.0, C_bulk_ppm - C_floor_ppm)
    C_sol = C_floor_ppm + np.minimum(C_mob, S_eq_ppm)
    return C_sol


def simulate_venting_drainage(t_years, S_vent, chi_vent=1.0, C_bulk0=500.0, C_floor=50.0):
    """
    Simulate time evolution of marker volatile depletion under surface venting.

    Parameters
    ----------
    t_years : numpy.ndarray
        Time in years.
    S_vent : float
        Surface venting volumetric sink rate [1/s].
    chi_vent : float, optional
        Volatile extraction efficiency factor.
    C_bulk0 : float, optional
        Initial bulk volatile concentration in ppmw.
    C_floor : float, optional
        Solid retention floor in ppmw.

    Returns
    -------
    tuple of numpy.ndarray
        (C_bulk, C_mob, C_vented_fraction)
    """
    t_sec = t_years * 365.25 * 86400.0
    C_mob0 = max(0.0, C_bulk0 - C_floor)
    decay = np.exp(-S_vent * chi_vent * t_sec)
    C_mob = C_mob0 * decay
    C_bulk = C_floor + C_mob
    drained_cum = C_mob0 * (1.0 - decay)
    vented_fraction = drained_cum / C_bulk0
    return C_bulk, C_mob, vented_fraction


def main():
    """Generate the four-panel volatile retention and venting benchmark figure."""
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "volatile_retention_benchmark.png")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(
        "Volatile Retention Floors & Low-Temperature Venting Drainage Coupling",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    # -------------------------------------------------------------------------
    # Panel (a): Retention Floor vs Temperature
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    T_arr = np.linspace(1000.0, 1800.0, 400)
    T_sol = 1400.0
    dT = 200.0

    laws = [
        ("constant_floor", "Constant Floor", "#1f77b4", "-"),
        ("linear_melt_blend", "Linear Melt Blend", "#ff7f0e", "--"),
        ("nams_exponential", "NAMs Exponential", "#2ca02c", "-"),
    ]

    for law_key, label, col, ls in laws:
        floor_h2o = compute_retention_floor(T_arr, 50.0, T_sol, dT, law=law_key)
        ax_a.plot(T_arr, floor_h2o, color=col, ls=ls, lw=2.2, label=f"H2O ({label})")

    floor_c_nams = compute_retention_floor(T_arr, 25.0, T_sol, dT, law="nams_exponential")
    ax_a.plot(T_arr, floor_c_nams, color="#9467bd", ls=":", lw=2.2, label="C (NAMs Exponential)")

    ax_a.axvline(T_sol, color="gray", ls="-.", lw=1.0, alpha=0.7, label=r"Solidus $T_{\mathrm{sol}} = 1400$ K")
    ax_a.set_xlabel("Temperature $T$ [K]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Retained Volatile Floor [ppmw]", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Temperature-Dependent Solid Retention Floor", fontsize=12, fontweight="bold")
    ax_a.set_xlim(1000, 1800)
    ax_a.set_ylim(-2, 55)
    ax_a.grid(True, ls=":", alpha=0.6)
    ax_a.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Decompression Exsolution with Vacuum Floor
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    P_arr = np.logspace(-5, 8, 500)  # 10^-5 Pa to 100 MPa
    C_bulk = 500.0
    C_floor = 50.0

    C_sol_with = compute_decompression_solubility(P_arr, C_bulk, C_floor, with_retention=True)
    C_sol_without = compute_decompression_solubility(P_arr, C_bulk, C_floor, with_retention=False)

    ax_b.plot(P_arr, C_sol_with, color="#2ca02c", lw=2.4, label="With Retention Floor (50 ppmw)")
    ax_b.plot(P_arr, C_sol_without, color="#d62728", ls="--", lw=2.0, label="Without Retention Floor (Classical)")
    ax_b.axhline(C_floor, color="#1f77b4", ls=":", lw=1.5, alpha=0.8, label=r"Vacuum Floor $C_{\mathrm{ret}}$ (50 ppmw)")
    ax_b.axhline(C_bulk, color="gray", ls="-.", lw=1.0, alpha=0.7, label=r"Bulk Inventory $C_{\mathrm{bulk}}$ (500 ppmw)")

    ax_b.set_xscale("log")
    ax_b.set_xlabel("Pore Fluid Pressure $P$ [Pa]", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Dissolved Volatile Concentration [ppmw]", fontsize=11, fontweight="bold")
    ax_b.set_title("(b) Decompression Exsolution & Vacuum Retention Limit", fontsize=12, fontweight="bold")
    ax_b.set_xlim(1.0e-5, 1.0e8)
    ax_b.set_ylim(0, 550)
    ax_b.grid(True, which="both", ls=":", alpha=0.6)
    ax_b.legend(loc="lower right", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Venting Drainage Kinetics & Floor Protection
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    t_kyr = np.linspace(0.0, 100.0, 500)  # 0 to 100 kyr
    s_vent_rates = [
        (1.0e-12, r"$S_{\mathrm{vent}} = 10^{-12}\ \mathrm{s}^{-1}$", "#1f77b4"),
        (3.0e-12, r"$S_{\mathrm{vent}} = 3\times 10^{-12}\ \mathrm{s}^{-1}$", "#ff7f0e"),
        (1.0e-11, r"$S_{\mathrm{vent}} = 10^{-11}\ \mathrm{s}^{-1}$", "#d62728"),
    ]

    for s_rate, s_label, col in s_vent_rates:
        c_b, c_m, _ = simulate_venting_drainage(t_kyr * 1000.0, s_rate, chi_vent=1.0, C_bulk0=500.0, C_floor=50.0)
        ax_c.plot(t_kyr, c_b, color=col, lw=2.2, label=f"Bulk: {s_label}")
        ax_c.plot(t_kyr, c_m, color=col, ls="--", lw=1.5, alpha=0.7)

    ax_c.axhline(50.0, color="black", ls=":", lw=1.8, label=r"Retention Floor $C_{\mathrm{ret}}$ (50 ppmw)")
    ax_c.set_xlabel("Venting Duration [kyr]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Marker Volatile Concentration [ppmw]", fontsize=11, fontweight="bold")
    ax_c.set_title(r"(c) Low-Temperature Vent Drainage Kinetics ($C_{\mathrm{bulk}}$ solid, $C_{\mathrm{mob}}$ dashed)", fontsize=11, fontweight="bold")
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0, 520)
    ax_c.grid(True, ls=":", alpha=0.6)
    ax_c.legend(loc="upper right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (d): Multi-Species Drained vs Retained Inventories
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    species_data = [
        ("H2O", 500.0, 50.0, "#1f77b4"),
        ("CO2 / C", 200.0, 25.0, "#2ca02c"),
        ("N2 / N", 50.0, 5.0, "#9467bd"),
        ("H2S / S", 300.0, 100.0, "#d62728"),
    ]

    s_ref = 3.0e-12
    t_eval = t_kyr * 1000.0

    for name, c0, c_fl, col in species_data:
        _, _, v_frac = simulate_venting_drainage(t_eval, s_ref, chi_vent=1.0, C_bulk0=c0, C_floor=c_fl)
        ax_d.plot(t_kyr, v_frac * 100.0, color=col, lw=2.2, label=rf"{name} ($C_0={c0:.0f}$, $C_{{\mathrm{{ret}}}}={c_fl:.0f}$)")

    ax_d.set_xlabel("Venting Duration [kyr]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel("Cumulative Vented Inventory [% of Initial]", fontsize=11, fontweight="bold")
    ax_d.set_title(r"(d) Multi-Species Vented Volatile Fraction ($S_{\mathrm{vent}} = 3\times 10^{-12}\ \mathrm{s}^{-1}$)", fontsize=11, fontweight="bold")
    ax_d.set_xlim(0, 100)
    ax_d.set_ylim(0, 100)
    ax_d.grid(True, ls=":", alpha=0.6)
    ax_d.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=300)
    plt.close()
    print(f"Generated benchmark figure: {out_fig}")


if __name__ == "__main__":
    main()
