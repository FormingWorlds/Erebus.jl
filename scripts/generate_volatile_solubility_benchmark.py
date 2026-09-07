#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for multi-species volatile solubility,
redox-dependent nitrogen chemistry, and organic devolatilization in Erebus.jl.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compute_iron_wustite_fO2(T_K, delta_IW=0.0):
    """
    Compute log10(fO2 [bar]) for iron-wüstite buffer.

    Parameters
    ----------
    T_K : float or numpy.ndarray
        Temperature in Kelvin.
    delta_IW : float, optional
        Offset relative to IW buffer in log10 units.

    Returns
    -------
    float or numpy.ndarray
        log10 of oxygen fugacity in bar.
    """
    return 6.541 - 28164.0 / T_K + delta_IW


def compute_water_solubility_melt(P_Pa, As=0.40):
    """
    Compute equilibrium water solubility in silicate melt.

    Parameters
    ----------
    P_Pa : float or numpy.ndarray
        Pore fluid pressure in Pascals.
    As : float, optional
        Burnham water solubility coefficient [wt% / MPa^0.5].

    Returns
    -------
    float or numpy.ndarray
        Dissolved water concentration in wt%.
    """
    P_pos = np.maximum(0.0, P_Pa)
    P_MPa = P_pos * 1.0e-6
    return As * np.sqrt(P_MPa)


def compute_nitrogen_solubility_melt(P_Pa, delta_IW, Kh=0.40, C_nitride=1.0e-3):
    """
    Compute nitrogen solubility partitioning in silicate melt.

    Parameters
    ----------
    P_Pa : float or numpy.ndarray
        Pore fluid pressure in Pascals.
    delta_IW : float or numpy.ndarray
        Oxygen fugacity offset relative to IW buffer.
    Kh : float, optional
        Henry coefficient for molecular N2 [ppm / bar].
    C_nitride : float, optional
        Chemical nitride capacity [wt% / bar^0.5].

    Returns
    -------
    tuple
        (total_ppm, physical_ppm, chemical_ppm)
    """
    f_N2 = np.maximum(0.0, P_Pa) * 1.0e-5
    physical_ppm = Kh * f_N2

    # Chemical nitride scaling with oxygen fugacity
    fO2_ratio = 10.0 ** delta_IW
    chemical_ppm = (C_nitride * 1.0e4) * np.sqrt(f_N2) * (fO2_ratio ** (-0.75))
    total_ppm = physical_ppm + chemical_ppm
    return total_ppm, physical_ppm, chemical_ppm


def compute_organic_nitrogen_yield(T_K, T_devol=550.0, delta_T=50.0):
    """
    Compute devolatilization yield of organic nitrogen.

    Parameters
    ----------
    T_K : float or numpy.ndarray
        Temperature in Kelvin.
    T_devol : float, optional
        Midpoint devolatilization temperature in Kelvin.
    delta_T : float, optional
        Transition width scale in Kelvin.

    Returns
    -------
    float or numpy.ndarray
        Devolatilized nitrogen fraction in [0, 1].
    """
    arg = (T_K - T_devol) / delta_T
    arg = np.clip(arg, -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-arg))


def main():
    """Generate multi-panel volatile solubility benchmark figure."""
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "volatile_solubility_benchmark.png")

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=300)
    plt.subplots_adjust(hspace=0.32, wspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Water Solubility vs Pore Pressure
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    P_MPa = np.linspace(0.0, 100.0, 500)
    P_Pa = P_MPa * 1.0e6

    As_values = [0.30, 0.40, 0.50]
    colors_a = ["#1f77b4", "#2ca02c", "#ff7f0e"]
    for As_val, col in zip(As_values, colors_a):
        w_H2O = compute_water_solubility_melt(P_Pa, As=As_val)
        ax_a.plot(
            P_MPa,
            w_H2O,
            color=col,
            lw=2.2,
            label=f"$A_s = {As_val:.2f}$ wt% / MPa$^{{0.5}}$",
        )

    ax_a.set_xlabel("Pore Fluid Pressure $P_f$ [MPa]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Dissolved $\\mathrm{H}_2\\mathrm{O}$ in Melt [wt%]", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Water Solubility vs Pressure (Square-Root Law)", fontsize=12, fontweight="bold")
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(0, 5.5)
    ax_a.grid(True, ls=":", alpha=0.6)
    ax_a.legend(loc="lower right", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Iron-Wüstite Oxygen Fugacity vs Temperature
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    T_range = np.linspace(800.0, 1800.0, 500)
    delta_IW_values = [-3.0, -2.0, -1.0, 0.0, +1.0, +2.0]
    cmap_b = plt.get_cmap("coolwarm", len(delta_IW_values))

    for idx, d_iw in enumerate(delta_IW_values):
        log_fO2 = compute_iron_wustite_fO2(T_range, delta_IW=d_iw)
        label_str = f"$\\Delta$IW = {d_iw:+.0f}" if d_iw != 0 else "IW Buffer ($\\Delta$IW = 0)"
        ax_b.plot(T_range, log_fO2, color=cmap_b(idx), lw=2.0, label=label_str)

    ax_b.set_xlabel("Melt Temperature $T$ [K]", fontsize=11, fontweight="bold")
    ax_b.set_ylabel(r"Oxygen Fugacity $\log_{10}(f_{\mathrm{O}_2}\ [\mathrm{bar}])$", fontsize=11, fontweight="bold")
    ax_b.set_title("(b) Iron-Wüstite Oxygen Fugacity vs Temperature", fontsize=12, fontweight="bold")
    ax_b.set_xlim(800, 1800)
    ax_b.grid(True, ls=":", alpha=0.6)
    ax_b.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Nitrogen Solubility Partitioning vs Redox State
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    d_IW_range = np.linspace(-4.0, 4.0, 500)
    P_fixed_Pa = 10.0e6  # 10 MPa = 100 bar
    T_fixed_K = 1500.0

    total_ppm, phys_ppm, chem_ppm = compute_nitrogen_solubility_melt(
        P_fixed_Pa, d_IW_range, Kh=0.40, C_nitride=1.0e-3
    )

    ax_c.plot(d_IW_range, total_ppm, color="#800080", lw=2.5, label="Total Nitrogen")
    ax_c.plot(d_IW_range, chem_ppm, color="#d62728", lw=2.0, ls="--", label="Chemical Nitride ($\\mathrm{N}^{3-}$)")
    ax_c.plot(d_IW_range, np.full_like(d_IW_range, phys_ppm), color="#1f77b4", lw=2.0, ls=":", label="Physical ($\\mathrm{N}_2$)")

    ax_c.axvline(0.0, color="gray", ls="-.", lw=1.0, alpha=0.7)
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"Redox State $\Delta$IW [log$_{10}$ units]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Nitrogen Concentration in Melt [ppm]", fontsize=11, fontweight="bold")
    ax_c.set_title("(c) Nitrogen Solubility vs Redox State ($P = 10$ MPa, $T = 1500$ K)", fontsize=12, fontweight="bold")
    ax_c.set_xlim(-4, 4)
    ax_c.set_ylim(1.0e-2, 5.0e4)
    ax_c.grid(True, which="both", ls=":", alpha=0.6)
    ax_c.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (d): Organic Nitrogen Devolatilization Yield
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    T_organic = np.linspace(300.0, 800.0, 500)
    T_devol_list = [500.0, 550.0, 600.0]
    colors_d = ["#2ca02c", "#1f77b4", "#d62728"]

    for T_mid, col in zip(T_devol_list, colors_d):
        yield_vals = compute_organic_nitrogen_yield(T_organic, T_devol=T_mid, delta_T=50.0)
        ax_d.plot(
            T_organic,
            yield_vals,
            color=col,
            lw=2.2,
            label=f"$T_{{\\mathrm{{devol}}}} = {T_mid:.0f}$ K",
        )

    ax_d.axhline(0.5, color="gray", ls="--", lw=1.0, alpha=0.7, label="Midpoint Yield ($0.5$)")
    ax_d.set_xlabel("Rock Temperature $T$ [K]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel("Devolatilization Yield $y(T)$ [-]", fontsize=11, fontweight="bold")
    ax_d.set_title("(d) Organic Nitrogen Devolatilization Kinetics", fontsize=12, fontweight="bold")
    ax_d.set_xlim(300, 800)
    ax_d.set_ylim(-0.02, 1.02)
    ax_d.grid(True, ls=":", alpha=0.6)
    ax_d.legend(loc="lower right", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=300)
    plt.close()
    print(f"Generated benchmark figure: {out_fig}")


if __name__ == "__main__":
    main()
