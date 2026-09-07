#!/usr/bin/env python3
"""
Generate diagnostic multi-panel figures for the comprehensive H-C-N-S volatile
solubility, speciation, graphite saturation, and SCSS ceilings in Erebus.jl.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def compute_iron_wustite_fO2(T_K, delta_IW=0.0):
    return 6.541 - 28164.0 / T_K + delta_IW


def compute_water_solubility_melt(P_Pa, law="burnham_dixon", As=0.40):
    P_pos = np.maximum(0.0, P_Pa)
    if law == "burnham_dixon":
        return As * np.sqrt(P_pos * 1.0e-6)
    p_bar = P_pos * 1.0e-5
    if law == "sossi_peridotite":
        return (524.0 * np.sqrt(p_bar)) * 1.0e-4
    elif law == "basalt_dixon":
        return (965.0 * np.sqrt(p_bar)) * 1.0e-4
    elif law == "newcombe_lunar":
        return (683.0 * np.sqrt(p_bar)) * 1.0e-4
    else:
        raise ValueError(f"Unknown water law: {law}")


def compute_h2_solubility_melt(p_H2_Pa, law="hirschmann2012"):
    p_pos = np.maximum(0.0, p_H2_Pa)
    p_bar = p_pos * 1.0e-5
    if law == "hirschmann2012":
        return 10.0 ** (1.10083602 + 0.52413928 * np.log10(np.maximum(p_bar, 1e-30)))
    elif law == "gaillard2003":
        return 0.163 * (p_bar ** 1.252)
    else:
        raise ValueError(f"Unknown H2 law: {law}")


def compute_co_solubility_melt(p_CO_Pa, p_total_Pa, law="armstrong2015"):
    p_co = np.maximum(0.0, p_CO_Pa)
    p_tot = np.maximum(0.0, p_total_Pa)
    p_co_bar = p_co * 1.0e-5
    p_tot_bar = p_tot * 1.0e-5
    if law == "armstrong2015":
        log_co = -0.738 + 0.876 * np.log10(np.maximum(p_co_bar, 1e-30)) - 5.44e-5 * p_tot_bar
        return 10.0 ** log_co
    elif law == "yoshioka2019_morb":
        co_wtp = 10.0 ** (-5.20 + 0.80 * np.log10(np.maximum(p_co_bar, 1e-30)))
        return co_wtp * 1.0e4 * (28.0101 / 12.011)
    else:
        raise ValueError(f"Unknown CO law: {law}")


def compute_ch4_solubility_melt(p_CH4_Pa, p_total_Pa):
    p_ch4_gpa = np.maximum(0.0, p_CH4_Pa) * 1.0e-9
    p_tot_gpa = np.maximum(0.0, p_total_Pa) * 1.0e-9
    return p_ch4_gpa * np.exp(4.93 - 1.93 * p_tot_gpa)


def compute_co2_solubility_melt(p_CO2_Pa, T_K):
    p_bar = np.maximum(0.0, p_CO2_Pa) * 1.0e-5
    x = 3.8e-7 * p_bar * np.exp(-23.0 * (p_bar - 1.0) / (83.15 * T_K))
    denom = np.maximum(36.6 - 44.0 * x, 1e-6)
    return 1.0e4 * (4400.0 * x) / denom


def compute_graphite_saturation_fugacity(T_K, log10_fO2):
    log_co = 5785.0 / T_K + 4.545 + 0.5 * log10_fO2
    log_co2 = 20590.0 / T_K - 0.043 + log10_fO2
    return 10.0 ** log_co, 10.0 ** log_co2


def compute_nitrogen_solubility_dasgupta(p_N2_Pa, p_tot_Pa, T_K, delta_IW, x_SiO2=0.56, x_Al2O3=0.11, x_TiO2=0.01):
    pN2_GPa = np.maximum(0.0, p_N2_Pa) * 1.0e-9
    ptot_GPa = np.maximum(0.0, p_tot_Pa) * 1.0e-9
    chem_exp = (5908.0 * np.sqrt(np.maximum(ptot_GPa, 1.0e-15))) / T_K - 1.6 * delta_IW
    chem_exp_clamped = np.clip(chem_exp, -100.0, 100.0)
    chem_ppm = np.sqrt(pN2_GPa) * np.exp(chem_exp_clamped)
    phys_prefactor = np.exp(4.67 + 7.11 * x_SiO2 - 13.06 * x_Al2O3 - 120.67 * x_TiO2)
    phys_ppm = pN2_GPa * phys_prefactor
    return phys_ppm + chem_ppm, phys_ppm, chem_ppm


def compute_sulfur_solubility_melt(p_S2_Pa, T_K, delta_IW, law="boulliung2023", sulfide_melt="basalt", include_sulfate=False, x_FeO=10.0):
    p_s2_bar = np.maximum(0.0, p_S2_Pa) * 1.0e-5
    log10_fO2 = compute_iron_wustite_fO2(T_K, delta_IW=delta_IW)
    fO2_bar = 10.0 ** log10_fO2

    if law == "boulliung2023":
        slope_s2 = 8045.7465 if sulfide_melt == "basalt" else 8921.0927
        logC_s2 = 0.225 - slope_s2 / T_K
        s_wtp = 10.0 ** (logC_s2 - 0.5 * (log10_fO2 - np.log10(np.maximum(p_s2_bar, 1e-30))))
        s_ppm = s_wtp * 1.0e4
        if include_sulfate:
            slope_s6 = 32333.5635 if sulfide_melt == "basalt" else 31586.2393
            logC_s6 = -12.948 + slope_s6 / T_K
            so4_wtp = 10.0 ** (logC_s6 + 0.5 * np.log10(np.maximum(p_s2_bar, 1e-30)) + 1.5 * log10_fO2)
            s_ppm += (so4_wtp * (32.065 / 96.06)) * 1.0e4
        return s_ppm
    elif law == "gaillard2022":
        ln_s = 13.8426 - 26476.0 / T_K + 0.124 * x_FeO + 0.5 * np.log(np.maximum(p_s2_bar / fO2_bar, 1e-30))
        return np.exp(np.clip(ln_s, -100.0, 100.0))
    else:
        raise ValueError(f"Unknown sulfur law: {law}")


def compute_scss(T_K, p_total_Pa, x_FeO=10.0):
    p_bar = np.maximum(0.0, p_total_Pa) * 1.0e-5
    fe_term = np.log(max(x_FeO, 0.1))
    ln_scss = 7.50 - 4500.0 / T_K + 0.90 * fe_term - 2.5e-4 * (p_bar / T_K)
    return np.exp(ln_scss)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, "hcns_solubility_benchmark.png")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
    plt.subplots_adjust(hspace=0.34, wspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Extended Water & H2 Solubility vs Pressure
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    P_bar = np.linspace(0.0, 200.0, 500)
    P_Pa = P_bar * 1.0e5

    w_bd = compute_water_solubility_melt(P_Pa, law="burnham_dixon", As=0.40)
    w_sossi = compute_water_solubility_melt(P_Pa, law="sossi_peridotite")
    w_dixon = compute_water_solubility_melt(P_Pa, law="basalt_dixon")
    w_newc = compute_water_solubility_melt(P_Pa, law="newcombe_lunar")

    ax_a.plot(P_bar, w_dixon, color="#1f77b4", lw=2.4, label="Basalt (Dixon 1995)")
    ax_a.plot(P_bar, w_newc, color="#2ca02c", lw=2.2, label="Lunar Basalt (Newcombe 2017)")
    ax_a.plot(P_bar, w_sossi, color="#ff7f0e", lw=2.2, label="Peridotite (Sossi 2023)")
    ax_a.plot(P_bar, w_bd, color="#9467bd", lw=1.8, ls="--", label="Standard Basalt (Burnham-Dixon)")

    ax_a.set_xlabel("Pore / Surface Pressure $P$ [bar]", fontsize=11, fontweight="bold")
    ax_a.set_ylabel(r"Dissolved $\mathrm{H}_2\mathrm{O}$ in Melt [wt%]", fontsize=11, fontweight="bold")
    ax_a.set_title("(a) Water Solubility Laws Across Melt Compositions", fontsize=12, fontweight="bold")
    ax_a.set_xlim(0, 200)
    ax_a.set_ylim(0, 2.0)
    ax_a.grid(True, ls=":", alpha=0.6)
    ax_a.legend(loc="upper left", fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Carbon Species Solubility & Graphite Ceiling
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    p_gas_bar = np.linspace(0.1, 100.0, 500)
    p_gas_Pa = p_gas_bar * 1.0e5
    p_tot_Pa = 100.0 * 1.0e5
    T_1500K = 1500.0

    co_ppm = compute_co_solubility_melt(p_gas_Pa, p_tot_Pa, law="armstrong2015")
    ch4_ppm = compute_ch4_solubility_melt(p_gas_Pa, p_tot_Pa)
    co2_ppm = compute_co2_solubility_melt(p_gas_Pa, T_1500K)

    # Graphite ceiling at IW-1
    f_co_gr, f_co2_gr = compute_graphite_saturation_fugacity(T_1500K, compute_iron_wustite_fO2(T_1500K, delta_IW=-1.0))
    co_gr_cap = compute_co_solubility_melt(f_co_gr * 1.0e5, p_tot_Pa, law="armstrong2015")
    co2_gr_cap = compute_co2_solubility_melt(f_co2_gr * 1.0e5, T_1500K)

    ax_b.plot(p_gas_bar, co2_ppm, color="#d62728", lw=2.4, label=r"$\mathrm{CO}_2$ Carbonate (Dixon 1995)")
    ax_b.plot(p_gas_bar, co_ppm, color="#ff7f0e", lw=2.2, label=r"$\mathrm{CO}$ Dissolved (Armstrong 2015)")
    ax_b.plot(p_gas_bar, ch4_ppm, color="#2ca02c", lw=2.0, label=r"$\mathrm{CH}_4$ Dissolved (Ardia 2013)")

    ax_b.axhline(co_gr_cap, color="#ff7f0e", ls=":", lw=1.5, label=r"Graphite Saturation $\mathrm{CO}$ cap ($\Delta$IW = -1)")
    ax_b.axhline(co2_gr_cap, color="#d62728", ls=":", lw=1.5, label=r"Graphite Saturation $\mathrm{CO}_2$ cap ($\Delta$IW = -1)")

    ax_b.set_yscale("log")
    ax_b.set_xlabel("Species Partial Pressure $p_i$ [bar]", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Dissolved Carbon Species in Melt [ppmw]", fontsize=11, fontweight="bold")
    ax_b.set_title(r"(b) Carbon Solubility Speciation & Graphite Saturation", fontsize=12, fontweight="bold")
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(1e-4, 5e3)
    ax_b.grid(True, which="both", ls=":", alpha=0.6)
    ax_b.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Nitrogen Solubility vs Redox State (Dasgupta vs Libourel)
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    d_IW_range = np.linspace(-4.0, 2.0, 500)
    p_N2_Pa = 10.0 * 1.0e5 # 10 bar N2
    p_tot_Pa = 100.0 * 1.0e5 # 100 bar total
    T_1600K = 1600.0

    tot_dasg, phys_dasg, chem_dasg = compute_nitrogen_solubility_dasgupta(
        p_N2_Pa, p_tot_Pa, T_1600K, d_IW_range, x_SiO2=0.56, x_Al2O3=0.11, x_TiO2=0.01
    )
    tot_lunar, _, chem_lunar = compute_nitrogen_solubility_dasgupta(
        p_N2_Pa, p_tot_Pa, T_1600K, d_IW_range, x_SiO2=0.48, x_Al2O3=0.16, x_TiO2=0.05
    )

    ax_c.plot(d_IW_range, tot_dasg, color="#800080", lw=2.5, label=r"Total N: Earth/Chondrite Mantle (Dasgupta 2022)")
    ax_c.plot(d_IW_range, tot_lunar, color="#e377c2", lw=2.0, ls="--", label=r"Total N: Lunar Basalt (Dasgupta 2022)")
    ax_c.plot(d_IW_range, chem_dasg, color="#d62728", lw=1.8, ls=":", label=r"Nitride $\mathrm{N}^{3-}$ Chemical Dissolution")
    ax_c.plot(d_IW_range, np.full_like(d_IW_range, phys_dasg), color="#1f77b4", lw=1.8, ls="-.", label=r"Molecular $\mathrm{N}_2$ Physical Dissolution")

    ax_c.axvline(0.0, color="gray", ls="-.", lw=1.0, alpha=0.7)
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"Redox Offset $\Delta$IW [$\log_{10}$ units]", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Nitrogen Concentration in Melt [ppmw]", fontsize=11, fontweight="bold")
    ax_c.set_title(r"(c) Compositional Nitrogen Solubility vs Redox State ($P = 10$ bar $\mathrm{N}_2$)", fontsize=12, fontweight="bold")
    ax_c.set_xlim(-4, 2)
    ax_c.set_ylim(1.0e-3, 5.0e3)
    ax_c.grid(True, which="both", ls=":", alpha=0.6)
    ax_c.legend(loc="upper right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (d): Sulfur Solubility & SCSS Saturation Ceiling
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    p_S2_Pa = 1.0 * 1.0e5 # 1 bar S2
    T_1500K = 1500.0
    p_tot_100bar = 100.0 * 1.0e5

    s_boul_pure = np.array([
        compute_sulfur_solubility_melt(p_S2_Pa, T_1500K, d, law="boulliung2023", sulfide_melt="basalt", include_sulfate=False)
        for d in d_IW_range
    ])
    s_boul_so4 = np.array([
        compute_sulfur_solubility_melt(p_S2_Pa, T_1500K, d, law="boulliung2023", sulfide_melt="basalt", include_sulfate=True)
        for d in d_IW_range
    ])
    s_gail = np.array([
        compute_sulfur_solubility_melt(p_S2_Pa, T_1500K, d, law="gaillard2022", x_FeO=10.0)
        for d in d_IW_range
    ])

    scss_10feo = compute_scss(T_1500K, p_tot_100bar, x_FeO=10.0)
    scss_20feo = compute_scss(T_1500K, p_tot_100bar, x_FeO=20.0)

    ax_d.plot(d_IW_range, s_boul_pure, color="#ff7f0e", lw=2.2, label=r"Boulliung & Wood (2023) Sulfide $\mathrm{S}^{2-}$")
    ax_d.plot(d_IW_range, s_boul_so4, color="#d62728", lw=2.4, ls="--", label=r"Boulliung & Wood (2023) Sulfide + Sulfate $\mathrm{SO}_4^{2-}$")
    ax_d.plot(d_IW_range, s_gail, color="#2ca02c", lw=2.0, ls="-.", label=r"Gaillard et al. (2022) Basalt")

    ax_d.axhline(scss_10feo, color="#7f7f7f", ls=":", lw=2.2, label=f"SCSS Ceiling (10 wt% FeO): {scss_10feo:.0f} ppmw")
    ax_d.axhline(scss_20feo, color="#bcbd22", ls=":", lw=2.0, label=f"SCSS Ceiling (20 wt% FeO): {scss_20feo:.0f} ppmw")

    # Shaded region where raw solubility exceeds SCSS (immiscible Fe-S matte precipitation)
    ax_d.axhspan(scss_10feo, 1e7, color="gold", alpha=0.12, label="Sulfide Liquid (Matte) Precipitation Zone")

    ax_d.axvline(0.0, color="gray", ls="-.", lw=1.0, alpha=0.7)
    ax_d.set_yscale("log")
    ax_d.set_xlabel(r"Redox Offset $\Delta$IW [$\log_{10}$ units]", fontsize=11, fontweight="bold")
    ax_d.set_ylabel("Sulfur Concentration in Melt [ppmw]", fontsize=11, fontweight="bold")
    ax_d.set_title(r"(d) Sulfur Solubility, Sulfate Transition & SCSS Ceiling ($T = 1500$ K)", fontsize=12, fontweight="bold")
    ax_d.set_xlim(-4, 2)
    ax_d.set_ylim(1.0, 1.0e6)
    ax_d.grid(True, which="both", ls=":", alpha=0.6)
    ax_d.legend(loc="upper right", fontsize=8.0, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=300)
    plt.close()
    print(f"Generated benchmark figure: {out_fig}")


if __name__ == "__main__":
    main()
