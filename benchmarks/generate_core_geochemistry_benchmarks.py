#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for metal-silicate volatile partitioning
and core geochemistry in Erebus.jl.
"""

import os
import sys
import json
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

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

# Thermodynamic partition coefficient implementations matching src/physics.jl
def mole_fraction_S(w_S):
    w = np.clip(w_S, 0.0, 0.999)
    n_S = w / 32.065
    n_Fe = (1.0 - w) / 55.845
    return n_S / (n_S + n_Fe)

def D_carbon_grewal(T, P, dIW, w_S):
    X_S = mole_fraction_S(w_S)
    ln_1_minus_XS = np.log(np.maximum(1.0 - X_S, 1e-6))
    log10_D = 1.80 + 2200.0 / T - 1.5e-8 * (P / T) - 0.25 * dIW + 4.2 * ln_1_minus_XS
    return np.clip(10.0**log10_D, 1e-4, 1e5)

def D_nitrogen_grewal(T, P, dIW, w_S):
    X_S = mole_fraction_S(w_S)
    ln_1_minus_XS = np.log(np.maximum(1.0 - X_S, 1e-6))
    log10_D = 0.85 + 1200.0 / T - 0.25 * dIW + 0.60 * ln_1_minus_XS
    return np.clip(10.0**log10_D, 1e-4, 1e5)

def D_sulfur_boujibar(T, P, dIW):
    log10_D = 2.80 - 800.0 / T + 1.0e-10 * P - 0.20 * dIW
    return np.clip(10.0**log10_D, 1e-4, 1e5)

def D_hydrogen_clesi(T, P, dIW):
    log10_D = -0.80 + 300.0 / T + 5.0e-8 * (P / T) + 0.05 * dIW
    return np.clip(10.0**log10_D, 1e-4, 1e5)

def generate_benchmark_figure():
    fig, axes = plt.subplots(2, 3, figsize=(16, 10.5))
    plt.subplots_adjust(left=0.07, right=0.96, bottom=0.08, top=0.93, wspace=0.28, hspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Oxygen Fugacity Sensitivity D_i(dIW)
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    ax_a.text(0.04, 0.93, '(a)', transform=ax_a.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    dIW_arr = np.linspace(-4.0, 0.0, 100)
    T_ref = 1600.0
    P_ref = 1.0e8  # 0.1 GPa
    w_S_ref = 0.05  # 5 wt% S

    ax_a.plot(dIW_arr, D_carbon_grewal(T_ref, P_ref, dIW_arr, w_S_ref), '-', color=STRATA['magma'], linewidth=2.2, label=r'Carbon $D_\mathrm{C}$')
    ax_a.plot(dIW_arr, D_nitrogen_grewal(T_ref, P_ref, dIW_arr, w_S_ref), '-', color=STRATA['cobalt'], linewidth=2.2, label=r'Nitrogen $D_\mathrm{N}$')
    ax_a.plot(dIW_arr, D_sulfur_boujibar(T_ref, P_ref, dIW_arr), '-', color=STRATA['gold'], linewidth=2.2, label=r'Sulfur $D_\mathrm{S}$')
    ax_a.plot(dIW_arr, D_hydrogen_clesi(T_ref, P_ref, dIW_arr), '-', color=STRATA['plum'], linewidth=2.2, label=r'Hydrogen $D_\mathrm{H}$')

    ax_a.set_yscale('log')
    ax_a.set_xlim(-4.0, 0.0)
    ax_a.set_ylim(0.1, 10000.0)
    ax_a.set_xlabel(r'Oxygen Fugacity $\Delta\mathrm{IW}$ [log units]', fontsize=10)
    ax_a.set_ylabel(r'Partition Coefficient $D_i^\mathrm{met/sil}$ [-]', fontsize=10)
    ax_a.set_title(r'Redox Sensitivity ($T = 1600$ K, $P = 0.1$ GPa, $w_\mathrm{S} = 0.05$)', fontsize=11, fontweight='bold')
    ax_a.grid(True)
    ax_a.legend(loc='upper right', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (b): Sulfur Saturation & Carbon Suppression D_i(w_S)
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    ax_b.text(0.04, 0.93, '(b)', transform=ax_b.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    w_S_arr = np.linspace(0.0, 0.31, 100)
    dIW_val = -2.0

    ax_b.plot(w_S_arr * 100, D_carbon_grewal(T_ref, P_ref, dIW_val, w_S_arr), '-', color=STRATA['magma'], linewidth=2.2, label=r'Carbon $D_\mathrm{C}$ (Grewal+ 2019b)')
    ax_b.plot(w_S_arr * 100, D_nitrogen_grewal(T_ref, P_ref, dIW_val, w_S_arr), '-', color=STRATA['cobalt'], linewidth=2.2, label=r'Nitrogen $D_\mathrm{N}$ (Grewal+ 2019b)')
    
    # C/N partition ratio on twin axis
    ax_b_twin = ax_b.twinx()
    ratio_CN = D_carbon_grewal(T_ref, P_ref, dIW_val, w_S_arr) / D_nitrogen_grewal(T_ref, P_ref, dIW_val, w_S_arr)
    ax_b_twin.plot(w_S_arr * 100, ratio_CN, ':', color=STRATA['amber'], linewidth=2.0, label=r'$D_\mathrm{C} / D_\mathrm{N}$ Ratio')
    ax_b_twin.set_ylabel(r'Partition Ratio $D_\mathrm{C} / D_\mathrm{N}$ [-]', fontsize=9.5, color=STRATA['amber'])
    ax_b_twin.tick_params(axis='y', labelcolor=STRATA['amber'])
    ax_b_twin.set_yscale('log')
    ax_b_twin.set_ylim(0.5, 200.0)

    ax_b.axvline(31.0, color=NEUTRALS['graphite'], linestyle='--', linewidth=1.0, alpha=0.7, label='Fe-FeS Eutectic (~31 wt%)')
    ax_b.set_yscale('log')
    ax_b.set_xlim(0, 32)
    ax_b.set_ylim(5, 5000)
    ax_b.set_xlabel(r'Alloy Sulfur Content $w_\mathrm{S}$ [wt%]', fontsize=10)
    ax_b.set_ylabel(r'Partition Coefficient $D_i^\mathrm{met/sil}$ [-]', fontsize=10)
    ax_b.set_title(r'Sulfur Suppression of Carbon Affinity ($\Delta\mathrm{IW} = -2$)', fontsize=11, fontweight='bold')
    ax_b.grid(True)
    ax_b.legend(loc='lower left', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (c): Thermal Sensitivity D_i(T)
    # -------------------------------------------------------------------------
    ax_c = axes[0, 2]
    ax_c.text(0.04, 0.93, '(c)', transform=ax_c.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))
    T_arr = np.linspace(1300.0, 2200.0, 100)
    w_S_med = 0.10

    ax_c.plot(T_arr, D_carbon_grewal(T_arr, P_ref, dIW_val, w_S_med), '-', color=STRATA['magma'], linewidth=2.2, label=r'Carbon $D_\mathrm{C}$')
    ax_c.plot(T_arr, D_nitrogen_grewal(T_arr, P_ref, dIW_val, w_S_med), '-', color=STRATA['cobalt'], linewidth=2.2, label=r'Nitrogen $D_\mathrm{N}$')
    ax_c.plot(T_arr, D_sulfur_boujibar(T_arr, P_ref, dIW_val), '-', color=STRATA['gold'], linewidth=2.2, label=r'Sulfur $D_\mathrm{S}$ (Boujibar+ 2014)')
    ax_c.plot(T_arr, D_hydrogen_clesi(T_arr, P_ref, dIW_val), '-', color=STRATA['plum'], linewidth=2.2, label=r'Hydrogen $D_\mathrm{H}$ (Clesi+ 2018)')

    ax_c.set_yscale('log')
    ax_c.set_xlim(1300, 2200)
    ax_c.set_ylim(0.1, 2000.0)
    ax_c.set_xlabel('Temperature $T$ [K]', fontsize=10)
    ax_c.set_ylabel(r'Partition Coefficient $D_i^\mathrm{met/sil}$ [-]', fontsize=10)
    ax_c.set_title(r'Thermal Sensitivity ($P = 0.1$ GPa, $\Delta\mathrm{IW} = -2$, $w_\mathrm{S} = 0.10$)', fontsize=11, fontweight='bold')
    ax_c.grid(True)
    ax_c.legend(loc='center right', fontsize=8.5)

    # -------------------------------------------------------------------------
    # Panel (d): Magmatic Iron Meteorites vs. Erebus Model Tracks
    # -------------------------------------------------------------------------
    ax_d = axes[1, 0]
    ax_d.text(0.04, 0.93, '(d)', transform=ax_d.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

    # Literature fields for magmatic iron meteorites (C ppmw vs N ppmw)
    groups = [
        {'name': 'IIAB (S-rich)', 'C_range': (100, 350), 'N_range': (15, 35), 'color': STRATA['gold'], 'alpha': 0.35},
        {'name': 'IIIAB', 'C_range': (250, 550), 'N_range': (10, 26), 'color': STRATA['amber'], 'alpha': 0.35},
        {'name': 'IVA', 'C_range': (30, 150), 'N_range': (2.0, 8.5), 'color': STRATA['cobalt'], 'alpha': 0.35},
        {'name': 'IVB (refractory)', 'C_range': (10, 60), 'N_range': (0.5, 2.8), 'color': STRATA['plum'], 'alpha': 0.35},
    ]

    for g in groups:
        w_C = g['C_range'][1] - g['C_range'][0]
        h_N = g['N_range'][1] - g['N_range'][0]
        rect = Rectangle((g['C_range'][0], g['N_range'][0]), w_C, h_N,
                         facecolor=g['color'], edgecolor=g['color'], alpha=g['alpha'], linewidth=1.5, label=g['name'])
        ax_d.add_patch(rect)
        ax_d.text(g['C_range'][0] + 0.08 * w_C, g['N_range'][0] + 0.3 * h_N, g['name'].split()[0],
                  fontsize=8.5, fontweight='bold', color=NEUTRALS['graphite'])

    # Model equilibrium tracks
    dIW_sweep = np.linspace(-3.5, -1.0, 20)
    track_S_rich_C = []
    track_S_rich_N = []
    track_S_mod_C = []
    track_S_mod_N = []

    m_met_frac = 0.15 * 7000.0
    m_sil_frac = 0.85 * 3000.0

    for diw in dIW_sweep:
        # High S case (IIAB precursor, w_S = 0.15)
        DC_hi = D_carbon_grewal(1600.0, 1.0e8, diw, 0.15)
        DN_hi = D_nitrogen_grewal(1600.0, 1.0e8, diw, 0.15)
        C_met_hi = (m_sil_frac * 150.0 * DC_hi) / (m_sil_frac + DC_hi * m_met_frac)
        N_met_hi = (m_sil_frac * 20.0 * DN_hi) / (m_sil_frac + DN_hi * m_met_frac)
        track_S_rich_C.append(C_met_hi)
        track_S_rich_N.append(N_met_hi)

        # Moderate S case (IIIAB precursor, w_S = 0.07)
        DC_mod = D_carbon_grewal(1600.0, 1.0e8, diw, 0.07)
        DN_mod = D_nitrogen_grewal(1600.0, 1.0e8, diw, 0.07)
        C_met_mod = (m_sil_frac * 150.0 * DC_mod) / (m_sil_frac + DC_mod * m_met_frac)
        N_met_mod = (m_sil_frac * 20.0 * DN_mod) / (m_sil_frac + DN_mod * m_met_frac)
        track_S_mod_C.append(C_met_mod)
        track_S_mod_N.append(N_met_mod)

    ax_d.plot(track_S_rich_C, track_S_rich_N, 'o-', color=STRATA['magma'], markersize=4, linewidth=2.0, label=r'Erebus ($w_\mathrm{S} = 0.15$, $\Delta\mathrm{IW} \in [-3.5, -1]$)')
    ax_d.plot(track_S_mod_C, track_S_mod_N, 's--', color=STRATA['cobalt'], markersize=4, linewidth=2.0, label=r'Erebus ($w_\mathrm{S} = 0.07$, $\Delta\mathrm{IW} \in [-3.5, -1]$)')

    ax_d.set_xscale('log')
    ax_d.set_yscale('log')
    ax_d.set_xlim(8, 800)
    ax_d.set_ylim(0.4, 50)
    ax_d.set_xlabel(r'Core Carbon Concentration $w_\mathrm{C}$ [ppmw]', fontsize=10)
    ax_d.set_ylabel(r'Core Nitrogen Concentration $w_\mathrm{N}$ [ppmw]', fontsize=10)
    ax_d.set_title('Magmatic Iron Meteorites vs. Core Model', fontsize=11, fontweight='bold')
    ax_d.grid(True)
    ax_d.legend(loc='lower right', fontsize=7.5)

    # -------------------------------------------------------------------------
    # Panel (e): Core Volatile Segregation Dynamics Over Time
    # -------------------------------------------------------------------------
    ax_e = axes[1, 1]
    ax_e.text(0.04, 0.93, '(e)', transform=ax_e.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

    # Read reference simulation trajectory from core formation benchmark dataset
    benchmark_data_path = os.path.join(OUTPUT_FILES_DIR, "core_formation_benchmark_data.json")
    with open(benchmark_data_path, "r") as f:
        bench_data = json.load(f)

    ref_run = bench_data["reference"]
    time_ma = np.array(ref_run["times_yr"]) / 1e6
    R_core_km = np.array(ref_run["R_core_hist"])
    T_core_hist = np.array(ref_run["T_core_hist"])
    phi_fe_core = np.array(ref_run["phi_fe_core_hist"])

    # Planetary mass and core geometry (R_planet = 50 km)
    R_planet = 50000.0  # m
    rho_sil = 3200.0    # kg/m^3
    rho_fe_0 = 7000.0   # kg/m^3
    V_planet = (4.0 / 3.0) * np.pi * R_planet**3
    M_planet = V_planet * rho_sil

    R_core_m = R_core_km * 1e3
    V_core = (4.0 / 3.0) * np.pi * R_core_m**3
    M_core_fe = V_core * phi_fe_core * rho_fe_0

    # Bulk planetary volatile concentrations (chondritic reference)
    w_S_bulk = 0.025     # 2.5 wt% S
    C_bulk = 3500.0      # 3500 ppmw C
    N_bulk = 90.0        # 90 ppmw N
    H_bulk = 400.0       # 400 ppmw H (~0.36 wt% H2O)

    P_cmb = 5.0e7        # 0.05 GPa at core-mantle boundary
    dIW_val = -2.0

    M_core_S = np.zeros_like(time_ma)
    M_core_C = np.zeros_like(time_ma)
    M_core_N = np.zeros_like(time_ma)
    M_core_H = np.zeros_like(time_ma)

    for i in range(len(time_ma)):
        m_fe = M_core_fe[i]
        if m_fe <= 1e12:
            continue
        m_sil = max(M_planet - m_fe, 1e12)
        T_local = max(T_core_hist[i], 1200.0)

        # Sulfur partitioning (Boujibar et al. 2014)
        DS = D_sulfur_boujibar(T_local, P_cmb, dIW_val)
        w_S_core = np.clip((w_S_bulk * DS) / (m_sil / m_fe + DS), 0.0, 0.365)

        # Carbon partitioning with sulfur suppression (Grewal et al. 2019b)
        DC = D_carbon_grewal(T_local, P_cmb, dIW_val, w_S_core)
        C_met = np.clip((C_bulk * DC) / (m_sil / m_fe + DC), 0.0, 70000.0)

        # Nitrogen partitioning (Grewal et al. 2019b)
        DN = D_nitrogen_grewal(T_local, P_cmb, dIW_val, w_S_core)
        N_met = np.clip((N_bulk * DN) / (m_sil / m_fe + DN), 0.0, 40000.0)

        # Hydrogen partitioning (Clesi et al. 2018)
        DH = D_hydrogen_clesi(T_local, P_cmb, dIW_val)
        H_met = np.clip((H_bulk * DH) / (m_sil / m_fe + DH), 0.0, 10000.0)

        M_core_S[i] = m_fe * w_S_core
        M_core_C[i] = m_fe * (C_met * 1e-6)
        M_core_N[i] = m_fe * (N_met * 1e-6)
        M_core_H[i] = m_fe * (H_met * 1e-6)

    ax_e.plot(time_ma, M_core_S / 1e16, '-', color=STRATA['gold'], linewidth=2.2, label=r'Core Sulfur ($M_\mathrm{S} / 10^{16}$ kg)')
    ax_e.plot(time_ma, M_core_C / 1e14, '-', color=STRATA['magma'], linewidth=2.2, label=r'Core Carbon ($M_\mathrm{C} / 10^{14}$ kg)')
    ax_e.plot(time_ma, M_core_H / 1e13, '-', color=STRATA['plum'], linewidth=2.2, label=r'Core Hydrogen ($M_\mathrm{H} / 10^{13}$ kg)')
    ax_e.plot(time_ma, M_core_N / 1e13, '-', color=STRATA['cobalt'], linewidth=2.2, label=r'Core Nitrogen ($M_\mathrm{N} / 10^{13}$ kg)')

    ax_e.axvspan(1.05, 1.6, color=NEUTRALS['cream'], alpha=0.5, label='Runaway Core Segregation')
    ax_e.set_xlim(0.5, 3.5)
    ax_e.set_ylim(0.0, 16.0)
    ax_e.set_xlabel('Time [Ma]', fontsize=10)
    ax_e.set_ylabel(r'Integrated Core Mass [$10^{13} - 10^{16}$ kg]', fontsize=10)
    ax_e.set_title('Core Volatile Delivery Timeline', fontsize=11, fontweight='bold')
    ax_e.grid(True)
    ax_e.legend(loc='center left', fontsize=8.0)

    # -------------------------------------------------------------------------
    # Panel (f): Planetary Reservoir Mass Partitioning Budgets
    # -------------------------------------------------------------------------
    ax_f = axes[1, 2]
    ax_f.text(0.04, 0.93, '(f)', transform=ax_f.transAxes, fontsize=12, fontweight='bold',
              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

    # Closed-system thermodynamic partition equilibrium at final differentiation conditions
    m_core_final = M_core_fe[-1]
    x_met = m_core_final / M_planet
    x_sil = 1.0 - x_met
    T_final = max(T_core_hist[-1], 1600.0)
    w_S_final = M_core_S[-1] / m_core_final

    DS_final = D_sulfur_boujibar(T_final, P_cmb, dIW_val)
    DC_final = D_carbon_grewal(T_final, P_cmb, dIW_val, w_S_final)
    DN_final = D_nitrogen_grewal(T_final, P_cmb, dIW_val, w_S_final)
    DH_final = D_hydrogen_clesi(T_final, P_cmb, dIW_val)

    D_map = {'Sulfur': DS_final, 'Carbon': DC_final, 'Nitrogen': DN_final, 'Hydrogen': DH_final}
    f_deg_map = {'Sulfur': 0.05, 'Carbon': 0.45, 'Nitrogen': 0.35, 'Hydrogen': 0.22}

    elements = ['Sulfur', 'Carbon', 'Nitrogen', 'Hydrogen']
    core_pct = []
    degassed_pct = []
    mantle_pct = []

    for el in elements:
        D = D_map[el]
        f_c = (x_met * D) / (x_sil + x_met * D)
        pct_c = f_c * 100.0
        pct_d = (100.0 - pct_c) * f_deg_map[el]
        pct_m = 100.0 - pct_c - pct_d
        core_pct.append(pct_c)
        degassed_pct.append(pct_d)
        mantle_pct.append(pct_m)

    core_pct = np.array(core_pct)
    degassed_pct = np.array(degassed_pct)
    mantle_pct = np.array(mantle_pct)

    y_pos = np.arange(len(elements))
    bar_height = 0.55

    ax_f.barh(y_pos, core_pct, height=bar_height, color=STRATA['gold'], edgecolor=NEUTRALS['mist'], label='Segregated Metallic Core')
    ax_f.barh(y_pos, mantle_pct, height=bar_height, left=core_pct, color=STRATA['cobalt'], edgecolor=NEUTRALS['mist'], label='Retained Silicate Mantle')
    ax_f.barh(y_pos, degassed_pct, height=bar_height, left=core_pct + mantle_pct, color=STRATA['magma'], edgecolor=NEUTRALS['mist'], label='Degassed & Vented Loss')

    for i in range(len(elements)):
        c = core_pct[i]
        m = mantle_pct[i]
        d = degassed_pct[i]
        if c > 8:
            ax_f.text(c / 2, y_pos[i], f'{c:.0f}%', ha='center', va='center', fontsize=8.5, fontweight='bold', color=NEUTRALS['graphite'] if c < 60 else 'white')
        if m > 8:
            ax_f.text(c + m / 2, y_pos[i], f'{m:.0f}%', ha='center', va='center', fontsize=8.5, fontweight='bold', color='white')
        if d > 8:
            ax_f.text(c + m + d / 2, y_pos[i], f'{d:.0f}%', ha='center', va='center', fontsize=8.5, fontweight='bold', color='white')

    ax_f.set_yticks(y_pos)
    ax_f.set_yticklabels(elements, fontsize=10, fontweight='bold')
    ax_f.set_xlim(0, 100)
    ax_f.set_xlabel('Elemental Mass Fraction [%]', fontsize=10)
    ax_f.set_title('Planetary Reservoir Allocation', fontsize=11, fontweight='bold')
    ax_f.grid(True, axis='x')
    ax_f.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=8.0)

    # Save outputs
    png_path_assets = os.path.join(ASSETS_DIR, "core_geochemistry_benchmark.png")
    pdf_path_assets = os.path.join(ASSETS_DIR, "core_geochemistry_benchmark.pdf")
    png_path_out = os.path.join(OUTPUT_FILES_DIR, "core_geochemistry_benchmark.png")
    pdf_path_out = os.path.join(OUTPUT_FILES_DIR, "core_geochemistry_benchmark.pdf")

    fig.savefig(png_path_assets, dpi=200)
    fig.savefig(pdf_path_assets)
    fig.savefig(png_path_out, dpi=200)
    fig.savefig(pdf_path_out)
    plt.close(fig)
    print(f"Successfully generated benchmark figure at {png_path_assets}")

if __name__ == "__main__":
    generate_benchmark_figure()
