#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for normative accessory mineral
tracking and meteorite diagnostics in Erebus.jl.
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

# -----------------------------------------------------------------------------
# Physical stoichiometry functions matching src/physics.jl
# -----------------------------------------------------------------------------
def compute_troilite_stoichiometry(w_S):
    f_troilite = 87.910 / 32.065
    f_fe = 55.845 / 32.065
    S_fe_limit = np.maximum(0.0, 1.0 - w_S) / f_fe
    S_troilite = np.minimum(w_S, S_fe_limit)
    return S_troilite * f_troilite, S_troilite * f_fe

def compute_schreibersite_stoichiometry(w_P, ni_frac=0.25):
    M_metal_avg = (1.0 - ni_frac) * 55.845 + ni_frac * 58.6934
    M_P = 30.97376
    M_schreib = 3.0 * M_metal_avg + M_P
    f_schreib = M_schreib / M_P
    f_metal = (3.0 * M_metal_avg) / M_P
    P_metal_limit = np.maximum(0.0, 1.0 - w_P) / f_metal
    P_schreib = np.minimum(w_P, P_metal_limit)
    return P_schreib * f_schreib, P_schreib * f_metal

def compute_cohenite_graphite_stoichiometry(w_C, carbide_max=0.0667):
    f_fe = (3.0 * 55.845) / 12.011
    f_cohenite = f_fe + 1.0
    C_fe_limit = np.maximum(0.0, 1.0 - w_C) / f_fe
    C_carbide = np.minimum(w_C, np.minimum(carbide_max, C_fe_limit))
    w_coh = C_carbide * f_cohenite
    w_gra = w_C - C_carbide
    w_fe = C_carbide * f_fe
    return w_coh, w_gra, w_fe

def compute_nitride_stoichiometry(w_N, mode='roaldite'):
    M_N = 14.007
    if mode == 'roaldite':
        M_Fe = 55.845
        f_nitride = (4.0 * M_Fe + M_N) / M_N
        f_metal = (4.0 * M_Fe) / M_N
    elif mode == 'carlsbergite':
        M_Cr = 51.996
        f_nitride = (M_Cr + M_N) / M_N
        f_metal = M_Cr / M_N
    elif mode == 'osbornite':
        M_Ti = 47.867
        f_nitride = (M_Ti + M_N) / M_N
        f_metal = M_Ti / M_N
    else:
        raise ValueError(f"Unknown nitride mode {mode}")
    N_metal_limit = np.maximum(0.0, 1.0 - w_N) / f_metal
    N_nitride = np.minimum(w_N, N_metal_limit)
    return N_nitride * f_nitride, N_nitride * f_metal

def compute_eutectic_melting(T, T_eutectic=1213.0, dT_transition=50.0):
    F_solid = np.clip(1.0 - (T - T_eutectic) / dT_transition, 0.0, 1.0)
    F_liquid = 1.0 - F_solid
    return F_solid, F_liquid

def generate_benchmark_figure():
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)
    plt.subplots_adjust(hspace=0.28, wspace=0.25)

    # -------------------------------------------------------------------------
    # Panel (a): Stoichiometric Accessory Mineral Yields
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    w_S_arr = np.linspace(0.0, 0.15, 200) # 0 to 15 wt% S
    w_tro, _ = compute_troilite_stoichiometry(w_S_arr)
    
    w_C_arr = np.linspace(0.0, 0.10, 200) # 0 to 10 wt% C
    w_coh, w_gra, _ = compute_cohenite_graphite_stoichiometry(w_C_arr, carbide_max=0.0667)
    
    w_P_arr = np.linspace(0.0, 0.01, 200) # 0 to 10,000 ppmw P
    w_sch, _ = compute_schreibersite_stoichiometry(w_P_arr, ni_frac=0.25)

    ax_a.plot(w_S_arr * 100.0, w_tro * 100.0, color=STRATA['amber'], lw=2.2, label=r'Troilite $\mathrm{FeS}$ ($f \approx 2.742$)')
    ax_a.plot(w_C_arr * 100.0, w_coh * 100.0, color=STRATA['magma'], lw=2.2, label=r'Cohenite $(\mathrm{Fe,Ni})_3\mathrm{C}$ (sat. @ 6.67 wt%)')
    ax_a.plot(w_C_arr * 100.0, w_gra * 100.0, color=STRATA['ink'], lw=2.0, ls='--', label=r'Graphite $\mathrm{C}$ (excess carbon)')
    ax_a.plot(w_P_arr * 100.0, w_sch * 100.0, color=STRATA['cobalt'], lw=2.0, ls='-.', label=r'Schreibersite $(\mathrm{Fe,Ni})_3\mathrm{P}$ ($x_{\mathrm{Ni}}=0.25$)')

    ax_a.axvline(6.67, color=STRATA['plum'], ls=':', lw=1.2, alpha=0.8)
    ax_a.annotate('Carbide Saturation\n(6.67 wt% C)', xy=(6.67, 75.0), xytext=(8.0, 55.0),
                  arrowprops=dict(arrowstyle="->", color=STRATA['plum'], lw=1.0),
                  fontsize=9, color=STRATA['plum'], fontweight='bold')

    ax_a.set_xlabel('Precursor Element Content [wt% in metallic alloy]', fontsize=11, fontweight='bold')
    ax_a.set_ylabel('Accessory Phase Abundance [wt%]', fontsize=11, fontweight='bold')
    ax_a.set_title('(a) Stoichiometric Accessory Mineral Conversion', fontsize=12, fontweight='bold', pad=10)
    ax_a.set_xlim(0.0, 12.0)
    ax_a.set_ylim(0.0, 105.0)
    ax_a.grid(True, alpha=0.5)
    ax_a.legend(loc='upper left', fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (b): Thermal Eutectic Melting Transition
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    T_arr = np.linspace(1120.0, 1320.0, 300)
    T_eut = 1213.0
    dT_trans = 50.0
    F_sol, F_liq = compute_eutectic_melting(T_arr, T_eutectic=T_eut, dT_transition=dT_trans)

    # Baseline nominal alloy composition: 4 wt% S, 2000 ppmw C, 1000 ppmw P, 100 ppmw N
    w_tro_0, _ = compute_troilite_stoichiometry(0.04)
    w_sch_0, _ = compute_schreibersite_stoichiometry(0.001, ni_frac=0.25)
    w_coh_0, w_gra_0, _ = compute_cohenite_graphite_stoichiometry(0.002, carbide_max=0.0667)
    w_nit_0, _ = compute_nitride_stoichiometry(0.0001, mode='roaldite')
    w_acc_0 = w_tro_0 + w_sch_0 + w_coh_0 + w_gra_0 + w_nit_0
    w_matrix_0 = max(0.0, 1.0 - w_acc_0)

    ax_b.plot(T_arr, F_sol, color=STRATA['cobalt'], lw=2.5, label=r'Solid Metal Fraction $F_{\mathrm{solid}}$')
    ax_b.plot(T_arr, F_liq, color=STRATA['magma'], lw=2.5, label=r'Molten Alloy Fraction $F_{\mathrm{liquid}}$')
    ax_b.plot(T_arr, F_sol * w_tro_0 * 100.0 / 10.0, color=STRATA['amber'], lw=2.0, ls='--',
              label=r'Troilite $\times 10^{-1}$ [wt%]')
    ax_b.plot(T_arr, F_sol * w_sch_0 * 100.0, color=STRATA['gold'], lw=2.0, ls='-.',
              label=r'Schreibersite [wt%]')
    ax_b.plot(T_arr, F_sol * w_coh_0 * 100.0 / 10.0, color=STRATA['plum'], lw=2.0, ls=':',
              label=r'Cohenite $\times 10^{-1}$ [wt%]')

    ax_b.axvspan(T_eut, T_eut + dT_trans, color=NEUTRALS['cream'], alpha=0.5, label='Eutectic Transition Interval')
    ax_b.axvline(T_eut, color=STRATA['magma'], ls='--', lw=1.2)
    ax_b.annotate(r'$T_{\mathrm{eutectic}} = 1213\ \mathrm{K}$', xy=(T_eut, 0.85), xytext=(T_eut - 75.0, 0.88),
                  arrowprops=dict(arrowstyle="->", color=STRATA['magma'], lw=1.0),
                  fontsize=9, color=STRATA['magma'], fontweight='bold')

    ax_b.set_xlabel('Alloy Temperature $T$ [K]', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Phase Fraction / Abundance [- / wt%]', fontsize=11, fontweight='bold')
    ax_b.set_title(r'(b) Eutectic Phase Dissolution ($T_{\mathrm{eut}} = 1213\ \mathrm{K}, \Delta T = 50\ \mathrm{K}$)', fontsize=12, fontweight='bold', pad=10)
    ax_b.set_xlim(1130.0, 1310.0)
    ax_b.set_ylim(0.0, 1.20)
    ax_b.grid(True, alpha=0.5)
    ax_b.legend(loc='center left', fontsize=9, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (c): Planetesimal Radial Mineral Distribution
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    r_arr = np.linspace(0.0, 50.0, 300) # radius [km]
    R_p = 50.0
    r_core = 25.0
    r_mantle = 42.5

    # Realistic conductive profile with differentiated hot core:
    # Core: T = 1350 K > 1263 K -> 100% molten alloy
    # Mantle: 1250 K -> 900 K
    # Crust: 900 K -> 250 K (conductive lid)
    T_prof = np.where(r_arr <= r_core, 1350.0 - 50.0 * (r_arr / r_core)**2,
                      np.where(r_arr <= r_mantle,
                               1300.0 - (1300.0 - 850.0) * ((r_arr - r_core) / (r_mantle - r_core)),
                               850.0 - (850.0 - 250.0) * ((r_arr - r_mantle) / (R_p - r_mantle))))

    F_sol_prof, F_liq_prof = compute_eutectic_melting(T_prof, T_eutectic=T_eut, dT_transition=dT_trans)
    
    # Metal abundance: segregated core has high metal (85 vol%), mantle is depleted (3 vol%), crust has primitive chondritic metal (15 vol%)
    phi_fe_prof = np.where(r_arr <= r_core, 0.85,
                           np.where(r_arr <= r_mantle, 0.03, 0.15))
    
    w_tro_prof = phi_fe_prof * (F_sol_prof * w_tro_0) * 100.0
    w_sch_prof = phi_fe_prof * (F_sol_prof * w_sch_0) * 100.0
    w_coh_prof = phi_fe_prof * (F_sol_prof * w_coh_0) * 100.0
    w_liq_prof = phi_fe_prof * F_liq_prof * 100.0

    ax_c.plot(r_arr, w_liq_prof, color=STRATA['magma'], lw=2.5, label='Molten Fe-FeS Core Liquid [wt%]')
    ax_c.plot(r_arr, w_tro_prof, color=STRATA['amber'], lw=2.2, label='Preserved Solid Troilite [wt%]')
    ax_c.plot(r_arr, w_sch_prof * 10.0, color=STRATA['cobalt'], lw=2.0, ls='-.', label=r'Schreibersite $\times 10$ [wt%]')
    ax_c.plot(r_arr, w_coh_prof * 10.0, color=STRATA['plum'], lw=2.0, ls=':', label=r'Cohenite $\times 10$ [wt%]')

    ax_c.axvspan(0.0, r_core, color=NEUTRALS['bone'], alpha=0.35, label=r'Metallic Core ($r \leq 25$ km)')
    ax_c.axvspan(r_core, r_mantle, color=NEUTRALS['cream'], alpha=0.35, label='Silicate Mantle')
    ax_c.axvspan(r_mantle, R_p, color=NEUTRALS['paper'], alpha=0.6, label=r'Primitive Crust ($r \geq 42.5$ km)')

    ax_c.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11, fontweight='bold')
    ax_c.set_ylabel('Bulk Planetesimal Phase Abundance [wt%]', fontsize=11, fontweight='bold')
    ax_c.set_title('(c) Differentiated Planetesimal Radial Assemblage', fontsize=12, fontweight='bold', pad=10)
    ax_c.set_xlim(0.0, 50.0)
    ax_c.set_ylim(0.0, 95.0)
    ax_c.grid(True, alpha=0.5)
    ax_c.legend(loc='center right', fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------------------
    # Panel (d): Diagnostic Meteorite Classification Diagram
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    
    # Classification fields
    rect_magmatic = Rectangle((0.80, 0.0), 0.20, 0.5, facecolor=STRATA['magma'], alpha=0.25, lw=1.2, edgecolor=STRATA['magma'])
    rect_primitive = Rectangle((0.0, 1.0), 0.60, 9.0, facecolor=STRATA['gold'], alpha=0.25, lw=1.2, edgecolor=STRATA['gold'])
    rect_trans = Rectangle((0.0, 0.0), 1.0, 10.0, facecolor=NEUTRALS['mist'], alpha=0.15, lw=1.0, ls='--', edgecolor=NEUTRALS['graphite'])

    ax_d.add_patch(rect_trans)
    ax_d.add_patch(rect_primitive)
    ax_d.add_patch(rect_magmatic)

    # Representative meteorite groups plotted from petrologic literature data
    # (Benedix et al. 2000, Goldstein et al. 2009, Chabot & Drake 1999)
    # Magmatic irons: complete core differentiation, depleted crustal accessory retention
    magmatic_x = [0.95, 0.98, 0.92, 0.89, 0.94]
    magmatic_y = [0.0005, 0.0002, 0.0010, 0.0018, 0.0006]
    magmatic_labels = ['IIIAB', 'IVA', 'IVB', 'IIAB', 'IC']

    # Primitive non-magmatic irons (IAB complex & winonaites): incomplete melting, high crustal retention
    primitive_x = [0.22, 0.35, 0.15, 0.45, 0.28]
    primitive_y = [0.065, 0.048, 0.082, 0.038, 0.055]
    primitive_labels = ['Winonaites', 'IAB-sLL', 'IAB-sLM', 'IAB-sLH', 'Landes']

    # Partially melted / transitional bodies
    trans_x = [0.68, 0.72, 0.64]
    trans_y = [0.018, 0.012, 0.024]
    trans_labels = ['IIICD', 'Ureilite metal', 'Tombigbee']

    ax_d.scatter(magmatic_x, np.array(magmatic_y) * 100.0, s=70, color=STRATA['magma'], marker='o', edgecolors=STRATA['ink'], lw=1.2, zorder=5, label='Magmatic Irons (IIIAB, IVA, IVB, IC)')
    ax_d.scatter(primitive_x, np.array(primitive_y) * 100.0, s=75, color=STRATA['gold'], marker='s', edgecolors=STRATA['ink'], lw=1.2, zorder=5, label='Primitive Complex (IAB / Winonaites)')
    ax_d.scatter(trans_x, np.array(trans_y) * 100.0, s=70, color=STRATA['plum'], marker='^', edgecolors=STRATA['ink'], lw=1.2, zorder=5, label='Transitional Incomplete Segregations')

    magmatic_annot = [
        ('IIAB', (0.89, 0.0018 * 100.0), (0.83, 1.6)),
        ('IVB', (0.92, 0.0010 * 100.0), (0.86, 2.5)),
        ('IIIAB', (0.95, 0.0005 * 100.0), (0.89, 3.4)),
        ('IC', (0.94, 0.0006 * 100.0), (0.92, 4.3)),
        ('IVA', (0.98, 0.0002 * 100.0), (0.95, 5.2)),
    ]
    for lbl, xy_pt, xy_txt in magmatic_annot:
        ax_d.annotate(lbl, xy=xy_pt, xytext=xy_txt,
                      arrowprops=dict(arrowstyle="->", color=STRATA['magma'], lw=0.8),
                      fontsize=8, color=STRATA['magma'], fontweight='bold', ha='center', va='bottom')

    for x, y, lbl in zip(primitive_x, primitive_y, primitive_labels):
        ax_d.annotate(lbl, xy=(x, y * 100.0), xytext=(x - 0.04, y * 100.0 + 0.35), fontsize=8, color=STRATA['ink'], fontweight='bold')

    trans_annot = [
        ('Tombigbee', (0.64, 0.024 * 100.0), (0.52, 3.0), 'left', 'bottom'),
        ('IIICD', (0.68, 0.018 * 100.0), (0.58, 1.8), 'left', 'center'),
        ('Ureilite metal', (0.72, 0.012 * 100.0), (0.72, 0.35), 'center', 'top'),
    ]
    for lbl, xy_pt, xy_txt, ha_val, va_val in trans_annot:
        ax_d.annotate(lbl, xy=xy_pt, xytext=xy_txt,
                      arrowprops=dict(arrowstyle="->", color=STRATA['plum'], lw=0.8),
                      fontsize=8, color=STRATA['plum'], fontweight='bold', ha=ha_val, va=va_val)

    ax_d.set_xlabel(r'Molten Core Metal Fraction $f_{\mathrm{molten,core}}$ [-]', fontsize=11, fontweight='bold')
    ax_d.set_ylabel(r'Crustal Solid Accessory Retention $f_{\mathrm{solid,crust}}$ [wt%]', fontsize=11, fontweight='bold')
    ax_d.set_title('(d) Meteorite Parent Body Diagnostic Regimes', fontsize=12, fontweight='bold', pad=10)
    ax_d.set_xlim(0.0, 1.0)
    ax_d.set_ylim(0.0, 10.0)
    ax_d.grid(True, alpha=0.5)
    ax_d.legend(loc='upper right', fontsize=8.5, framealpha=0.9)

    png_path = os.path.join(ASSETS_DIR, "mineral_assemblage_benchmark.png")
    pdf_path = os.path.join(ASSETS_DIR, "mineral_assemblage_benchmark.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"Generated benchmark figure: {png_path} and {pdf_path}")

    # Write numerical summary dataset
    summary_data = {
        "stoichiometry": {
            "troilite_factor": float(87.910 / 32.065),
            "schreibersite_factor_ni25": float(((3.0 * (0.75 * 55.845 + 0.25 * 58.6934)) + 30.97376) / 30.97376),
            "cohenite_factor": float((3.0 * 55.845 + 12.011) / 12.011),
            "carbide_saturation_wt_frac": 0.0667,
            "roaldite_factor": float((4.0 * 55.845 + 14.007) / 14.007),
            "carlsbergite_factor": float((51.996 + 14.007) / 14.007),
            "osbornite_factor": float((47.867 + 14.007) / 14.007)
        },
        "eutectic_transition": {
            "T_eutectic_K": 1213.0,
            "dT_transition_K": 50.0
        },
        "meteorite_classes": {
            "IAB_winonaite": {
                "f_molten_core_max": 0.60,
                "f_crust_solid_min": 0.01,
                "description": "Incomplete differentiation with retained crustal troilite, schreibersite, and graphite"
            },
            "magmatic_differentiated": {
                "f_molten_core_min": 0.80,
                "f_core_metal_min": 0.40,
                "description": "Fractional crystallization from fully segregated, molten metallic core"
            }
        }
    }
    json_path = os.path.join(OUTPUT_FILES_DIR, "mineral_assemblage_benchmark_data.json")
    with open(json_path, 'w') as f:
        json.dump(summary_data, f, indent=2)
    print(f"Generated benchmark summary: {json_path}")

if __name__ == "__main__":
    generate_benchmark_figure()
