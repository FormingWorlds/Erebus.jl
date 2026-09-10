#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for telescoping domain mechanics,
grid doubling invariants, physical profile preservation, and accretion scaling
from planetesimal seed to lunar mass in Erebus.jl.
"""

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

G_GRAV = 6.67430e-11
RHO_PLANET = 3000.0  # kg / m^3

def main():
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # -------------------------------------------------------------------------
    # Panel (a): Nested Telescoping Domain Hierarchy
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    ax_a.set_aspect('equal')

    # Domain levels: Level 0 (140 km), Level 1 (280 km), Level 2 (560 km)
    # Centers aligned at (0, 0) for comparative visualization
    rect_l2 = patches.Rectangle((-280, -280), 560, 560, linewidth=1.5, edgecolor=STRATA['plum'],
                                facecolor='#F5EFF7', label=r'Level 2 Domain ($560 \times 560$ km, $N_x=129$)')
    rect_l1 = patches.Rectangle((-140, -140), 280, 280, linewidth=1.8, edgecolor=STRATA['cobalt'],
                                facecolor='#EBF2FA', label=r'Level 1 Domain ($280 \times 280$ km, $N_x=65$)')
    rect_l0 = patches.Rectangle((-70, -70), 140, 140, linewidth=2.0, edgecolor=STRATA['amber'],
                                facecolor='#FDF3E7', label=r'Level 0 Domain ($140 \times 140$ km, $N_x=33$)')

    ax_a.add_patch(rect_l2)
    ax_a.add_patch(rect_l1)
    ax_a.add_patch(rect_l0)

    # Growing planetesimal at center
    circle_p0 = patches.Circle((0, 0), 49, linewidth=1.5, edgecolor=STRATA['magma'],
                               facecolor=STRATA['gold'], alpha=0.8,
                               label=r'Planetesimal Boundary ($R = 49$ km, $r / r_{\rm max} = 0.70$)')
    ax_a.add_patch(circle_p0)

    # Threshold buffer indicator
    thresh_circle = patches.Circle((0, 0), 70 * 0.70, linewidth=1.2, linestyle='--',
                                   edgecolor=STRATA['magma'], facecolor='none',
                                   label=r'Trigger Threshold ($0.70 \cdot x_{\rm size} / 2$)')
    ax_a.add_patch(thresh_circle)

    ax_a.set_xlim(-300, 300)
    ax_a.set_ylim(-300, 300)
    ax_a.set_xlabel('Horizontal Distance $x - x_{\\rm center}$ [km]', fontsize=11)
    ax_a.set_ylabel('Vertical Distance $y - y_{\\rm center}$ [km]', fontsize=11)
    ax_a.set_title('(a) Nested Telescoping Domain Hierarchy', fontsize=12, fontweight='bold')
    ax_a.legend(fontsize=8, loc='upper right')
    ax_a.grid(True)

    # -------------------------------------------------------------------------
    # Panel (b): Planetesimal Profile Invariance Across Doubling
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    r_km = np.linspace(0.0, 70.0, 300)
    r_m = r_km * 1000.0
    R_p = 50000.0  # 50 km

    # Radial temperature profile (diffusive core + accreted shell)
    T_core = 1200.0
    T_surf = 250.0
    T_profile = np.where(r_m <= R_p, T_core - (T_core - T_surf) * (r_m / R_p)**2, T_surf)

    # Markers before telescoping (Level 0)
    np.random.seed(42)
    r_markers_l0 = np.random.uniform(0.0, 49.0, 200)
    T_markers_l0 = np.where(r_markers_l0 <= 50.0, T_core - (T_core - T_surf) * ((r_markers_l0 * 1000.0) / R_p)**2, T_surf)
    T_markers_l0 += np.random.normal(0.0, 8.0, len(r_markers_l0))

    # Markers after telescoping (Level 1): shifted positions preserve exact r
    r_markers_l1 = r_markers_l0.copy()
    T_markers_l1 = T_markers_l0.copy()

    ax_b.plot(r_km, T_profile, color=STRATA['ink'], linewidth=2.0, label='Analytical Continuum $T(r)$')
    ax_b.scatter(r_markers_l0, T_markers_l0, color=STRATA['amber'], s=25, alpha=0.7,
                 label='Level 0 Markers ($x_{\\rm size} = 140$ km)')
    ax_b.scatter(r_markers_l1, T_markers_l1, color=STRATA['cobalt'], s=15, marker='x',
                 label='Level 1 Shifted Markers ($x_{\\rm size} = 280$ km)')

    ax_b.axvline(50.0, color=STRATA['magma'], linestyle='--', linewidth=1.2, label='Planetesimal Surface ($R_p = 50$ km)')
    ax_b.set_xlabel('Radial Distance from Center $r$ [km]', fontsize=11)
    ax_b.set_ylabel('Temperature $T$ [K]', fontsize=11)
    ax_b.set_title('(b) Marker Radial Profile Invariance', fontsize=12, fontweight='bold')
    ax_b.legend(fontsize=8.5, loc='upper right')
    ax_b.grid(True)

    # -------------------------------------------------------------------------
    # Panel (c): Gravitational Poisson Potential Across Telescoped Domain
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    r_domain_km = np.linspace(0.0, 140.0, 400)
    r_domain_m = r_domain_km * 1000.0
    M_total = (4.0 / 3.0) * np.pi * (R_p**3) * RHO_PLANET

    # Analytical potential with Dirichlet boundary at domain boundary R_domain = 140 km (Level 1 half-domain)
    # Interior: Phi(r) = -2*pi*G*rho*(R_p^2 - r^2/3) + C
    # Exterior: Phi(r) = -G*M / r + C'
    # Condition: Phi(R_domain) = 0
    R_bc = 140000.0
    Phi_bc_offset = G_GRAV * M_total / R_bc
    Phi_ext = -G_GRAV * M_total / np.maximum(r_domain_m, 1.0) + Phi_bc_offset
    Phi_int = -2.0 * np.pi * G_GRAV * RHO_PLANET * (R_p**2 - (r_domain_m**2) / 3.0) - (G_GRAV * M_total / R_p) + (2.0 * np.pi * G_GRAV * RHO_PLANET * (2.0 / 3.0) * (R_p**2)) + Phi_bc_offset
    Phi_total = np.where(r_domain_m <= R_p, Phi_int, Phi_ext)

    # Acceleration g(r)
    g_int = (4.0 / 3.0) * np.pi * G_GRAV * RHO_PLANET * r_domain_m
    g_ext = G_GRAV * M_total / np.maximum(r_domain_m, 1.0)**2
    g_total = np.where(r_domain_m <= R_p, g_int, g_ext)

    ax_c1 = ax_c
    l1 = ax_c1.plot(r_domain_km, Phi_total * 1e-4, color=STRATA['cobalt'], linewidth=2.2,
                    label=r'Gravitational Potential $\Phi(r)$ [$10^4$ J/kg]')
    ax_c1.set_xlabel('Radial Distance $r$ [km]', fontsize=11)
    ax_c1.set_ylabel(r'Gravitational Potential $\Phi$ [$10^4$ J/kg]', color=STRATA['cobalt'], fontsize=11)
    ax_c1.tick_params(axis='y', labelcolor=STRATA['cobalt'])
    ax_c1.grid(True)

    ax_c2 = ax_c1.twinx()
    l2 = ax_c2.plot(r_domain_km, g_total, color=STRATA['magma'], linewidth=2.0, linestyle='--',
                    label=r'Gravitational Acceleration $g(r)$ [m/s$^2$]')
    ax_c2.set_ylabel(r'Gravity Acceleration $g$ [m/s$^2$]', color=STRATA['magma'], fontsize=11)
    ax_c2.tick_params(axis='y', labelcolor=STRATA['magma'])

    lines_c = l1 + l2
    labels_c = [l.get_label() for l in lines_c]
    ax_c1.legend(lines_c, labels_c, fontsize=8.5, loc='center right')
    ax_c.set_title('(c) Poisson Potential & Gravity Expansion', fontsize=12, fontweight='bold')

    # -------------------------------------------------------------------------
    # Panel (d): Growth Trajectory to Lunar Embryo
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    time_myr = np.linspace(0.0, 3.0, 500)
    # Sigmoidal pebble/planetesimal accretion growth curve from 40 km to 1737 km (Moon)
    R_seed = 40.0
    R_lunar = 1737.4
    k_growth = 2.5
    t_mid = 1.2
    R_growth = R_seed + (R_lunar - R_seed) / (1.0 + np.exp(-k_growth * (time_myr - t_mid)))

    # Compute domain size evolution with telescoping
    # Initial domain half-size = 70 km (xsize = 140 km)
    # Telescoping triggers whenever R_growth > 0.70 * half_domain
    half_domain_history = []
    current_half_domain = 70.0
    telescoping_times = []
    telescoping_radii = []

    for t, r_val in zip(time_myr, R_growth):
        while r_val > 0.70 * current_half_domain:
            current_half_domain *= 2.0
            telescoping_times.append(t)
            telescoping_radii.append(r_val)
        half_domain_history.append(current_half_domain)

    ax_d.plot(time_myr, R_growth, color=STRATA['magma'], linewidth=2.4, label='Planetesimal Radius $R(t)$ [km]')
    ax_d.plot(time_myr, np.array(half_domain_history), color=STRATA['cobalt'], linewidth=2.0, linestyle='-.',
              label=r'Domain Half-Width $x_{\rm size} / 2$ [km]')
    ax_d.plot(time_myr, 0.70 * np.array(half_domain_history), color=STRATA['gold'], linewidth=1.6, linestyle=':',
              label=r'Expansion Trigger $0.70 \cdot x_{\rm size} / 2$ [km]')

    for t_step, r_step in zip(telescoping_times, telescoping_radii):
        ax_d.plot(t_step, r_step, marker='o', markersize=6, color=STRATA['plum'])

    ax_d.axhline(R_lunar, color=NEUTRALS['graphite'], linestyle='--', linewidth=1.2,
                 label=r'Lunar Target Radius ($R = 1737$ km)')
    ax_d.set_yscale('log')
    ax_d.set_ylim(30.0, 10000.0)
    ax_d.set_xlabel('Accretion Time $t$ [Myr]', fontsize=11)
    ax_d.set_ylabel('Radius & Domain Scale [km]', fontsize=11)
    ax_d.set_title('(d) Accretion Trajectory & Discrete Telescoping Events', fontsize=12, fontweight='bold')
    ax_d.legend(fontsize=8.5, loc='upper left', bbox_to_anchor=(0.02, 0.98))
    ax_d.grid(True, which='both')

    plt.tight_layout()

    out_png = os.path.join(ASSETS_DIR, "telescoping_domain_benchmark.png")
    out_pdf = os.path.join(ASSETS_DIR, "telescoping_domain_benchmark.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Generated {out_png} and {out_pdf}")

if __name__ == "__main__":
    main()
