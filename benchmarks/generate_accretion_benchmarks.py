#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for planetesimal accretion engine,
pebble accretion regimes, Safronov focusing, impact heating, and onion-shell thermal
structure in Erebus.jl.
"""

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

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

# Constants
G_GRAV = 6.67430e-11
M_SUN = 1.98847e30
AU_METERS = 1.495978707e11
SEC_PER_YEAR = 3.15576e7
K_BOLTZMANN = 1.380649e-23
M_PROTON = 1.67262192e-27

# -----------------------------------------------------------------------------
# Physics functions
# -----------------------------------------------------------------------------
def compute_keplerian_frequency(a_m, M_star=M_SUN):
    return np.sqrt(G_GRAV * M_star / (a_m**3))

def compute_sound_speed(T_gas, gamma=1.0, mu_gas=2.34):
    return np.sqrt(gamma * K_BOLTZMANN * T_gas / (mu_gas * M_PROTON))

def compute_pebble_accretion(M, M_star, a_m, Sigma_peb, St, c_s, alpha_turb):
    Omega_K = compute_keplerian_frequency(a_m, M_star)
    v_K = np.sqrt(G_GRAV * M_star / a_m)
    H_gas = c_s / Omega_K
    H_peb = H_gas * np.sqrt(alpha_turb / (alpha_turb + St))
    R_B = G_GRAV * M / (c_s**2)
    R_H = a_m * ((M / (3.0 * M_star))**(1.0 / 3.0))

    eta_disk = 1.5 * (c_s / v_K)**2
    v_rel = eta_disk * v_K
    rho_peb = Sigma_peb / (np.sqrt(2.0 * np.pi) * H_peb)

    # Bondi regime: capture radius bounded by Hill radius (Lambrechts & Johansen 2012 eq. 6)
    r_acc_bondi = np.minimum(R_H, 2.0 * np.sqrt((St / Omega_K) * G_GRAV * M / np.maximum(v_rel, 1.0e-6)))
    dM_dt_bondi_2d = 2.0 * r_acc_bondi * Sigma_peb * v_rel
    dM_dt_bondi_3d = np.pi * (r_acc_bondi**2) * rho_peb * v_rel
    dM_dt_bondi = np.minimum(dM_dt_bondi_2d, dM_dt_bondi_3d)

    # Hill regime
    r_acc_hill = R_H * (St**(1.0 / 3.0))
    v_H = Omega_K * R_H
    dM_dt_hill_2d = 2.0 * r_acc_hill * Sigma_peb * v_H
    dM_dt_hill_3d = np.pi * (r_acc_hill**2) * rho_peb * v_H
    dM_dt_hill = np.minimum(dM_dt_hill_2d, dM_dt_hill_3d)

    M_trans = np.sqrt(1.0 / 3.0) * (v_rel**3) / (G_GRAV * Omega_K)
    dM_dt_auto = dM_dt_bondi if M < M_trans else dM_dt_hill
    return dM_dt_bondi, dM_dt_hill, dM_dt_auto, R_H, R_B, H_peb

def compute_safronov_accretion(M, R, Sigma_pl, v_disp, Omega_K):
    v_esc_sq = 2.0 * G_GRAV * M / R
    Theta = v_esc_sq / (2.0 * (v_disp**2))
    F_g = 1.0 + 2.0 * Theta
    return np.pi * (R**2) * Sigma_pl * Omega_K * F_g

def compute_impact_heating(M, R, cp=1000.0, h_impact=0.5, v_inf=0.0):
    u_acc = (G_GRAV * M / R) + 0.5 * (v_inf**2)
    return h_impact * u_acc / cp

def main():
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plt.subplots_adjust(hspace=0.32, wspace=0.28)

    # -------------------------------------------------------------------------
    # Panel (a): Pebble Accretion Regimes & Transition Mass
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    a_m = 2.5 * AU_METERS
    T_gas = 140.0
    c_s = compute_sound_speed(T_gas)
    Omega_K = compute_keplerian_frequency(a_m)
    Sigma_peb = 30.0  # kg/m^2
    St = 0.05
    alpha_turb = 1.0e-3

    M_arr = np.logspace(16, 24, 300)
    dM_bondi = []
    dM_hill = []
    dM_auto = []
    for m in M_arr:
        b, h, a, _, _, _ = compute_pebble_accretion(m, M_SUN, a_m, Sigma_peb, St, c_s, alpha_turb)
        dM_bondi.append(b)
        dM_hill.append(h)
        dM_auto.append(a)

    v_K = np.sqrt(G_GRAV * M_SUN / a_m)
    v_rel = 1.5 * (c_s / v_K)**2 * v_K
    M_trans = np.sqrt(1.0 / 3.0) * (v_rel**3) / (G_GRAV * Omega_K)

    ax_a.loglog(M_arr, dM_bondi, label=r'Bondi Pebble Accretion ($\dot{M}_{\rm B}$)',
                color=STRATA['cobalt'], linestyle='--', linewidth=1.8)
    ax_a.loglog(M_arr, dM_hill, label=r'Hill Pebble Accretion ($\dot{M}_{\rm H}$)',
                color=STRATA['amber'], linestyle='--', linewidth=1.8)
    ax_a.loglog(M_arr, dM_auto, label=r'Auto Transition (:pebble_auto)',
                color=STRATA['magma'], linewidth=2.4)
    ax_a.axvline(M_trans, color=NEUTRALS['graphite'], linestyle=':', linewidth=1.5,
                 label=r'Transition Mass $M_{\rm trans}$')

    ax_a.set_xlabel('Planetesimal Mass $M$ [kg]', fontsize=11)
    ax_a.set_ylabel(r'Accretion Rate $\dot{M}$ [kg / s]', fontsize=11)
    ax_a.set_title('(a) Pebble Accretion Regimes & Bondi-Hill Transition', fontsize=12, fontweight='bold')
    ax_a.legend(fontsize=9, loc='upper left')
    ax_a.grid(True, which='both')

    # -------------------------------------------------------------------------
    # Panel (b): Safronov Gravitational Focusing
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    R_arr = np.linspace(10000.0, 100000.0, 200)
    rho_bulk = 3000.0
    Sigma_pl = 100.0
    v_disps = [50.0, 100.0, 200.0, 500.0]
    colors = [STRATA['cobalt'], STRATA['gold'], STRATA['amber'], STRATA['magma']]

    for vd, col in zip(v_disps, colors):
        rates = []
        for r in R_arr:
            m = (4.0 / 3.0) * np.pi * (r**3) * rho_bulk
            rates.append(compute_safronov_accretion(m, r, Sigma_pl, vd, Omega_K))
        ax_b.plot(R_arr / 1000.0, np.array(rates) / 1.0e12,
                  label=rf'$\sigma_v = {vd:.0f}$ m/s', color=col, linewidth=2.0)

    # Geometric limit (no focusing, Theta -> 0)
    geom_rates = [np.pi * (r**2) * Sigma_pl * Omega_K / 1.0e12 for r in R_arr]
    ax_b.plot(R_arr / 1000.0, geom_rates, label='Geometric limit ($F_g = 1$)',
              color=NEUTRALS['graphite'], linestyle=':', linewidth=1.6)

    ax_b.set_xlabel('Planetesimal Radius $R$ [km]', fontsize=11)
    ax_b.set_ylabel(r'Accretion Rate $\dot{M}$ [$10^{12}$ kg / s]', fontsize=11)
    ax_b.set_title('(b) Safronov Gravitational Focusing Swarm', fontsize=12, fontweight='bold')
    ax_b.legend(fontsize=9, loc='upper left')
    ax_b.grid(True)

    # -------------------------------------------------------------------------
    # Panel (c): Impact Heating Temperature Rise
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    R_vals_km = np.linspace(10.0, 150.0, 200)
    h_efficiencies = [0.2, 0.5, 0.8, 1.0]
    h_colors = [STRATA['cobalt'], STRATA['gold'], STRATA['amber'], STRATA['magma']]

    for h_eff, col in zip(h_efficiencies, h_colors):
        dT_vals = []
        for r_km in R_vals_km:
            r = r_km * 1000.0
            m = (4.0 / 3.0) * np.pi * (r**3) * rho_bulk
            dT_vals.append(compute_impact_heating(m, r, cp=1000.0, h_impact=h_eff, v_inf=200.0))
        ax_c.plot(R_vals_km, dT_vals, label=rf'$h_{{\rm impact}} = {h_eff:.1f}$',
                  color=col, linewidth=2.0)

    ax_c.set_xlabel('Planetesimal Radius $R$ [km]', fontsize=11)
    ax_c.set_ylabel(r'Impact Temperature Rise $\Delta T_{\rm impact}$ [K]', fontsize=11)
    ax_c.set_title('(c) Accretion Impact Heating Thermodynamics', fontsize=12, fontweight='bold')
    ax_c.legend(fontsize=9, loc='upper left')
    ax_c.grid(True)

    # -------------------------------------------------------------------------
    # Panel (d): Onion-Shell Thermal & Clock Structure
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    # Radial coordinate from center to 60 km
    r_coords_km = np.linspace(0.0, 60.0, 300)
    # Radiogenic clock: core accreted at t=0, outer shells up to 2.5 Myr
    t_acc_myr = np.where(r_coords_km <= 20.0, 0.0,
                         (r_coords_km - 20.0) / 40.0 * 2.5)
    tau_al_myr = 0.717
    Q_norm = np.exp(-t_acc_myr / tau_al_myr)

    # Volatile water fraction step across snowline (suppose snowline crossed at 45 km)
    H2O_wtpct = np.where(r_coords_km < 45.0, 0.1, 10.0)

    ax_d1 = ax_d
    line1 = ax_d1.plot(r_coords_km, Q_norm, color=STRATA['magma'], linewidth=2.2,
                       label=r'Normalized Radiogenic Power $Q(t_{\rm acc}) / Q_0$')
    ax_d1.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11)
    ax_d1.set_ylabel(r'Radiogenic Power $Q_{\rm rad} / Q_0$ [-]', color=STRATA['magma'], fontsize=11)
    ax_d1.tick_params(axis='y', labelcolor=STRATA['magma'])
    ax_d1.grid(True)

    ax_d2 = ax_d1.twinx()
    line2 = ax_d2.plot(r_coords_km, H2O_wtpct, color=STRATA['cobalt'], linewidth=2.0,
                       linestyle='--', label=r'Accreted $\rm H_2O$ Content [wt%]')
    ax_d2.set_ylabel(r'Water Content [wt%]', color=STRATA['cobalt'], fontsize=11)
    ax_d2.tick_params(axis='y', labelcolor=STRATA['cobalt'])

    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax_d1.legend(lines, labels, fontsize=8.5, loc='center left')
    ax_d.set_title('(d) Onion-Shell Radiogenic Clock & Snowline Volatiles', fontsize=12, fontweight='bold')

    plt.tight_layout()

    out_png = os.path.join(ASSETS_DIR, "planetesimal_accretion_benchmark.png")
    out_pdf = os.path.join(ASSETS_DIR, "planetesimal_accretion_benchmark.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Generated {out_png} and {out_pdf}")

if __name__ == "__main__":
    main()
