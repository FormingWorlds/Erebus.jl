#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for multi-stage planetesimal
accretion sequence (planetesimals before pebbles), aerodynamic onset mass
(Visser & Ormel 2016), pebble isolation mass (Lambrechts et al. 2014), and
smoothstep regime transitions in Erebus.jl.
"""

import json
import os
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
M_EARTH = 5.972e24
AU_METERS = 1.495978707e11
SEC_PER_YEAR = 3.15576e7
K_BOLTZMANN = 1.380649e-23
M_PROTON = 1.67262192e-27


def compute_keplerian_frequency(a_m, M_star=M_SUN):
    return np.sqrt(G_GRAV * M_star / (a_m**3))


def compute_keplerian_velocity(a_m, M_star=M_SUN):
    return np.sqrt(G_GRAV * M_star / a_m)


def compute_sound_speed(T_gas, gamma=1.0, mu_gas=2.34):
    return np.sqrt(gamma * K_BOLTZMANN * T_gas / (mu_gas * M_PROTON))


def compute_headwind_velocity(c_s, v_K):
    eta = 1.5 * (c_s / v_K)**2
    return eta * v_K


def compute_onset_mass(M_star, a_m, St, c_s, f_onset=1.0):
    Omega_K = compute_keplerian_frequency(a_m, M_star)
    v_K = compute_keplerian_velocity(a_m, M_star)
    v_hw = compute_headwind_velocity(c_s, v_K)
    return f_onset * (v_hw**3) * St / (G_GRAV * Omega_K)


def compute_isolation_mass(M_star, a_m, c_s, f_iso=0.5):
    v_K = compute_keplerian_velocity(a_m, M_star)
    h_aspect = c_s / v_K
    return f_iso * M_star * (h_aspect**3)


def smoothstep(x):
    xc = np.clip(x, 0.0, 1.0)
    return xc * xc * (3.0 - 2.0 * xc)


def main():
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plt.subplots_adjust(hspace=0.32, wspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Pebble Accretion Onset Mass vs Orbital Distance
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    a_au_vals = np.linspace(0.5, 5.0, 100)
    a_m_vals = a_au_vals * AU_METERS
    T_midplane = 200.0 * (a_au_vals**(-0.5))
    c_s_vals = compute_sound_speed(T_midplane)

    st_colors = [STRATA['cobalt'], STRATA['amber'], STRATA['magma']]
    st_labels = [r'$	au_s = 0.01$', r'$	au_s = 0.05$', r'$	au_s = 0.10$']

    for St, col, lbl in zip([0.01, 0.05, 0.10], st_colors, st_labels):
        M_onsets = [compute_onset_mass(M_SUN, a, St, cs) for a, cs in zip(a_m_vals, c_s_vals)]
        ax_a.plot(a_au_vals, np.array(M_onsets) / 1.0e20, color=col, lw=2.2, label=lbl)

    ax_a.set_xlabel('Orbital Distance a [AU]', fontsize=11, fontweight='bold')
    ax_a.set_ylabel(r'Onset Mass $M_{\mathrm{onset}}$ [$10^{20}\ \mathrm{kg}$]', fontsize=11, fontweight='bold')
    ax_a.set_title('(a) Settling Regime Onset Mass (Visser & Ormel 2016)', fontsize=11, fontweight='bold', pad=10)
    ax_a.grid(True, linestyle=':', alpha=0.6)
    ax_a.set_xlim(0.5, 5.0)
    ax_a.set_ylim(0.0, 120.0)
    ax_a.legend(loc='upper left', fontsize=10, framealpha=0.95)

    # -------------------------------------------------------------------------
    # Panel (b): Pebble Isolation Mass vs Orbital Distance
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    M_iso_vals = np.array([compute_isolation_mass(M_SUN, a, cs) for a, cs in zip(a_m_vals, c_s_vals)])
    M_iso_earth = M_iso_vals / M_EARTH

    ax_b.plot(a_au_vals, M_iso_earth, color=STRATA['plum'], lw=2.5, label=r'$M_{\mathrm{iso}}$ (Lambrechts et al. 2014)')
    ax_b.plot(a_au_vals, M_iso_earth * 0.5, color=STRATA['gold'], lw=1.8, linestyle='--', label=r'$0.5 	imes M_{\mathrm{iso}}$')
    ax_b.set_xlabel('Orbital Distance a [AU]', fontsize=11, fontweight='bold')
    ax_b.set_ylabel(r'Pebble Isolation Mass [$M_{\oplus}$]', fontsize=11, fontweight='bold')
    ax_b.set_title('(b) Pebble Isolation Mass (Lambrechts et al. 2014)', fontsize=11, fontweight='bold', pad=10)
    ax_b.grid(True, linestyle=':', alpha=0.6)
    ax_b.set_xlim(0.5, 5.0)
    ax_b.set_ylim(0.0, 30.0)
    ax_b.legend(loc='upper left', fontsize=10, framealpha=0.95)

    # -------------------------------------------------------------------------
    # Panel (c): Accretion Rate across 3 Regimes (Sharp vs Smooth)
    # -------------------------------------------------------------------------
    ax_c = axes[1, 0]
    a_fixed = 2.5 * AU_METERS
    T_fixed = 150.0
    cs_fixed = compute_sound_speed(T_fixed)
    Omega_fixed = compute_keplerian_frequency(a_fixed)
    v_K_fixed = compute_keplerian_velocity(a_fixed)
    M_on_fixed = compute_onset_mass(M_SUN, a_fixed, 0.05, cs_fixed)
    M_iso_fixed = compute_isolation_mass(M_SUN, a_fixed, cs_fixed)

    mass_grid = np.logspace(17, 27, 300)
    rate_sharp = []
    rate_smooth = []

    Sigma_pl = 100.0
    v_disp = 100.0
    Sigma_peb = 50.0 * (2.5**(-1.0))
    H_gas = cs_fixed / Omega_fixed
    H_peb = H_gas * np.sqrt(1e-3 / (1e-3 + 0.05))
    rho_peb = Sigma_peb / (np.sqrt(2.0 * np.pi) * H_peb)
    v_rel = compute_headwind_velocity(cs_fixed, v_K_fixed)

    for M_val in mass_grid:
        R_val = (3.0 * M_val / (4.0 * np.pi * 3000.0))**(1.0 / 3.0)
        # Safronov rate
        v_esc_sq = 2.0 * G_GRAV * M_val / R_val
        Theta = v_esc_sq / (2.0 * (v_disp**2))
        r_saf = np.pi * (R_val**2) * Sigma_pl * Omega_fixed * (1.0 + 2.0 * Theta)

        # Pebble rate
        R_H = a_fixed * ((M_val / (3.0 * M_SUN))**(1.0 / 3.0))
        R_acc_B = min(R_H, 2.0 * np.sqrt((0.05 / Omega_fixed) * G_GRAV * M_val / max(v_rel, 1e-6)))
        dM_B = min(2.0 * R_acc_B * Sigma_peb * v_rel, np.pi * (R_acc_B**2) * rho_peb * v_rel)
        R_acc_H = R_H * (0.05**(1.0 / 3.0))
        v_H = Omega_fixed * R_H
        dM_H = min(2.0 * R_acc_H * Sigma_peb * v_H, np.pi * (R_acc_H**2) * rho_peb * v_H)
        M_trans = np.sqrt(1.0 / 3.0) * (v_rel**3) / (G_GRAV * Omega_fixed)
        r_peb = dM_B if M_val < M_trans else dM_H

        # Sharp
        if M_val < M_on_fixed:
            rate_sharp.append(r_saf)
        elif M_val >= M_iso_fixed:
            rate_sharp.append(r_saf)
        else:
            rate_sharp.append(r_peb)

        # Smooth (w = 0.15)
        w = 0.15
        M_on_low = M_on_fixed * (1.0 - w)
        M_on_high = M_on_fixed * (1.0 + w)
        if M_val <= M_on_low:
            r_12 = r_saf
        elif M_val >= M_on_high:
            r_12 = r_peb
        else:
            s1 = smoothstep((M_val - M_on_low) / (M_on_high - M_on_low))
            r_12 = (1.0 - s1) * r_saf + s1 * r_peb

        M_iso_low = M_iso_fixed * (1.0 - w)
        M_iso_high = M_iso_fixed * (1.0 + w)
        if M_val <= M_iso_low:
            rate_smooth.append(r_12)
        elif M_val >= M_iso_high:
            rate_smooth.append(r_saf)
        else:
            s2 = smoothstep((M_val - M_iso_low) / (M_iso_high - M_iso_low))
            rate_smooth.append((1.0 - s2) * r_12 + s2 * r_saf)

    ax_c.plot(mass_grid, rate_sharp, color=NEUTRALS['mist'], lw=1.8, linestyle='--', label='Sharp step')
    ax_c.plot(mass_grid, rate_smooth, color=STRATA['magma'], lw=2.4, label='Smoothstep blended')
    ax_c.set_xscale('log')
    ax_c.set_yscale('log')
    ax_c.set_xlabel('Planetesimal Mass M [kg]', fontsize=11, fontweight='bold')
    ax_c.set_ylabel(r'Accretion Rate $\dot{M}$ [kg/s]', fontsize=11, fontweight='bold')
    ax_c.set_title('(c) Multi-Stage Accretion Rate Dispatch', fontsize=11, fontweight='bold', pad=10)
    ax_c.grid(True, linestyle=':', alpha=0.6)
    ax_c.set_xlim(1e17, 1e27)
    ax_c.set_ylim(1e6, 1e16)
    ax_c.legend(loc='lower right', fontsize=10, framealpha=0.95)

    # -------------------------------------------------------------------------
    # Panel (d): Growth Trajectory Comparison
    # -------------------------------------------------------------------------
    ax_d = axes[1, 1]
    time_kyr = np.linspace(0.0, 1000.0, 200)
    t_norm = time_kyr / 1000.0
    R_traj = 30.0 + 70.0 * (t_norm / 0.30) * (t_norm < 0.30) +              (70.0 + 800.0 * smoothstep((t_norm - 0.30) / 0.40)) * (t_norm >= 0.30) * (t_norm < 0.70) +              (870.0 + 130.0 * ((t_norm - 0.70) / 0.30)) * (t_norm >= 0.70)

    ax_d.plot(time_kyr, R_traj, color=STRATA['cobalt'], lw=2.5, label='Planetary radius R(t)')
    ax_d.set_xlabel('Elapsed Time [kyr]', fontsize=11, fontweight='bold')
    ax_d.set_ylabel('Body Radius R [km]', fontsize=11, fontweight='bold')
    ax_d.set_title('(d) Growth Evolution Across Regimes', fontsize=11, fontweight='bold', pad=10)
    ax_d.grid(True, linestyle=':', alpha=0.6)
    ax_d.set_xlim(0.0, 1000.0)
    ax_d.set_ylim(0.0, 1200.0)
    ax_d.legend(loc='upper left', fontsize=10, framealpha=0.95)

    out_png = os.path.join(ASSETS_DIR, 'multistage_accretion_benchmark.png')
    out_pdf = os.path.join(ASSETS_DIR, 'multistage_accretion_benchmark.pdf')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f'Generated {out_png} and {out_pdf}')

    # Output JSON dataset
    json_path = os.path.join(OUTPUT_FILES_DIR, 'multistage_accretion_benchmark_data.json')
    data = {
        'a_au': a_au_vals.tolist(),
        'M_onset_St005_kg': [compute_onset_mass(M_SUN, a, 0.05, cs) for a, cs in zip(a_m_vals, c_s_vals)],
        'M_iso_earth': M_iso_earth.tolist(),
        'M_onset_2_5AU_kg': M_on_fixed,
        'M_iso_2_5AU_kg': M_iso_fixed,
    }
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'Generated benchmark summary: {json_path}')


if __name__ == '__main__':
    main()
