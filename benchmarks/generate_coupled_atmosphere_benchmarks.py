#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for coupled 1D atmosphere,
disk gas envelope capture, Ormel et al. (2015) recycling limit, Guillot (2010)
semi-grey radiative equilibrium, and Zahnle & Kasting (1986) crossover escape in Erebus.jl.
"""

import json
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import simpson

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

# Physical constants
G = 6.6743e-11 # m^3/(kg s^2)
KB = 1.380649e-23 # J/K
SIGMA_SB = 5.670374419e-8 # W/(m^2 K^4)
AU_METERS = 1.495978707e11 # m
M_SUN_KG = 1.98847e30 # kg
AMU_KG = 1.66053906660e-27 # kg

SPECIES_MASSES = {
    'H2': 2.016 * AMU_KG,
    'CH4': 16.04 * AMU_KG,
    'NH3': 17.03 * AMU_KG,
    'H2O': 18.015 * AMU_KG,
    'CO': 28.01 * AMU_KG,
    'N2': 28.0134 * AMU_KG,
    'CO2': 44.01 * AMU_KG,
    'SO2': 64.066 * AMU_KG,
}


def compute_guillot_profile(tau_arr, T_int, T_eqm, gamma):
    """Guillot (2010) semi-grey temperature profile."""
    sqrt3 = np.sqrt(3.0)
    inv_gam_sqrt3 = 1.0 / (gamma * sqrt3)
    gam_over_sqrt3 = gamma / sqrt3
    exp_term = np.exp(-gamma * tau_arr * sqrt3)
    bracket = 2.0 / 3.0 + inv_gam_sqrt3 + (gam_over_sqrt3 - inv_gam_sqrt3) * exp_term
    T4 = (3.0 / 4.0) * (T_int**4) * (tau_arr + 2.0 / 3.0) + (3.0 / 4.0) * (T_eqm**4) * bracket
    return np.maximum(0.0, T4)**0.25


def compute_envelope_mass(R_planet, rho_bulk=3000.0, rho_disk=1.0e-9, c_s=300.0, a_orb=AU_METERS, f_rec=0.10):
    """Compute bound isothermal envelope mass with Ormel recycling limit."""
    M_planet = (4.0 * np.pi / 3.0) * (R_planet**3) * rho_bulk
    R_Bondi = G * M_planet / (c_s**2)
    R_Hill = a_orb * (M_planet / (3.0 * M_SUN_KG))**(1.0 / 3.0)
    R_cap = min(R_Bondi, R_Hill)

    if R_cap <= R_planet:
        return 0.0, 0.0

    r_grid = np.linspace(R_planet, R_cap, 128)
    psi = (G * M_planet / (c_s**2)) * (1.0 / r_grid - 1.0 / R_cap)
    rho_env = rho_disk * np.exp(np.minimum(50.0, np.maximum(0.0, psi)))
    integrand = 4.0 * np.pi * (r_grid**2) * rho_env
    M_iso = simpson(integrand, x=r_grid)

    M_rec = f_rec * (4.0 * np.pi / 3.0) * (R_cap**3) * rho_disk
    return min(M_iso, M_rec), M_rec


def main():
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(hspace=0.32, wspace=0.30)

    # -------------------------------------------------------------------------
    # Panel (a): Guillot (2010) Semi-Grey Radiative Equilibrium Profile
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    tau_vals = np.logspace(-2, 2, 200)
    T_int = 100.0 # K
    T_eqm = 250.0 # K
    gamma_list = [0.01, 0.1, 1.0, 5.0]
    colors_gamma = [STRATA['cobalt'], STRATA['plum'], STRATA['amber'], STRATA['gold']]

    for gamma, col in zip(gamma_list, colors_gamma):
        T_prof = compute_guillot_profile(tau_vals, T_int, T_eqm, gamma)
        ax.plot(T_prof, tau_vals, label=rf"$\gamma = {gamma}$", color=col, lw=2.2)

    # Skin temperature mark (factor 2^-0.25)
    T_skin_limit = (0.5 * (T_int**4) + 0.5 * (1.0 + 0.1) * (T_eqm**4))**0.25
    ax.axvline(T_eqm, color=NEUTRALS['graphite'], ls='--', lw=1.2, label=r"$T_{\mathrm{eqm}} = 250\,$K")

    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlim(180, 550)
    ax.set_ylim(100, 0.01)
    ax.set_xlabel("Equilibrium Temperature $T$ [K]", fontsize=11, fontweight='bold')
    ax.set_ylabel(r"Optical Depth $\tau_{\mathrm{LW}}$", fontsize=11, fontweight='bold')
    ax.set_title("(a) Guillot (2010) Radiative Equilibrium", fontsize=12, fontweight='bold', pad=10)
    ax.grid(True)
    ax.legend(loc='lower left', fontsize=9.5)

    # -------------------------------------------------------------------------
    # Panel (b): Disk Envelope Capture & Recycling Limit (Ormel et al. 2015)
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    R_planets = np.linspace(100e3, 2000e3, 100) # 100 km to 2000 km
    M_envs = []
    M_recs = []

    for Rp in R_planets:
        M_e, M_r = compute_envelope_mass(Rp)
        M_envs.append(M_e)
        M_recs.append(M_r)

    M_envs = np.array(M_envs)
    M_recs = np.array(M_recs)
    R_km = R_planets / 1e3

    ax.plot(R_km, M_recs, label=r"Ormel et al. (2015) Recycling Limit", color=STRATA['magma'], lw=2.2, ls='--')
    ax.plot(R_km, M_envs, label=r"Captured Bound Envelope $M_{\mathrm{env}}^*$", color=STRATA['cobalt'], lw=2.4)

    # Annotations for transition points
    idx_500 = np.argmin(np.abs(R_km - 500.0))
    idx_1500 = np.argmin(np.abs(R_km - 1500.0))
    ax.scatter([R_km[idx_500], R_km[idx_1500]], [M_envs[idx_500], M_envs[idx_1500]],
               color=STRATA['amber'], s=45, zorder=5)
    ax.text(R_km[idx_500] + 50, M_envs[idx_500] * 2.0, "Vesta-scale", fontsize=9.5, color=NEUTRALS['graphite'])
    ax.text(R_km[idx_1500] - 350, M_envs[idx_1500] * 1.5, "Lunar embryo", fontsize=9.5, color=NEUTRALS['graphite'])

    ax.set_yscale('log')
    ax.set_ylim(1e12, 1e20)
    ax.set_xlim(100, 2000)
    ax.set_xlabel(r"Planetesimal Radius $R_{\mathrm{planet}}$ [km]", fontsize=11, fontweight='bold')
    ax.set_ylabel(r"Envelope Mass $M_{\mathrm{env}}$ [kg]", fontsize=11, fontweight='bold')
    ax.set_title("(b) Disk Envelope Capture & Recycling Limit", fontsize=12, fontweight='bold', pad=10)
    ax.grid(True)
    ax.legend(loc='upper left', fontsize=9.5)

    # -------------------------------------------------------------------------
    # Panel (c): Thermal Blanketing & Effective Heat Transfer Coefficient
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    tau_arr = np.linspace(0.0, 20.0, 150)
    T_surf = 350.0
    T_amb = 200.0
    emissivity = 0.90
    h_bare = 4.0 * emissivity * SIGMA_SB * (((T_surf + T_amb) / 2.0)**3)
    h_eff = h_bare / (1.0 + 0.75 * tau_arr)

    ax.plot(tau_arr, h_eff, label=r"Effective Radiative HTC $h_{\mathrm{rad,eff}}$", color=STRATA['amber'], lw=2.4)
    ax.axhline(h_bare, color=NEUTRALS['graphite'], ls=':', lw=1.2, label=rf"Bare Surface Limit ($h = {h_bare:.2f}\,$W/m$^2$K)")

    # Secondary y-axis for relative thermal insulation
    ax_sec = ax.twinx()
    insulation_factor = 1.0 + 0.75 * tau_arr
    ax_sec.plot(tau_arr, insulation_factor, color=STRATA['plum'], lw=1.8, ls='-.', label="Insulation Factor $(1 + 3\\tau/4)$")
    ax_sec.set_ylabel(r"Thermal Blanket Impedance $R_{\mathrm{blanket}}$", fontsize=11, fontweight='bold', color=STRATA['plum'])
    ax_sec.tick_params(axis='y', labelcolor=STRATA['plum'])

    ax.set_xlim(0, 20)
    ax.set_ylim(0, h_bare * 1.1)
    ax.set_xlabel(r"Atmospheric Longwave Optical Depth $\tau_{\mathrm{LW}}$", fontsize=11, fontweight='bold')
    ax.set_ylabel(r"$h_{\mathrm{rad,eff}}$ [W/(m$^2$ K)]", fontsize=11, fontweight='bold')
    ax.set_title("(c) Greenhouse Thermal Blanketing", fontsize=12, fontweight='bold', pad=10)
    ax.grid(True)

    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax_sec.get_legend_handles_labels()
    ax.legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right', fontsize=9.0)

    # -------------------------------------------------------------------------
    # Panel (d): Zahnle-Kasting (1986) Multi-Species Crossover Escape
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    Phi_H2 = np.logspace(16, 21, 200) # molecules / (m^2 s)
    T_exo = 250.0 # K
    M_planet = 1.0e23 # kg
    R_planet = 1.5e6 # m
    g_surf = G * M_planet / (R_planet**2)
    b_diff = 1.0e21
    m_carrier = SPECIES_MASSES['H2']

    # Crossover mass
    m_c = m_carrier + (KB * T_exo * Phi_H2) / (b_diff * g_surf * 1.0)
    m_c_amu = m_c / AMU_KG

    # Fractionation efficiencies for different volatile species
    species_plot = [
        ('CH4', STRATA['cobalt'], '-'),
        ('H2O', STRATA['plum'], '-'),
        ('CO', STRATA['amber'], '-'),
        ('CO2', STRATA['magma'], '-'),
        ('SO2', STRATA['ink'], '-')
    ]

    for sp, col, ls in species_plot:
        m_sp = SPECIES_MASSES[sp]
        x_drag = np.maximum(0.0, 1.0 - (m_sp - m_carrier) / np.maximum(1e-30, m_c - m_carrier))
        ax.plot(Phi_H2, x_drag, label=rf"{sp} ({m_sp/AMU_KG:.1f} Da)", color=col, lw=2.0, ls=ls)

    ax.set_xscale('log')
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(1e16, 1e21)
    ax.set_xlabel(r"Carrier Hydrogen Escape Flux $\Phi_{\mathrm{H}_2}$ [$\mathrm{m}^{-2}\mathrm{s}^{-1}$]", fontsize=11, fontweight='bold')
    ax.set_ylabel(r"Crossover Drag Efficiency $x_j$", fontsize=11, fontweight='bold')
    ax.set_title("(d) Zahnle-Kasting Hydrodynamic Crossover Drag", fontsize=12, fontweight='bold', pad=10)
    ax.grid(True)
    ax.legend(loc='lower right', fontsize=9.0)

    # Save figure
    png_path = os.path.join(ASSETS_DIR, "coupled_atmosphere_benchmark.png")
    pdf_path = os.path.join(ASSETS_DIR, "coupled_atmosphere_benchmark.pdf")
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Generated benchmark figure:\n  {png_path}\n  {pdf_path}")

    # Output validation dataset
    val_data = {
        "guillot_profiles": {
            f"gamma_{gamma}": compute_guillot_profile(np.array([0.01, 0.1, 1.0, 10.0]), T_int, T_eqm, gamma).tolist()
            for gamma in gamma_list
        },
        "envelope_capture": {
            "R_km": [200.0, 500.0, 1000.0, 1500.0, 2000.0],
            "M_env_bound": [compute_envelope_mass(r * 1e3)[0] for r in [200.0, 500.0, 1000.0, 1500.0, 2000.0]],
            "M_rec_limit": [compute_envelope_mass(r * 1e3)[1] for r in [200.0, 500.0, 1000.0, 1500.0, 2000.0]],
        },
        "crossover_escape": {
            "Phi_H2": [1e17, 1e18, 1e19, 1e20],
            "m_c_amu": ((m_carrier + (KB * T_exo * np.array([1e17, 1e18, 1e19, 1e20])) / (b_diff * g_surf)) / AMU_KG).tolist()
        }
    }
    json_path = os.path.join(OUTPUT_FILES_DIR, "coupled_atmosphere_validation.json")
    with open(json_path, "w") as f:
        json.dump(val_data, f, indent=2)
    print(f"Generated validation dataset:\n  {json_path}")


if __name__ == "__main__":
    main()
