#!/usr/bin/env python3
"""
Generate benchmark figure and validation dataset for hydrothermal subgrid convection
parameterization in Erebus.jl.
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

# -----------------------------------------------------------------------------
# Physics functions matching src/physics.jl
# -----------------------------------------------------------------------------
def compute_porous_rayleigh_darcy(rho_f, cp_f, g, alpha_f, K, dT, H, mu_f, k_cond):
    if dT <= 0.0 or K <= 0.0 or H <= 0.0 or g <= 0.0:
        return 0.0
    num = (rho_f**2) * cp_f * g * alpha_f * K * dT * H
    den = mu_f * k_cond
    return num / den

def compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, dT, H, mu_f, k_f):
    if dT <= 0.0 or H <= 0.0 or g <= 0.0:
        return 0.0
    num = (rho_f**2) * cp_f * g * alpha_f * dT * (H**3)
    den = mu_f * k_f
    return num / den

def compute_hydrothermal_nusselt(Ra_m, Ra, phi, phi_start=0.30, phi_end=0.70,
                                 Ra_m_crit=4.0*np.pi**2, Ra_crit=1100.0,
                                 c_porous=1.0, c_free=0.088):
    Nu_porous = 1.0 if Ra_m <= Ra_m_crit else 1.0 + c_porous * (Ra_m / Ra_m_crit - 1.0)
    Nu_free = 1.0 if Ra <= Ra_crit else np.maximum(1.0, c_free * np.cbrt(Ra))

    if phi <= phi_start:
        return Nu_porous
    elif phi >= phi_end:
        return Nu_free
    else:
        xi = (phi - phi_start) / (phi_end - phi_start)
        w_phi = xi * xi * (3.0 - 2.0 * xi)
        log_Nu = (1.0 - w_phi) * np.log10(Nu_porous) + w_phi * np.log10(Nu_free)
        return 10.0**log_Nu

def compute_effective_hydrothermal_conductivity(k_cond, Nu, Pe_cell=0.0, Pe_crit=2.0,
                                               resolution_weighting=True, picard_damping=1.0,
                                               k_floor=1.0e-3, k_cutoff=1.0e6, k_prev=None):
    if k_prev is None:
        k_prev = k_cond
    k_raw = Nu * k_cond
    k_target = np.maximum(k_cond, np.clip(k_raw, k_floor, k_cutoff))

    if resolution_weighting and Pe_cell > 0.0:
        w_res = np.clip(Pe_cell / Pe_crit, 0.0, 1.0)
        k_res = (1.0 - w_res) * k_target + w_res * k_cond
    else:
        k_res = k_target

    k_damped = (1.0 - picard_damping) * k_prev + picard_damping * k_res
    return np.maximum(k_cond, np.clip(k_damped, k_floor, k_cutoff))


def main():
    # Standard planetesimal aquifer baseline properties
    rho_f = 1000.0
    cp_f = 4184.0
    g = 0.5
    alpha_f = 2.0e-4
    H = 10000.0
    mu_f = 1.0e-3
    k_cond = 2.5
    k_f = 0.6
    Ra_m_crit = 4.0 * np.pi**2
    Ra_crit = 1100.0

    fig, axs = plt.subplots(2, 2, figsize=(13, 10))

    # -------------------------------------------------------------------------
    # Panel (a): Ra_m vs Permeability K
    # -------------------------------------------------------------------------
    ax_a = axs[0, 0]
    K_vals = np.logspace(-15, -11, 200)
    dT_targets = [10.0, 25.0, 50.0, 100.0]
    colors_a = [STRATA['cobalt'], STRATA['gold'], STRATA['amber'], STRATA['magma']]

    for dT_val, col in zip(dT_targets, colors_a):
        Ra_m_arr = [compute_porous_rayleigh_darcy(rho_f, cp_f, g, alpha_f, K, dT_val, H, mu_f, k_cond)
                    for K in K_vals]
        ax_a.loglog(K_vals, Ra_m_arr, label=f"ΔT = {int(dT_val)} K", color=col, lw=2.0)

    ax_a.axhline(Ra_m_crit, color=NEUTRALS['graphite'], linestyle='--', lw=1.5,
                 label=r"Onset threshold $Ra_{m,\mathrm{crit}} = 4\pi^2$")
    ax_a.set_xlabel(r"Aquifer permeability $K$ [$\mathrm{m}^2$]", fontsize=11, fontweight='medium')
    ax_a.set_ylabel(r"Porous Rayleigh-Darcy number $Ra_m$ [-]", fontsize=11, fontweight='medium')
    ax_a.set_title("(a) Porous Convection Onset & Scaling", fontsize=12, fontweight='bold', pad=10)
    ax_a.grid(True)
    ax_a.legend(loc='lower right', fontsize=9)

    # -------------------------------------------------------------------------
    # Panel (b): Blended Nusselt Number vs Porosity phi
    # -------------------------------------------------------------------------
    ax_b = axs[0, 1]
    phi_arr = np.linspace(0.0, 1.0, 300)
    permeabilities = [1.0e-13, 5.0e-13, 2.0e-12]
    styles_b = ['-', '--', '-.']

    dT_fixed = 50.0
    Ra_free_fixed = compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, dT_fixed, H, mu_f, k_f)

    for K_perm, st, col in zip(permeabilities, styles_b, [STRATA['cobalt'], STRATA['amber'], STRATA['magma']]):
        Ra_m_fixed = compute_porous_rayleigh_darcy(rho_f, cp_f, g, alpha_f, K_perm, dT_fixed, H, mu_f, k_cond)
        Nu_vals = [compute_hydrothermal_nusselt(Ra_m_fixed, Ra_free_fixed, p) for p in phi_arr]
        ax_b.semilogy(phi_arr, Nu_vals, label=f"$K = {K_perm:.0e}$ m²", color=col, linestyle=st, lw=2.0)

    # Transition zone highlight
    ax_b.axvspan(0.30, 0.70, color=NEUTRALS['bone'], alpha=0.35,
                 label="Smoothstep transition zone $[0.3, 0.7]$")
    ax_b.axvline(0.30, color=STRATA['plum'], linestyle=':', lw=1.2)
    ax_b.axvline(0.70, color=STRATA['plum'], linestyle=':', lw=1.2)

    ax_b.text(0.15, 2.5, "Porous\nDarcy", color=STRATA['cobalt'], fontsize=10, ha='center', fontweight='bold')
    ax_b.text(0.50, 15.0, "Smoothstep\nBlending", color=STRATA['plum'], fontsize=10, ha='center', fontweight='bold')
    ax_b.text(0.85, 200.0, "Free Fluid\nBoundary Layer", color=STRATA['amber'], fontsize=10, ha='center', fontweight='bold')

    ax_b.set_xlabel(r"Porosity $\phi$ [-]", fontsize=11, fontweight='medium')
    ax_b.set_ylabel(r"Effective Nusselt number $Nu$ [-]", fontsize=11, fontweight='medium')
    ax_b.set_title("(b) Nusselt Scaling Across Porosity Regimes", fontsize=12, fontweight='bold', pad=10)
    ax_b.grid(True)
    ax_b.legend(loc='lower right', fontsize=9)

    # -------------------------------------------------------------------------
    # Panel (c): Effective Conductivity Enhancement Across phi - dT
    # -------------------------------------------------------------------------
    ax_c = axs[1, 0]
    phi_grid = np.linspace(0.05, 0.95, 100)
    dT_grid = np.linspace(5.0, 100.0, 100)
    PHI, DT = np.meshgrid(phi_grid, dT_grid)
    KEFF_RATIO = np.zeros_like(PHI)

    K_ref = 5.0e-13
    for i in range(len(dT_grid)):
        for j in range(len(phi_grid)):
            dt_cur = DT[i, j]
            phi_cur = PHI[i, j]
            ram_cur = compute_porous_rayleigh_darcy(rho_f, cp_f, g, alpha_f, K_ref, dt_cur, H, mu_f, k_cond)
            ra_cur = compute_free_fluid_rayleigh(rho_f, cp_f, g, alpha_f, dt_cur, H, mu_f, k_f)
            nu_cur = compute_hydrothermal_nusselt(ram_cur, ra_cur, phi_cur)
            # Quadratic temperature ramp across dT_min matching apply_hydrothermal_convection_closure
            w_T = np.clip(dt_cur / 5.0, 0.0, 1.0)**2
            nu_ramped = 10.0**(w_T * np.log10(nu_cur)) if nu_cur > 1.0 else 1.0
            keff_cur = compute_effective_hydrothermal_conductivity(k_cond, nu_ramped)
            KEFF_RATIO[i, j] = keff_cur / k_cond

    levels = np.logspace(0, 3, 10)
    cs = ax_c.contourf(PHI, DT, KEFF_RATIO, levels=levels, norm=mpl.colors.LogNorm(vmin=1.0, vmax=1000.0),
                       cmap='magma', alpha=0.85)
    cbar = fig.colorbar(cs, ax=ax_c)
    cbar.set_label(r"Enhancement factor $k_{\mathrm{eff}} / k_{\mathrm{cond}}$ [-]", fontsize=10)

    ax_c.axvline(0.30, color='#FFFFFF', linestyle='--', lw=1.2)
    ax_c.axvline(0.70, color='#FFFFFF', linestyle='--', lw=1.2)
    ax_c.set_xlabel(r"Porosity $\phi$ [-]", fontsize=11, fontweight='medium')
    ax_c.set_ylabel(r"Temperature contrast $\Delta T$ [K]", fontsize=11, fontweight='medium')
    ax_c.set_title(r"(c) Conductivity Ratio $k_{\mathrm{eff}} / k_{\mathrm{cond}}$ ($K=5\times 10^{-13}\,\mathrm{m}^2$)",
                   fontsize=12, fontweight='bold', pad=10)

    # -------------------------------------------------------------------------
    # Panel (d): Cell-Péclet Resolution Damping
    # -------------------------------------------------------------------------
    ax_d = axs[1, 1]
    Pe_vals = np.linspace(0.0, 3.0, 200)
    Nu_targets = [2.0, 5.0, 10.0, 25.0, 50.0]
    colors_d = [STRATA['cobalt'], STRATA['gold'], STRATA['amber'], STRATA['magma'], STRATA['plum']]

    for Nu_t, col in zip(Nu_targets, colors_d):
        keff_ratio_pe = [compute_effective_hydrothermal_conductivity(k_cond, Nu_t, Pe_cell=pe, Pe_crit=2.0) / k_cond
                         for pe in Pe_vals]
        ax_d.plot(Pe_vals, keff_ratio_pe, label=f"$Nu = {int(Nu_t)}$", color=col, lw=2.0)

    ax_d.axvline(2.0, color=NEUTRALS['graphite'], linestyle='--', lw=1.5,
                 label=r"Resolved grid threshold $Pe_{\mathrm{crit}} = 2.0$")
    ax_d.axhline(1.0, color=NEUTRALS['mist'], linestyle=':', lw=1.2)
    ax_d.set_xlabel(r"Cell-Péclet number $Pe_{\mathrm{cell}} = v_{\mathrm{Darcy}} \Delta x / \kappa_f$ [-]",
                    fontsize=11, fontweight='medium')
    ax_d.set_ylabel(r"Realized enhancement factor $k_{\mathrm{eff}} / k_{\mathrm{cond}}$ [-]",
                    fontsize=11, fontweight='medium')
    ax_d.set_title("(d) Grid-Resolution Damping (Prevent Double Counting)",
                   fontsize=12, fontweight='bold', pad=10)
    ax_d.grid(True)
    ax_d.legend(loc='upper right', fontsize=9)

    plt.tight_layout()

    out_png = os.path.join(ASSETS_DIR, "hydrothermal_convection_benchmark.png")
    out_pdf = os.path.join(ASSETS_DIR, "hydrothermal_convection_benchmark.pdf")
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Rendered benchmark plots to:\n  {out_png}\n  {out_pdf}")

if __name__ == "__main__":
    main()
