#!/usr/bin/env python3
"""
Generate verification benchmark figure for Darcy thermal buoyancy in Erebus.jl.

Produces docs/src/assets/darcy_buoyancy_verification.png (300 DPI) and .svg.
Reconciled with exact code constants:
  tmfluidphase = 273.0 K
  alphafluidm = 5.0e-5 K⁻¹ (baseline)
  rhoH2Ofluid = 1000.0 kg/m³
  rhoH2Ofluidice = 917.0 kg/m³
"""

import os
import numpy as np
import matplotlib.pyplot as plt

rho0_water = 1000.0  # kg/m³
rho0_ice = 917.0     # kg/m³
T_melt = 273.0       # K (tmfluidphase in Erebus constants)

alphas = [5.0e-5, 2.0e-4, 5.0e-4]  # 1/K
alpha_labels = [r'$\alpha_f = 5\times 10^{-5}\ \mathrm{K}^{-1}$ (code baseline)',
                r'$\alpha_f = 2\times 10^{-4}\ \mathrm{K}^{-1}$ (ambient water)',
                r'$\alpha_f = 5\times 10^{-4}\ \mathrm{K}^{-1}$ (hydrothermal)']
colors = ['#1f77b4', '#d62728', '#2ca02c']

eta_f = 1.0e-3  # Pa s
g = 9.81        # m/s²
k_phis = [1.0e-12, 1.0e-13, 1.0e-14]  # m²

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# -------------------------------------------------------------
# Panel (a): Density vs Temperature with Ice-Water Phase Transition
# -------------------------------------------------------------
T_sub = np.linspace(240, T_melt, 100)
T_sup = np.linspace(T_melt, 700, 300)

# Sub-freezing ice density
ax1.plot(T_sub, np.full_like(T_sub, rho0_ice), color='black', lw=2.0, label=r'Ice phase ($\rho_{\mathrm{ice}} = 917\ \mathrm{kg/m}^3$)')
ax1.plot([T_melt, T_melt], [rho0_ice, rho0_water], color='black', ls=':', lw=1.5)

for alpha, label, color in zip(alphas, alpha_labels, colors):
    # Liquid water thermal expansion: rho_f(T) = rho0 * max(0.1, 1 - alpha*(T - T_melt))
    rho_liquid = rho0_water * np.maximum(0.1, 1.0 - alpha * (T_sup - T_melt))
    ax1.plot(T_sup, rho_liquid, label=label, color=color, lw=2.0)

ax1.axvline(T_melt, color='gray', ls='--', lw=1.0, alpha=0.7, label=r'$T_{\mathrm{melt}} = 273.0\ \mathrm{K}$')
ax1.axhline(rho0_water, color='black', ls=':', lw=1.0, alpha=0.3)

ax1.set_xlabel('Temperature $T$ [K]', fontsize=11)
ax1.set_ylabel(r'Fluid Density $\rho_f(T)$ [kg m$^{-3}$]', fontsize=11)
ax1.set_title('(a) Pore Fluid Thermal Expansion Equation of State', fontsize=11, fontweight='bold')
ax1.set_xlim(240, 700)
ax1.set_ylim(750, 1050)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='lower left', fontsize=8.5, framealpha=0.9)

# -------------------------------------------------------------
# Panel (b): Vertical Buoyant Darcy Velocity vs Delta T
# -------------------------------------------------------------
delta_T = np.linspace(0, 200, 200)
alpha_baseline = 5.0e-5

k_styles = ['-', '--', '-.']
k_labels = [r'$k_\phi = 10^{-12}\ \mathrm{m}^2$',
            r'$k_\phi = 10^{-13}\ \mathrm{m}^2$',
            r'$k_\phi = 10^{-14}\ \mathrm{m}^2$']

for k, k_style, k_label in zip(k_phis, k_styles, k_labels):
    # Buoyant Darcy velocity: q_y = (k_phi / eta_f) * rho_0 * alpha_f * Delta_T * g
    q_darcy = (k / eta_f) * (rho0_water * alpha_baseline * delta_T * g)  # m/s
    q_darcy_m_yr = q_darcy * (365.25 * 86400.0)                         # m/yr

    ax2.plot(delta_T, q_darcy_m_yr, label=k_label, color='#1f77b4', ls=k_style, lw=2.0)

ax2.set_xlabel(r'Temperature Contrast $\Delta T = T - T_{\mathrm{melt}}$ [K]', fontsize=11)
ax2.set_ylabel(r'Buoyant Darcy Discharge Velocity $|q_{yD}|$ [m yr$^{-1}$]', fontsize=11)
ax2.set_title(r'(b) Darcy Buoyancy Scaling ($\alpha_f = 5\times 10^{-5}\ \mathrm{K}^{-1}$)', fontsize=11, fontweight='bold')
ax2.set_xlim(0, 200)
ax2.set_ylim(bottom=0)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=9, framealpha=0.9)

plt.tight_layout()

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'src', 'assets')
os.makedirs(out_dir, exist_ok=True)

png_path = os.path.join(out_dir, 'darcy_buoyancy_verification.png')
svg_path = os.path.join(out_dir, 'darcy_buoyancy_verification.svg')

plt.savefig(png_path, dpi=300)
plt.savefig(svg_path)
plt.close()

print(f"Generated {png_path} and {svg_path}")
