#!/usr/bin/env python3
"""
Generate verification benchmark figure for temperature-dependent fluid viscosity in Erebus.jl.

Produces docs/src/assets/fluid_viscosity_temperature.png (300 DPI) and .svg.
Compares Arrhenius parametrization against experimental liquid water viscosity (IAPWS standard).
"""

import os
import numpy as np
import matplotlib.pyplot as plt

T_melt = 273.0        # K
T0 = 293.15           # K (20 °C)
eta0 = 1.0e-3         # Pa s
eta_ice = 1.0e12      # Pa s
R_gas = 8.314462618   # J/(mol K)

# IAPWS / NIST experimental data for saturated liquid water
# [T (K), eta (Pa s)]
nist_T = np.array([273.15, 293.15, 323.15, 373.15, 473.15, 573.15])
nist_eta = np.array([1.792e-3, 1.002e-3, 5.47e-4, 2.82e-4, 1.35e-4, 8.6e-5])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# -------------------------------------------------------------
# Panel (a): Viscosity vs Temperature (Semi-log)
# -------------------------------------------------------------
T_liquid = np.linspace(T_melt, 650, 300)

Ea_values = [15.0e3, 10.0e3, 20.0e3]
Ea_labels = [r'$E_a = 15\ \mathrm{kJ/mol}$ (Erebus default)',
             r'$E_a = 10\ \mathrm{kJ/mol}$',
             r'$E_a = 20\ \mathrm{kJ/mol}$']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for Ea, label, color in zip(Ea_values, Ea_labels, colors):
    log_ratio = (Ea / R_gas) * (1.0 / T_liquid - 1.0 / T0)
    eta = eta0 * np.exp(log_ratio)
    eta = np.clip(eta, 1.0e-5, 1.0e12)
    ax1.plot(T_liquid, eta, label=label, color=color, lw=2.0)

# Constant mode
ax1.plot(T_liquid, np.full_like(T_liquid, eta0), 'k--', lw=1.5, label=r'Constant mode ($\eta_0 = 10^{-3}\ \mathrm{Pa\ s}$)')

# Experimental NIST data points
ax1.plot(nist_T, nist_eta, 'ro', markersize=6, label='IAPWS / NIST water data', zorder=5)

ax1.set_yscale('log')
ax1.set_xlabel('Temperature $T$ [K]', fontsize=11)
ax1.set_ylabel(r'Dynamic Fluid Viscosity $\eta_f(T)$ [Pa s]', fontsize=11)
ax1.set_title('(a) Fluid Viscosity vs Temperature', fontsize=11, fontweight='bold')
ax1.set_xlim(270, 650)
ax1.set_ylim(5.0e-5, 3.0e-3)
ax1.grid(True, which='both', alpha=0.3)
ax1.legend(loc='upper right', fontsize=8.5, framealpha=0.9)

# -------------------------------------------------------------
# Panel (b): Hydrothermal Darcy Mobility Enhancement
# -------------------------------------------------------------
# Relative mobility ratio: (k_phi / eta_f) / (k_phi / eta0) = eta0 / eta_f
for Ea, label, color in zip(Ea_values, Ea_labels, colors):
    log_ratio = (Ea / R_gas) * (1.0 / T_liquid - 1.0 / T0)
    eta = eta0 * np.exp(log_ratio)
    mobility_enhancement = eta0 / np.clip(eta, 1.0e-5, 1.0e12)
    ax2.plot(T_liquid, mobility_enhancement, label=label, color=color, lw=2.0)

ax2.plot(nist_T, eta0 / nist_eta, 'ro', markersize=6, label='IAPWS / NIST data', zorder=5)
ax2.axhline(1.0, color='black', ls='--', lw=1.0, alpha=0.5)

ax2.set_xlabel('Temperature $T$ [K]', fontsize=11)
ax2.set_ylabel(r'Mobility Ratio $\eta_0 / \eta_f(T)$ [-]', fontsize=11)
ax2.set_title(r'(b) Darcy Flux Enhancement Factor', fontsize=11, fontweight='bold')
ax2.set_xlim(270, 650)
ax2.set_ylim(0, 45)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper left', fontsize=8.5, framealpha=0.9)

plt.tight_layout()

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'src', 'assets')
os.makedirs(out_dir, exist_ok=True)

png_path = os.path.join(out_dir, 'fluid_viscosity_temperature.png')
svg_path = os.path.join(out_dir, 'fluid_viscosity_temperature.svg')

plt.savefig(png_path, dpi=300)
plt.savefig(svg_path)
plt.close()

print(f"Generated {png_path} and {svg_path}")
