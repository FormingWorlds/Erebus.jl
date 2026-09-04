#!/usr/bin/env python3
"""
Generate verification benchmark figure for dynamic hydrofracturing permeability enhancement in Erebus.jl.

Produces docs/src/assets/hydrofracture_verification.png (300 DPI) and .svg.
Displays effective permeability response across compressive, tensile, and hydrofracture regimes.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

sigma_t = 10.0e6      # 10 MPa rock tensile strength
kappa_frac = 1.0e3    # Default multiplier
kmax = 1.0e-9         # 10^-9 m² ceiling

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

# -------------------------------------------------------------
# Panel (a): Permeability vs Terzaghi Effective Stress
# -------------------------------------------------------------
# Peff = Pt - Pf in MPa (compressive > 0, tensile < 0)
Peff_arr_MPa = np.linspace(-50.0, 20.0, 500)
Peff_arr = Peff_arr_MPa * 1.0e6

k0_values = [1.0e-16, 1.0e-15, 1.0e-14]
k0_labels = [r'$k_0 = 10^{-16}\ \mathrm{m}^2$',
             r'$k_0 = 10^{-15}\ \mathrm{m}^2$',
             r'$k_0 = 10^{-14}\ \mathrm{m}^2$']
colors = ['#1f77b4', '#2ca02c', '#d62728']

for k0, label, color in zip(k0_values, k0_labels, colors):
    overpressure = np.maximum(0.0, -Peff_arr - sigma_t)
    norm_op = overpressure / sigma_t
    factor = 1.0 + kappa_frac * (norm_op ** 1.0)
    keff = np.clip(k0 * factor, k0, kmax)
    ax1.plot(Peff_arr_MPa, keff, label=label, color=color, lw=2.0)

# Mark failure boundary
ax1.axvline(-sigma_t * 1e-6, color='black', ls='--', lw=1.2, alpha=0.8, label=r'Tensile failure ($-\sigma_t = -10\ \mathrm{MPa}$)')
ax1.axvspan(-50, -sigma_t * 1e-6, alpha=0.1, color='red', label='Hydrofracture regime')
ax1.axvspan(-sigma_t * 1e-6, 20, alpha=0.08, color='green', label='Intact matrix regime')

ax1.set_yscale('log')
ax1.set_xlabel(r'Effective Stress $P_{\mathrm{eff}} = P_t - P_f$ [MPa]', fontsize=11)
ax1.set_ylabel(r'Effective Permeability $k_{\phi}^{\mathrm{eff}}$ [$\mathrm{m}^2$]', fontsize=11)
ax1.set_title('(a) Permeability vs Effective Stress', fontsize=11, fontweight='bold')
ax1.set_xlim(-50, 20)
ax1.set_ylim(1.0e-17, 1.0e-8)
ax1.grid(True, which='both', alpha=0.3)
ax1.legend(loc='lower left', fontsize=8.0, framealpha=0.9)

# -------------------------------------------------------------
# Panel (b): Permeability Enhancement vs Normalized Overpressure
# -------------------------------------------------------------
norm_op_range = np.linspace(0.0, 5.0, 300)

gamma_values = [0.5, 1.0, 2.0]
gamma_labels = [r'$\gamma = 0.5$ (sub-linear)',
                r'$\gamma = 1.0$ (linear, default)',
                r'$\gamma = 2.0$ (quadratic)']
gamma_colors = ['#9467bd', '#1f77b4', '#ff7f0e']

for g, label, color in zip(gamma_values, gamma_labels, gamma_colors):
    enhancement = 1.0 + kappa_frac * (norm_op_range ** g)
    ax2.plot(norm_op_range, enhancement, label=label, color=color, lw=2.0)

ax2.set_yscale('log')
ax2.set_xlabel(r'Normalized Overpressure $[P_f - (P_t + \sigma_t)] / \sigma_t$ [-]', fontsize=11)
ax2.set_ylabel(r'Enhancement Factor $k_{\mathrm{eff}} / k_0$ [-]', fontsize=11)
ax2.set_title(r'(b) Permeability Multiplier Scaling ($\kappa_{\mathrm{frac}} = 10^3$)', fontsize=11, fontweight='bold')
ax2.set_xlim(0, 5)
ax2.set_ylim(1, 2.0e4)
ax2.grid(True, which='both', alpha=0.3)
ax2.legend(loc='lower right', fontsize=8.5, framealpha=0.9)

plt.tight_layout()

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'src', 'assets')
os.makedirs(out_dir, exist_ok=True)

png_path = os.path.join(out_dir, 'hydrofracture_verification.png')
svg_path = os.path.join(out_dir, 'hydrofracture_verification.svg')

plt.savefig(png_path, dpi=300)
plt.savefig(svg_path)
plt.close()

print(f"Generated {png_path} and {svg_path}")
