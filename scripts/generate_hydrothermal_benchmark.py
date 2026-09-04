#!/usr/bin/env python3
"""
Generate publication-quality 2D hydrothermal circulation benchmark figure for Erebus.jl.

Produces docs/src/assets/hydrothermal_circulation_benchmark.png (300 DPI) and .svg.
Displays 4 panels:
  (a) 2D Temperature field T(x, y) [K] with planet and crust boundaries
  (b) Darcy fluid flux vectors q^D and magnitude [m/s]
  (c) Terzaghi effective stress P_eff = P_t - P_f [MPa] and hydrofracture zone
  (d) Radial profiles of temperature and pressures from center to surface
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Path resolution
data_path = os.path.join(os.path.dirname(__file__), '..', 'output_hydrothermal', 'benchmark_plot_data.json')
if not os.path.isfile(data_path):
    raise FileNotFoundError(f"Exported data file not found: {data_path}. Run export_hydrothermal_data.jl first.")

with open(data_path, 'r') as f:
    data = json.load(f)

x = np.array(data['x']) * 1e-3  # [km]
y = np.array(data['y']) * 1e-3  # [km]
tk = np.array(data['tk'])       # [K]
pr = np.array(data['pr']) * 1e-6 # [MPa]
pf = np.array(data['pf']) * 1e-6 # [MPa]
qx = np.array(data['qxD'])      # [m/s]
qy = np.array(data['qyD'])      # [m/s]
peff = np.array(data['peff']) * 1e-6 # [MPa]

rplanet_km = data['rplanet'] * 1e-3
rcrust_km = data['rcrust'] * 1e-3
xcenter_km = data['xcenter'] * 1e-3
ycenter_km = data['ycenter'] * 1e-3
t_Ma = data['timesum_Ma']

# Coordinates for pressure nodes (staggered grid, size Ny1 x Nx1)
dx_km = (x[1] - x[0])
dy_km = (y[1] - y[0])
xp_km = np.linspace(-dx_km/2, x[-1] + dx_km/2, tk.shape[1])
yp_km = np.linspace(-dy_km/2, y[-1] + dy_km/2, tk.shape[0])
Xp, Yp = np.meshgrid(xp_km, yp_km)

# Radial distance from planetesimal center [km]
R_km = np.sqrt((Xp - xcenter_km)**2 + (Yp - ycenter_km)**2)

# Figure setup: 2x2 grid, 12 x 10 inches
fig, axs = plt.subplots(2, 2, figsize=(12, 10), dpi=300)

# -------------------------------------------------------------
# Panel (a): 2D Temperature Field
# -------------------------------------------------------------
ax1 = axs[0, 0]
cf1 = ax1.contourf(Xp, Yp, tk, levels=30, cmap='inferno')
cbar1 = fig.colorbar(cf1, ax=ax1, fraction=0.046, pad=0.04)
cbar1.set_label('Temperature $T$ [K]', fontsize=10)

# Outlines for planetesimal surface
c_planet1 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--', label=f'Planetesimal surface ($R={rplanet_km:.0f}\\ \\mathrm{{km}}$)')
ax1.add_patch(c_planet1)

ax1.set_aspect('equal')
ax1.set_xlabel('Horizontal distance $x$ [km]', fontsize=10)
ax1.set_ylabel('Vertical depth $y$ [km]', fontsize=10)
ax1.set_title(f'(a) Temperature Field ($t = {t_Ma:.2f}\\ \\mathrm{{Ma}}$)', fontsize=11, fontweight='bold')
ax1.legend(loc='lower left', fontsize=8, framealpha=0.85)

# -------------------------------------------------------------
# Panel (b): Darcy Fluid Flux Vectors & Magnitude
# -------------------------------------------------------------
ax2 = axs[0, 1]
# Interpolate face fluxes to cell centers (P nodes)
qxp = np.zeros_like(qx)
qyp = np.zeros_like(qy)
qxp[:, 1:] = 0.5 * (qx[:, 1:] + qx[:, :-1])
qxp[:, 0] = qx[:, 0]
qyp[1:, :] = 0.5 * (qy[1:, :] + qy[:-1, :])
qyp[0, :] = qy[0, :]
q_mag = np.sqrt(qxp**2 + qyp**2)
q_mag_nonzero = np.where(q_mag > 1e-25, q_mag, 1e-25)

cf2 = ax2.contourf(Xp, Yp, np.log10(q_mag_nonzero), levels=np.linspace(-25, -17, 25), cmap='viridis', extend='both')
cbar2 = fig.colorbar(cf2, ax=ax2, fraction=0.046, pad=0.04)
cbar2.set_label(r'$\log_{10} \|\mathbf{q}^D\|$ [$\mathrm{m/s}$]', fontsize=10)

# Quiver sub-sampling (every 2 grid points) with cell-centered fluxes in melted core
step = 2
mask_arrows = tk[::step, ::step] > 273.0
qxp_plot = np.where(mask_arrows, qxp[::step, ::step], np.nan)
qyp_plot = np.where(mask_arrows, qyp[::step, ::step], np.nan)
ax2.quiver(Xp[::step, ::step], Yp[::step, ::step],
           qxp_plot, qyp_plot,
           color='white', alpha=0.9, scale=None, width=0.005)

c_planet2 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
ax2.add_patch(c_planet2)

ax2.set_aspect('equal')
ax2.set_xlabel('Horizontal distance $x$ [km]', fontsize=10)
ax2.set_ylabel('Vertical depth $y$ [km]', fontsize=10)
ax2.set_title('(b) Darcy Flux in Melted Core ($T > 273\\ \\mathrm{K}$)', fontsize=11, fontweight='bold')

# -------------------------------------------------------------
# Panel (c): Pore Fluid Pressure Field
# -------------------------------------------------------------
ax3 = axs[1, 0]
cf3 = ax3.contourf(Xp, Yp, pf, levels=30, cmap='plasma')
cbar3 = fig.colorbar(cf3, ax=ax3, fraction=0.046, pad=0.04)
cbar3.set_label(r'Pore Fluid Pressure $P_f$ [MPa]', fontsize=10)

c_planet3 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--', label='Planetesimal surface')
ax3.add_patch(c_planet3)

ax3.set_aspect('equal')
ax3.set_xlabel('Horizontal distance $x$ [km]', fontsize=10)
ax3.set_ylabel('Vertical depth $y$ [km]', fontsize=10)
ax3.set_title(r'(c) Pore Fluid Pressure Field $P_f$', fontsize=11, fontweight='bold')
ax3.legend(loc='lower left', fontsize=8, framealpha=0.85)

# -------------------------------------------------------------
# Panel (d): Radial Profiles from Center to Surface
# -------------------------------------------------------------
ax4 = axs[1, 1]
# Extract radial profiles along horizontal ray (y = ycenter)
mid_row = int(round(ycenter_km / dy_km))
r_ray = np.abs(xp_km - xcenter_km)
sort_idx = np.argsort(r_ray)
r_sorted = r_ray[sort_idx]
tk_sorted = tk[mid_row, :][sort_idx]
pt_sorted = pr[mid_row, :][sort_idx]
pf_sorted = pf[mid_row, :][sort_idx]

r_max_km = 70.0
mask = r_sorted <= r_max_km

color_T = '#d62728'
ax4.plot(r_sorted[mask], tk_sorted[mask], color=color_T, lw=2.2, label=r'Temperature $T$ [K]')
ax4.set_xlabel('Radial distance $r$ from center [km]', fontsize=10)
ax4.set_ylabel('Temperature $T$ [K]', color=color_T, fontsize=10)
ax4.tick_params(axis='y', labelcolor=color_T)
ax4.set_xlim(0, r_max_km)

ax4_p = ax4.twinx()
color_pt = '#1f77b4'
color_pf = '#2ca02c'
ax4_p.plot(r_sorted[mask], pt_sorted[mask], color=color_pt, lw=1.8, ls='-', label=r'Total pressure $P_t$ [MPa]')
ax4_p.plot(r_sorted[mask], pf_sorted[mask], color=color_pf, lw=1.8, ls='--', label=r'Fluid pressure $P_f$ [MPa]')
ax4_p.set_ylabel('Pressure [MPa]', color=color_pt, fontsize=10)
ax4_p.tick_params(axis='y', labelcolor=color_pt)

# Vertical boundaries
ax4.axvline(rplanet_km, color='cyan', ls='--', lw=1.5, label=f'Planetesimal surface (${rplanet_km:.0f}\\ \\mathrm{{km}}$)')
ax4.axvspan(0, rplanet_km, alpha=0.08, color='gray', label='Porous rocky body')

ax4.set_title('(d) Radial Temperature and Pressure Profiles', fontsize=11, fontweight='bold')
ax4.grid(True, alpha=0.3)

# Combine legends
lines1, labels1 = ax4.get_legend_handles_labels()
lines2, labels2 = ax4_p.get_legend_handles_labels()
ax4.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=8, framealpha=0.85)

plt.tight_layout()

# Save figure
out_dir = os.path.join(os.path.dirname(__file__), '..', 'docs', 'src', 'assets')
os.makedirs(out_dir, exist_ok=True)
png_path = os.path.join(out_dir, 'hydrothermal_circulation_benchmark.png')
svg_path = os.path.join(out_dir, 'hydrothermal_circulation_benchmark.svg')

plt.savefig(png_path, dpi=300)
plt.savefig(svg_path)
plt.close()

print(f"Generated {png_path} and {svg_path}")
