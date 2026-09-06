#!/usr/bin/env python3
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys

if len(sys.argv) > 1 and sys.argv[1].strip():
    data_path = sys.argv[1]
else:
    data_path = os.path.join(os.path.dirname(__file__), "..", "output_hydrothermal_reaction_on_32", "reaction_plot_data.json")

if not os.path.isfile(data_path):
    raise FileNotFoundError(f"Data file not found: {data_path}. Run export_reaction_data.jl first.")

with open(data_path, 'r') as f:
    data = json.load(f)

time_Ma = np.array(data['time_Ma'])
mean_T = np.array(data['mean_T'])
max_T = np.array(data['max_T'])
mean_XW = np.array(data['mean_XW'])
max_XW = np.array(data['max_XW'])
mean_q = np.array(data['mean_q'])
water_solid = np.array(data['water_solid'])
water_fluid = np.array(data['water_fluid'])

tk = np.array(data['tk'])
pf = np.array(data['pf']) * 1e-6 # MPa
XWS = np.array(data['XWS'])
DQPF = np.array(data['DQPF'])
DHP = np.array(data['DHP'])

x = np.array(data['x']) * 1e-3
y = np.array(data['y']) * 1e-3
rplanet_km = data['rplanet'] * 1e-3
xcenter_km = data['xcenter'] * 1e-3
ycenter_km = data['ycenter'] * 1e-3

dx_km = (x[1] - x[0]) if len(x) > 1 else 1.0
dy_km = (y[1] - y[0]) if len(y) > 1 else 1.0
xp_km = np.linspace(-dx_km/2, x[-1] + dx_km/2, tk.shape[1]) - xcenter_km
yp_km = np.linspace(-dy_km/2, y[-1] + dy_km/2, tk.shape[0]) - ycenter_km
Xp, Yp = np.meshgrid(xp_km, yp_km)

# Okabe-Ito colorblind friendly palette
C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_RED = "#D55E00"
C_SKY = "#56B4E9"
C_PURPLE = "#CC79A7"
C_BLACK = "#000000"

fig = plt.figure(figsize=(15, 10), dpi=300)

# 1. Time series: Temperatures
ax1 = plt.subplot(231)
ax1.plot(time_Ma, mean_T, color=C_BLUE, ls='-', lw=2, label='Mean T')
ax1.plot(time_Ma, max_T, color=C_RED, ls='--', lw=2, label='Max T')
ax1.axhline(1416.0, color=C_BLACK, ls=':', lw=1.5, label='Rock Solidus (1416 K)')
ax1.axhline(273.0, color=C_SKY, ls='-.', lw=1.5, label='Ice Melting (273 K)')
ax1.set_xlabel('Time [Ma]')
ax1.set_ylabel('Temperature [K]')
ax1.set_ylim(100.0, 1500.0)
ax1.set_title('(a) Temperature Evolution')
ax1.legend(loc='best', fontsize=8)
ax1.grid(True, alpha=0.3)

# 2. Time series: Water Budget
ax2 = plt.subplot(232)
ax2.plot(time_Ma, water_solid, color=C_GREEN, ls='-', lw=2, label='Solid (Mineral)')
ax2.plot(time_Ma, water_fluid, color=C_SKY, ls='--', lw=2, label='Fluid (Pores)')
ax2.plot(time_Ma, water_solid + water_fluid, color=C_BLACK, ls=':', lw=1.5, label='Total')
ax2.set_xlabel('Time [Ma]')
ax2.set_ylabel('Water Mass [kg]')
ax2.set_title('(b) Global Water Budget')
ax2.legend(loc='best', fontsize=8)
ax2.grid(True, alpha=0.3)

# 3. Time series: Mean XW and Darcy flux
ax3 = plt.subplot(233)
ax3.plot(time_Ma, mean_XW, color=C_PURPLE, ls='-', lw=2, label='Mean $X_W$')
ax3.set_xlabel('Time [Ma]')
ax3.set_ylabel('Mean $X_W$', color=C_PURPLE)
ax3.tick_params(axis='y', labelcolor=C_PURPLE)
ax3_twin = ax3.twinx()
ax3_twin.plot(time_Ma, mean_q, color=C_ORANGE, ls='--', lw=2, label='Mean |q|')
ax3_twin.set_ylabel('Mean Darcy flux [m/s]', color=C_ORANGE)
ax3_twin.tick_params(axis='y', labelcolor=C_ORANGE)
ax3.set_title('(c) Reaction & Circulation')
ax3.grid(True, alpha=0.3)

# Helper for 2D fields
def format_2d_ax(ax, title, has_ylabel=True):
    ax.set_facecolor('white')
    ax.add_patch(Circle((0.0, 0.0), rplanet_km, fill=False, edgecolor='black', lw=1.2, ls='-'))
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.set_xlabel("x [km]")
    if has_ylabel:
        ax.set_ylabel("y [km]")
    ax.set_xlim(xp_km[0], xp_km[-1])
    ax.set_ylim(yp_km[0], yp_km[-1])

# 4. 2D Field: Reaction Extent XW
ax4 = plt.subplot(234)
xw_levels = np.linspace(0.0, 1.0, 31)
cf4 = ax4.contourf(Xp, Yp, XWS, levels=xw_levels, cmap='cividis_r', vmin=0.0, vmax=1.0)
clip_c4 = Circle((0.0, 0.0), rplanet_km, transform=ax4.transData)
try:
    cf4.set_clip_path(clip_c4)
except Exception:
    for col in cf4.collections:
        col.set_clip_path(clip_c4)
fig.colorbar(cf4, ax=ax4, label='$X_W$', ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
format_2d_ax(ax4, '(d) Hydration Extent ($X_W$)', has_ylabel=True)

# 5. 2D Field: DQPF
ax5 = plt.subplot(235)
dqpf_limit = 1.0e-14
dqpf_levels = np.linspace(-dqpf_limit, dqpf_limit, 41)
cf5 = ax5.contourf(Xp, Yp, DQPF, levels=dqpf_levels, cmap='PuOr_r', vmin=-dqpf_limit, vmax=dqpf_limit, extend='both')
clip_c5 = Circle((0.0, 0.0), rplanet_km, transform=ax5.transData)
try:
    cf5.set_clip_path(clip_c5)
except Exception:
    for col in cf5.collections:
        col.set_clip_path(clip_c5)
fig.colorbar(cf5, ax=ax5, label='DQPF [1/s]', ticks=[-1.0e-14, -5.0e-15, 0.0, 5.0e-15, 1.0e-14], format='%.1e')
format_2d_ax(ax5, '(e) Fluid Source Term (DQPF)', has_ylabel=False)

# 6. 2D Field: DHP
ax6 = plt.subplot(236)
dhp_limit = 1.0e-8
dhp_levels = np.linspace(-dhp_limit, dhp_limit, 41)
cf6 = ax6.contourf(Xp, Yp, DHP, levels=dhp_levels, cmap='PuOr_r', vmin=-dhp_limit, vmax=dhp_limit, extend='both')
clip_c6 = Circle((0.0, 0.0), rplanet_km, transform=ax6.transData)
try:
    cf6.set_clip_path(clip_c6)
except Exception:
    for col in cf6.collections:
        col.set_clip_path(clip_c6)
fig.colorbar(cf6, ax=ax6, label='DHP [W/m³]', ticks=[-1.0e-8, -5.0e-9, 0.0, 5.0e-9, 1.0e-8], format='%.1e')
format_2d_ax(ax6, '(f) Latent Heat (DHP)', has_ylabel=False)

plt.tight_layout()

# Save
out_dir = os.path.dirname(data_path)
png_path = os.path.join(out_dir, 'reaction_benchmark_summary.png')
plt.savefig(png_path, dpi=300)
print(f"Generated {png_path}")

# Also save to docs assets if this is the 128x128 benchmark run
if "128" in out_dir:
    docs_asset = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'docs', 'src', 'assets', 'hydrothermal_reaction_128.png')
    plt.savefig(docs_asset, dpi=300)
    print(f"Saved benchmark summary asset to {docs_asset}")
