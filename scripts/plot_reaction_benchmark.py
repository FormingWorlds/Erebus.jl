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
xp_km = np.linspace(-dx_km/2, x[-1] + dx_km/2, tk.shape[1])
yp_km = np.linspace(-dy_km/2, y[-1] + dy_km/2, tk.shape[0])
Xp, Yp = np.meshgrid(xp_km, yp_km)

fig = plt.figure(figsize=(15, 10), dpi=300)

# 1. Time series: Temperatures
ax1 = plt.subplot(231)
ax1.plot(time_Ma, mean_T, 'b-', label='Mean T')
ax1.plot(time_Ma, max_T, 'r--', label='Max T')
ax1.set_xlabel('Time [Ma]')
ax1.set_ylabel('Temperature [K]')
ax1.set_title('(a) Temperature Evolution')
ax1.legend()

# 2. Time series: Water Budget
ax2 = plt.subplot(232)
ax2.plot(time_Ma, water_solid, 'g-', label='Solid (Mineral)')
ax2.plot(time_Ma, water_fluid, 'c--', label='Fluid (Pores)')
ax2.plot(time_Ma, water_solid + water_fluid, 'k:', label='Total')
ax2.set_xlabel('Time [Ma]')
ax2.set_ylabel('Water Mass [kg]')
ax2.set_title('(b) Global Water Budget')
ax2.legend()

# 3. Time series: Mean XW and Darcy flux
ax3 = plt.subplot(233)
ax3.plot(time_Ma, mean_XW, 'purple', label='Mean XW')
ax3.set_xlabel('Time [Ma]')
ax3.set_ylabel('Mean X_W', color='purple')
ax3.tick_params(axis='y', labelcolor='purple')
ax3_twin = ax3.twinx()
ax3_twin.plot(time_Ma, mean_q, 'orange', ls='--', label='Mean |q|')
ax3_twin.set_ylabel('Mean Darcy flux [m/s]', color='orange')
ax3_twin.tick_params(axis='y', labelcolor='orange')
ax3.set_title('(c) Reaction & Circulation')

# 4. 2D Field: Reaction Extent XW
ax4 = plt.subplot(234)
cf4 = ax4.contourf(Xp, Yp, XWS, levels=30, cmap='YlGnBu')
fig.colorbar(cf4, ax=ax4, label='$X_W$')
c4 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
ax4.add_patch(c4)
ax4.set_aspect('equal')
ax4.set_title(f'(d) Hydration Extent ($X_W$)')

# 5. 2D Field: DQPF
ax5 = plt.subplot(235)
cf5 = ax5.contourf(Xp, Yp, DQPF, levels=30, cmap='RdBu_r')
fig.colorbar(cf5, ax=ax5, label='DQPF [1/s]')
c5 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
ax5.add_patch(c5)
ax5.set_aspect('equal')
ax5.set_title('(e) Fluid Source Term (DQPF)')

# 6. 2D Field: DHP
ax6 = plt.subplot(236)
cf6 = ax6.contourf(Xp, Yp, DHP, levels=30, cmap='coolwarm')
fig.colorbar(cf6, ax=ax6, label='DHP [W/m³]')
c6 = Circle((xcenter_km, ycenter_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
ax6.add_patch(c6)
ax6.set_aspect('equal')
ax6.set_title('(f) Latent Heat (DHP)')

plt.tight_layout()

# Save
out_dir = os.path.dirname(data_path)
png_path = os.path.join(out_dir, 'reaction_benchmark_summary.png')
plt.savefig(png_path, dpi=300)
print(f"Generated {png_path}")
