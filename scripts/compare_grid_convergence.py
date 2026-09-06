#!/usr/bin/env python3
"""
Compare 2D hydrothermal reaction benchmark across spatial grid resolutions (32x32 vs 128x128).
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path_32 = os.path.join(repo_root, "output_hydrothermal_reaction_on_32", "reaction_plot_data.json")
path_128 = os.path.join(repo_root, "output_hydrothermal_reaction_on_128", "reaction_plot_data.json")

with open(path_32) as f:
    d32 = json.load(f)
with open(path_128) as f:
    d128 = json.load(f)

t32 = np.array(d32["time_Ma"])
t128 = np.array(d128["time_Ma"])

fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=300)

# 1. Temperature comparison
ax = axes[0, 0]
ax.plot(t128, d128["max_T"], 'r-', lw=2, label='Max T (128x128)')
ax.plot(t32, d32["max_T"], 'r--', lw=1.5, label='Max T (32x32)')
ax.plot(t128, d128["mean_T"], 'b-', lw=2, label='Mean T (128x128)')
ax.plot(t32, d32["mean_T"], 'b--', lw=1.5, label='Mean T (32x32)')
ax.set_xlabel('Time [Ma]')
ax.set_ylabel('Temperature [K]')
ax.set_title('(a) Thermal Evolution')
ax.legend(loc='best', frameon=True)
ax.grid(True, alpha=0.3)

# 2. Hydration extent
ax = axes[0, 1]
ax.plot(t128, d128["max_XW"], color='darkgreen', ls='-', lw=2, label='Max $X_W$ (128x128)')
ax.plot(t32, d32["max_XW"], color='lightgreen', ls='--', lw=1.5, label='Max $X_W$ (32x32)')
ax.plot(t128, d128["mean_XW"], color='purple', ls='-', lw=2, label='Mean $X_W$ (128x128)')
ax.plot(t32, d32["mean_XW"], color='orchid', ls='--', lw=1.5, label='Mean $X_W$ (32x32)')
ax.set_xlabel('Time [Ma]')
ax.set_ylabel('Melt/Hydration Fraction $X_W$')
ax.set_title('(b) Hydration Reaction Extent')
ax.legend(loc='best', frameon=True)
ax.grid(True, alpha=0.3)

# 3. Water budget
ax = axes[1, 0]
tot128 = np.array(d128["water_solid"]) + np.array(d128["water_fluid"])
tot32 = np.array(d32["water_solid"]) + np.array(d32["water_fluid"])
ax.plot(t128, d128["water_solid"], 'g-', lw=2, label='Mineral Water (128x128)')
ax.plot(t32, d32["water_solid"], 'g--', lw=1.5, label='Mineral Water (32x32)')
ax.plot(t128, d128["water_fluid"], 'c-', lw=2, label='Pore Fluid (128x128)')
ax.plot(t32, d32["water_fluid"], 'c--', lw=1.5, label='Pore Fluid (32x32)')
ax.plot(t128, tot128, 'k-', lw=1.5, label='Total Inventory (128x128)')
ax.plot(t32, tot32, 'k:', lw=1.5, label='Total Inventory (32x32)')
ax.set_xlabel('Time [Ma]')
ax.set_ylabel('Water Mass [kg]')
ax.set_title('(c) Water Mass Partition')
ax.legend(loc='best', frameon=True)
ax.grid(True, alpha=0.3)

# 4. Darcy flux / Circulation intensity
ax = axes[1, 1]
ax.plot(t128, d128["mean_q"], color='darkorange', ls='-', lw=2, label='Mean |q| (128x128)')
ax.plot(t32, d32["mean_q"], color='gold', ls='--', lw=1.5, label='Mean |q| (32x32)')
ax.set_xlabel('Time [Ma]')
ax.set_ylabel('Darcy Flux Magnitude [m/s]')
ax.set_yscale('log')
ax.set_title('(d) Hydrothermal Circulation Rate')
ax.legend(loc='best', frameon=True)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_fig = os.path.join(repo_root, "docs", "src", "assets", "hydrothermal_grid_convergence.png")
plt.savefig(out_fig, dpi=300)
print(f"Generated grid convergence comparison plot: {out_fig}")
