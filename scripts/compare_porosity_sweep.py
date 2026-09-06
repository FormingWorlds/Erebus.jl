#!/usr/bin/env python3
"""
Compare hydrothermal hydration regimes across initial porosity sweep (phi = 0.20, 0.35, 0.50).
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dirs = [
    ("phi = 0.20", os.path.join(repo_root, "output_sweep_phi20")),
    ("phi = 0.35", os.path.join(repo_root, "output_sweep_phi35")),
    ("phi = 0.50", os.path.join(repo_root, "output_sweep_phi50"))
]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=300)

for idx, (label, dpath) in enumerate(dirs):
    json_file = os.path.join(dpath, "reaction_plot_data.json")
    if not os.path.isfile(json_file):
        print(f"Skipping {label}: {json_file} not ready yet.")
        continue

    with open(json_file) as f:
        data = json.load(f)

    t = np.array(data["time_Ma"])
    col = colors[idx]

    # Panel (a): Peak & Mean Temperature
    ax = axes[0, 0]
    ax.plot(t, data["max_T"], color=col, ls='-', lw=2, label=f"{label} (Max)")
    ax.plot(t, data["mean_T"], color=col, ls='--', lw=1.2, label=f"{label} (Mean)")

    # Panel (b): Hydration Extent
    ax = axes[0, 1]
    ax.plot(t, data["mean_XW"], color=col, ls='-', lw=2, label=label)

    # Panel (c): Solid Bound Water Fraction of Total Water
    ax = axes[1, 0]
    ws = np.array(data["water_solid"])
    wf = np.array(data["water_fluid"])
    wtot = ws + wf
    frac_bound = ws / np.maximum(1e-10, wtot)
    ax.plot(t, frac_bound * 100.0, color=col, ls='-', lw=2, label=label)

    # Panel (d): Darcy flux (circulation rate)
    ax = axes[1, 1]
    ax.plot(t, data["mean_q"], color=col, ls='-', lw=2, label=label)

# Configure axes
ax = axes[0, 0]
ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Temperature [K]")
ax.set_title("(a) Thermal Evolution")
ax.legend(loc="best", frameon=True, fontsize=9)
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Mean Wet-Solid Fraction $\\bar{X}_W$")
ax.set_title("(b) Global Hydration Extent")
ax.legend(loc="best", frameon=True)
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Mineral-Bound Water [% of total]")
ax.set_title("(c) Water Sequestration Efficiency")
ax.legend(loc="best", frameon=True)
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Mean Darcy Flux [m/s]")
ax.set_yscale("log")
ax.set_title("(d) Hydrothermal Circulation Vigor")
ax.legend(loc="best", frameon=True)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_fig = os.path.join(repo_root, "docs", "src", "assets", "hydrothermal_porosity_sweep.png")
plt.savefig(out_fig, dpi=300)
print(f"Generated porosity sweep comparison plot: {out_fig}")
