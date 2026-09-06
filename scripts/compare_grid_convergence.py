#!/usr/bin/env python3
"""
Compare 2D hydrothermal reaction benchmark across spatial grid resolutions (32x32, 64x64, 128x128).
"""
import os
import json
import numpy as np
import matplotlib.pyplot as plt

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

configs = [
    {"label": "32x32", "dir": "output_hydrothermal_reaction_on_32", "color": "#0072B2", "ls": "-", "lw": 1.8},
    {"label": "64x64", "dir": "output_hydrothermal_reaction_on_64", "color": "#009E73", "ls": "--", "lw": 1.8},
    {"label": "128x128", "dir": "output_hydrothermal_reaction_on_128", "color": "#D55E00", "ls": "-.", "lw": 2.0},
]

datasets = []
for cfg in configs:
    json_path = os.path.join(repo_root, cfg["dir"], "reaction_plot_data.json")
    if os.path.isfile(json_path):
        with open(json_path) as f:
            d = json.load(f)
        d["config"] = cfg
        d["time_Ma"] = np.array(d["time_Ma"])
        datasets.append(d)
        print(f"Loaded {cfg['label']} with {len(d['time_Ma'])} points (time range: {d['time_Ma'][0]:.2f} - {d['time_Ma'][-1]:.2f} Ma)")
    else:
        print(f"Skipping {cfg['label']}: {json_path} not found")

if not datasets:
    raise FileNotFoundError("No reaction_plot_data.json files found for grid convergence comparison.")

fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=300)

# (a) Thermal evolution
ax = axes[0, 0]
for d in datasets:
    cfg = d["config"]
    t = d["time_Ma"]
    ax.plot(t, d["max_T"], color=cfg["color"], ls=cfg["ls"], lw=cfg["lw"], label=f"Max T ({cfg['label']})")
    ax.plot(t, d["mean_T"], color=cfg["color"], ls=":", lw=1.2, alpha=0.8, label=f"Mean T ({cfg['label']})")

ax.axhline(1416.0, color="#000000", ls=":", lw=1.5, label="Rock Solidus (1416 K)")
ax.axhline(273.0, color="#56B4E9", ls="-.", lw=1.5, label="Ice Melting (273 K)")
ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Temperature [K]")
ax.set_ylim(100.0, 1500.0)
ax.set_title("(a) Thermal Evolution")
ax.legend(loc="upper right", frameon=True, fontsize=7.5, ncol=2)
ax.grid(True, alpha=0.3)

# (b) Hydration extent
ax = axes[0, 1]
for d in datasets:
    cfg = d["config"]
    t = d["time_Ma"]
    ax.plot(t, d["mean_XW"], color=cfg["color"], ls=cfg["ls"], lw=cfg["lw"], label=f"Mean $X_W$ ({cfg['label']})")
    ax.plot(t, d["max_XW"], color=cfg["color"], ls=":", lw=1.2, alpha=0.7, label=f"Max $X_W$ ({cfg['label']})")

ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Hydration Fraction $X_W$")
ax.set_ylim(-0.02, 1.05)
ax.set_title("(b) Hydration Reaction Extent")
ax.legend(loc="lower right", frameon=True, fontsize=7.5, ncol=2)
ax.grid(True, alpha=0.3)

# (c) Water mass partition: Mineral and Total Water
ax = axes[1, 0]
for d in datasets:
    cfg = d["config"]
    t = d["time_Ma"]
    ax.plot(t, d["water_solid"], color=cfg["color"], ls=cfg["ls"], lw=cfg["lw"], label=f"Mineral Water ({cfg['label']})")
    tot = np.array(d["water_solid"]) + np.array(d["water_fluid"])
    ax.plot(t, tot, color=cfg["color"], ls=":", lw=1.2, alpha=0.8, label=f"Total Water ({cfg['label']})")

ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Water Mass [kg]")
ax.set_title("(c) Water Mass Partition")
ax.legend(loc="lower right", frameon=True, fontsize=7.0, ncol=2)
ax.grid(True, alpha=0.3)

# (d) Hydrothermal Circulation Rate (Darcy flux)
ax = axes[1, 1]
for d in datasets:
    cfg = d["config"]
    t = d["time_Ma"]
    ax.plot(t, d["mean_q"], color=cfg["color"], ls=cfg["ls"], lw=cfg["lw"], label=f"Mean |q| ({cfg['label']})")

ax.set_xlabel("Time [Ma]")
ax.set_ylabel("Darcy Flux Magnitude [m/s]")
ax.set_yscale("log")
ax.set_title("(d) Hydrothermal Circulation Rate")
ax.legend(loc="upper right", frameon=True, fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_fig = os.path.join(repo_root, "docs", "src", "assets", "hydrothermal_grid_convergence.png")
plt.savefig(out_fig, dpi=300)
print(f"Generated grid convergence comparison plot: {out_fig}")
