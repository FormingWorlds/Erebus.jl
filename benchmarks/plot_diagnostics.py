#!/usr/bin/env python3
"""
Diagnostic plotting script for Erebus.jl Soft Turbulence Model and Magma Ocean Benchmarks.
Applies Interra visual identity guidelines (palette, fonts, clean layout).
"""

import json
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Interra palette
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

# Load benchmark data
data_path = os.path.join(os.path.dirname(__file__), "..", "output_files", "soft_turbulence_benchmark_data.json")
if not os.path.exists(data_path):
    print(f"Error: {data_path} not found. Run generate_soft_turbulence_benchmarks.jl first.")
    sys.exit(1)

with open(data_path, "r") as f:
    data = json.load(f)

out_dir = os.path.join(os.path.dirname(__file__), "..", "output_files")
os.makedirs(out_dir, exist_ok=True)

Fm = np.array(data["Fm_vals"])
dT = np.array(data["dT_vals"])
w_T = np.array(data["w_T_vals"])

# -----------------------------------------------------------------------------
# Figure 1: Conductivity Regularization & Derivative Comparison
# -----------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

# Panel (a): keff vs Fm
ax1.text(0.04, 0.93, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

# Plot i2elvis reference (discontinuous step at Fm = 0.40)
ax1.plot(Fm, data["k_i2elvis_curves"]["100"], '--', color=NEUTRALS['mist'], linewidth=2.0,
         label=r'Raw i2elvis ($F_m=0.40$ step)')

# Plot Erebus smooth regularized curves for different eta_fluid
colors_eta = {'10': STRATA['gold'], '100': STRATA['amber'], '1000': STRATA['magma']}
labels_eta = {
    '10': r'Erebus ($\eta_\mathrm{fluid} = 10\ \mathrm{Pa\cdot s}$)',
    '100': r'Erebus ($\eta_\mathrm{fluid} = 100\ \mathrm{Pa\cdot s}$, baseline)',
    '1000': r'Erebus ($\eta_\mathrm{fluid} = 1000\ \mathrm{Pa\cdot s}$)'
}

for key in ['10', '100', '1000']:
    ax1.plot(Fm, data["k_erebus_curves"][key], '-', color=colors_eta[key], linewidth=2.2, label=labels_eta[key])

# Shaded transition zone
ax1.axvspan(0.30, 0.50, color=STRATA['gold'], alpha=0.12, label=r'Transition window $[0.30, 0.50]$')

ax1.set_yscale('log')
ax1.set_xlim(0.0, 1.0)
ax1.set_ylim(1.0, 5e6)
ax1.set_xlabel(r'Silicate Melt Fraction $F_m$ [-]', fontsize=11)
ax1.set_ylabel(r'Effective Thermal Conductivity $k_\mathrm{eff}$ [$\mathrm{W/(m\cdot K)}$]', fontsize=11)
ax1.set_title('Thermal Conductivity Scaling', fontsize=12, fontweight='bold', pad=10)
ax1.grid(True)
ax1.legend(loc='lower right', fontsize=8.5)

# Panel (b): Logarithmic derivative d(log10 k)/dFm (smoothness proof)
ax2.text(0.04, 0.93, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

for key in ['10', '100', '1000']:
    ax2.plot(Fm, data["dk_dFm_erebus"][key], '-', color=colors_eta[key], linewidth=2.2, label=labels_eta[key])

# Mark i2elvis delta singularity
ax2.axvline(0.40, color=STRATA['ink'], linestyle=':', linewidth=2.0, label=r'i2elvis Dirac $\delta$-spike at $0.40$')

ax2.set_xlim(0.20, 0.60)
ax2.set_ylim(-2, 45)
ax2.set_xlabel(r'Silicate Melt Fraction $F_m$ [-]', fontsize=11)
ax2.set_ylabel(r'$\mathrm{d}(\log_{10} k_\mathrm{eff}) / \mathrm{d}F_m$ [-]', fontsize=11)
ax2.set_title('Logarithmic Gradient Smoothness', fontsize=12, fontweight='bold', pad=10)
ax2.grid(True)
ax2.legend(loc='upper right', fontsize=8.5)

# Panel (c): Surface temperature difference gating
ax3.text(0.04, 0.93, '(c)', transform=ax3.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

ax3.plot(dT, w_T, '-', color=STRATA['cobalt'], linewidth=2.4, label=r'Weighting factor $w_T = [\mathrm{clamp}(\Delta T / \Delta T_\mathrm{min}, 0, 1)]^2$')
ax3.axvline(10.0, color=STRATA['amber'], linestyle='--', linewidth=1.5, label=r'$\Delta T_\mathrm{min} = 10\ \mathrm{K}$ anchor')
ax3.axhspan(0, 1, color=STRATA['cobalt'], alpha=0.08)

ax3.set_xlim(0.0, 30.0)
ax3.set_ylim(-0.05, 1.05)
ax3.set_xlabel(r'Temperature Difference $\Delta T = T - T_\mathrm{surface}$ [K]', fontsize=11)
ax3.set_ylabel(r'Thermal Boundary Weight $w_T$ [-]', fontsize=11)
ax3.set_title('Surface Boundary Weighting', fontsize=12, fontweight='bold', pad=10)
ax3.grid(True)
ax3.legend(loc='lower right', fontsize=8.5)

fig.savefig(os.path.join(out_dir, "diagnostic_conductivity_regularization.pdf"))
fig.savefig(os.path.join(out_dir, "diagnostic_conductivity_regularization.png"), dpi=300)
plt.close(fig)
print("Saved diagnostic_conductivity_regularization.[pdf,png]")

# -----------------------------------------------------------------------------
# Figure 2: Planetesimal Cooling and Magma Ocean Freezing Comparison
# -----------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

res_off = data["res_off"]
res_on = data["res_on"]

r_km = np.array(res_off["r"]) / 1000.0 # to km
t_off = np.array(res_off["times_yr"])
t_on = np.array(res_on["times_yr"])

# Panel (a): Core Temperature Evolution T_core(t)
ax1.text(0.04, 0.93, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

ax1.plot(t_off / 1000.0, res_off["T_core_hist"], '-', color=STRATA['cobalt'], linewidth=2.2, label='Conduction only (Soft Turb. OFF)')
ax1.plot(t_on / 1000.0, res_on["T_core_hist"], '-', color=STRATA['magma'], linewidth=2.2, label='Regularized Soft Turb. ON')
ax1.axhline(1800.0, color=STRATA['gold'], linestyle=':', linewidth=1.5, label='Liquidus $T_l = 1800$ K')
ax1.axhline(1400.0, color=STRATA['ink'], linestyle=':', linewidth=1.5, label='Solidus $T_s = 1400$ K')

ax1.set_xlim(0, 50.0)
ax1.set_ylim(1300, 1900)
ax1.set_xlabel('Time [kyr]', fontsize=11)
ax1.set_ylabel(r'Central Temperature $T_\mathrm{core}$ [K]', fontsize=11)
ax1.set_title('Core Thermal Quenching', fontsize=12, fontweight='bold', pad=10)
ax1.grid(True)
ax1.legend(loc='lower left', fontsize=8.5)

# Panel (b): Radial Temperature Profiles T(r) at Selected Times
ax2.text(0.04, 0.93, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

times_to_plot = ["0.0", "5000.0", "15000.0", "50000.0"]
time_colors = [STRATA['gold'], STRATA['amber'], STRATA['magma'], STRATA['plum']]

for t_key, col in zip(times_to_plot, time_colors):
    t_kyr = float(t_key) / 1000.0
    if t_key in res_off["snapshot_T"]:
        ax2.plot(r_km, res_off["snapshot_T"][t_key], ':', color=col, linewidth=1.6)
    if t_key in res_on["snapshot_T"]:
        ax2.plot(r_km, res_on["snapshot_T"][t_key], '-', color=col, linewidth=2.2, label=f'$t = {t_kyr:.0f}$ kyr')

# Dummy lines for legend
ax2.plot([], [], '-', color=NEUTRALS['graphite'], linewidth=2.0, label='Soft Turb. ON')
ax2.plot([], [], ':', color=NEUTRALS['graphite'], linewidth=1.6, label='Soft Turb. OFF')

ax2.set_xlim(0, 50)
ax2.set_ylim(200, 1900)
ax2.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11)
ax2.set_ylabel('Temperature $T(r)$ [K]', fontsize=11)
ax2.set_title('Radial Temperature Structure', fontsize=12, fontweight='bold', pad=10)
ax2.grid(True)
ax2.legend(loc='lower left', fontsize=8.5)

# Panel (c): Magma Ocean Solidification Front r_melt(t)
ax3.text(0.04, 0.93, '(c)', transform=ax3.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

ax3.plot(t_off / 1000.0, np.array(res_off["melt_radius_hist"]) / 1000.0, '-', color=STRATA['cobalt'], linewidth=2.2, label='Conduction only (Soft Turb. OFF)')
ax3.plot(t_on / 1000.0, np.array(res_on["melt_radius_hist"]) / 1000.0, '-', color=STRATA['magma'], linewidth=2.2, label='Regularized Soft Turb. ON')

ax3.set_xlim(0, 50.0)
ax3.set_ylim(0, 35.0)
ax3.set_xlabel('Time [kyr]', fontsize=11)
ax3.set_ylabel(r'Magma Ocean Radius ($F_m \geq 0.40$) [km]', fontsize=11)
ax3.set_title('Magma Ocean Crystallization Front', fontsize=12, fontweight='bold', pad=10)
ax3.grid(True)
ax3.legend(loc='upper right', fontsize=8.5)

fig.savefig(os.path.join(out_dir, "diagnostic_planetesimal_cooling.pdf"))
fig.savefig(os.path.join(out_dir, "diagnostic_planetesimal_cooling.png"), dpi=300)
plt.close(fig)
print("Saved diagnostic_planetesimal_cooling.[pdf,png]")

# -----------------------------------------------------------------------------
# Figure 3: Numerical Resolution Convergence & Boundary Invariant Verification
# -----------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

res_50 = data["res_grid50"]
res_100 = data["res_on"]
res_200 = data["res_grid200"]

r_50_km = np.array(res_50["r"]) / 1000.0
r_100_km = np.array(res_100["r"]) / 1000.0
r_200_km = np.array(res_200["r"]) / 1000.0

# Panel (a): Radial Temperature at t = 15 kyr across Resolutions
ax1.text(0.04, 0.93, '(a)', transform=ax1.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

ax1.plot(r_50_km, res_50["snapshot_T"]["15000.0"], '--', color=STRATA['cobalt'], linewidth=2.0, label=r'$N_r = 50$ cells ($\Delta r = 1.0$ km)')
ax1.plot(r_100_km, res_100["snapshot_T"]["15000.0"], '-', color=STRATA['amber'], linewidth=2.0, label=r'$N_r = 100$ cells ($\Delta r = 0.5$ km)')
ax1.plot(r_200_km, res_200["snapshot_T"]["15000.0"], ':', color=STRATA['magma'], linewidth=2.2, label=r'$N_r = 200$ cells ($\Delta r = 0.25$ km)')

ax1.set_xlim(0, 50)
ax1.set_ylim(250, 1900)
ax1.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11)
ax1.set_ylabel('Temperature $T(r)$ [K] at $t = 15$ kyr', fontsize=11)
ax1.set_title('Resolution Convergence: $T(r)$ Profile', fontsize=12, fontweight='bold', pad=10)
ax1.grid(True)
ax1.legend(loc='lower left', fontsize=8.5)

# Panel (b): Melt Fraction Profile Fm(r) at t = 15 kyr
ax2.text(0.04, 0.93, '(b)', transform=ax2.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

ax2.plot(r_50_km, res_50["snapshot_Fm"]["15000.0"], '--', color=STRATA['cobalt'], linewidth=2.0, label=r'$N_r = 50$')
ax2.plot(r_100_km, res_100["snapshot_Fm"]["15000.0"], '-', color=STRATA['amber'], linewidth=2.0, label=r'$N_r = 100$')
ax2.plot(r_200_km, res_200["snapshot_Fm"]["15000.0"], ':', color=STRATA['magma'], linewidth=2.2, label=r'$N_r = 200$')
ax2.axhline(0.40, color=STRATA['ink'], linestyle=':', linewidth=1.5, label=r'Rheological Transition $F_\mathrm{crit} = 0.40$')

ax2.set_xlim(0, 40)
ax2.set_ylim(-0.05, 1.05)
ax2.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11)
ax2.set_ylabel('Melt Fraction $F_m(r)$ [-] at $t = 15$ kyr', fontsize=11)
ax2.set_title('Resolution Convergence: Melt Distribution', fontsize=12, fontweight='bold', pad=10)
ax2.grid(True)
ax2.legend(loc='lower left', fontsize=8.5)

# Panel (c): Convergence Error Profile Relative to Nr = 200
ax3.text(0.04, 0.93, '(c)', transform=ax3.transAxes, fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=NEUTRALS['mist'], alpha=0.9))

# Interpolate onto Nr=50 and Nr=100 grid points
T_ref = np.interp(r_50_km, r_200_km, res_200["snapshot_T"]["15000.0"])
err_50 = np.abs(np.array(res_50["snapshot_T"]["15000.0"]) - T_ref)

T_ref_100 = np.interp(r_100_km, r_200_km, res_200["snapshot_T"]["15000.0"])
err_100 = np.abs(np.array(res_100["snapshot_T"]["15000.0"]) - T_ref_100)

ax3.plot(r_50_km, err_50, '-', color=STRATA['cobalt'], linewidth=2.0, label=r'$|T_{50} - T_{200}|$ ($\mathrm{max} = ' + f'{np.max(err_50):.1f}' + r'\ \mathrm{K}$)')
ax3.plot(r_100_km, err_100, '-', color=STRATA['amber'], linewidth=2.0, label=r'$|T_{100} - T_{200}|$ ($\mathrm{max} = ' + f'{np.max(err_100):.1f}' + r'\ \mathrm{K}$)')

ax3.set_xlim(0, 50)
ax3.set_ylim(0, max(np.max(err_50)*1.1, 10.0))
ax3.set_xlabel('Planetesimal Radius $r$ [km]', fontsize=11)
ax3.set_ylabel(r'Absolute Temperature Error $|\Delta T|$ [K]', fontsize=11)
ax3.set_title('Grid Resolution Error Convergence', fontsize=12, fontweight='bold', pad=10)
ax3.grid(True)
ax3.legend(loc='upper right', fontsize=8.5)

fig.savefig(os.path.join(out_dir, "diagnostic_resolution_convergence.pdf"))
fig.savefig(os.path.join(out_dir, "diagnostic_resolution_convergence.png"), dpi=300)
plt.close(fig)
print("Saved diagnostic_resolution_convergence.[pdf,png]")

