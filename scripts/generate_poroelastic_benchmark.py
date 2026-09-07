#!/usr/bin/env python3
"""
Generate verification benchmark figure for poroelastic constitutive limits in Erebus.jl.

Produces docs/src/assets/poroelastic_verification.png (300 DPI) and .svg.
Reconciled with exact formulas from src/physics.jl:
  compute_drained_compressibility(beta_phi, phi, beta_s)
  compute_biot_willis_coefficient(beta_d, beta_s)
  compute_skempton_coefficient(beta_d, phi, beta_s, beta_f)
"""

import os
import numpy as np
import matplotlib.pyplot as plt

G_p = 1.0e10  # Elastic shear modulus [Pa]
phi_min = 1.0e-4
phi_max = 0.999


def compute_drained_compressibility(beta_phi, phi, beta_s):
    bphi = np.maximum(beta_phi, 0.0)
    bsolid = max(beta_s, 0.0)
    phi_eff = np.clip(phi, phi_min, phi_max)
    return (bphi + bsolid) / (1.0 - phi_eff)


def compute_biot_willis_coefficient(beta_d, beta_s):
    if beta_s <= 0.0:
        return np.ones_like(beta_d) if isinstance(beta_d, np.ndarray) else 1.0
    ratio = beta_s / np.maximum(beta_d, 1.0e-30)
    return np.clip(1.0 - ratio, 0.0, 1.0)


def compute_skempton_coefficient(beta_d, phi, beta_s, beta_f):
    bsolid = max(beta_s, 0.0)
    bfluid = np.maximum(beta_f, 0.0)
    phi_eff = np.clip(phi, phi_min, phi_max)
    num = beta_d - bsolid
    denom = num + phi_eff * (bfluid - bsolid)
    safe_denom = np.where(denom <= 0.0, 1.0, denom)
    ratio = np.where(denom <= 0.0, 1.0, num / safe_denom)
    return np.clip(ratio, 0.0, 1.0)


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(repo_root, "docs", "src", "assets")
    os.makedirs(out_dir, exist_ok=True)
    png_path = os.path.join(out_dir, "poroelastic_verification.png")
    svg_path = os.path.join(out_dir, "poroelastic_verification.svg")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), dpi=300)

    # -------------------------------------------------------------
    # Panel (a): Biot-Willis Coupling vs Porosity
    # -------------------------------------------------------------
    phi_arr = np.linspace(0.001, 0.50, 400)
    beta_phi_arr = phi_arr / G_p

    beta_s_values = [0.0, 1.0e-11, 2.5e-11, 5.0e-11]
    labels_a = [
        r"$\beta_s = 0$ (Incompressible grains, $K_{\mathrm{BW}} \equiv 1$)",
        r"$\beta_s = 1.0\times 10^{-11}\ \mathrm{Pa}^{-1}$",
        r"$\beta_s = 2.5\times 10^{-11}\ \mathrm{Pa}^{-1}$",
        r"$\beta_s = 5.0\times 10^{-11}\ \mathrm{Pa}^{-1}$",
    ]
    colors_a = ["#1b9e77", "#2b5c8f", "#7570b3", "#e7298a"]

    for beta_s, label, col in zip(beta_s_values, labels_a, colors_a):
        beta_d = compute_drained_compressibility(beta_phi_arr, phi_arr, beta_s)
        kbw = compute_biot_willis_coefficient(beta_d, beta_s)
        ax1.plot(phi_arr, kbw, label=label, color=col, lw=2.2)

    ax1.set_xlabel(r"Porosity $\phi$ [-]", fontsize=11, fontweight="bold")
    ax1.set_ylabel(r"Biot-Willis coefficient $K_{\mathrm{BW}}$ [-]", fontsize=11, fontweight="bold")
    ax1.set_title("(a) Biot-Willis Coupling vs Porosity", fontsize=11, fontweight="bold")
    ax1.set_xlim(0.0, 0.50)
    ax1.set_ylim(0.40, 1.05)
    ax1.grid(True, ls=":", alpha=0.6)
    ax1.legend(loc="lower right", fontsize=8.5, framealpha=0.9)

    # -------------------------------------------------------------
    # Panel (b): Undrained Pore Pressure Response (B)
    # -------------------------------------------------------------
    beta_f_arr = np.logspace(-11, -8, 400)  # Pa^-1
    beta_s_ref = 2.5e-11  # Pa^-1

    phi_samples = [0.01, 0.05, 0.10, 0.30]
    colors_b = ["#1b9e77", "#1f78b4", "#7570b3", "#e31a1c"]

    for phi_val, col in zip(phi_samples, colors_b):
        b_phi = phi_val / G_p
        b_d = compute_drained_compressibility(b_phi, phi_val, beta_s_ref)
        B_val = compute_skempton_coefficient(b_d, phi_val, beta_s_ref, beta_f_arr)
        ax2.plot(beta_f_arr, B_val, label=rf"$\phi = {phi_val:.2f}$", color=col, lw=2.2)

    # Reference water compressibility
    beta_f_water = 4.0e-10
    ax2.axvline(
        beta_f_water,
        color="gray",
        ls="--",
        lw=1.8,
        label=r"Water ($\beta_f = 4\times 10^{-10}\ \mathrm{Pa}^{-1}$)",
    )

    ax2.set_xscale("log")
    ax2.set_xlabel(r"Pore fluid compressibility $\beta_f$ [$\mathrm{Pa}^{-1}$]", fontsize=11, fontweight="bold")
    ax2.set_ylabel(r"Skempton coefficient $B$ [-]", fontsize=11, fontweight="bold")
    ax2.set_title(r"(b) Undrained Pore Pressure Response ($B$)", fontsize=11, fontweight="bold")
    ax2.set_xlim(1.0e-11, 1.0e-8)
    ax2.set_ylim(0.0, 1.05)
    ax2.grid(True, which="both", ls=":", alpha=0.6)
    ax2.legend(loc="lower left", fontsize=8.5, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.savefig(svg_path)
    plt.close()
    print(f"Generated {png_path} and {svg_path}")


if __name__ == "__main__":
    main()
