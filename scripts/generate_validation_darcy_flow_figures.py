#!/usr/bin/env python3
"""Plot the spherical Darcy flow validation figures in a light and a dark version.

The profile figure compares the uniform-permeability case with the exact solution, the
jump figure does the same for the case with a low-permeability core.

Usage: python3 scripts/generate_validation_darcy_flow_figures.py <data_dir> [assets_dir] [preview_dir]

The data folder comes from scripts/export_validation_darcy_flow_data.jl.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

THEMES = {
    "light": {
        "fig": "#FFFFFF",
        "ax": "#F8FAFC",
        "fg": "#1E293B",
        "grid": "#E2E8F0",
        "ref": "#0072B2",
        "cont": "#CC79A7",
        "c65": "#D55E00",
        "c33": "#B8860B",
    },
    "dark": {
        "fig": "#0F172A",
        "ax": "#1E293B",
        "fg": "#F1F5F9",
        "grid": "#334155",
        "ref": "#56B4E9",
        "cont": "#E7A6C9",
        "c65": "#FF8C42",
        "c33": "#F0C419",
    },
}
R_PLANET_M = 50_000.0
NX_PROFILE = 65
PA_PER_MPA = 1.0e6
PA_PER_GPA = 1.0e9
A_JUMP_OVER_R = 0.5


def load(path):
    """Read a CSV file with a header into a named structured array."""
    return np.genfromtxt(path, delimiter=",", names=True)


def style(fig, axes, t):
    """Apply the colours of theme `t` to a figure and its axes."""
    fig.patch.set_facecolor(t["fig"])
    for ax in axes:
        ax.set_facecolor(t["ax"])
        ax.tick_params(colors=t["fg"])
        ax.xaxis.label.set_color(t["fg"])
        ax.yaxis.label.set_color(t["fg"])
        ax.title.set_color(t["fg"])
        ax.grid(color=t["grid"], linewidth=0.8, which="both")
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color(t["fg"])


def legend(ax, t, **kwargs):
    """Draw a legend with the colours of theme `t`."""
    leg = ax.legend(frameon=True, fontsize=9, **kwargs)
    leg.get_frame().set_facecolor(t["ax"])
    leg.get_frame().set_edgecolor(t["grid"])
    for text in leg.get_texts():
        text.set_color(t["fg"])


def figure(data, t):
    """Plot the excess pressure profiles and the profile error against the grid size."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    curves = data["curves"]
    x = curves["r_m"] / R_PLANET_M
    a.plot(
        x,
        curves["dp_exact_on_Pa"] / PA_PER_MPA,
        color=t["ref"],
        lw=2.2,
        label="Exact, sphere (d = 3)",
    )
    a.plot(
        x,
        curves["dp_exact_off_Pa"] / PA_PER_MPA,
        color=t["cont"],
        lw=2.0,
        ls="--",
        label="Exact, disk (d = 2)",
    )
    for key, col, mk, label in (
        ("on", t["c65"], "o", "Erebus, spherical_flow = true"),
        ("off", t["c33"], "D", "Erebus, spherical_flow = false"),
    ):
        d = data[key]
        a.plot(
            d["r_m"] / R_PLANET_M,
            d["dp_Pa"] / PA_PER_MPA,
            ls="none",
            marker=mk,
            ms=2.2,
            color=col,
            alpha=0.7,
            label=label,
            rasterized=True,
        )
    a.set_xlabel("Radius r / R")
    a.set_ylabel("Excess fluid pressure p - P_vent (MPa)")
    a.set_title(f"(a) Pressure profile, Nx = {NX_PROFILE}")
    legend(a, t, loc="lower left")

    conv = data["conv"]
    nx = conv["nx"]
    # The two errors nearly coincide: hollow wide markers below, filled small ones on top.
    b.loglog(
        nx,
        conv["l2_off"],
        color=t["c33"],
        marker="D",
        ms=11,
        mfc="none",
        mew=1.6,
        lw=2.6,
        label="spherical_flow = false",
    )
    b.loglog(
        nx,
        conv["l2_on"],
        color=t["c65"],
        marker="o",
        ms=5,
        lw=1.2,
        label="spherical_flow = true",
    )
    b.loglog(
        nx,
        conv["l2_on"][0] * nx[0] / nx,
        color=t["fg"],
        lw=1.2,
        ls=":",
        label="First order, 1/Nx",
    )
    b.set_xticks(nx, [str(int(n)) for n in nx])
    b.minorticks_off()
    b.set_xlabel("Grid nodes per side, Nx")
    b.set_ylabel("L2 error of p - P_vent, r <= 0.9 R")
    b.set_title("(b) Profile error against grid size")
    legend(b, t, loc="lower left")
    style(fig, (a, b), t)
    return fig


def jump_figure(data, t):
    """Plot the jump-case pressure profiles and the profile error against the grid size.

    Parameters
    ----------
    data : dict
        Structured arrays ``curves``, ``conv``, ``off`` and ``on`` of the jump case.
    t : dict
        Colours of the theme.

    Returns
    -------
    matplotlib.figure.Figure
        The figure with panel (a) the profiles in GPa and panel (b) the error.
    """
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    curves = data["curves"]
    x = curves["r_m"] / R_PLANET_M
    a.plot(
        x,
        curves["dp_exact_on_Pa"] / PA_PER_GPA,
        color=t["ref"],
        lw=2.2,
        label="Exact, sphere",
    )
    a.plot(
        x,
        curves["dp_exact_off_Pa"] / PA_PER_GPA,
        color=t["cont"],
        lw=2.0,
        ls="--",
        label="Exact, disk",
    )
    for key, col, mk, label in (
        ("on", t["c65"], "o", "Erebus, sphere"),
        ("off", t["c33"], "D", "Erebus, disk"),
    ):
        d = data[key]
        a.plot(
            d["r_m"] / R_PLANET_M,
            d["dp_Pa"] / PA_PER_GPA,
            ls="none",
            marker=mk,
            ms=2.2,
            color=col,
            alpha=0.7,
            label=label,
            rasterized=True,
        )
    a.axvline(A_JUMP_OVER_R, color=t["fg"], lw=1.2, ls=":", label="Jump, r = R/2")
    a.set_xlabel("Radius r / R")
    a.set_ylabel("Excess fluid pressure p - P_vent (GPa)")
    a.set_title(f"(a) Pressure profile, Nx = {NX_PROFILE}")
    legend(a, t, loc="upper right")

    conv = data["conv"]
    nx = conv["nx"]
    b.loglog(
        nx,
        conv["l2_off"],
        color=t["c33"],
        marker="D",
        ms=11,
        mfc="none",
        mew=1.6,
        lw=2.6,
        label="spherical_flow = false",
    )
    b.loglog(
        nx,
        conv["l2_on"],
        color=t["c65"],
        marker="o",
        ms=5,
        lw=1.2,
        label="spherical_flow = true",
    )
    b.set_xticks(nx, [str(int(n)) for n in nx])
    b.minorticks_off()
    b.set_xlabel("Grid nodes per side, Nx")
    b.set_ylabel("L2 error of p - P_vent, r <= 0.9 R")
    b.set_title("(b) Profile error against grid size")
    legend(b, t, loc="lower left")
    style(fig, (a, b), t)
    return fig


def main(data_dir, assets_dir, preview_dir=None):
    """Read the data folder and write both figures in both themes.

    Parameters
    ----------
    data_dir : str or Path
        Folder with the CSV files.
    assets_dir : str or Path
        Folder for the SVG figures.
    preview_dir : str or Path, optional
        Folder for PNG previews.

    Both output folders are created when they are missing.
    """
    data_dir, assets_dir = Path(data_dir), Path(assets_dir)
    assets_dir.mkdir(parents=True, exist_ok=True)
    if preview_dir is not None:
        Path(preview_dir).mkdir(parents=True, exist_ok=True)
    data = {
        "curves": load(data_dir / "reference_curves.csv"),
        "conv": load(data_dir / "convergence.csv"),
        "off": load(data_dir / f"profile_off_{NX_PROFILE}.csv"),
        "on": load(data_dir / f"profile_on_{NX_PROFILE}.csv"),
    }
    jump = {
        "curves": load(data_dir / "jump_reference_curves.csv"),
        "conv": load(data_dir / "jump_convergence.csv"),
        "off": load(data_dir / f"jump_profile_off_{NX_PROFILE}.csv"),
        "on": load(data_dir / f"jump_profile_on_{NX_PROFILE}.csv"),
    }
    for name, t in THEMES.items():
        fig = figure(data, t)
        fig.savefig(
            assets_dir / f"validation_darcy_flow_profile_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            fig.savefig(
                Path(preview_dir) / f"profile_{name}.png",
                dpi=90,
                facecolor=fig.get_facecolor(),
            )
        plt.close(fig)
        fig = jump_figure(jump, t)
        fig.savefig(
            assets_dir / f"validation_darcy_flow_jump_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            fig.savefig(
                Path(preview_dir) / f"jump_{name}.png",
                dpi=90,
                facecolor=fig.get_facecolor(),
            )
        plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(
        sys.argv[1],
        sys.argv[2] if len(sys.argv) > 2 else "docs/src/assets",
        sys.argv[3] if len(sys.argv) > 3 else None,
    )
