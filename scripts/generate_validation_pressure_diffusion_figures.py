#!/usr/bin/env python3
"""Plot the spherical pressure diffusion validation figure in a light and a dark version.

Panel (a) compares the fluid overpressure at the grid nodes after 5 and 20 steps with the
backward-Euler series of the sphere and the disk; panel (b) shows the centre and profile
errors after 20 steps against the grid size with a first-order guide.

Usage: python3 scripts/generate_validation_pressure_diffusion_figures.py <data_dir> [assets_dir] [preview_dir]

The data folder comes from scripts/export_validation_pressure_diffusion_data.jl.
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
NX_PROFILE = 65
STEPS = (5, 20)


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
    kwargs.setdefault("fontsize", 9)
    leg = ax.legend(frameon=True, **kwargs)
    leg.get_frame().set_facecolor(t["ax"])
    leg.get_frame().set_edgecolor(t["grid"])
    for text in leg.get_texts():
        text.set_color(t["fg"])


def profile_panel(a, data, t):
    """Plot the fluid overpressure against radius after 5 and 20 steps.

    Parameters
    ----------
    a : matplotlib.axes.Axes
        Panel to draw on.
    data : dict
        Structured arrays ``curves``, ``off`` and ``on``.
    t : dict
        Colours of the theme.
    """
    curves = data["curves"]
    x = curves["x"]
    for m in STEPS:
        a.plot(
            x,
            curves[f"p_be_m{m}_on"],
            color=t["ref"],
            lw=2.2,
            label="Series, sphere (d = 3)" if m == STEPS[0] else None,
        )
        a.plot(
            x,
            curves[f"p_be_m{m}_off"],
            color=t["cont"],
            lw=2.0,
            ls="--",
            label="Series, disk (d = 2)" if m == STEPS[0] else None,
        )
    for key, col, mk, label in (
        ("on", t["c65"], "o", "Erebus, spherical_flow = true"),
        ("off", t["c33"], "D", "Erebus, spherical_flow = false"),
    ):
        d = data[key]
        order = np.argsort(d["x"])
        for m in STEPS:
            a.plot(
                d["x"][order],
                d[f"p_m{m}"][order],
                ls="none",
                marker=mk,
                ms=2.2,
                color=col,
                alpha=0.7,
                label=label if m == STEPS[0] else None,
                rasterized=True,
            )
    a.text(0.03, 0.9, "m = 5", transform=a.transAxes, color=t["fg"], fontsize=9)
    a.text(0.03, 0.42, "m = 20", transform=a.transAxes, color=t["fg"], fontsize=9)
    a.set_xlabel("Radius r / R")
    a.set_ylabel("Fluid overpressure (Pf - P_amb) / P0")
    a.set_title(f"(a) Pressure after 5 and 20 steps, Nx = {NX_PROFILE}")
    a.set_ylim(-0.05, 1.45)
    legend(a, t, loc="upper right", fontsize=8)


def error_panel(b, data, t):
    """Plot the centre and profile errors after 20 steps against the grid size.

    Parameters
    ----------
    b : matplotlib.axes.Axes
        Panel to draw on.
    data : dict
        Structured array ``conv``.
    t : dict
        Colours of the theme.
    """
    conv = data["conv"]
    nx = conv["nx"]
    for key, col, mk, label in (
        ("off", t["c33"], "D", "L2 error, spherical_flow = false"),
        ("on", t["c65"], "o", "L2 error, spherical_flow = true"),
    ):
        b.loglog(nx, conv[f"l2_{key}"], color=col, marker=mk, ms=6, lw=1.6, label=label)
    for key, col, mk, label in (
        ("off", t["c33"], "s", "Centre error, spherical_flow = false"),
        ("on", t["c65"], "^", "Centre error, spherical_flow = true"),
    ):
        b.loglog(
            nx,
            conv[f"centre_{key}"],
            color=col,
            marker=mk,
            ms=6,
            mfc="none",
            lw=1.2,
            ls="--",
            label=label,
        )
    b.loglog(
        nx,
        0.5 * conv["centre_off"][0] * nx[0] / nx,
        color=t["fg"],
        lw=1.2,
        ls=":",
        label="First order, 1/Nx",
    )
    b.set_xticks(nx, [str(int(n)) for n in nx])
    b.minorticks_off()
    b.set_ylim(3.0e-3, 1.5)
    b.set_xlabel("Grid nodes per side, Nx")
    b.set_ylabel("Error of (Pf - P_amb) / P0 after 20 steps")
    b.set_title("(b) Errors against grid size")
    legend(b, t, loc="upper right", fontsize=8)


def figure(data, t):
    """Plot the pressure profiles and the errors against the grid size.

    Parameters
    ----------
    data : dict
        Structured arrays ``curves``, ``conv``, ``off`` and ``on``.
    t : dict
        Colours of the theme.

    Returns
    -------
    matplotlib.figure.Figure
        The figure with panel (a) the profiles and panel (b) the errors.
    """
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    profile_panel(a, data, t)
    error_panel(b, data, t)
    style(fig, (a, b), t)
    return fig


def main(data_dir, assets_dir, preview_dir=None):
    """Read the data folder and write the figure in both themes.

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
    for name, t in THEMES.items():
        fig = figure(data, t)
        fig.savefig(
            assets_dir / f"validation_pressure_diffusion_profile_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            fig.savefig(
                Path(preview_dir) / f"pressure_diffusion_profile_{name}.png",
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
