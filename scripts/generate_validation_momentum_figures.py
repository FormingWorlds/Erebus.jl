#!/usr/bin/env python3
"""Plot the spherical momentum validation figure in a light and a dark version.

Panel (a) compares the azimuthal velocity of case A on the axes through the centre with
the exact profiles with and without the spherical metric terms; panel (b) shows the
errors of the record against the grid size with a second-order guide.

Usage: python3 scripts/generate_validation_momentum_figures.py <data_dir> [assets_dir] [preview_dir]

The data folder comes from scripts/export_validation_momentum_data.jl.
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
        "cb": "#009E73",
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
        "cb": "#3DDC97",
    },
}
NX_PROFILE = 65


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
    """Plot the azimuthal velocity of case A against radius with and without the option.

    Parameters
    ----------
    a : matplotlib.axes.Axes
        Panel to draw on.
    data : dict
        Structured arrays ``curves``, ``on`` and ``off``.
    t : dict
        Colours of the theme.
    """
    curves = data["curves"]
    a.plot(
        curves["rho"], curves["u_on"], color=t["ref"], lw=2.2, label="Exact, option on"
    )
    a.plot(
        curves["rho"],
        curves["u_off"],
        color=t["cont"],
        lw=2.0,
        ls="--",
        label="Exact, option off",
    )
    for key, col, mk, label in (
        ("on", t["c65"], "o", "Erebus, spherical_flow = true"),
        ("off", t["c33"], "D", "Erebus, spherical_flow = false"),
    ):
        d = data[key]
        a.plot(d["rho"], d["u"], ls="none", marker=mk, ms=3.0, color=col, label=label)
    a.axhline(0.0, color=t["fg"], lw=0.6)
    a.set_xlabel("Radius r / R")
    a.set_ylabel("Azimuthal velocity u_theta eta / (F0 R^2)")
    a.set_title(f"(a) Case A velocity on the axes, Nx = {NX_PROFILE}")
    a.set_xlim(0.0, 1.0)
    legend(a, t, loc="upper right", fontsize=8)


def error_panel(b, data, t):
    """Plot the errors of the record against the grid size.

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
    rows = (
        ("l2_a_on", t["c65"], "o", "-", "Case A L2, option on"),
        ("l2_a_off", t["c33"], "D", "-", "Case A L2, option off"),
        ("l2_b", t["cb"], "s", "-", "Case B L2"),
        ("s33_b", t["cb"], "v", "--", "Case B S33 L2"),
        ("centre_fit_a_on", t["ref"], "^", "--", "Case A fitted centre, option on"),
        ("centre_node_a_on", t["cont"], "x", ":", "Case A centre node, option on"),
    )
    for key, col, mk, ls, label in rows:
        b.loglog(nx, conv[key], color=col, marker=mk, ms=6, lw=1.4, ls=ls, label=label)
    b.loglog(
        nx,
        0.5 * conv["l2_a_on"][0] * (nx[0] / nx) ** 2,
        color=t["fg"],
        lw=1.2,
        ls=":",
        label="Second order, 1/Nx^2",
    )
    b.set_xticks(nx, [str(int(n)) for n in nx])
    b.minorticks_off()
    b.set_ylim(3.0e-4, 0.3)
    b.set_xlabel("Grid nodes per side, Nx")
    b.set_ylabel("Relative error")
    b.set_title("(b) Errors against grid size")
    legend(b, t, loc="upper right", fontsize=7.5)


def figure(data, t):
    """Plot the case A profiles and the errors against the grid size.

    Parameters
    ----------
    data : dict
        Structured arrays ``curves``, ``conv``, ``on`` and ``off``.
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
        "on": load(data_dir / f"profile_on_{NX_PROFILE}.csv"),
        "off": load(data_dir / f"profile_off_{NX_PROFILE}.csv"),
    }
    for name, t in THEMES.items():
        fig = figure(data, t)
        fig.savefig(
            assets_dir / f"validation_momentum_profile_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            fig.savefig(
                Path(preview_dir) / f"momentum_profile_{name}.png",
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
