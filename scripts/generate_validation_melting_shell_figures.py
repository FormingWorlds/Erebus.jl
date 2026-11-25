#!/usr/bin/env python3
"""Plot the melting shell validation figure in a light and a dark version.

Usage: python3 scripts/generate_validation_melting_shell_figures.py <run_dir> [assets_dir] [preview_dir]

The run folder comes from scripts/run_validation_melting_shell.jl.
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
        "band": "#CBD5E1",
        "ref": "#0072B2",
        "c65": "#D55E00",
        "c129": "#009E73",
    },
    "dark": {
        "fig": "#0F172A",
        "ax": "#1E293B",
        "fg": "#F1F5F9",
        "grid": "#334155",
        "band": "#475569",
        "ref": "#3B9BD6",
        "c65": "#E8742A",
        "c129": "#21A985",
    },
}
R_KM = 1.0e-3
T_SOL_K = 1400.0
T_LIQ_K = 1800.0
CENTRE_ATOL_NX65_K = 20.0
GRID_MARKERS = {65: "D", 129: "o"}
MUTANTS = (("nx65_L0", "L = 0"), ("nx65", "L"), ("nx65_L1.5", "1.5 L"))


def load(path):
    """Read a CSV file with a header into a named structured array."""
    return np.atleast_1d(
        np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )


def result(run_dir, tag):
    """Return the result row of `tag` as a dict, from header.csv and result_<tag>.csv."""
    header = (run_dir / "header.csv").read_text().strip().split(",")
    row = (run_dir / f"result_{tag}.csv").read_text().strip().split(",")
    return dict(zip(header, row))


def style(fig, axes, t):
    """Apply the colours of theme `t` to a figure and its axes."""
    fig.patch.set_facecolor(t["fig"])
    for ax in axes:
        ax.set_facecolor(t["ax"])
        ax.tick_params(colors=t["fg"], which="both")
        ax.xaxis.label.set_color(t["fg"])
        ax.yaxis.label.set_color(t["fg"])
        ax.title.set_color(t["fg"])
        ax.grid(color=t["grid"], linewidth=0.8)
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


def figure(run_dir, t):
    """Plot the peak temperature profiles and the centre temperature of the mutants.

    Parameters
    ----------
    run_dir : Path
        Folder with profile_1d.csv, profile_nx65.csv, profile_nx129.csv and the result
        rows of the runs.
    t : dict
        Theme colours.

    Returns
    -------
    matplotlib.figure.Figure
        The two-panel figure.
    """
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(11.0, 4.4), constrained_layout=True, width_ratios=(1.5, 1.0)
    )
    ref = load(run_dir / "profile_1d.csv")
    a.plot(ref["r_m"] * R_KM, ref["tmax_K"], color=t["ref"], lw=2.2, label="1-D model")
    for nx in (65, 129):
        p = load(run_dir / f"profile_nx{nx}.csv")
        a.plot(
            p["r_m"] * R_KM,
            p["tmax_K"],
            ls="none",
            marker=GRID_MARKERS[nx],
            ms=4.5 if nx == 65 else 3.5,
            color=t[f"c{nx}"],
            label=f"Erebus, Nx = {nx}",
        )
    # Labels at the right end, above their lines, clear of the data (below the solidus
    # there) and of the melt radius ticks
    for temp, name in ((T_SOL_K, "solidus"), (T_LIQ_K, "liquidus")):
        a.axhline(temp, color=t["fg"], lw=1.1, ls=":")
        a.annotate(
            f"{name} {temp:.0f} K",
            (19.8, temp),
            xytext=(0, 4),
            textcoords="offset points",
            ha="right",
            va="bottom",
            color=t["fg"],
            fontsize=8.5,
        )
    marks = [(t["ref"], "-", float(result(run_dir, "nx65")["r_melt_ref"]))]
    marks += [
        (t[f"c{nx}"], "--", float(result(run_dir, f"nx{nx}")["r_melt"])) for nx in (65, 129)
    ]
    for color, ls, r_melt in marks:
        a.axvline(r_melt * R_KM, ymin=0.9, ymax=1.0, color=color, ls=ls, lw=2.0)
    a.set_xlim(0.0, 20.0)
    a.set_ylim(250.0, 1880.0)
    a.set_xlabel("Radius r (km)")
    a.set_ylabel(r"Peak temperature $T_\mathrm{max}$ (K)")
    a.set_title("(a) Peak temperature; ticks: melt radius")
    legend(a, t, loc="lower left")

    centre_ref = float(result(run_dir, "nx65")["T_centre_ref"])
    b.axhspan(
        centre_ref - CENTRE_ATOL_NX65_K,
        centre_ref + CENTRE_ATOL_NX65_K,
        color=t["band"],
        alpha=0.6,
        lw=0,
        label="1-D value ± 20 K (pass band)",
    )
    b.axhline(centre_ref, color=t["ref"], lw=2.0, label="1-D model at L")
    x = np.arange(len(MUTANTS))
    y = [float(result(run_dir, tag)["T_centre"]) for tag, _ in MUTANTS]
    b.plot(
        x,
        y,
        ls="none",
        marker="D",
        ms=8,
        color=t["c65"],
        label="Erebus, Nx = 65",
    )
    for xi, yi in zip(x, y):
        b.annotate(
            f"{yi - centre_ref:+.1f} K",
            (xi, yi),
            xytext=(10, -3),
            textcoords="offset points",
            color=t["fg"],
            fontsize=9,
        )
    b.set_xticks(x, [name for _, name in MUTANTS])
    b.set_xlim(-0.4, len(MUTANTS) - 0.2)
    b.set_xlabel("Latent heat of the Erebus run (L = 4e5 J/kg)")
    b.set_ylabel(r"Centre peak temperature (K)")
    b.set_ylim(1570.0, 1790.0)
    b.set_title("(b) Centre temperature, Nx = 65")
    legend(b, t, loc="lower left")
    style(fig, (a, b), t)
    return fig


def main(run_dir, assets_dir, preview_dir=None):
    """Read the run folder and write the figure in both themes."""
    run_dir, assets_dir = Path(run_dir), Path(assets_dir)
    for name, t in THEMES.items():
        fig = figure(run_dir, t)
        fig.savefig(
            assets_dir / f"validation_melting_shell_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            fig.savefig(
                Path(preview_dir) / f"melting_shell_{name}.png",
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
