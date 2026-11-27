#!/usr/bin/env python3
"""Plot the magma-ocean conductivity validation figure in light and dark themes.

Usage: python3 scripts/generate_validation_magma_ocean_conductivity_figures.py <run_dir> [assets_dir] [preview_dir]

The run folder comes from scripts/run_validation_magma_ocean_conductivity.jl.
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
        "muted": "#64748B",
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
        "muted": "#94A3B8",
        "grid": "#334155",
        "band": "#475569",
        "ref": "#3B9BD6",
        "c65": "#E8742A",
        "c129": "#21A985",
    },
}
R_KM = 1.0e-3
T_SOL_K = 1400.0
T_MUSH_K = 1520.0
T_LIQ_K = 1800.0
CENTRE_ATOL_NX65_K = 8.0
GRID_MARKERS = {65: "D", 129: "o"}
MUTANTS = (
    ("nx65_exp02", r"$\beta = 0.2$"),
    ("nx65", r"$\beta = 1/3$"),
    ("nx65_exp05", r"$\beta = 0.5$"),
)


def load(path):
    """Read a CSV file with a header into a structured numpy array."""
    return np.atleast_1d(
        np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )


def result(run_dir, tag):
    """Return the result row of `tag` as a dict from header.csv and result_<tag>.csv."""
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
    """Draw a legend styled with theme `t`."""
    leg = ax.legend(frameon=True, fontsize=8.8, **kwargs)
    leg.get_frame().set_facecolor(t["ax"])
    leg.get_frame().set_edgecolor(t["grid"])
    for text in leg.get_texts():
        text.set_color(t["fg"])


def figure(run_dir, t):
    """Plot peak temperature profiles and centre temperature across exponents.

    Parameters
    ----------
    run_dir : Path
        Folder containing profile_1d.csv, profile_nx65.csv, profile_nx129.csv and
        result_<tag>.csv files.
    t : dict
        Theme colour palette.

    Returns
    -------
    matplotlib.figure.Figure
        Two-panel validation figure.
    """
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(11.2, 4.5), constrained_layout=True, width_ratios=(1.45, 1.0)
    )
    ref = load(run_dir / "profile_1d.csv")
    a.plot(
        ref["r_m"] * R_KM,
        ref["tmax_off_K"],
        color=t["muted"],
        lw=1.7,
        ls="--",
        label="1-D, soft turb. off",
    )
    a.plot(
        ref["r_m"] * R_KM,
        ref["tmax_K"],
        color=t["ref"],
        lw=2.2,
        label="1-D model",
    )
    for nx in (65, 129):
        p = load(run_dir / f"profile_nx{nx}.csv")
        a.plot(
            p["r_m"] * R_KM,
            p["tmax_K"],
            ls="none",
            marker=GRID_MARKERS[nx],
            ms=4.4 if nx == 65 else 3.4,
            color=t[f"c{nx}"],
            label=f"Erebus, Nx = {nx}",
        )
    for temp, name, x_pos, dy, ha, va in (
        (T_SOL_K, "solidus", 0.35, -4, "left", "top"),
        (T_MUSH_K, "mush onset", 0.35, -4, "left", "top"),
        (T_LIQ_K, "liquidus", 19.8, 5, "right", "bottom"),
    ):
        a.plot([0.0, 20.0], [temp, temp], color=t["fg"], lw=1.0, ls=":")
        a.annotate(
            f"{name} {temp:.0f} K",
            (x_pos, temp),
            xytext=(0, dy),
            textcoords="offset points",
            ha=ha,
            va=va,
            color=t["fg"],
            fontsize=8.2,
        )
    marks = [(t["ref"], "-", float(result(run_dir, "nx65")["r_mush_ref"]))]
    marks += [
        (t[f"c{nx}"], "--", float(result(run_dir, f"nx{nx}")["r_mush"]))
        for nx in (65, 129)
    ]
    for color, ls, r_mush in marks:
        a.plot(
            [r_mush * R_KM, r_mush * R_KM],
            [1870.0, 1980.0],
            color=color,
            ls=ls,
            lw=2.0,
        )
    a.set_xlim(0.0, 20.0)
    a.set_ylim(250.0, 1980.0)
    a.set_xlabel("Radius r (km)")
    a.set_ylabel(r"Peak temperature $T_\mathrm{max}$ (K)")
    a.set_title("(a) Peak temperature; ticks: mush radius")
    legend(a, t, loc="lower left")

    centre_ref = float(result(run_dir, "nx65")["T_centre_ref"])
    b.axhspan(
        centre_ref - CENTRE_ATOL_NX65_K,
        centre_ref + CENTRE_ATOL_NX65_K,
        color=t["band"],
        alpha=0.6,
        lw=0,
        label="1-D value ± 8 K (pass band)",
    )
    b.plot(
        [-0.45, 2.58],
        [centre_ref, centre_ref],
        color=t["ref"],
        lw=2.0,
        label=r"1-D model ($\beta = 1/3$)",
    )
    x = np.arange(len(MUTANTS))
    y = [float(result(run_dir, tag)["T_centre"]) for tag, _ in MUTANTS]
    for idx, (xi, yi) in enumerate(zip(x, y)):
        b.plot(
            [xi],
            [yi],
            ls="none",
            marker="D",
            ms=8,
            color=t["c65"],
            label="Erebus, Nx = 65" if idx == 0 else "_nolegend_",
        )
        offset = (10, 4) if yi > centre_ref else (10, -11)
        b.annotate(
            f"{yi - centre_ref:+.1f} K",
            (xi, yi),
            xytext=offset,
            textcoords="offset points",
            color=t["fg"],
            fontsize=8.8,
        )
    b.set_xticks(x, [name for _, name in MUTANTS])
    b.set_xlim(-0.45, 2.58)
    b.set_ylim(1558.0, 1608.0)
    b.set_xlabel(r"Soft-turbulence exponent $\beta$")
    b.set_ylabel("Centre peak temperature (K)")
    b.set_title("(b) Exponent sensitivity, Nx = 65")
    legend(b, t, loc="upper right")
    style(fig, (a, b), t)
    return fig


def main(run_dir, assets_dir, preview_dir=None):
    """Read the run folder and save light and dark SVG figures."""
    run_dir, assets_dir = Path(run_dir), Path(assets_dir)
    assets_dir.mkdir(parents=True, exist_ok=True)
    for name, t in THEMES.items():
        fig = figure(run_dir, t)
        fig.savefig(
            assets_dir / f"validation_magma_ocean_conductivity_{name}.svg",
            facecolor=fig.get_facecolor(),
            dpi=200,
        )
        if preview_dir is not None:
            Path(preview_dir).mkdir(parents=True, exist_ok=True)
            fig.savefig(
                Path(preview_dir) / f"magma_ocean_conductivity_{name}.png",
                dpi=110,
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
