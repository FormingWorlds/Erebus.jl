#!/usr/bin/env python3
"""Plot the spherical conduction validation figures in a light and a dark version.

Usage: python3 scripts/generate_validation_conduction_figures.py <data_dir> [assets_dir] [preview_dir]

The data folder comes from scripts/export_validation_conduction_data.jl.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

THEMES = {
    "light": {"fig": "#FFFFFF", "ax": "#F8FAFC", "fg": "#1E293B", "grid": "#E2E8F0",
              "ref": "#0072B2", "cont": "#CC79A7", "fv": "#009E73",
              "c65": "#D55E00", "c33": "#B8860B"},
    "dark": {"fig": "#0F172A", "ax": "#1E293B", "fg": "#F1F5F9", "grid": "#334155",
             "ref": "#56B4E9", "cont": "#E7A6C9", "fv": "#35C9A1",
             "c65": "#FF8C42", "c33": "#F0C419"},
}
R_KM = 1.0e-3


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


def node_panels(data, case, ref_col, ref_label, t):
    """Plot a temperature profile panel and a difference panel for the grids Nx = 33 and 65."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    curves = data["curves"]
    a.plot(curves["r_m"] * R_KM, curves[ref_col], color=t["ref"], lw=2.2, label=ref_label)
    if case == "transient":
        a.plot(curves["r_m"] * R_KM, curves["dT_cont_K"], color=t["cont"], lw=1.8, ls=":",
               label="Continuous-time solution")
        a.plot(data["fv"]["r_m"][::25] * R_KM, data["fv"]["dT_K"][::25], ls="none", marker="s",
               ms=4, mfc="none", color=t["fv"], label="1-D finite-volume check")
    for nx, col, mk in ((33, t["c33"], "D"), (65, t["c65"], "o")):
        d = data[f"{case}_{nx}"]
        a.plot(d["r_m"] * R_KM, d["dT_sim_K"], ls="none", marker=mk, ms=2.2, color=col, alpha=0.7,
               label=f"Erebus, Nx = {nx}", rasterized=True)
        b.plot(d["r_m"] * R_KM, d["dT_sim_K"] - d["dT_ref_K"], ls="none", marker=mk, ms=2.2,
               color=col, alpha=0.7, label=f"Nx = {nx}", rasterized=True)
    b.axhline(0.0, color=t["ref"], lw=1.5)
    a.set_xlabel("Radius r (km)")
    a.set_ylabel("Temperature rise above the surface (K)")
    a.set_title("(a) Temperature profile")
    b.set_xlabel("Radius r (km)")
    b.set_ylabel("Erebus minus reference (K)")
    b.set_title("(b) Difference from the reference")
    legend(a, t, loc="lower left")
    legend(b, t, loc="best", markerscale=2.5)
    style(fig, (a, b), t)
    return fig


def guard_panels(data, t):
    """Plot the centre rise and L2 error of the marker seeds against the regression bands."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    seeds = data["seeds"]
    x = np.arange(len(seeds))
    labels = [str(int(s)) for s in seeds["seed"]]
    heat = data["heat"]
    panels = (
        (a, "dT_center_K", 0.005 * seeds["dT_center_K"][0], "Centre temperature rise (K)",
         "(a) Centre temperature rise, Nx = 65"),
        (b, "l2_rel", 0.005, "L2 profile error", "(b) L2 profile error, Nx = 65"),
    )
    for ax, col, half, ylabel, title in panels:
        rec = seeds[col][0]
        ax.axhspan(rec - half, rec + half, color=t["ref"], alpha=0.22, label="Regression band")
        ax.axhline(rec, color=t["ref"], lw=1.6, label="Recorded value (seed 42)")
        ax.plot(x, seeds[col], ls="none", marker="o", ms=6, color=t["c65"], label="Marker seeds")
        ax.axhline(float(heat[col]), color=t["cont"], lw=1.8, ls="--",
                   label="Heating raised by 2 %")
        ax.set_xticks(x, labels)
        ax.set_xlabel("Marker seed (solver.seed)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    legend(a, t, loc="upper left", bbox_to_anchor=(0.0, 0.9))
    legend(b, t, loc="center", bbox_to_anchor=(0.5, 0.25))
    style(fig, (a, b), t)
    return fig


def main(data_dir, assets_dir, preview_dir=None):
    """Read the data folder and write each figure in both themes.

    Parameters
    ----------
    data_dir : str or Path
        Folder with the CSV files.
    assets_dir : str or Path
        Folder for the SVG figures.
    preview_dir : str or Path, optional
        Folder for PNG previews.
    """
    data_dir, assets_dir = Path(data_dir), Path(assets_dir)
    data = {
        "curves": load(data_dir / "reference_curves.csv"),
        "fv": load(data_dir / "fv_1d.csv"),
        "seeds": load(data_dir / "seeds.csv"),
        "heat": load(data_dir / "heating_plus_2pct.csv"),
    }
    for case in ("steady", "transient"):
        for nx in (33, 65):
            data[f"{case}_{nx}"] = load(data_dir / f"{case}_{nx}.csv")
    for name, t in THEMES.items():
        figs = {
            "steady": node_panels(data, "steady", "dT_steady_K", "Exact steady solution", t),
            "transient": node_panels(data, "transient", "dT_be_K",
                                     "Backward-Euler reference", t),
            "bands": guard_panels(data, t),
        }
        for key, fig in figs.items():
            fig.savefig(assets_dir / f"validation_conduction_{key}_{name}.svg",
                        facecolor=fig.get_facecolor(), dpi=200)
            if preview_dir is not None:
                fig.savefig(Path(preview_dir) / f"{key}_{name}.png", dpi=90,
                            facecolor=fig.get_facecolor())
            plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(
        sys.argv[1],
        sys.argv[2] if len(sys.argv) > 2 else "docs/src/assets",
        sys.argv[3] if len(sys.argv) > 3 else None,
    )
