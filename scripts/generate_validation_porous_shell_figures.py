#!/usr/bin/env python3
"""Plot the porous shell validation figures in a light and a dark version.

Usage: python3 scripts/generate_validation_porous_shell_figures.py <data_dir> [assets_dir] [preview_dir]

The data folder comes from scripts/export_validation_porous_shell_data.jl.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

THEMES = {
    "light": {
        "fig": "#FFFFFF",
        "ax": "#F8FAFC",
        "fg": "#1E293B",
        "grid": "#E2E8F0",
        "ref": "#0072B2",
        "c65": "#D55E00",
        "c129": "#009E73",
        "air": "#CC79A7",
        "phi010": "#CC79A7",
        "phi025": "#0072B2",
        "phi050": "#D55E00",
        "phi075": "#009E73",
    },
    "dark": {
        "fig": "#0F172A",
        "ax": "#1E293B",
        "fg": "#F1F5F9",
        "grid": "#334155",
        "ref": "#56B4E9",
        "c65": "#FF8C42",
        "c129": "#35C9A1",
        "air": "#E7A6C9",
        "phi010": "#E7A6C9",
        "phi025": "#56B4E9",
        "phi050": "#FF8C42",
        "phi075": "#35C9A1",
    },
}
R_KM = 1.0e-3
R_PLANET_M = 20_000.0
L_BOX_M = 2.4 * R_PLANET_M
T_SINTER_K = 700.0
RUN_LIMIT_MYR = 20.0
ZOOM_KM = 8.0
DEPTH_REL_BOUND = 0.15
MAIN_CASE = {65: "s65_25_025", 129: "s129_25_025"}
GRID_MARKERS = {65: "D", 129: "o"}
PHI_STYLES = {
    0.1: ("phi010", "v"),
    0.25: ("phi025", "o"),
    0.5: ("phi050", "s"),
    0.75: ("phi075", "^"),
}


def load(path):
    """Read a CSV file with a header into a named structured array.

    Parameters
    ----------
    path : Path
        CSV file.

    Returns
    -------
    numpy.ndarray
        Structured array with at least one dimension; string columns stay strings.
    """
    return np.atleast_1d(
        np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    )


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


def suite_row(summary, tag):
    """Return the summary row of `tag`, or None when the run is missing."""
    rows = summary[summary["tag"] == tag]
    return rows[0] if len(rows) else None


def profile_panels(data, t):
    """Plot the peak temperature and porosity of the main case against radius.

    Parameters
    ----------
    data : dict
        Loaded data; uses ``summary`` and the ``profile_<tag>`` entries of `MAIN_CASE`.
    t : dict
        Theme colours.

    Returns
    -------
    matplotlib.figure.Figure
        The two-panel figure.
    """
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    grids = [nx for nx, tag in MAIN_CASE.items() if f"profile_{tag}" in data]
    ref = data[f"profile_{MAIN_CASE[max(grids)]}"]
    r0_km = R_PLANET_M * R_KM - ZOOM_KM
    for ax, col, ref_col in ((a, "tmax_K", "tmax_ref_K"), (b, "phi", "phi_ref")):
        ax.plot(
            ref["r_m"] * R_KM,
            ref[ref_col],
            color=t["ref"],
            lw=2.2,
            label="1-D reference",
        )
        for nx in grids:
            p = data[f"profile_{MAIN_CASE[nx]}"]
            ax.plot(
                p["r_m"] * R_KM,
                p[col],
                ls="none",
                marker=GRID_MARKERS[nx],
                ms=4.5,
                color=t[f"c{nx}"],
                label=f"Erebus, Nx = {nx}",
            )
        ax.set_xlim(r0_km, R_PLANET_M * R_KM)
        ax.set_xlabel("Radius r (km)")
    a.axhline(
        T_SINTER_K,
        color=t["fg"],
        lw=1.2,
        ls=":",
        label=f"No sintering below {T_SINTER_K:.0f} K",
    )
    b.set_ylim(-0.01, 0.32)
    edges = [
        (t["ref"], "-", suite_row(data["summary"], MAIN_CASE[max(grids)])["d_ref"])
    ]
    edges += [
        (t[f"c{nx}"], "--", suite_row(data["summary"], MAIN_CASE[nx])["d_shell"])
        for nx in grids
    ]
    for ax in (a, b):
        for color, ls, depth in edges:
            ax.axvline(
                (R_PLANET_M - depth) * R_KM,
                ymin=0.88,
                ymax=1.0,
                color=color,
                ls=ls,
                lw=2.0,
            )
    a.set_ylabel(r"Peak temperature $T_\mathrm{max}$ (K)")
    a.set_title(f"(a) Peak temperature, outer {ZOOM_KM:.0f} km")
    b.set_ylabel(r"Porosity $\phi$")
    b.set_title(f"(b) Final porosity, outer {ZOOM_KM:.0f} km")
    handles, labels = b.get_legend_handles_labels()
    handles.append(Line2D([], [], color=t["fg"], lw=2.0, ls="--"))
    labels.append(r"Shell edge $R_p - d$ (ticks)")
    legend(a, t, loc="lower left")
    legend(b, t, handles=handles, labels=labels, loc="lower right")
    style(fig, (a, b), t)
    return fig


def depth_panel(data, t):
    """Plot the Erebus shell depth against the 1-D shell depth for the suite runs.

    Parameters
    ----------
    data : dict
        Loaded data; uses ``summary``.
    t : dict
        Theme colours.

    Returns
    -------
    matplotlib.figure.Figure
        The one-panel figure.
    """
    fig, ax = plt.subplots(figsize=(5.6, 4.8), constrained_layout=True)
    s = data["summary"]
    top = 1.2 * max(np.max(s["d_ref"]), np.max(s["d_shell"])) * R_KM + 0.5
    x = np.linspace(0.0, top, 200)
    ax.plot(x, x, color=t["fg"], lw=1.4, label="1:1, equal depths")
    for nx, ls in ((65, "--"), (129, ":")):
        dx_km = L_BOX_M / (nx - 1) * R_KM
        half = np.maximum(dx_km, DEPTH_REL_BOUND * x)
        ax.plot(
            x,
            x + half,
            color=t[f"c{nx}"],
            ls=ls,
            lw=1.3,
            label=f"Bound, Nx = {nx}: max(Δx, 15 %)",
        )
        ax.plot(x, x - half, color=t[f"c{nx}"], ls=ls, lw=1.3)
    air = s["k_air"] > 0.0
    for nx in (65, 129):
        rows = s[(s["Nx"] == nx) & ~air]
        if len(rows):
            ax.plot(
                rows["d_ref"] * R_KM,
                rows["d_shell"] * R_KM,
                ls="none",
                marker=GRID_MARKERS[nx],
                ms=8,
                color=t[f"c{nx}"],
                label=f"Erebus, Nx = {nx}",
            )
    if np.any(air):
        rows = s[air]
        ax.plot(
            rows["d_ref"] * R_KM,
            rows["d_shell"] * R_KM,
            ls="none",
            marker="s",
            ms=11,
            mfc="none",
            mew=2.0,
            color=t["air"],
            label="air k = 3000 W/(m K)",
        )
    ax.set_xlim(0.0, top)
    ax.set_ylim(0.0, top)
    ax.set_aspect("equal")
    ax.set_xlabel(r"1-D shell depth $d_\mathrm{ref}$ (km)")
    ax.set_ylabel(r"Erebus shell depth $d$ (km)")
    ax.set_title("Shell depth, Erebus against the 1-D reference")
    legend(ax, t, loc="upper left")
    style(fig, (ax,), t)
    return fig


def table_panel(data, t):
    """Plot the 1-D shell depth over the planet radius against the formation time.

    Parameters
    ----------
    data : dict
        Loaded data; uses ``table``.
    t : dict
        Theme colours.

    Returns
    -------
    matplotlib.figure.Figure
        The one-panel figure.
    """
    fig, ax = plt.subplots(figsize=(8.4, 4.6), constrained_layout=True)
    tab = data["table"]
    handles, labels = [], []
    for phi, (key, mk) in PHI_STYLES.items():
        rows = np.sort(tab[np.isclose(tab["phi_init"], phi)], order="t_form_myr")
        if not len(rows):
            continue
        y = rows["d_shell_m"] / R_PLANET_M
        ax.plot(rows["t_form_myr"], y, color=t[key], lw=2.0)
        stopped = np.isclose(rows["t_stop_myr"], RUN_LIMIT_MYR)
        for sel, face in ((~stopped, t[key]), (stopped, t["ax"])):
            ax.plot(
                rows["t_form_myr"][sel],
                y[sel],
                ls="none",
                marker=mk,
                ms=8,
                mew=1.8,
                color=t[key],
                mfc=face,
            )
        handles.append(Line2D([], [], color=t[key], lw=2.0, marker=mk, ms=8))
        labels.append(rf"$\phi_\mathrm{{init}}$ = {phi:g}")
    four_dx = 4.0 * L_BOX_M / 128 / R_PLANET_M
    ax.axhline(four_dx, color=t["fg"], lw=1.3, ls="--")
    handles.append(Line2D([], [], color=t["fg"], lw=1.3, ls="--"))
    labels.append(r"4 $\Delta x$ at Nx = 129")
    handles.append(
        Line2D([], [], ls="none", marker="o", ms=8, mew=1.8, color=t["fg"], mfc="none")
    )
    labels.append(f"Open: stopped at the\n{RUN_LIMIT_MYR:.0f} Myr run limit")
    ax.set_yscale("log")
    ax.set_xlabel(r"Formation time $t_\mathrm{form}$ (Myr)")
    ax.set_ylabel(r"1-D shell depth over radius $d_\mathrm{shell} / R_p$")
    ax.set_title(r"1-D shell depth, $R_p$ = 20 km")
    legend(
        ax,
        t,
        handles=handles,
        labels=labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
    )
    style(fig, (ax,), t)
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

    Raises
    ------
    FileNotFoundError
        If the summary, the table or both main case profiles are missing.
    """
    data_dir, assets_dir = Path(data_dir), Path(assets_dir)
    data = {
        "summary": load(data_dir / "suite_summary.csv"),
        "table": load(data_dir / "table_1d.csv"),
    }
    for tag in MAIN_CASE.values():
        if (data_dir / f"profile_{tag}.csv").is_file():
            data[f"profile_{tag}"] = load(data_dir / f"profile_{tag}.csv")
    if not any(f"profile_{tag}" in data for tag in MAIN_CASE.values()):
        raise FileNotFoundError(f"no main case profile in {data_dir}")
    used = {
        "profiles": ", ".join(
            tag for tag in MAIN_CASE.values() if f"profile_{tag}" in data
        ),
        "depths": ", ".join(data["summary"]["tag"]),
        "table": f"{len(data['table'])} rows of table_1d.csv",
    }
    for name, t in THEMES.items():
        figs = {
            "profiles": profile_panels(data, t),
            "depths": depth_panel(data, t),
            "table": table_panel(data, t),
        }
        for key, fig in figs.items():
            fig.savefig(
                assets_dir / f"validation_porous_shell_{key}_{name}.svg",
                facecolor=fig.get_facecolor(),
                dpi=200,
            )
            if preview_dir is not None:
                fig.savefig(
                    Path(preview_dir) / f"{key}_{name}.png",
                    dpi=90,
                    facecolor=fig.get_facecolor(),
                )
            plt.close(fig)
    for key, rows in used.items():
        print(f"{key}: {rows}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(
        sys.argv[1],
        sys.argv[2] if len(sys.argv) > 2 else "docs/src/assets",
        sys.argv[3] if len(sys.argv) > 3 else None,
    )
