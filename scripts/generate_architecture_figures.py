#!/usr/bin/env python3
"""Generate architecture and planetary evolution documentation figures for Erebus.jl.

Generates two standalone SVG figures in docs/src/assets/:
- Figure 1: erebus_execution_pipeline.svg (Execution spine and equation cards)
- Figure 2: erebus_planetesimal_stages.svg (Three-stage planetesimal evolution cutaway)

Installation:
    pip install git+https://github.com/FormingWorlds/fwl-visual-language.git
    pip install git+https://github.com/FormingWorlds/fwl-mpl.git

Usage:
    python3 scripts/generate_architecture_figures.py
    python3 scripts/generate_architecture_figures.py --fig 1
    python3 scripts/generate_architecture_figures.py --fig 2
    python3 scripts/generate_architecture_figures.py --fonttype path
    python3 scripts/generate_architecture_figures.py --fonttype none
"""

import argparse
from pathlib import Path

import fwl_mpl
import fwl_visual_language as fvl
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from fwl_visual_language._web import mix, text_safe
from fwl_visual_language.documenter import _tokens
from matplotlib import patheffects as pe
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Wedge

REPO_ROOT = Path(__file__).resolve().parent.parent
ASSETS_DIR = REPO_ROOT / "docs" / "src" / "assets"

fwl_mpl.use("light")

TOK = _tokens("light")
PAPER_ROLES = fvl.roles("paper")
INK_RAMP = fvl.ramp("ink")

BG = PAPER_ROLES["background"]
SURFACE = PAPER_ROLES["surface"]
SUB = INK_RAMP["100"]
LINE = PAPER_ROLES["line"]
INK_300 = INK_RAMP["300"]
AXIS = PAPER_ROLES["axis"]
MUTED = PAPER_ROLES["muted"]
TEXT = PAPER_ROLES["text"]
INK = PAPER_ROLES["heading"]


def C(name: str) -> dict:
    """Return fwl-visual-language tokens for a named palette hue.

    Parameters
    ----------
    name : str
        Token palette hue name (e.g., 'lagoon', 'fern').

    Returns
    -------
    dict
        Color dictionary with accent, solid, text, wash, tint, mid, and border keys.
    """
    acc = TOK["hues"][name]
    return {
        "accent": acc,
        "solid": text_safe(acc, INK, [BG]),
        "text": TOK["safe"][name],
        "wash": mix(BG, acc, 0.10),
        "tint": mix(BG, acc, 0.24),
        "mid": mix(BG, acc, 0.46),
        "border": acc,
    }


C_LAGOON = C("lagoon")
C_SIGNAL = C("signal")
C_FERN = C("fern")
C_MAGMA = C("magma")
C_GARNET = C("garnet")
C_IRIS = C("iris")
C_OCHRE = C("ochre")

C_INK = {
    "accent": MUTED,
    "solid": INK_RAMP["600"],
    "text": INK_RAMP["700"],
    "wash": SURFACE,
    "tint": SUB,
    "mid": LINE,
    "border": AXIS,
}


def new_equal_canvas(width_in: float, height_in: float, w_units: float, h_units: float):
    """Create a 1:1 equal-aspect canvas so geometry is strictly isotropic.

    Parameters
    ----------
    width_in : float
        Figure width in inches.
    height_in : float
        Figure height in inches.
    w_units : float
        Coordinate width limit in arbitrary unit space.
    h_units : float
        Coordinate height limit in arbitrary unit space.

    Returns
    -------
    tuple
        Tuple of (matplotlib.figure.Figure, matplotlib.axes.Axes).
    """
    fig = plt.figure(figsize=(width_in, height_in), dpi=200, facecolor=BG, layout="none")
    ax = fig.add_axes([0, 0, 1, 1], facecolor=BG)
    ax.set_xlim(0, w_units)
    ax.set_ylim(h_units, 0)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig, ax


def draw_box(ax, x, y, w, h, facecolor=BG, edgecolor=LINE, lw=2.0, radius=14, linestyle="-", alpha=1.0, zorder=2):
    """Draw a rounded bounding box patch onto axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left coordinate.
    y : float
        Top coordinate.
    w : float
        Width of box.
    h : float
        Height of box.
    facecolor : str, optional
        Fill color, by default BG.
    edgecolor : str, optional
        Border color, by default LINE.
    lw : float, optional
        Line width, by default 2.0.
    radius : float, optional
        Corner rounding radius, by default 14.
    linestyle : str, optional
        Stroke style, by default "-".
    alpha : float, optional
        Opacity fraction, by default 1.0.
    zorder : int, optional
        Layer ordering, by default 2.

    Returns
    -------
    matplotlib.patches.FancyBboxPatch
        The added patch object.
    """
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        linestyle=linestyle,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def draw_pill(ax, cx, cy, w, h, text, facecolor=BG, edgecolor=INK, text_color=TEXT, fontsize=12.0, weight="medium", radius=10, lw=1.8, zorder=4):
    """Draw a rounded label pill with centered text.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    cx : float
        Center x coordinate.
    cy : float
        Center y coordinate.
    w : float
        Pill width.
    h : float
        Pill height.
    text : str
        Label string.
    facecolor : str, optional
        Background fill color, by default BG.
    edgecolor : str, optional
        Border outline color, by default INK.
    text_color : str, optional
        Text foreground color, by default TEXT.
    fontsize : float, optional
        Font point size, by default 12.0.
    weight : str, optional
        Font weight, by default 'medium'.
    radius : float, optional
        Corner radius, by default 10.
    lw : float, optional
        Border stroke width, by default 1.8.
    zorder : int, optional
        Layer ordering, by default 4.

    Returns
    -------
    None
    """
    draw_box(ax, cx - w / 2.0, cy - h / 2.0, w, h, facecolor=facecolor, edgecolor=edgecolor, lw=lw, radius=radius, zorder=zorder)
    ax.text(cx, cy + 1, text, fontsize=fontsize, fontweight=weight, color=text_color, ha="center", va="center", zorder=zorder + 1)


def draw_arrow(ax, x1, y1, x2, y2, color=INK, lw=2.0, style="-|>", scale=14, rad=0.0, ls="-", zorder=5):
    """Draw a directed arrow patch between two points.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x1 : float
        Starting x coordinate.
    y1 : float
        Starting y coordinate.
    x2 : float
        Target x coordinate.
    y2 : float
        Target y coordinate.
    color : str, optional
        Stroke and arrowhead color, by default INK.
    lw : float, optional
        Line width, by default 2.0.
    style : str, optional
        Arrowhead style, by default '-|>'.
    scale : float, optional
        Mutation scale for head size, by default 14.
    rad : float, optional
        Curvature radius, by default 0.0.
    ls : str, optional
        Line stroke style, by default '-'.
    zorder : int, optional
        Layer ordering, by default 5.

    Returns
    -------
    matplotlib.patches.FancyArrowPatch
        The added arrow patch.
    """
    conn = f"arc3,rad={rad}" if rad != 0.0 else "arc3,rad=0.0"
    arr = FancyArrowPatch(
        (x1, y1), (x2, y2),
        connectionstyle=conn,
        arrowstyle=style,
        mutation_scale=scale,
        linewidth=lw,
        linestyle=ls,
        color=color,
        zorder=zorder,
    )
    ax.add_patch(arr)
    return arr


def draw_block_arrow(ax, x1, y1, x2, y2, facecolor=SUB, edgecolor=INK, width=12.0, head_width=24.0, head_len=14.0, lw=1.6, zorder=6):
    """Draw a 2D polygonal block arrow for flux vectors.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x1 : float
        Base center x coordinate.
    y1 : float
        Base center y coordinate.
    x2 : float
        Tip x coordinate.
    y2 : float
        Tip y coordinate.
    facecolor : str, optional
        Fill color, by default SUB.
    edgecolor : str, optional
        Border line color, by default INK.
    width : float, optional
        Shaft thickness, by default 12.0.
    head_width : float, optional
        Arrowhead base width, by default 24.0.
    head_len : float, optional
        Arrowhead longitudinal length, by default 14.0.
    lw : float, optional
        Outline stroke width, by default 1.6.
    zorder : int, optional
        Layer ordering, by default 6.

    Returns
    -------
    matplotlib.patches.Polygon or None
        Polygon patch or None if zero length.
    """
    dx = x2 - x1
    dy = y2 - y1
    length = np.hypot(dx, dy)
    if length < 1e-6:
        return
    ux, uy = dx / length, dy / length
    nx, ny = -uy, ux
    shaft_len = max(0.0, length - head_len)
    hw = width / 2.0
    hhw = head_width / 2.0

    bx, by = x1 + ux * shaft_len, y1 + uy * shaft_len
    pts = [
        (x1 + nx * hw, y1 + ny * hw),
        (bx + nx * hw, by + ny * hw),
        (bx + nx * hhw, by + ny * hhw),
        (x2, y2),
        (bx - nx * hhw, by - ny * hhw),
        (bx - nx * hw, by - ny * hw),
        (x1 - nx * hw, y1 - ny * hw),
    ]
    poly = Polygon(pts, closed=True, facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, zorder=zorder)
    ax.add_patch(poly)


# =============================================================================
# Figure 1: Execution Pipeline & Equation Cards (1A-2)
# =============================================================================

W1, H1 = 1300.0, 1040.0
CX1 = 650.0
NW1, NH1 = 264.0, 42.0
SL1, SR1 = CX1 - NW1 / 2, CX1 + NW1 / 2

# Vertical coordinates for Figure 1 nodes (shifted up by 50 units since title banner is removed)
DY1 = -50.0
Y1 = {
    "acc": 190.0 + DY1, "heat": 254.0 + DY1, "p2m": 318.0 + DY1,
    "hyd": 420.0 + DY1, "sd": 480.0 + DY1, "th": 540.0 + DY1,
    "mk": 654.0 + DY1, "vent": 718.0 + DY1, "adv": 782.0 + DY1, "out": 846.0 + DY1,
}
CONT1 = (480.0, 140.0 + DY1, 340.0, 852.0)
ITER1 = (506.0, 380.0 + DY1, 288.0, 242.0)
DIA_Y1, DIA_HW1, DIA_HH1 = 936.0 + DY1, 74.0, 32.0
LOOP_X1, RETRY_X1, TOP_Y1 = 494.0, 806.0, 158.0 + DY1
FS_NODE1, FS_TITLE1, FS_BADGE1, FS_EQ1, FS_PILL1 = 14.5, 16.0, 11.5, 14.0, 13.0

LEFT_X1, RIGHT_X1, BW1, BH1 = 24.0, 856.0, 420.0, 170.0
LE1, RE1 = 462.0, 838.0

LAYOUT1 = {
    "acc": (LEFT_X1, 120.0 + DY1, [[(SL1, Y1["acc"] + NH1 / 2), (LEFT_X1 + BW1, Y1["acc"] + NH1 / 2)]]),
    "heat": (LEFT_X1, 308.0 + DY1, [[(SL1, Y1["heat"] + NH1 / 2), (LE1, Y1["heat"] + NH1 / 2), (LE1, 372.0 + DY1), (LEFT_X1 + BW1, 372.0 + DY1)]]),
    "hyd": (LEFT_X1, 496.0 + DY1, [[(SL1, Y1["hyd"] + NH1 / 2), (LE1, Y1["hyd"] + NH1 / 2), (LE1, 560.0 + DY1), (LEFT_X1 + BW1, 560.0 + DY1)]]),
    "vent": (LEFT_X1, 684.0 + DY1, [[(SL1, Y1["vent"] + NH1 / 2), (LEFT_X1 + BW1, Y1["vent"] + NH1 / 2)]]),
    "p2m": (RIGHT_X1, 120.0 + DY1, [[(SR1, Y1["p2m"] + NH1 / 2), (RE1, Y1["p2m"] + NH1 / 2), (RE1, 250.0 + DY1), (RIGHT_X1, 250.0 + DY1)]]),
    "sd": (RIGHT_X1, 308.0 + DY1, [[(SR1, Y1["sd"] + NH1 / 2), (RE1, Y1["sd"] + NH1 / 2), (RE1, 440.0 + DY1), (RIGHT_X1, 440.0 + DY1)]]),
    "th": (RIGHT_X1, 496.0 + DY1, [[(SR1, Y1["th"] + NH1 / 2), (RIGHT_X1, Y1["th"] + NH1 / 2)]]),
    "mk": (RIGHT_X1, 684.0 + DY1, [
        [(SR1, Y1["mk"] + NH1 / 2), (RE1, Y1["mk"] + NH1 / 2), (RE1, 730.0 + DY1), (RIGHT_X1, 730.0 + DY1)],
        [(SR1, Y1["adv"] + NH1 / 2), (RIGHT_X1, Y1["adv"] + NH1 / 2)],
    ]),
}

COLS1 = {
    "acc": C_LAGOON,
    "heat": C_OCHRE,
    "hyd": C_FERN,
    "vent": C_IRIS,
    "p2m": C_SIGNAL,
    "sd": C_MAGMA,
    "th": C_GARNET,
    "mk": C_INK,
}


def line_halo(lw: float):
    """Return path effects for white outline halo behind strokes.

    Parameters
    ----------
    lw : float
        Line width of core stroke.

    Returns
    -------
    list
        List of matplotlib path effect objects.
    """
    return [pe.Stroke(linewidth=lw + 7, foreground=BG), pe.Normal()]


def draw_polyline(ax, pts, color, lw=2.2, ls="-", z=6, bridge=True):
    """Draw a polyline across coordinates with optional background bridge halo.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    pts : list of tuple
        Sequence of (x, y) coordinates.
    color : str
        Stroke color.
    lw : float, optional
        Line width, by default 2.2.
    ls : str or tuple, optional
        Line style, by default "-".
    z : int, optional
        Layer ordering, by default 6.
    bridge : bool, optional
        Whether to apply background halo effect, by default True.

    Returns
    -------
    None
    """
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="butt",
            path_effects=line_halo(lw) if bridge else None, zorder=z)


def render_figure_1(out_path: Path):
    """Render Figure 1: Erebus execution pipeline and equation cards.

    Parameters
    ----------
    out_path : pathlib.Path
        Target SVG output file path.

    Returns
    -------
    None
    """
    fig, ax = new_equal_canvas(15.34, 15.34 * H1 / W1, W1, H1)

    # Spine execution boxes
    def node(key, label, col):
        y = Y1[key]
        draw_box(ax, SL1, y, NW1, NH1, facecolor=col["wash"], edgecolor=col["accent"], lw=2.0, radius=10, zorder=8)
        ax.text(CX1, y + NH1 / 2 + 1, label, fontsize=FS_NODE1, fontweight="medium", color=col["text"], ha="center", va="center", zorder=9)

    def down(y1, y2):
        draw_arrow(ax, CX1, y1, CX1, y2, color=INK, lw=2.0, scale=13, zorder=8)

    draw_box(ax, SL1, 80.0 + DY1, NW1, 44, facecolor=BG, edgecolor=INK, lw=2.0, radius=10, zorder=8)
    ax.text(CX1, 103.0 + DY1, "Setup & Restart", fontsize=FS_NODE1, fontweight="medium", color=INK, ha="center", va="center", zorder=9)
    down(124.0 + DY1, Y1["acc"])

    x, y, w, h = CONT1
    draw_box(ax, x, y, w, h, facecolor=SURFACE, edgecolor=INK, lw=2.2, radius=18, zorder=1)
    ax.text(x + w - 14, y + h - 18, "simulation_loop", fontsize=11.5, color=MUTED, ha="right", va="center", fontfamily="Geist Mono", zorder=3)

    node("acc", "Ambient & Accretion", C_LAGOON)
    down(Y1["acc"] + NH1, Y1["heat"])
    node("heat", "Heating & Pyrolysis", C_OCHRE)
    down(Y1["heat"] + NH1, Y1["p2m"])
    node("p2m", "P2M, Melt & Gravity", C_SIGNAL)
    down(Y1["p2m"] + NH1, ITER1[1])

    ix, iy, iw, ih = ITER1
    draw_box(ax, ix, iy, iw, ih, facecolor=C_MAGMA["wash"], edgecolor=C_MAGMA["accent"], lw=2.0, radius=14, zorder=2)
    ax.text(CX1, iy + 20, r"Outer iteration ($\leq n_\mathrm{iter}$)", fontsize=13.0, fontweight="medium", color=C_MAGMA["text"], ha="center", va="center", zorder=4)
    node("hyd", "Hydration Reaction", C_FERN)
    down(Y1["hyd"] + NH1, Y1["sd"])
    node("sd", "Stokes-Darcy", C_MAGMA)
    down(Y1["sd"] + NH1, Y1["th"])
    node("th", "Thermal & Segregation", C_GARNET)
    ax.text(CX1, 604.0 + DY1, r"not converged: restore, $\Delta t \leftarrow \Delta t/2$", fontsize=12.0, color=C_MAGMA["text"], fontstyle="italic", ha="center", va="center", zorder=4)
    down(iy + ih, Y1["mk"])

    node("mk", "Marker Stress & Diffusion", C_INK)
    down(Y1["mk"] + NH1, Y1["vent"])
    node("vent", "Venting & Atmosphere", C_IRIS)
    down(Y1["vent"] + NH1, Y1["adv"])
    node("adv", "Advect & Telescope", C_INK)
    down(Y1["adv"] + NH1, Y1["out"])

    yout = Y1["out"]
    pts = [(SL1 + 14, yout), (SR1 + 6, yout), (SR1 - 14, yout + NH1), (SL1 - 6, yout + NH1)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=C_GARNET["wash"], edgecolor=C_GARNET["accent"], lw=2.0, zorder=8))
    ax.text(CX1, yout + NH1 / 2 + 1, "Diagnostics & Output", fontsize=FS_NODE1, fontweight="medium", color=C_GARNET["text"], ha="center", va="center", zorder=9)
    down(yout + NH1, DIA_Y1 - DIA_HH1)

    d = [(CX1, DIA_Y1 - DIA_HH1), (CX1 + DIA_HW1, DIA_Y1), (CX1, DIA_Y1 + DIA_HH1), (CX1 - DIA_HW1, DIA_Y1)]
    ax.add_patch(Polygon(d, closed=True, facecolor=BG, edgecolor=INK, lw=2.0, zorder=8))
    ax.text(CX1, DIA_Y1 + 1, r"$t > t_\mathrm{end}$ ?", fontsize=FS_NODE1, fontweight="medium", color=INK, ha="center", va="center", zorder=9)

    # Next step loop-back: diamond -> left rail -> top -> into Ambient & Accretion
    loop_col = INK_RAMP["700"]
    draw_polyline(ax, [(CX1 - DIA_HW1, DIA_Y1), (LOOP_X1, DIA_Y1), (LOOP_X1, TOP_Y1), (CX1 - 70, TOP_Y1), (CX1 - 70, Y1["acc"] - 12)],
                  loop_col, lw=2.2, z=5, bridge=False)
    draw_arrow(ax, CX1 - 70, Y1["acc"] - 14, CX1 - 70, Y1["acc"], color=loop_col, lw=2.2, scale=13, zorder=8)
    ax.text(537, DIA_Y1 - 15, "No", fontsize=13.0, fontweight="medium", color=loop_col, ha="center", va="center", zorder=9)
    ax.text(537, DIA_Y1 + 17, r"$t \leftarrow t+\Delta t$", fontsize=12.5, fontweight="medium", color=loop_col, ha="center", va="center", zorder=9)

    # Non-convergence retry: iteration box -> right rail -> top -> into Ambient & Accretion
    rc = C_MAGMA["accent"]
    draw_polyline(ax, [(ix + iw, 604.0 + DY1), (RETRY_X1, 604.0 + DY1), (RETRY_X1, TOP_Y1), (CX1 + 70, TOP_Y1), (CX1 + 70, Y1["acc"] - 12)],
                  rc, lw=2.2, ls=(0, (5, 3)), z=5, bridge=False)
    draw_arrow(ax, CX1 + 70, Y1["acc"] - 14, CX1 + 70, Y1["acc"], color=rc, lw=2.2, scale=13, zorder=8)
    ax.text((CX1 + 70 + RETRY_X1) / 2, TOP_Y1 + 17, "retry", fontsize=12.0, fontweight="medium", color=C_MAGMA["text"], ha="center", va="center", zorder=9)

    draw_arrow(ax, CX1, DIA_Y1 + DIA_HH1, CX1, 1022.0 + DY1, color=C_FERN["solid"], lw=2.4, scale=13, zorder=8)
    ax.text(CX1 + 14, 1006.0 + DY1, "Yes", fontsize=13.0, fontweight="medium", color=C_FERN["text"], ha="left", va="center", zorder=9)
    draw_box(ax, CX1 - 92, 1022.0 + DY1, 184, 44, facecolor=BG, edgecolor=INK, lw=2.0, radius=22, zorder=8)
    ax.text(CX1, 1044.0 + DY1, "End of Run", fontsize=FS_NODE1, fontweight="medium", color=INK, ha="center", va="center", zorder=9)

    # Connectors between cards and spine
    for key, (_, _, paths) in LAYOUT1.items():
        for p in paths:
            draw_polyline(ax, p, COLS1[key]["accent"], lw=2.2, z=6)

    # Legend at bottom left
    y0 = 1010.0 + DY1
    draw_box(ax, 24, y0, 420, 72, facecolor=BG, edgecolor=AXIS, lw=1.6, radius=12, zorder=2)
    draw_box(ax, 40, y0 + 18, 82, 36, facecolor=BG, edgecolor=INK, lw=1.6, radius=8, zorder=3)
    ax.text(81, y0 + 37, "Step", fontsize=12.5, fontweight="medium", color=INK, ha="center", va="center", zorder=4)
    draw_polyline(ax, [(140, y0 + 36), (228, y0 + 36)], C_MAGMA["accent"], lw=2.2, ls=(0, (5, 3)), z=3, bridge=False)
    ax.text(184, y0 + 22, "retry", fontsize=11.5, fontweight="medium", color=C_MAGMA["text"], ha="center", va="center", zorder=4)
    ax.add_patch(Polygon([(296, y0 + 14), (342, y0 + 36), (296, y0 + 58), (250, y0 + 36)], closed=True,
                         facecolor=BG, edgecolor=INK, lw=1.6, zorder=3))
    ax.text(296, y0 + 37, "Check", fontsize=11.5, fontweight="medium", color=INK, ha="center", va="center", zorder=4)
    ax.add_patch(Polygon([(368, y0 + 18), (430, y0 + 18), (420, y0 + 54), (358, y0 + 54)], closed=True,
                         facecolor=C_GARNET["wash"], edgecolor=C_GARNET["accent"], lw=1.6, zorder=3))
    ax.text(394, y0 + 37, "I/O", fontsize=12.5, fontweight="medium", color=C_GARNET["text"], ha="center", va="center", zorder=4)

    # 8 Equation Cards
    cards1 = {
        "acc": ("Accretion & Ambient", "accretion.jl",
                [r"$w_\mathrm{3D} = 2r$", r"$\int w_\mathrm{3D}\,\mathrm{d}A = \frac{4}{3}\pi R^3$"],
                ["Pebble Accretion", "Sticky Air"]),
        "heat": ("Heating & Pyrolysis", "radioactivity.jl",
                 [r"$Q_i = \frac{f_i r_i E_i}{\tau_i} e^{-t/\tau_i}$", r"$i \in \{^{26}\mathrm{Al},\, ^{60}\mathrm{Fe}\}$"],
                 ["Radiogenic Decay", "C, N, H Pyrolysis"]),
        "hyd": ("Hydration Reaction", "serpentinization.jl",
                [r"$\mathrm{W_{sil}} \rightleftharpoons \mathrm{D_{sil}} + \mathrm{H_2O}$",
                 r"$\Delta G = \Delta H - T\Delta S + P_\mathrm{f}\Delta V + RT\ln(X_\mathrm{D}/X_\mathrm{W})$"],
                ["Hydration", "Dehydration"]),
        "vent": ("Venting & Atmosphere", "atmosphere.jl",
                 [r"$T^4(\tau) = \frac{3}{4}T_\mathrm{int}^4\left(\frac{2}{3}+\tau\right) + \frac{3}{4}T_\mathrm{eq}^4\,f_\gamma(\tau)$",
                  r"$f_\gamma(\tau) = \frac{2}{3} + \frac{1}{\gamma\sqrt{3}} + \left(\frac{\gamma}{\sqrt{3}} - \frac{1}{\gamma\sqrt{3}}\right) e^{-\gamma\tau\sqrt{3}}$"],
                 ["H, C, N, S Degassing", "Jeans, XUV Escape"]),
        "p2m": ("P2M, Melt & Gravity", "gravity.jl",
                [r"$\nabla^2\Phi = \frac{8}{3}\pi G\rho$", r"$F_\mathrm{m} = (T - T_\mathrm{s})\,/\,(T_\mathrm{l} - T_\mathrm{s})$"],
                ["Marker to Grid", r"Surface $\sigma T^4$"]),
        "sd": ("Stokes-Darcy", "stokes_darcy.jl",
               [r"$\nabla\cdot\boldsymbol{\sigma}' - \nabla P + \rho\mathbf{g} = 0$",
                r"$\nabla\cdot\mathbf{v}_\mathrm{s} = -D\ln(1-\phi)/Dt$"],
               [r"Darcy $\mathbf{q}_\mathrm{D}$", "Plastic Yielding"]),
        "th": ("Thermal & Segregation", "thermal.jl",
               [r"$\rho c_p\, \partial T/\partial t = \nabla\cdot(k_\mathrm{eff}\nabla T)$",
                r"$+\,H_\mathrm{r} + H_\mathrm{s} + H_\mathrm{a} + Q_\mathrm{seg} + H_\mathrm{L}$"],
               ["Metal & Melt Seg.", r"Soft-Turb. $k_\mathrm{eff}$"]),
        "mk": ("Markers & Domain", "advect_markers.jl",
               [r"$\mathrm{d}\mathbf{x}_m/\mathrm{d}t = \mathbf{v}(\mathbf{x}_m)$", "4th-order Runge-Kutta"],
               ["Subgrid Diffusion", "Domain Doubling"]),
    }

    for key, (x_card, y_card, _) in LAYOUT1.items():
        col = COLS1[key]
        title, mod, eqs, pills = cards1[key]
        draw_box(ax, x_card, y_card, BW1, BH1, facecolor=col["wash"], edgecolor=col["accent"], lw=2.4, radius=16, zorder=2)
        ax.text(x_card + 18, y_card + 26, title, fontsize=FS_TITLE1, fontweight="semibold", color=col["text"], ha="left", va="center", zorder=4)
        bw = len(mod) * 8.6 + 20
        draw_box(ax, x_card + BW1 - bw - 14, y_card + 12, bw, 28, facecolor=BG, edgecolor=col["border"], lw=1.2, radius=6, zorder=3)
        ax.text(x_card + BW1 - bw / 2 - 14, y_card + 26, mod, fontsize=FS_BADGE1, color=MUTED, fontfamily="Geist Mono", ha="center", va="center", zorder=4)
        draw_box(ax, x_card + 14, y_card + 48, BW1 - 28, 70, facecolor=BG, edgecolor=col["border"], lw=1.5, radius=10, zorder=3)
        for j, eq in enumerate(eqs):
            ax.text(x_card + BW1 / 2, y_card + 68 + j * 30, eq, fontsize=FS_EQ1, fontweight="normal", color=INK, ha="center", va="center", zorder=4)
        pw = (BW1 - 28 - 12) / 2
        for i, lbl in enumerate(pills):
            px = x_card + 14 + pw / 2 + i * (pw + 12)
            draw_pill(ax, px, y_card + 144, pw, 36, lbl, facecolor=BG, edgecolor=col["accent"], text_color=col["text"], fontsize=FS_PILL1)

    fig.savefig(out_path, format="svg", bbox_inches=None)
    plt.close(fig)


# =============================================================================
# Figure 2: Planetesimal Evolution in Three Stages (2A-1)
# =============================================================================

FS_TITLE2, FS_PILL2, FS_LAYER2, FS_FLUX2 = 15.5, 13.0, 13.0, 15.0


def P(cx: float, cy: float, r: float, th: float):
    """Compute Cartesian point from polar coordinates.

    Parameters
    ----------
    cx : float
        Center x coordinate.
    cy : float
        Center y coordinate.
    r : float
        Radius from center.
    th : float
        Angle in degrees (270 degrees is top).

    Returns
    -------
    tuple of float
        (x, y) coordinates.
    """
    t = np.deg2rad(th)
    return cx + r * np.cos(t), cy + r * np.sin(t)


def tangent_rot(th: float) -> float:
    """Compute tangent rotation angle for text along circular arcs.

    Parameters
    ----------
    th : float
        Angle in degrees.

    Returns
    -------
    float
        Readable rotation angle in degrees.
    """
    rot = (270.0 - th) % 180.0
    return rot - 180.0 if rot > 90.0 else rot


def polar_label(ax, c, r, th, s, fs, color, bg, tangent=True, **kw):
    """Render polar label text with background stroke halo.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    c : tuple of float
        Center coordinates (cx, cy).
    r : float
        Radius from center.
    th : float
        Angle in degrees.
    s : str
        Label string.
    fs : float
        Font size.
    color : str
        Text color.
    bg : str or None
        Halo background color.
    tangent : bool, optional
        Whether to align text along arc tangent, by default True.
    **kw : dict
        Additional keyword arguments for ax.text.

    Returns
    -------
    None
    """
    x, y = P(*c, r, th)
    rot = tangent_rot(th) if tangent else 0.0
    effects = [pe.Stroke(linewidth=4, foreground=bg), pe.Normal()] if bg is not None else None
    ax.text(x, y, s, fontsize=fs, fontweight="medium", color=color, ha="center", va="center",
            rotation=rot, rotation_mode="anchor", zorder=12,
            path_effects=effects, **kw)


def radial_block(ax, c, th, r1, r2, col_face, col_edge, w=14, hw=28, hl=16, z=11):
    """Draw a radial block arrow pointing along an angle from r1 to r2.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    c : tuple of float
        Center coordinates (cx, cy).
    th : float
        Angle in degrees.
    r1 : float
        Starting radius.
    r2 : float
        Ending radius.
    col_face : str
        Fill color.
    col_edge : str
        Outline color.
    w : float, optional
        Shaft width, by default 14.
    hw : float, optional
        Head width, by default 28.
    hl : float, optional
        Head length, by default 16.
    z : int, optional
        Layer ordering, by default 11.

    Returns
    -------
    None
    """
    x1, y1 = P(*c, r1, th)
    x2, y2 = P(*c, r2, th)
    draw_block_arrow(ax, x1, y1, x2, y2, facecolor=col_face, edgecolor=col_edge,
                     width=w, head_width=hw, head_len=hl, lw=1.6, zorder=z)


def dots(ax, c, pts, face, r=6.5, z=10):
    """Draw circular marker dots at polar coordinates.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    c : tuple of float
        Center coordinates (cx, cy).
    pts : list of tuple
        List of (radius, angle_deg) tuples.
    face : str
        Fill color.
    r : float, optional
        Circle radius, by default 6.5.
    z : int, optional
        Layer ordering, by default 10.

    Returns
    -------
    None
    """
    for rr, th in pts:
        ax.add_patch(Circle(P(*c, rr, th), r, facecolor=face, edgecolor=INK, lw=1.5, zorder=z))


def convection(ax, c, r, th, size=19, color=INK, z=10):
    """Draw convection circulation swirl arrows.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    c : tuple of float
        Center coordinates (cx, cy).
    r : float
        Radius from center.
    th : float
        Angle in degrees.
    size : float, optional
        Horizontal swirl radius, by default 19.
    color : str, optional
        Arrow color, by default INK.
    z : int, optional
        Layer ordering, by default 10.

    Returns
    -------
    None
    """
    x, y = P(*c, r, th)
    draw_arrow(ax, x - size, y + 4, x + size, y + 4, color=color, lw=2.2, rad=0.9, scale=13, zorder=z)
    draw_arrow(ax, x + size, y - 4, x - size, y - 4, color=color, lw=2.2, rad=0.9, scale=13, zorder=z)


def group_box(ax, x, y, w, h, col, title, pills):
    """Draw a process grouping box containing a title and two pills.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left coordinate.
    y : float
        Top coordinate.
    w : float
        Width.
    h : float
        Height.
    col : dict
        Color token dictionary.
    title : str
        Group title.
    pills : list of str
        List of 2 pill label strings.

    Returns
    -------
    None
    """
    draw_box(ax, x, y, w, h, facecolor=col["wash"], edgecolor=col["accent"], lw=2.4, radius=16, zorder=3)
    ax.text(x + w / 2, y + 27, title, fontsize=FS_TITLE2, fontweight="semibold", color=col["text"], ha="center", va="center", zorder=4)
    gap = 10
    pw = (w - 24 - gap) / 2
    for i, s in enumerate(pills):
        draw_pill(ax, x + 12 + pw / 2 + i * (pw + gap), y + h - 32, pw, 38, s, facecolor=BG,
                  edgecolor=col["accent"], text_color=col["text"], fontsize=FS_PILL2, radius=8, lw=1.6)


def legend2(ax, x, y, w=380, h=156):
    """Render legend box for Figure 2.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left coordinate.
    y : float
        Top coordinate.
    w : float, optional
        Box width, by default 380.
    h : float, optional
        Box height, by default 156.

    Returns
    -------
    None
    """
    draw_box(ax, x, y, w, h, facecolor=BG, edgecolor=AXIS, lw=1.6, radius=14, zorder=3)
    rows = [y + 28, y + 62, y + 96, y + 130]
    for dx, dy in [(-15, 3), (0, -4), (15, 3)]:
        ax.add_patch(Circle((x + 48 + dx, rows[0] + dy), 6.5, facecolor=BG, edgecolor=INK, lw=1.5, zorder=4))
    for dx, dy in [(-13, -2), (3, 4), (16, -3)]:
        ax.add_patch(Circle((x + 48 + dx, rows[1] + dy), 6.5, facecolor=C_GARNET["accent"], edgecolor=INK, lw=1.4, zorder=4))
    draw_box(ax, x + 20, rows[2] - 13, 56, 26, facecolor=SURFACE, edgecolor=AXIS, lw=1.6, radius=7, zorder=4)
    draw_block_arrow(ax, x + 20, rows[3], x + 78, rows[3], facecolor=INK_300, edgecolor=INK,
                     width=12, head_width=24, head_len=15, zorder=4)
    for yy, s in zip(rows, ["Volatile gas (H, C, N, S)", "Fe-FeS metal", "Process group", "Energy or mass flux"]):
        ax.text(x + 98, yy + 1, s, fontsize=13.5, color=TEXT, ha="left", va="center", zorder=4)


def render_figure_2(out_path: Path):
    """Render Figure 2: Planetesimal evolution across three physical stages.

    Parameters
    ----------
    out_path : pathlib.Path
        Target SVG output file path.

    Returns
    -------
    None
    """
    W2, H2 = 1400.0, 960.0
    fig, ax = new_equal_canvas(15.4, 15.4 * H2 / W2, W2, H2)

    # Shift planetary center and boxes up by 30 units since title banner is removed
    DY2 = -30.0
    c = (700.0, 500.0 + DY2)
    R_HALO, R_OUT, R_MID, R_IN = 262.0, 232.0, 180.0, 112.0

    ax.add_patch(Circle(c, R_HALO, facecolor=SURFACE, edgecolor=LINE, lw=1.8, zorder=1))

    # Sector I: cold accretion (150..270)
    sig = C_SIGNAL
    ax.add_patch(Wedge(c, R_OUT, 150, 270, width=R_OUT - R_MID, facecolor=sig["mid"], edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_MID, 150, 270, facecolor=sig["wash"], edgecolor=INK, lw=1.6, zorder=2))
    polar_label(ax, c, (R_OUT + R_MID) / 2, 222, "Accreting Shell", FS_LAYER2, sig["text"], sig["mid"])
    polar_label(ax, c, 112, 212, "Porous Ice + Rock", FS_LAYER2, sig["text"], sig["wash"], tangent=False)
    lag = C_LAGOON
    for th in (196, 248):
        radial_block(ax, c, th, 318, R_OUT + 4, lag["mid"], lag["accent"])
    polar_label(ax, c, 330, 261, r"$\dot{M}_\mathrm{acc}$", FS_FLUX2, lag["text"], BG, tangent=False)

    # Sector II: hydrothermal alteration (270..390)
    fern, och = C_FERN, C_OCHRE
    ax.add_patch(Wedge(c, R_OUT, 270, 390, width=R_OUT - R_MID, facecolor=fern["mid"], edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_MID, 270, 390, width=R_MID - R_IN, facecolor=fern["tint"], edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_IN, 270, 390, facecolor=och["tint"], edgecolor=INK, lw=1.6, zorder=2))
    polar_label(ax, c, (R_OUT + R_MID) / 2, 295, "Crust", FS_LAYER2, fern["text"], fern["mid"])
    polar_label(ax, c, (R_MID + R_IN) / 2, 292, "Serpentinization", FS_LAYER2, fern["text"], fern["tint"])
    polar_label(ax, c, 60, 330, "Dehydration", FS_LAYER2 - 1.0, och["text"], och["tint"], tangent=False)
    # Darcy percolation vector confined within permeable solid crust
    radial_block(ax, c, 336, R_IN + 20, R_OUT - 12, fern["tint"], fern["text"])
    polar_label(ax, c, (R_IN + 20 + R_OUT - 12) / 2, 347, r"$\mathbf{q}_\mathrm{D}$", FS_FLUX2, fern["text"], fern["tint"], tangent=False)
    convection(ax, c, (R_MID + R_IN) / 2, 370, color=fern["text"])

    # Sector III: magma ocean and core (30..150)
    mag, gar, iris = C_MAGMA, C_GARNET, C_IRIS
    R_CORE, R_MO, R_CR = 100.0, 204.0, 232.0
    ax.add_patch(Wedge(c, R_HALO, 30, 150, width=R_HALO - R_CR, facecolor=iris["tint"], edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_CR, 30, 150, width=R_CR - R_MO, facecolor=SUB, edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_MO, 30, 150, width=R_MO - R_CORE, facecolor=mag["mid"], edgecolor=INK, lw=1.6, zorder=2))
    ax.add_patch(Wedge(c, R_CORE, 30, 150, facecolor=gar["accent"], edgecolor=INK, lw=1.8, zorder=2))
    polar_label(ax, c, 52, 90, "Fe-FeS Core", FS_LAYER2, BG, gar["accent"], tangent=False)
    polar_label(ax, c, 152, 90, "Magma Ocean", FS_LAYER2 + 0.5, mag["text"], mag["mid"], tangent=False)
    polar_label(ax, c, (R_CR + R_MO) / 2, 134, "Crust", 11.5, MUTED, SUB)
    polar_label(ax, c, (R_HALO + R_CR) / 2, 136, "Atmosphere", 11.5, iris["text"], iris["tint"])
    convection(ax, c, 152, 126, color=mag["text"])
    dots(ax, c, [(118, 62), (132, 70), (120, 78)], gar["accent"])
    # Metal settling drift-flux vector descending toward core
    radial_block(ax, c, 50, 150, R_CORE + 2, gar["tint"], gar["accent"])
    polar_label(ax, c, 178, 38, r"$\mathbf{J}_\mathrm{Fe}$", FS_FLUX2 - 1, gar["text"], mag["mid"], tangent=False)
    dots(ax, c, [(186, 100), (182, 82)], BG)
    dots(ax, c, [(247, 84), (247, 98), (247, 112)], BG, r=6.0)
    # Surface radiative heat flux radiating into space
    radial_block(ax, c, 62, R_CR, 282, mag["tint"], mag["text"])
    polar_label(ax, c, 304, 54, r"$F_\mathrm{rad}$", FS_FLUX2, mag["text"], BG, tangent=False)

    for th in (30, 150, 270):
        x2, y2 = P(*c, R_HALO, th)
        ax.plot([c[0], x2], [c[1], y2], color=INK, lw=2.6, zorder=8)

    for x_pill, y_pill, col, s in [(470, 112.0 + DY2, sig, "I. Cold Accretion"),
                                  (930, 112.0 + DY2, fern, "II. Hydrothermal Alteration"),
                                  (700, 812.0 + DY2, mag, "III. Magma Ocean & Core")]:
        draw_pill(ax, x_pill, y_pill, 330, 42, s, facecolor=col["wash"], edgecolor=col["accent"],
                  text_color=col["text"], fontsize=14.0, radius=14)

    # 8 Process group boxes
    L, R, BW, BH = 24, 1051, 325, 112
    group_box(ax, L, 175.0 + DY2, BW, BH, lag, "Accretion & Disk", ["Pebbles", r"Disk $T_\mathrm{amb}$"])
    group_box(ax, L, 335.0 + DY2, BW, BH, sig, "Porous Ice & Rock", ["Accreted Ices", "Pore Sealing"])
    group_box(ax, L, 495.0 + DY2, BW, BH, gar, "Core Formation", ["Segregation", r"Partition $D_i$"])
    group_box(ax, L, 655.0 + DY2, BW, BH, iris, "Volatiles & Redox", ["H, C, N, S", "IW / QFM"])
    group_box(ax, R, 175.0 + DY2, BW, BH, fern, "Hydrothermal", [r"Darcy $\mathbf{q}_\mathrm{D}$", "Hydrofracture"])
    group_box(ax, R, 335.0 + DY2, BW, BH, och, "Heating & Pyrolysis", [r"$^{26}$Al, $^{60}$Fe", "C, N, H Pyrolysis"])
    group_box(ax, R, 495.0 + DY2, BW, BH, mag, "Magma Ocean", [r"$F_\mathrm{m} > F_\mathrm{crit}$", r"Soft-Turb. $k_\mathrm{eff}$"])
    group_box(ax, R, 655.0 + DY2, BW, BH, iris, "Atmosphere & Loss", [r"Guillot $T(\tau)$", "Jeans, XUV"])

    legend2(ax, 995, 820.0 + DY2)

    fig.savefig(out_path, format="svg", bbox_inches=None)
    plt.close(fig)


def main():
    """Parse command-line arguments and generate architecture figures.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(description="Generate Erebus.jl documentation architecture figures.")
    parser.add_argument("--fig", choices=["1", "2", "all"], default="all", help="Figure to generate (1, 2, or all)")
    parser.add_argument("--fonttype", choices=["path", "none"], default="path", help="SVG fonttype ('path' or 'none')")
    parser.add_argument("--outdir", type=Path, default=ASSETS_DIR, help="Output directory")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams["svg.fonttype"] = args.fonttype

    if args.fig in ("1", "all"):
        f1_path = args.outdir / "erebus_execution_pipeline.svg"
        print(f"Generating Figure 1 -> {f1_path} (fonttype={args.fonttype})...")
        render_figure_1(f1_path)
        sz = f1_path.stat().st_size
        print(f"Figure 1 generated: {sz} bytes ({sz / 1024:.1f} KB)")

    if args.fig in ("2", "all"):
        f2_path = args.outdir / "erebus_planetesimal_stages.svg"
        print(f"Generating Figure 2 -> {f2_path} (fonttype={args.fonttype})...")
        render_figure_2(f2_path)
        sz = f2_path.stat().st_size
        print(f"Figure 2 generated: {sz} bytes ({sz / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
