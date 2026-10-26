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
from matplotlib import patheffects as pe
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Wedge

REPO_ROOT = Path(__file__).resolve().parent.parent
ASSETS_DIR = REPO_ROOT / "docs" / "src" / "assets"


def hex_to_rgb(hex_str: str) -> tuple[float, float, float]:
    """Convert hex colour string to sRGB tuple with components in [0, 1].

    Parameters
    ----------
    hex_str : str
        Hexadecimal color string (e.g., '#ffffff').

    Returns
    -------
    tuple of float
        (r, g, b) color tuple.
    """
    h = hex_str.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))


def rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    """Convert sRGB tuple with components in [0, 1] to hex colour string.

    Parameters
    ----------
    rgb : tuple of float
        (r, g, b) color tuple.

    Returns
    -------
    str
        Hexadecimal color string.
    """
    return f"#{int(round(rgb[0] * 255)):02x}{int(round(rgb[1] * 255)):02x}{int(round(rgb[2] * 255)):02x}"


def mix(a: str, b: str, t: float) -> str:
    """Return the colour a fraction t of the way from a to b in sRGB.

    Parameters
    ----------
    a : str
        Starting hex color.
    b : str
        Ending hex color.
    t : float
        Interpolation weight in [0, 1].

    Returns
    -------
    str
        Interpolated hex color.
    """
    rgb_a = hex_to_rgb(a)
    rgb_b = hex_to_rgb(b)
    mixed = tuple(x + (y - x) * t for x, y in zip(rgb_a, rgb_b, strict=True))
    return rgb_to_hex(mixed)


def text_safe(
    color: str, ink: str, backgrounds: tuple[str, ...], min_contrast: float = 4.5
) -> str:
    """Shift color toward ink until contrast ratio meets threshold on all backgrounds.

    Parameters
    ----------
    color : str
        Initial foreground color.
    ink : str
        Target high-contrast ink color.
    backgrounds : tuple of str
        Collection of background colors to test against.
    min_contrast : float, optional
        Minimum acceptable contrast ratio, by default 4.5.

    Returns
    -------
    str
        Contrast-safe hex color.
    """
    for step in range(101):
        candidate = mix(color, ink, step / 100.0)
        if all(fvl.contrast(candidate, bg) >= min_contrast for bg in backgrounds):
            return candidate
    return ink


def resolve_tokens(mode: str = "light") -> dict:
    """Resolve colour tokens from fwl_visual_language without private imports.

    Parameters
    ----------
    mode : str, optional
        Palette mode ('light' or 'dark'), by default 'light'.

    Returns
    -------
    dict
        Resolved design token dictionary.
    """
    r = fvl.roles("paper" if mode == "light" else "screen")
    bg, ink = r["background"], r["text"]
    sub = mix(bg, ink, 0.045 if mode == "light" else 0.06)
    hover = mix(bg, ink, 0.09)
    raw_palette = dict(fvl.palette("paper" if mode == "light" else "screen"))
    hues = {n: c for n, c in raw_palette.items() if n != "ink"}
    washes = {n: mix(bg, c, 0.10) for n, c in hues.items()}
    surfaces = (bg, sub, hover, *washes.values())
    safe = {n: text_safe(c, ink, surfaces) for n, c in hues.items()}
    return {
        "bg": bg,
        "sub": sub,
        "ink": ink,
        "muted": text_safe(r["muted"], ink, surfaces),
        "line": r["line"],
        "hover": hover,
        "link": text_safe(r["link"], ink, surfaces),
        "hues": hues,
        "washes": washes,
        "surfaces": surfaces,
        "safe": safe,
    }


def get_palette_hues(
    tok: dict,
    ink_ramp: dict,
    bg: str,
    ink: str,
    muted: str,
    surface: str,
    sub: str,
    line: str,
    axis: str,
) -> dict:
    """Return color dictionaries for named palette hues.

    Parameters
    ----------
    tok : dict
        Resolved tokens.
    ink_ramp : dict
        Ink ramp tokens.
    bg : str
        Background color.
    ink : str
        Text color.
    muted : str
        Muted color.
    surface : str
        Surface color.
    sub : str
        Subdued color.
    line : str
        Line color.
    axis : str
        Axis color.

    Returns
    -------
    dict
        Mapping of hue names to color dictionaries.
    """

    def c_hue(name: str) -> dict:
        acc = tok["hues"][name]
        return {
            "accent": acc,
            "solid": text_safe(acc, ink, [bg]),
            "text": tok["safe"][name],
            "wash": mix(bg, acc, 0.10),
            "tint": mix(bg, acc, 0.24),
            "mid": mix(bg, acc, 0.46),
            "border": acc,
        }

    c_ink = {
        "accent": muted,
        "solid": ink_ramp["600"],
        "text": ink_ramp["700"],
        "wash": surface,
        "tint": sub,
        "mid": line,
        "border": axis,
    }

    return {
        "lagoon": c_hue("lagoon"),
        "signal": c_hue("signal"),
        "fern": c_hue("fern"),
        "magma": c_hue("magma"),
        "garnet": c_hue("garnet"),
        "iris": c_hue("iris"),
        "ochre": c_hue("ochre"),
        "ink": c_ink,
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
        Canvas user coordinate system width.
    h_units : float
        Canvas user coordinate system height.

    Returns
    -------
    tuple
        (fig, ax) matplotlib figure and axes.
    """
    fig = plt.figure(figsize=(width_in, height_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w_units)
    ax.set_ylim(h_units, 0)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return fig, ax


def draw_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    facecolor: str,
    edgecolor: str,
    lw: float = 1.6,
    radius: float = 10.0,
    zorder: int = 3,
):
    """Draw a rounded rectangular box patch.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left coordinate.
    y : float
        Top coordinate.
    w : float
        Box width.
    h : float
        Box height.
    facecolor : str
        Fill color.
    edgecolor : str
        Outline stroke color.
    lw : float, optional
        Line width, by default 1.6.
    radius : float, optional
        Corner rounding radius, by default 10.0.
    zorder : int, optional
        Layer ordering, by default 3.

    Returns
    -------
    matplotlib.patches.FancyBboxPatch
        The added patch.
    """
    p = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        mutation_aspect=1.0,
        zorder=zorder,
    )
    ax.add_patch(p)
    return p


def draw_pill(
    ax,
    cx: float,
    cy: float,
    w: float,
    h: float,
    text: str,
    facecolor: str,
    edgecolor: str,
    text_color: str,
    fontsize: float = 12.0,
    radius: float = 8.0,
    lw: float = 1.4,
    zorder: int = 4,
):
    """Draw a pill badge with centered label text.

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
        Pill label string.
    facecolor : str
        Fill color.
    edgecolor : str
        Outline stroke color.
    text_color : str
        Text label color.
    fontsize : float, optional
        Font size, by default 12.0.
    radius : float, optional
        Corner radius, by default 8.0.
    lw : float, optional
        Line width, by default 1.4.
    zorder : int, optional
        Layer ordering, by default 4.

    Returns
    -------
    None
    """
    draw_box(
        ax,
        cx - w / 2,
        cy - h / 2,
        w,
        h,
        facecolor=facecolor,
        edgecolor=edgecolor,
        lw=lw,
        radius=radius,
        zorder=zorder,
    )
    ax.text(
        cx,
        cy + 1,
        text,
        fontsize=fontsize,
        fontweight="medium",
        color=text_color,
        ha="center",
        va="center",
        zorder=zorder + 1,
    )


def draw_arrow(
    ax,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    color: str,
    lw: float = 2.0,
    rad: float = 0.0,
    scale: float = 12.0,
    zorder: int = 5,
):
    """Draw a curved or straight line arrow with head.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x1 : float
        Start x coordinate.
    y1 : float
        Start y coordinate.
    x2 : float
        End x coordinate.
    y2 : float
        End y coordinate.
    color : str
        Stroke and head color.
    lw : float, optional
        Line width, by default 2.0.
    rad : float, optional
        Curvature radius, by default 0.0.
    scale : float, optional
        Arrowhead scale factor, by default 12.0.
    zorder : int, optional
        Layer ordering, by default 5.

    Returns
    -------
    matplotlib.patches.FancyArrowPatch
        The added patch.
    """
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle=f"-|>,head_length={scale * 0.7:.1f},head_width={scale * 0.45:.1f}",
        connectionstyle=f"arc3,rad={rad}",
        color=color,
        linewidth=lw,
        zorder=zorder,
    )
    ax.add_patch(arrow)
    return arrow


def draw_block_arrow(
    ax,
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    facecolor: str,
    edgecolor: str,
    width: float = 12.0,
    head_width: float = 24.0,
    head_len: float = 14.0,
    lw: float = 1.6,
    zorder: int = 6,
):
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
    facecolor : str
        Fill color.
    edgecolor : str
        Border line color.
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
        return None
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
    poly = Polygon(
        pts,
        closed=True,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=lw,
        zorder=zorder,
    )
    ax.add_patch(poly)
    return poly


def line_halo(lw: float, bg: str):
    """Return path effects for white outline halo behind strokes.

    Parameters
    ----------
    lw : float
        Line width of core stroke.
    bg : str
        Background color.

    Returns
    -------
    list
        List of matplotlib path effect objects.
    """
    return [pe.Stroke(linewidth=lw + 7, foreground=bg), pe.Normal()]


def draw_polyline(ax, pts, color, bg, lw=2.2, ls="-", z=6, bridge=True):
    """Draw a polyline across coordinates with optional background bridge halo.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    pts : list of tuple
        Sequence of (x, y) coordinates.
    color : str
        Stroke color.
    bg : str
        Background color for halo.
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
    ax.plot(
        xs,
        ys,
        color=color,
        lw=lw,
        ls=ls,
        solid_capstyle="butt",
        path_effects=line_halo(lw, bg) if bridge else None,
        zorder=z,
    )


# =============================================================================
# Figure 1: Execution Pipeline & Equation Cards
# =============================================================================


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
    fwl_mpl.use("light")
    tok = resolve_tokens("light")
    paper_roles = fvl.roles("paper")
    ink_ramp = fvl.ramp("ink")

    bg = paper_roles["background"]
    surface = paper_roles["surface"]
    sub = ink_ramp["100"]
    line = paper_roles["line"]
    axis = paper_roles["axis"]
    muted = paper_roles["muted"]
    ink = paper_roles["heading"]

    pal = get_palette_hues(tok, ink_ramp, bg, ink, muted, surface, sub, line, axis)

    w1, h1 = 1300.0, 1040.0
    cx1 = 650.0
    nw1, nh1 = 264.0, 40.0
    sl1, sr1 = cx1 - nw1 / 2, cx1 + nw1 / 2

    # Vertical coordinates for Figure 1 nodes
    y1 = {
        "acc": 140.0,
        "heat": 202.0,
        "p2m": 264.0,
        "hyd": 358.0,
        "sd": 418.0,
        "th": 478.0,
        "mk": 594.0,
        "vent": 656.0,
        "adv": 718.0,
        "out": 816.0,
    }
    cont1 = (480.0, 92.0, 340.0, 856.0)
    iter1 = (504.0, 324.0, 292.0, 238.0)
    dia_y1, dia_hw1, dia_hh1 = 896.0, 92.0, 28.0
    loop_x1, retry_x1, top_y1 = 494.0, 810.0, 110.0
    fs_node1, fs_title1, fs_badge1, fs_eq1, fs_pill1 = 13.5, 14.5, 11.5, 13.0, 12.0

    left_x1, right_x1, bw1, bh1 = 24.0, 856.0, 420.0, 168.0
    le1, re1 = 462.0, 838.0

    layout1 = {
        "acc": (
            left_x1,
            70.0,
            [[(sl1, y1["acc"] + nh1 / 2), (left_x1 + bw1, y1["acc"] + nh1 / 2)]],
        ),
        "heat": (
            left_x1,
            258.0,
            [
                [
                    (sl1, y1["heat"] + nh1 / 2),
                    (le1, y1["heat"] + nh1 / 2),
                    (le1, 322.0),
                    (left_x1 + bw1, 322.0),
                ]
            ],
        ),
        "hyd": (
            left_x1,
            446.0,
            [
                [
                    (sl1, y1["hyd"] + nh1 / 2),
                    (le1, y1["hyd"] + nh1 / 2),
                    (le1, 510.0),
                    (left_x1 + bw1, 510.0),
                ]
            ],
        ),
        "vent": (
            left_x1,
            634.0,
            [[(sl1, y1["vent"] + nh1 / 2), (left_x1 + bw1, y1["vent"] + nh1 / 2)]],
        ),
        "p2m": (
            right_x1,
            70.0,
            [
                [
                    (sr1, y1["p2m"] + nh1 / 2),
                    (re1, y1["p2m"] + nh1 / 2),
                    (re1, 200.0),
                    (right_x1, 200.0),
                ]
            ],
        ),
        "sd": (
            right_x1,
            258.0,
            [
                [
                    (sr1, y1["sd"] + nh1 / 2),
                    (re1, y1["sd"] + nh1 / 2),
                    (re1, 390.0),
                    (right_x1, 390.0),
                ]
            ],
        ),
        "th": (
            right_x1,
            446.0,
            [[(sr1, y1["th"] + nh1 / 2), (right_x1, y1["th"] + nh1 / 2)]],
        ),
        "mk": (
            right_x1,
            634.0,
            [
                [
                    (sr1, y1["mk"] + nh1 / 2),
                    (re1, y1["mk"] + nh1 / 2),
                    (re1, 680.0),
                    (right_x1, 680.0),
                ],
                [(sr1, y1["adv"] + nh1 / 2), (right_x1, y1["adv"] + nh1 / 2)],
            ],
        ),
    }

    cols1 = {
        "acc": pal["lagoon"],
        "heat": pal["ochre"],
        "hyd": pal["fern"],
        "vent": pal["iris"],
        "p2m": pal["signal"],
        "sd": pal["magma"],
        "th": pal["garnet"],
        "mk": pal["ink"],
    }

    fig, ax = new_equal_canvas(15.34, 15.34 * h1 / w1, w1, h1)

    # Spine execution boxes
    def node(key, label, col):
        y = y1[key]
        draw_box(
            ax,
            sl1,
            y,
            nw1,
            nh1,
            facecolor=col["wash"],
            edgecolor=col["accent"],
            lw=2.0,
            radius=10,
            zorder=8,
        )
        ax.text(
            cx1,
            y + nh1 / 2 + 1,
            label,
            fontsize=fs_node1,
            fontweight="medium",
            color=col["text"],
            ha="center",
            va="center",
            zorder=9,
        )

    def down(y_start, y_end):
        draw_arrow(ax, cx1, y_start, cx1, y_end, color=ink, lw=2.0, scale=13, zorder=8)

    draw_box(
        ax, sl1, 30.0, nw1, 42, facecolor=bg, edgecolor=ink, lw=2.0, radius=10, zorder=8
    )
    ax.text(
        cx1,
        52.0,
        "Setup & Restart",
        fontsize=fs_node1,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=9,
    )
    down(72.0, y1["acc"])

    x, y, w, h = cont1
    draw_box(
        ax, x, y, w, h, facecolor=surface, edgecolor=ink, lw=2.2, radius=18, zorder=1
    )
    ax.text(
        x + w - 14,
        y + h - 18,
        "simulation_loop",
        fontsize=11.5,
        color=muted,
        ha="right",
        va="center",
        fontfamily="Geist Mono",
        zorder=3,
    )

    node("acc", "1. Ambient & Accretion", pal["lagoon"])
    down(y1["acc"] + nh1, y1["heat"])
    node("heat", "2. Radiogenic Heating", pal["ochre"])
    down(y1["heat"] + nh1, y1["p2m"])
    node("p2m", "3. Marker Mapping & Gravity", pal["signal"])
    down(y1["p2m"] + nh1, iter1[1])

    ix, iy, iw, ih = iter1
    draw_box(
        ax,
        ix,
        iy,
        iw,
        ih,
        facecolor=pal["magma"]["wash"],
        edgecolor=pal["magma"]["accent"],
        lw=2.0,
        radius=14,
        zorder=2,
    )
    ax.text(
        cx1,
        iy + 16,
        r"Outer iteration ($\leq \mathrm{max\_iters}$)",
        fontsize=12.5,
        fontweight="semibold",
        color=pal["magma"]["text"],
        ha="center",
        va="center",
        zorder=4,
    )

    node("hyd", "4. Hydration Reactions", pal["fern"])
    down(y1["hyd"] + nh1, y1["sd"])
    node("sd", "5. Hydro-Mechanical Stokes-Darcy", pal["magma"])
    down(y1["sd"] + nh1, y1["th"])
    node("th", "6. Thermal & Segregation", pal["garnet"])

    # Inner iteration loopback: Stage 6 back to Stage 4
    inner_loop_col = pal["fern"]["solid"]
    draw_polyline(
        ax,
        [
            (cx1 - 70, y1["th"] + nh1),
            (cx1 - 70, y1["th"] + nh1 + 10),
            (ix + 16, y1["th"] + nh1 + 10),
            (ix + 16, y1["hyd"] + nh1 / 2),
            (sl1, y1["hyd"] + nh1 / 2),
        ],
        inner_loop_col,
        bg,
        lw=2.0,
        z=5,
        bridge=False,
    )
    draw_arrow(
        ax,
        ix + 16,
        y1["hyd"] + nh1 / 2,
        sl1,
        y1["hyd"] + nh1 / 2,
        color=inner_loop_col,
        lw=2.0,
        scale=11,
        zorder=8,
    )
    ax.text(
        ix + 20,
        y1["sd"] + nh1 / 2,
        r"iter $< \mathrm{max\_iters}$",
        fontsize=10.5,
        color=inner_loop_col,
        rotation=90,
        ha="center",
        va="center",
        zorder=6,
    )

    down(iy + ih, y1["mk"])

    node("mk", "7. Marker Stress & Diffusion", pal["ink"])
    down(y1["mk"] + nh1, y1["vent"])
    node("vent", "8. Fluid Venting & Atmosphere", pal["iris"])
    down(y1["vent"] + nh1, y1["adv"])
    node("adv", "9. Advect, Replenish & Telescope", pal["ink"])

    # Time advance: t <- t + dt between Stage 9 and Stage 10
    down(y1["adv"] + nh1, 770.0)
    draw_box(
        ax,
        cx1 - 65,
        770.0,
        130,
        26,
        facecolor=bg,
        edgecolor=axis,
        lw=1.6,
        radius=6,
        zorder=8,
    )
    ax.text(
        cx1,
        783.0,
        r"$t \leftarrow t + \Delta t$",
        fontsize=12.0,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=9,
    )
    down(796.0, y1["out"])

    yout = y1["out"]
    pts = [
        (sl1 + 14, yout),
        (sr1 + 6, yout),
        (sr1 - 14, yout + nh1),
        (sl1 - 6, yout + nh1),
    ]
    ax.add_patch(
        Polygon(
            pts,
            closed=True,
            facecolor=pal["garnet"]["wash"],
            edgecolor=pal["garnet"]["accent"],
            lw=2.0,
            zorder=8,
        )
    )
    ax.text(
        cx1,
        yout + nh1 / 2 + 1,
        "10. Diagnostics & Output",
        fontsize=fs_node1,
        fontweight="medium",
        color=pal["garnet"]["text"],
        ha="center",
        va="center",
        zorder=9,
    )
    down(yout + nh1, dia_y1 - dia_hh1)

    # Exit check diamond
    d = [
        (cx1, dia_y1 - dia_hh1),
        (cx1 + dia_hw1, dia_y1),
        (cx1, dia_y1 + dia_hh1),
        (cx1 - dia_hw1, dia_y1),
    ]
    ax.add_patch(Polygon(d, closed=True, facecolor=bg, edgecolor=ink, lw=2.0, zorder=8))
    ax.text(
        cx1,
        dia_y1 + 1,
        r"$t > t_\mathrm{end} \text{ or step} = n_\mathrm{steps} \text{ ?}$",
        fontsize=11.5,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=9,
    )

    # Next step loop-back: diamond -> left rail -> top -> into Ambient & Accretion
    loop_col = ink_ramp["700"]
    draw_polyline(
        ax,
        [
            (cx1 - dia_hw1, dia_y1),
            (loop_x1, dia_y1),
            (loop_x1, top_y1),
            (cx1 - 70, top_y1),
            (cx1 - 70, y1["acc"] - 12),
        ],
        loop_col,
        bg,
        lw=2.2,
        z=5,
        bridge=False,
    )
    draw_arrow(
        ax,
        cx1 - 70,
        y1["acc"] - 14,
        cx1 - 70,
        y1["acc"],
        color=loop_col,
        lw=2.2,
        scale=13,
        zorder=8,
    )
    ax.text(
        537,
        dia_y1 - 15,
        "No",
        fontsize=13.0,
        fontweight="medium",
        color=loop_col,
        ha="center",
        va="center",
        zorder=9,
    )

    # Non-convergence retry: starts at Stokes-Darcy (Stage 5) -> right rail -> top -> into Ambient & Accretion
    rc = pal["magma"]["accent"]
    draw_polyline(
        ax,
        [
            (sr1, y1["sd"] + nh1 / 2),
            (retry_x1, y1["sd"] + nh1 / 2),
            (retry_x1, top_y1),
            (cx1 + 70, top_y1),
            (cx1 + 70, y1["acc"] - 12),
        ],
        rc,
        bg,
        lw=2.2,
        ls=(0, (5, 3)),
        z=5,
        bridge=False,
    )
    draw_arrow(
        ax,
        cx1 + 70,
        y1["acc"] - 14,
        cx1 + 70,
        y1["acc"],
        color=rc,
        lw=2.2,
        scale=13,
        zorder=8,
    )
    ax.text(
        (cx1 + 70 + retry_x1) / 2,
        top_y1 + 17,
        "retry",
        fontsize=12.0,
        fontweight="medium",
        color=pal["magma"]["text"],
        ha="center",
        va="center",
        zorder=9,
    )
    ax.text(
        retry_x1 - 10,
        y1["sd"] + nh1 / 2 - 12,
        r"plastic not converged: restore, $\Delta t \leftarrow \min(\Delta t/2, \Delta t_\mathrm{next})$",
        fontsize=10.5,
        color=pal["magma"]["text"],
        fontstyle="italic",
        ha="right",
        va="bottom",
        zorder=9,
    )
    ax.text(
        retry_x1 - 10,
        y1["sd"] + nh1 / 2 + 12,
        r"(abort after solver.max_dt_reductions)",
        fontsize=10.0,
        color=muted,
        ha="right",
        va="top",
        zorder=9,
    )

    draw_arrow(
        ax,
        cx1,
        dia_y1 + dia_hh1,
        cx1,
        960.0,
        color=pal["fern"]["solid"],
        lw=2.4,
        scale=13,
        zorder=8,
    )
    ax.text(
        cx1 + 14,
        944.0,
        "Yes",
        fontsize=13.0,
        fontweight="medium",
        color=pal["fern"]["text"],
        ha="left",
        va="center",
        zorder=9,
    )
    draw_box(
        ax,
        cx1 - 92,
        960.0,
        184,
        44,
        facecolor=bg,
        edgecolor=ink,
        lw=2.0,
        radius=22,
        zorder=8,
    )
    ax.text(
        cx1,
        982.0,
        "End of Run",
        fontsize=fs_node1,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=9,
    )

    # Connectors between cards and spine
    for key, (_, _, paths) in layout1.items():
        for p in paths:
            draw_polyline(ax, p, cols1[key]["accent"], bg, lw=2.2, z=6)

    # Legend at bottom left
    y0 = 840.0
    draw_box(
        ax, 24, y0, 420, 74, facecolor=bg, edgecolor=axis, lw=1.6, radius=12, zorder=2
    )
    draw_box(
        ax, 40, y0 + 20, 78, 34, facecolor=bg, edgecolor=ink, lw=1.6, radius=8, zorder=3
    )
    ax.text(
        79,
        y0 + 38,
        "Stage",
        fontsize=12.0,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=4,
    )

    draw_polyline(
        ax,
        [(132, y0 + 37), (210, y0 + 37)],
        pal["magma"]["accent"],
        bg,
        lw=2.2,
        ls=(0, (5, 3)),
        z=3,
        bridge=False,
    )
    ax.text(
        171,
        y0 + 23,
        "retry",
        fontsize=11.5,
        fontweight="medium",
        color=pal["magma"]["text"],
        ha="center",
        va="center",
        zorder=4,
    )

    ax.add_patch(
        Polygon(
            [(266, y0 + 17), (308, y0 + 37), (266, y0 + 57), (224, y0 + 37)],
            closed=True,
            facecolor=bg,
            edgecolor=ink,
            lw=1.6,
            zorder=3,
        )
    )
    ax.text(
        266,
        y0 + 38,
        "Check",
        fontsize=11.0,
        fontweight="medium",
        color=ink,
        ha="center",
        va="center",
        zorder=4,
    )

    ax.add_patch(
        Polygon(
            [(330, y0 + 20), (386, y0 + 20), (378, y0 + 54), (322, y0 + 54)],
            closed=True,
            facecolor=pal["garnet"]["wash"],
            edgecolor=pal["garnet"]["accent"],
            lw=1.6,
            zorder=3,
        )
    )
    ax.text(
        352,
        y0 + 38,
        "I/O",
        fontsize=12.0,
        fontweight="medium",
        color=pal["garnet"]["text"],
        ha="center",
        va="center",
        zorder=4,
    )

    draw_polyline(
        ax,
        [(394, y0 + 37), (430, y0 + 37)],
        pal["fern"]["accent"],
        bg,
        lw=2.2,
        z=3,
        bridge=False,
    )

    # 8 Equation Cards
    cards1 = {
        "acc": (
            "1. Ambient & Accretion",
            "accretion.jl, geometry.jl",
            [
                r"$w_\mathrm{3D} = 2r \quad (\text{out-of-plane weight})$",
                r"$\int w_\mathrm{3D}\,\mathrm{d}A = \frac{4}{3}\pi R^3$",
            ],
            ["Pebble Accretion", "Sticky Air"],
        ),
        "heat": (
            "2. Radiogenic Heating",
            "radioactivity.jl",
            [
                r"$Q_i = \frac{f_i r_i E_i}{\tau_{1/2, i}/\ln 2} e^{-t\ln 2/\tau_{1/2, i}}$",
                r"$i \in \{^{26}\mathrm{Al},\, ^{60}\mathrm{Fe}\}$",
            ],
            ["Radiogenic Decay", "C, N, H Pyrolysis"],
        ),
        "hyd": (
            "4. Hydration Reactions",
            "serpentinization.jl",
            [
                r"$\mathrm{W_{sil}} \rightleftharpoons \mathrm{D_{sil}} + \mathrm{H_2O} \quad(\Delta t_\mathrm{r}\text{ relax})$",
                r"$\Delta G = \Delta H - T\Delta S + P_\mathrm{f}\Delta V + R_\mathrm{gas}T\ln(X_\mathrm{D}/X_\mathrm{W})$",
            ],
            ["Hydration", "Dehydration"],
        ),
        "vent": (
            "8. Fluid Venting & Atmosphere",
            "atmosphere.jl",
            [
                r"$T^4(\tau_\mathrm{IR}) = \frac{3}{4}T_\mathrm{int}^4\left(\frac{2}{3}+\tau_\mathrm{IR}\right) + \frac{3}{4}T_\mathrm{eq}^4\,f_\gamma(\tau_\mathrm{IR})$",
                r"$S_\mathrm{vent} = \max(0, \phi - \phi_\mathrm{min})/\Delta t$",
            ],
            ["H, C, N, S Degassing", "Jeans, XUV Escape"],
        ),
        "p2m": (
            "3. Marker Mapping & Gravity",
            "melting_convection.jl, gravity.jl",
            [
                r"$\nabla^2\Phi = \frac{8}{3}\pi G\rho \quad (F_\mathrm{grav})$",
                r"$F_\mathrm{m} = (T - T_\mathrm{s})\,/\,(T_\mathrm{l} - T_\mathrm{s})$",
            ],
            ["Marker to Grid", r"Surface $\sigma_\mathrm{SB} T^4$"],
        ),
        "sd": (
            "5. Hydro-Mechanical Stokes-Darcy",
            "stokes_darcy.jl",
            [
                r"$\nabla\cdot\mathbf{v}_\mathrm{s} + \frac{P_\mathrm{t}-P_\mathrm{f}}{\eta_\phi(1-\phi)} + \frac{\beta_\mathrm{d} K_\mathrm{BW}}{\Delta t}(\dots) = 0$",
                r"$\mathbf{q}_\mathrm{D} = -\frac{k}{\eta_\mathrm{f}}(\nabla P_\mathrm{f} - \rho_\mathrm{f}\mathbf{g})$",
            ],
            [r"Darcy $\mathbf{q}_\mathrm{D}$", "Plastic Yielding"],
        ),
        "th": (
            "6. Thermal & Segregation",
            "thermal.jl, segregation.jl",
            [
                r"$\rho c_{p,\mathrm{app}}\,\frac{\mathrm{D}T}{\mathrm{D}t} = \nabla\cdot(k_\mathrm{eff}\nabla T) + H_\mathrm{r} + H_\mathrm{s} + H_\mathrm{a} + H_\mathrm{seg} + H_\mathrm{L}$",
                r"$\mathbf{J}_\mathrm{Fe} = -\frac{k_\mathrm{Fe}}{\eta_\mathrm{Fe}}(\nabla P_\mathrm{t} - \rho_\mathrm{Fe}\mathbf{g})$",
            ],
            ["Metal & Melt Seg.", r"Soft-Turb. $k_\mathrm{eff}$"],
        ),
        "mk": (
            "7 & 9. Markers, Diffusion & Doubling",
            "particles.jl, advect_markers.jl",
            [
                r"$\mathrm{d}\mathbf{x}_m/\mathrm{d}t = \mathbf{v}(\mathbf{x}_m) \quad(\text{RK4})$",
                r"$\Delta \sigma_m = (\sigma_\mathrm{grid} - \sigma_m)[1 - e^{-d_s \Delta t/\tau_\mathrm{M}}]$",
            ],
            ["Subgrid Diffusion", "Domain Doubling"],
        ),
    }

    for key, (x_card, y_card, _) in layout1.items():
        col = cols1[key]
        title, mod, eqs, pills = cards1[key]
        draw_box(
            ax,
            x_card,
            y_card,
            bw1,
            bh1,
            facecolor=col["wash"],
            edgecolor=col["accent"],
            lw=2.4,
            radius=16,
            zorder=2,
        )
        ax.text(
            x_card + 18,
            y_card + 26,
            title,
            fontsize=fs_title1,
            fontweight="semibold",
            color=col["text"],
            ha="left",
            va="center",
            zorder=4,
        )
        bw = len(mod) * 7.5 + 18
        draw_box(
            ax,
            x_card + bw1 - bw - 14,
            y_card + 12,
            bw,
            28,
            facecolor=bg,
            edgecolor=col["border"],
            lw=1.2,
            radius=6,
            zorder=3,
        )
        ax.text(
            x_card + bw1 - bw / 2 - 14,
            y_card + 26,
            mod,
            fontsize=fs_badge1,
            color=muted,
            fontfamily="Geist Mono",
            ha="center",
            va="center",
            zorder=4,
        )
        draw_box(
            ax,
            x_card + 14,
            y_card + 48,
            bw1 - 28,
            68,
            facecolor=bg,
            edgecolor=col["border"],
            lw=1.5,
            radius=10,
            zorder=3,
        )
        for j, eq in enumerate(eqs):
            ax.text(
                x_card + bw1 / 2,
                y_card + 67 + j * 29,
                eq,
                fontsize=fs_eq1,
                fontweight="normal",
                color=ink,
                ha="center",
                va="center",
                zorder=4,
            )
        pw = (bw1 - 28 - 12) / 2
        for i, lbl in enumerate(pills):
            px = x_card + 14 + pw / 2 + i * (pw + 12)
            draw_pill(
                ax,
                px,
                y_card + 142,
                pw,
                34,
                lbl,
                facecolor=bg,
                edgecolor=col["accent"],
                text_color=col["text"],
                fontsize=fs_pill1,
            )

    fig.savefig(out_path, format="svg", bbox_inches=None, metadata={"Date": None})
    plt.close(fig)


# =============================================================================
# Figure 2: Planetesimal Evolution in Three Stages
# =============================================================================


def p_polar(cx: float, cy: float, r: float, th: float):
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
    x, y = p_polar(*c, r, th)
    rot = tangent_rot(th) if tangent else 0.0
    effects = (
        [pe.Stroke(linewidth=4, foreground=bg), pe.Normal()] if bg is not None else None
    )
    ax.text(
        x,
        y,
        s,
        fontsize=fs,
        fontweight="medium",
        color=color,
        ha="center",
        va="center",
        rotation=rot,
        rotation_mode="anchor",
        zorder=12,
        path_effects=effects,
        **kw,
    )


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
    x1, y1 = p_polar(*c, r1, th)
    x2, y2 = p_polar(*c, r2, th)
    draw_block_arrow(
        ax,
        x1,
        y1,
        x2,
        y2,
        facecolor=col_face,
        edgecolor=col_edge,
        width=w,
        head_width=hw,
        head_len=hl,
        lw=1.6,
        zorder=z,
    )


def dots(ax, c, pts, face, ink, r=6.5, z=10):
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
    ink : str
        Outline stroke color.
    r : float, optional
        Circle radius, by default 6.5.
    z : int, optional
        Layer ordering, by default 10.

    Returns
    -------
    None
    """
    for rr, th in pts:
        ax.add_patch(
            Circle(
                p_polar(*c, rr, th), r, facecolor=face, edgecolor=ink, lw=1.5, zorder=z
            )
        )


def convection(ax, c, r, th, ink, size=19, color=None, z=10):
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
    ink : str
        Default arrow color.
    size : float, optional
        Horizontal swirl radius, by default 19.
    color : str, optional
        Arrow color override.
    z : int, optional
        Layer ordering, by default 10.

    Returns
    -------
    None
    """
    col = color if color is not None else ink
    x, y = p_polar(*c, r, th)
    draw_arrow(
        ax,
        x - size,
        y + 4,
        x + size,
        y + 4,
        color=col,
        lw=2.2,
        rad=0.9,
        scale=13,
        zorder=z,
    )
    draw_arrow(
        ax,
        x + size,
        y - 4,
        x - size,
        y - 4,
        color=col,
        lw=2.2,
        rad=0.9,
        scale=13,
        zorder=z,
    )


def group_box(ax, x, y, w, h, col, bg, title, pills, fs_title=15.0, fs_pill=13.0):
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
    bg : str
        Background color.
    title : str
        Group title.
    pills : list of str
        List of 2 pill label strings.
    fs_title : float, optional
        Title font size, by default 15.0.
    fs_pill : float, optional
        Pill font size, by default 13.0.

    Returns
    -------
    None
    """
    draw_box(
        ax,
        x,
        y,
        w,
        h,
        facecolor=col["wash"],
        edgecolor=col["accent"],
        lw=2.4,
        radius=16,
        zorder=3,
    )
    ax.text(
        x + w / 2,
        y + 27,
        title,
        fontsize=fs_title,
        fontweight="semibold",
        color=col["text"],
        ha="center",
        va="center",
        zorder=4,
    )
    gap = 10
    pw = (w - 24 - gap) / 2
    for i, s in enumerate(pills):
        draw_pill(
            ax,
            x + 12 + pw / 2 + i * (pw + gap),
            y + h - 30,
            pw,
            36,
            s,
            facecolor=bg,
            edgecolor=col["accent"],
            text_color=col["text"],
            fontsize=fs_pill,
            radius=8,
            lw=1.6,
        )


def legend2(
    ax, x, y, bg, axis, ink, text_col, surface, ink_300, garnet_acc, w=380, h=180
):
    """Render legend box for Figure 2.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x : float
        Left coordinate.
    y : float
        Top coordinate.
    bg : str
        Background color.
    axis : str
        Axis border color.
    ink : str
        Ink stroke color.
    text_col : str
        Text color.
    surface : str
        Surface color.
    ink_300 : str
        Muted ink fill color.
    garnet_acc : str
        Garnet accent color.
    w : float, optional
        Box width, by default 380.
    h : float, optional
        Box height, by default 180.

    Returns
    -------
    None
    """
    draw_box(ax, x, y, w, h, facecolor=bg, edgecolor=axis, lw=1.6, radius=14, zorder=3)
    rows = [y + 26, y + 58, y + 90, y + 122, y + 154]
    for dx, dy in [(-15, 3), (0, -4), (15, 3)]:
        ax.add_patch(
            Circle(
                (x + 48 + dx, rows[0] + dy),
                6.5,
                facecolor=bg,
                edgecolor=ink,
                lw=1.5,
                zorder=4,
            )
        )
    for dx, dy in [(-13, -2), (3, 4), (16, -3)]:
        ax.add_patch(
            Circle(
                (x + 48 + dx, rows[1] + dy),
                6.5,
                facecolor=garnet_acc,
                edgecolor=ink,
                lw=1.4,
                zorder=4,
            )
        )
    draw_box(
        ax,
        x + 20,
        rows[2] - 12,
        56,
        24,
        facecolor=surface,
        edgecolor=axis,
        lw=1.6,
        radius=7,
        zorder=4,
    )
    draw_block_arrow(
        ax,
        x + 20,
        rows[3],
        x + 78,
        rows[3],
        facecolor=ink_300,
        edgecolor=ink,
        width=12,
        head_width=24,
        head_len=15,
        zorder=4,
    )
    # Swirl arrows for convection
    draw_arrow(
        ax,
        x + 28,
        rows[4] + 3,
        x + 68,
        rows[4] + 3,
        color=ink,
        lw=2.0,
        rad=0.9,
        scale=11,
        zorder=4,
    )
    draw_arrow(
        ax,
        x + 68,
        rows[4] - 3,
        x + 28,
        rows[4] - 3,
        color=ink,
        lw=2.0,
        rad=0.9,
        scale=11,
        zorder=4,
    )

    labels = [
        "Volatile gas (H, C, N, S)",
        "Fe-FeS metal alloy",
        "Process group",
        "Energy or mass flux",
        "Convection circulation",
    ]
    for yy, s in zip(rows, labels):
        ax.text(
            x + 98,
            yy + 1,
            s,
            fontsize=13.0,
            color=text_col,
            ha="left",
            va="center",
            zorder=4,
        )


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
    fwl_mpl.use("light")
    tok = resolve_tokens("light")
    paper_roles = fvl.roles("paper")
    ink_ramp = fvl.ramp("ink")

    bg = paper_roles["background"]
    surface = paper_roles["surface"]
    sub = ink_ramp["100"]
    line = paper_roles["line"]
    ink_300 = ink_ramp["300"]
    axis = paper_roles["axis"]
    muted = paper_roles["muted"]
    text_col = paper_roles["text"]
    ink = paper_roles["heading"]

    pal = get_palette_hues(tok, ink_ramp, bg, ink, muted, surface, sub, line, axis)

    w2, h2 = 1400.0, 960.0
    fig, ax = new_equal_canvas(15.4, 15.4 * h2 / w2, w2, h2)

    c = (700.0, 470.0)
    r_halo, r_out, r_mid, r_in = 262.0, 232.0, 180.0, 112.0
    fs_layer2, fs_flux2 = 13.0, 14.5

    ax.add_patch(Circle(c, r_halo, facecolor=surface, edgecolor=line, lw=1.8, zorder=1))

    # Sector I: cold accretion (150..270)
    sig = pal["signal"]
    ax.add_patch(
        Wedge(
            c,
            r_out,
            150,
            270,
            width=r_out - r_mid,
            facecolor=sig["mid"],
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(
            c, r_mid, 150, 270, facecolor=sig["wash"], edgecolor=ink, lw=1.6, zorder=2
        )
    )
    polar_label(
        ax,
        c,
        (r_out + r_mid) / 2,
        222,
        "Accreting Shell",
        fs_layer2,
        sig["text"],
        sig["mid"],
    )
    polar_label(
        ax,
        c,
        112,
        212,
        "Porous Ice + Rock",
        fs_layer2,
        sig["text"],
        sig["wash"],
        tangent=False,
    )
    lag = pal["lagoon"]
    for th in (196, 248):
        radial_block(ax, c, th, 318, r_out + 4, lag["mid"], lag["accent"])
    polar_label(
        ax,
        c,
        278,
        248,
        r"$\dot{M}_\mathrm{acc}$",
        fs_flux2,
        lag["text"],
        bg,
        tangent=False,
    )

    # Sector II: hydrothermal alteration (270..390)
    fern, och = pal["fern"], pal["ochre"]
    ax.add_patch(
        Wedge(
            c,
            r_out,
            270,
            390,
            width=r_out - r_mid,
            facecolor=fern["mid"],
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(
            c,
            r_mid,
            270,
            390,
            width=r_mid - r_in,
            facecolor=fern["tint"],
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(c, r_in, 270, 390, facecolor=och["tint"], edgecolor=ink, lw=1.6, zorder=2)
    )
    polar_label(
        ax, c, (r_out + r_mid) / 2, 295, "Crust", fs_layer2, fern["text"], fern["mid"]
    )
    polar_label(
        ax,
        c,
        (r_mid + r_in) / 2,
        312,
        "Serpentinization",
        fs_layer2,
        fern["text"],
        fern["tint"],
    )
    polar_label(
        ax,
        c,
        60,
        330,
        "Dehydration",
        fs_layer2 - 1.0,
        och["text"],
        och["tint"],
        tangent=False,
    )
    # Darcy percolation vector confined within permeable solid crust
    radial_block(ax, c, 342, r_in + 20, r_out - 12, fern["tint"], fern["text"])
    polar_label(
        ax,
        c,
        (r_in + 20 + r_out - 12) / 2,
        354,
        r"$\mathbf{q}_\mathrm{D}$",
        fs_flux2,
        fern["solid"],
        bg,
        tangent=False,
    )
    convection(ax, c, (r_mid + r_in) / 2, 372, ink, color=fern["text"])

    # Sector III: magma ocean and core (30..150)
    mag, gar, iris = pal["magma"], pal["garnet"], pal["iris"]
    r_core, r_mo, r_cr = 100.0, 204.0, 232.0
    ax.add_patch(
        Wedge(
            c,
            r_halo,
            30,
            150,
            width=r_halo - r_cr,
            facecolor=iris["tint"],
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(
            c,
            r_cr,
            30,
            150,
            width=r_cr - r_mo,
            facecolor=sub,
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(
            c,
            r_mo,
            30,
            150,
            width=r_mo - r_core,
            facecolor=mag["mid"],
            edgecolor=ink,
            lw=1.6,
            zorder=2,
        )
    )
    ax.add_patch(
        Wedge(
            c, r_core, 30, 150, facecolor=gar["accent"], edgecolor=ink, lw=1.8, zorder=2
        )
    )
    polar_label(
        ax, c, 52, 90, "Fe-FeS Core", fs_layer2, bg, gar["accent"], tangent=False
    )
    polar_label(
        ax,
        c,
        152,
        90,
        "Magma Ocean",
        fs_layer2 + 0.5,
        mag["text"],
        mag["mid"],
        tangent=False,
    )
    polar_label(ax, c, (r_cr + r_mo) / 2, 134, "Solid Lid", 12.5, ink, bg)
    polar_label(ax, c, (r_halo + r_cr) / 2, 136, "Atmosphere", 12.5, iris["text"], bg)
    convection(ax, c, 152, 126, ink, color=mag["text"])
    dots(ax, c, [(118, 62), (132, 70), (120, 78)], gar["accent"], ink)
    # Metal settling drift-flux vector descending toward core
    radial_block(ax, c, 50, 150, r_core + 2, gar["tint"], gar["accent"])
    polar_label(
        ax,
        c,
        178,
        38,
        r"$\mathbf{J}_\mathrm{Fe}$",
        fs_flux2 - 1,
        gar["text"],
        mag["mid"],
        tangent=False,
    )
    dots(ax, c, [(186, 100), (182, 82)], bg, ink)
    dots(ax, c, [(247, 84), (247, 98), (247, 112)], bg, ink, r=6.0)
    # Surface radiative heat flux radiating into space
    radial_block(ax, c, 62, r_cr, 282, mag["tint"], mag["text"])
    polar_label(
        ax, c, 304, 54, r"$F_\mathrm{rad}$", fs_flux2, mag["text"], bg, tangent=False
    )

    for th in (30, 150, 270):
        x2, y2 = p_polar(*c, r_halo, th)
        ax.plot([c[0], x2], [c[1], y2], color=ink, lw=2.6, zorder=8)

    for x_pill, y_pill, col, s in [
        (470, 82.0, sig, "I. Cold Accretion"),
        (930, 82.0, fern, "II. Hydrothermal Alteration"),
        (700, 782.0, mag, "III. Magma Ocean & Core"),
    ]:
        draw_pill(
            ax,
            x_pill,
            y_pill,
            330,
            42,
            s,
            facecolor=col["wash"],
            edgecolor=col["accent"],
            text_color=col["solid"],
            fontsize=14.0,
            radius=14,
            lw=2.0,
        )

    # 8 Process group boxes with Stage tags
    left_col, right_col, box_w, box_h = 24, 1051, 325, 112
    group_box(
        ax,
        left_col,
        145.0,
        box_w,
        box_h,
        lag,
        bg,
        "[Stage I] Accretion & Disk",
        ["Pebbles", r"Disk $T_\mathrm{amb}$"],
    )
    group_box(
        ax,
        left_col,
        305.0,
        box_w,
        box_h,
        sig,
        bg,
        "[Stage I] Porous Ice & Rock",
        ["Accreted Ices", "Pore Sealing"],
    )
    group_box(
        ax,
        left_col,
        465.0,
        box_w,
        box_h,
        gar,
        bg,
        "[Stage III] Core Formation",
        ["Segregation", r"Partition $D_i$"],
    )
    group_box(
        ax,
        left_col,
        625.0,
        box_w,
        box_h,
        iris,
        bg,
        "[Stage III] Volatiles & Redox",
        ["H, C, N, S", "IW / QFM"],
    )
    group_box(
        ax,
        right_col,
        145.0,
        box_w,
        box_h,
        fern,
        bg,
        "[Stage II] Hydrothermal",
        [r"Darcy $\mathbf{q}_\mathrm{D}$", "Hydrofracture"],
    )
    group_box(
        ax,
        right_col,
        305.0,
        box_w,
        box_h,
        och,
        bg,
        "[Stage II] Heating & Pyrolysis",
        [r"$^{26}$Al, $^{60}$Fe", "C, N, H Pyrolysis"],
    )
    group_box(
        ax,
        right_col,
        465.0,
        box_w,
        box_h,
        mag,
        bg,
        "[Stage III] Magma Ocean",
        [r"$F_\mathrm{m} > F_\mathrm{crit}$", r"Soft-Turb. $k_\mathrm{eff}$"],
    )
    group_box(
        ax,
        right_col,
        625.0,
        box_w,
        box_h,
        iris,
        bg,
        "[Stage III] Atmosphere & Loss",
        [r"Guillot $T(\tau_\mathrm{IR})$", "Jeans, XUV"],
    )

    legend2(ax, 995, 765.0, bg, axis, ink, text_col, surface, ink_300, gar["accent"])

    fig.savefig(out_path, format="svg", bbox_inches=None, metadata={"Date": None})
    plt.close(fig)


def main():
    """Parse command-line arguments and generate architecture figures.

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description="Generate Erebus.jl documentation architecture figures."
    )
    parser.add_argument(
        "--fig",
        choices=["1", "2", "all"],
        default="all",
        help="Figure to generate (1, 2, or all)",
    )
    parser.add_argument(
        "--fonttype",
        choices=["path", "none"],
        default="path",
        help="SVG fonttype ('path' or 'none')",
    )
    parser.add_argument(
        "--outdir", type=Path, default=ASSETS_DIR, help="Output directory"
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams["svg.fonttype"] = args.fonttype
    mpl.rcParams["svg.hashsalt"] = "erebus_reproducible_hash_2026"

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
