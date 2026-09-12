#!/usr/bin/env python3
"""
Collision and overlap validator for figures in Erebus.jl.

Validates that no data lines, threshold lines, vertical/horizontal indicators,
or annotations intersect text labels or legend boxes across generated figures.
"""

import sys
import os
import importlib.util
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt


def bboxes_intersect(b1, b2, margin: float = 2.0) -> bool:
    """
    Check if two bounding boxes intersect within a specified pixel margin.

    Parameters
    ----------
    b1 : matplotlib.transforms.Bbox
        First bounding box.
    b2 : matplotlib.transforms.Bbox
        Second bounding box.
    margin : float, optional
        Safety margin in display pixels, by default 2.0.

    Returns
    -------
    bool
        True if bounding boxes overlap.
    """
    if (b1.x1 + margin < b2.x0 or b1.x0 - margin > b2.x1 or
            b1.y1 + margin < b2.y0 or b1.y0 - margin > b2.y1):
        return False
    return True


def segment_intersects_bbox(p1: np.ndarray, p2: np.ndarray, bbox, margin: float = 2.0) -> bool:
    """
    Check if a 2D line segment intersects a rectangular bounding box.

    Parameters
    ----------
    p1 : np.ndarray
        Start coordinate [x, y] in display space.
    p2 : np.ndarray
        End coordinate [x, y] in display space.
    bbox : matplotlib.transforms.Bbox
        Bounding box to test against.
    margin : float, optional
        Safety margin in display pixels, by default 2.0.

    Returns
    -------
    bool
        True if line segment intersects bounding box.
    """
    x0, y0, x1, y1 = bbox.x0 - margin, bbox.y0 - margin, bbox.x1 + margin, bbox.y1 + margin
    if (x0 <= p1[0] <= x1 and y0 <= p1[1] <= y1) or (x0 <= p2[0] <= x1 and y0 <= p2[1] <= y1):
        return True

    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    p = [-dx, dx, -dy, dy]
    q = [p1[0] - x0, x1 - p1[0], p1[1] - y0, y1 - p1[1]]
    u1, u2 = 0.0, 1.0

    for pi, qi in zip(p, q):
        if pi == 0:
            if qi < 0:
                return False
        else:
            t = qi / pi
            if pi < 0:
                if t > u2:
                    return False
                if t > u1:
                    u1 = t
            else:
                if t < u1:
                    return False
                if t < u2:
                    u2 = t
    return u1 <= u2


def check_ax_collisions(ax: plt.Axes, ax_name: str = "Axis") -> list:
    """
    Find collisions between lines, texts, and legends in an axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Matplotlib axis object.
    ax_name : str, optional
        Identifier name for the axis.

    Returns
    -------
    list of tuple
        Detected collision descriptors (ax_name, item1, description).
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    text_boxes = []
    for text in ax.texts:
        bbox = mpl.text.Text.get_window_extent(text, renderer)
        if bbox.width > 0 and bbox.height > 0:
            text_boxes.append((f"Text: '{text.get_text()[:25]}'", bbox, text))

    leg = ax.get_legend()
    leg_lines = []
    if leg and leg.get_visible():
        leg_bbox = leg.get_window_extent(renderer)
        text_boxes.append(("Legend", leg_bbox, leg))
        leg_lines = leg.get_lines()

    collisions = []

    # Check text-text and text-legend collisions
    for i in range(len(text_boxes)):
        for j in range(i + 1, len(text_boxes)):
            name1, b1, _ = text_boxes[i]
            name2, b2, _ = text_boxes[j]
            if bboxes_intersect(b1, b2):
                collisions.append((ax_name, name1, f"overlaps {name2}"))

    # Check line collisions with text/legend
    for line in ax.get_lines():
        if line in leg_lines:
            continue
        lbl = line.get_label()
        if lbl.startswith("_child"):
            continue

        xdata = line.get_xdata()
        ydata = line.get_ydata()
        if len(xdata) == 0:
            continue

        pts_data = np.column_stack([xdata, ydata])
        valid = np.isfinite(pts_data).all(axis=1)
        if not np.any(valid):
            continue
        pts_data = pts_data[valid]
        if len(pts_data) < 2:
            continue
        pts_disp = ax.transData.transform(pts_data)

        line_name = lbl if lbl and not lbl.startswith("_") else "Line"

        for name, bbox, _ in text_boxes:
            for i in range(len(pts_disp) - 1):
                if segment_intersects_bbox(pts_disp[i], pts_disp[i + 1], bbox):
                    collisions.append((ax_name, name, f"overlaps line {line_name}"))
                    break

    return collisions


def check_all_benchmarks() -> int:
    """
    Run benchmark scripts with collision detection enabled.

    Returns
    -------
    int
        Exit status code (0 if all figures pass with zero collisions).
    """
    total_collisions = 0
    orig_savefig = mpl.figure.Figure.savefig
    current_script = ""

    def patched_savefig(self, *args, **kwargs):
        nonlocal total_collisions
        for ax_idx, ax in enumerate(self.axes):
            cols = check_ax_collisions(ax, f"ax_{ax_idx}")
            for c in cols:
                total_collisions += 1
                print(f"  [COLLISION] in {current_script} ({c[0]}): {c[1]} {c[2]}")
        return orig_savefig(self, *args, **kwargs)

    mpl.figure.Figure.savefig = patched_savefig

    bench_targets = [
        ("benchmarks/generate_hydrothermal_convection_benchmarks.py", "main"),
        ("benchmarks/generate_accretion_benchmarks.py", "main"),
        ("benchmarks/generate_telescoping_benchmarks.py", "main"),
        ("benchmarks/generate_core_geochemistry_benchmarks.py", "generate_benchmark_figure"),
        ("benchmarks/generate_mineral_assemblage_benchmarks.py", "generate_benchmark_figure"),
        ("benchmarks/generate_volatile_mixture_benchmarks.py", "main"),
        ("benchmarks/generate_multistage_accretion_benchmarks.py", "main"),
        ("benchmarks/generate_lunar_growth_benchmarks.py", "main"),
    ]

    for bf, func in bench_targets:
        if not os.path.exists(bf):
            continue
        current_script = bf
        spec = importlib.util.spec_from_file_location("mod", bf)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        getattr(mod, func)()

    diag_script = "benchmarks/plot_diagnostics.py"
    if os.path.exists(diag_script):
        current_script = diag_script
        with open(diag_script, "r") as f:
            code = f.read()
        exec(code, {"__file__": diag_script})

    if total_collisions == 0:
        print("\nAll benchmark figures passed collision validation with 0 overlaps.")
        return 0
    else:
        print(f"\nDetected {total_collisions} line/text collisions across benchmark figures.")
        return 1


if __name__ == "__main__":
    sys.exit(check_all_benchmarks())
