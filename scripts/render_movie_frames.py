#!/usr/bin/env python3
"""
Render movie frames and encode an MP4 video of the hydrothermal reaction benchmark.
"""
import os
import glob
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from multiprocessing import Pool
import subprocess

import sys

run_tag = sys.argv[1] if len(sys.argv) > 1 else "128x128"
out_name = sys.argv[2] if len(sys.argv) > 2 else "hydrothermal_reaction_128"
frame_dir = sys.argv[3] if len(sys.argv) > 3 else "/tmp/erebus_movie_frames"
png_dir = f"/tmp/erebus_movie_pngs_{out_name}"
os.makedirs(png_dir, exist_ok=True)

# Read metadata
meta_path = os.path.join(frame_dir, "meta.bin")
with open(meta_path, "rb") as f:
    nx, ny = struct.unpack("ii", f.read(8))
    rplanet, xc, yc = struct.unpack("fff", f.read(12))
    x = np.frombuffer(f.read(nx * 4), dtype=np.float32)
    y = np.frombuffer(f.read(ny * 4), dtype=np.float32)

rplanet_km = rplanet * 1e-3
xc_km = xc * 1e-3
yc_km = yc * 1e-3
xp_km = np.linspace(x[0]*1e-3, x[-1]*1e-3, nx) - xc_km
yp_km = np.linspace(y[0]*1e-3, y[-1]*1e-3, ny) - yc_km
Xp, Yp = np.meshgrid(xp_km, yp_km)

def render_frame(fpath):
    idx_str = os.path.basename(fpath).replace("frame_", "").replace(".bin", "")
    png_path = os.path.join(png_dir, f"frame_{idx_str}.png")

    with open(fpath, "rb") as f:
        timesum = struct.unpack("f", f.read(4))[0]
        time_Ma = timesum / (365.25 * 86400 * 1e6)
        n_elem = ny * nx
        tk = np.frombuffer(f.read(n_elem * 4), dtype=np.float32).reshape(ny, nx).T
        XWS = np.frombuffer(f.read(n_elem * 4), dtype=np.float32).reshape(ny, nx).T
        DQPF = np.frombuffer(f.read(n_elem * 4), dtype=np.float32).reshape(ny, nx).T
        DHP = np.frombuffer(f.read(n_elem * 4), dtype=np.float32).reshape(ny, nx).T

    # Fixed scales and normalizations
    t_min, t_max = 170.0, 1200.0
    t_levels = np.linspace(t_min, t_max, 41)
    t_norm = matplotlib.colors.Normalize(vmin=t_min, vmax=t_max)
    t_ticks = [200, 400, 600, 800, 1000, 1200]

    xw_min, xw_max = 0.0, 1.0
    xw_levels = np.linspace(xw_min, xw_max, 31)
    xw_norm = matplotlib.colors.Normalize(vmin=xw_min, vmax=xw_max)
    xw_ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

    dqpf_limit = 1.0e-14
    dqpf_levels = np.linspace(-dqpf_limit, dqpf_limit, 41)
    dqpf_norm = matplotlib.colors.Normalize(vmin=-dqpf_limit, vmax=dqpf_limit)
    dqpf_ticks = [-1.0e-14, -5.0e-15, 0.0, 5.0e-15, 1.0e-14]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=150)
    fig.suptitle(f"Erebus 2D Hydrothermal Benchmark ({run_tag})  |  Time = {time_Ma:5.2f} Ma", fontsize=15, fontweight='bold')

    def format_movie_ax(ax, title, has_ylabel=True):
        ax.set_facecolor('white')
        ax.add_patch(Circle((0.0, 0.0), rplanet_km, fill=False, edgecolor='black', lw=1.2, ls='-'))
        ax.set_aspect('equal')
        ax.set_title(title)
        ax.set_xlabel("x [km]")
        if has_ylabel:
            ax.set_ylabel("y [km]")
        ax.set_xlim(xp_km[0], xp_km[-1])
        ax.set_ylim(yp_km[0], yp_km[-1])

    # Panel 1: Temperature (viridis: perceptually uniform and colorblind friendly)
    ax = axes[0]
    cf = ax.contourf(Xp, Yp, tk, levels=t_levels, norm=t_norm, cmap='viridis', extend='both')
    clip_c0 = Circle((0.0, 0.0), rplanet_km, transform=ax.transData)
    try:
        cf.set_clip_path(clip_c0)
    except Exception:
        for col in cf.collections:
            col.set_clip_path(clip_c0)
    cb = fig.colorbar(cf, ax=ax, ticks=t_ticks, label="Temperature [K]", fraction=0.046, pad=0.04)
    cb.ax.set_ylim(t_min, t_max)
    format_movie_ax(ax, "Temperature $T$", has_ylabel=True)

    # Panel 2: Hydration extent XW (cividis_r: inverted so hydration is blue)
    ax = axes[1]
    cf = ax.contourf(Xp, Yp, XWS, levels=xw_levels, norm=xw_norm, cmap='cividis_r', extend='both')
    clip_c1 = Circle((0.0, 0.0), rplanet_km, transform=ax.transData)
    try:
        cf.set_clip_path(clip_c1)
    except Exception:
        for col in cf.collections:
            col.set_clip_path(clip_c1)
    cb = fig.colorbar(cf, ax=ax, ticks=xw_ticks, label="Hydration Extent $X_W$", fraction=0.046, pad=0.04)
    cb.ax.set_ylim(xw_min, xw_max)
    format_movie_ax(ax, "Hydration Reaction Front ($X_W$)", has_ylabel=False)

    # Panel 3: Fluid source term DQPF (PuOr_r: diverging colorblind friendly)
    ax = axes[2]
    cf = ax.contourf(Xp, Yp, DQPF, levels=dqpf_levels, norm=dqpf_norm, cmap='PuOr_r', extend='both')
    clip_c2 = Circle((0.0, 0.0), rplanet_km, transform=ax.transData)
    try:
        cf.set_clip_path(clip_c2)
    except Exception:
        for col in cf.collections:
            col.set_clip_path(clip_c2)
    cb = fig.colorbar(cf, ax=ax, ticks=dqpf_ticks, label="DQPF [1/s] (Fluid Exchange)", fraction=0.046, pad=0.04, format='%.1e')
    cb.ax.set_ylim(-dqpf_limit, dqpf_limit)
    format_movie_ax(ax, "Fluid Mass Exchange (DQPF)", has_ylabel=False)

    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.12, top=0.88, wspace=0.35)
    plt.savefig(png_path)
    plt.close(fig)
    return png_path

if __name__ == '__main__':
    files = sorted(glob.glob(os.path.join(frame_dir, "frame_*.bin")))
    print(f"Rendering {len(files)} frames using multiprocessing...")
    with Pool() as pool:
        results = pool.map(render_frame, files)

    print(f"Rendered {len(results)} PNG frames in {png_dir}")

    # Encode MP4 video with ffmpeg
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    mp4_out = os.path.join(repo_root, "docs", "src", "assets", f"{out_name}.mp4")
    gif_out = os.path.join(repo_root, "docs", "src", "assets", f"{out_name}.gif")
    env = os.environ.copy()
    env["DYLD_FALLBACK_LIBRARY_PATH"] = "/opt/homebrew/Cellar/x265/4.2/lib"
    cmd_mp4 = [
        "ffmpeg", "-y",
        "-r", "20",
        "-i", os.path.join(png_dir, "frame_%04d.png"),
        "-vcodec", "libx264",
        "-crf", "22",
        "-pix_fmt", "yuv420p",
        mp4_out
    ]
    subprocess.run(cmd_mp4, env=env, check=True)
    print(f"Successfully generated MP4 video: {mp4_out}")

    cmd_gif = [
        "ffmpeg", "-y",
        "-i", mp4_out,
        "-vf", "fps=12,scale=1200:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse",
        gif_out
    ]
    subprocess.run(cmd_gif, env=env, check=True)
    print(f"Successfully generated GIF: {gif_out}")
