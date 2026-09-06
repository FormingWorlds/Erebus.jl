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

frame_dir = "/tmp/erebus_movie_frames"
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
xp_km = np.linspace(x[0]*1e-3, x[-1]*1e-3, nx)
yp_km = np.linspace(y[0]*1e-3, y[-1]*1e-3, ny)
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

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=150)
    fig.suptitle(f"Erebus 2D Hydrothermal Benchmark ({run_tag})  |  Time = {time_Ma:5.2f} Ma", fontsize=15, fontweight='bold')

    # Panel 1: Temperature
    ax = axes[0]
    cf = ax.contourf(Xp, Yp, tk, levels=35, cmap='inferno', vmin=170, vmax=4300)
    fig.colorbar(cf, ax=ax, label="Temperature [K]")
    c = Circle((xc_km, yc_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
    ax.add_patch(c)
    ax.set_aspect('equal')
    ax.set_title("Temperature $T$")
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")

    # Panel 2: Hydration extent XW
    ax = axes[1]
    cf = ax.contourf(Xp, Yp, XWS, levels=30, cmap='YlGnBu', vmin=0.0, vmax=1.0)
    fig.colorbar(cf, ax=ax, label="Hydration Extent $X_W$")
    c = Circle((xc_km, yc_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
    ax.add_patch(c)
    ax.set_aspect('equal')
    ax.set_title("Hydration Reaction Front ($X_W$)")
    ax.set_xlabel("x [km]")

    # Panel 3: Fluid source term DQPF
    ax = axes[2]
    q_max = max(1e-12, float(np.percentile(np.abs(DQPF), 99.5)))
    cf = ax.contourf(Xp, Yp, DQPF, levels=30, cmap='RdBu_r', vmin=-q_max, vmax=q_max)
    fig.colorbar(cf, ax=ax, label="DQPF [1/s] (Fluid Exchange)")
    c = Circle((xc_km, yc_km), rplanet_km, fill=False, edgecolor='cyan', lw=1.5, ls='--')
    ax.add_patch(c)
    ax.set_aspect('equal')
    ax.set_title("Fluid Mass Exchange (DQPF)")
    ax.set_xlabel("x [km]")

    plt.tight_layout()
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
