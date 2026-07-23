#!/usr/bin/env python3
"""
PlotSC.py — Animate Square Cylinder (SC) solver output.

Reads Probe_1_Map.csv (pressure), Probe_2u_Field.csv, Probe_2v_Field.csv
from a result directory and writes Animation_SC.mp4.

Usage:
    cd /path/to/HMT
    source .venv/bin/activate
    python ioPlots/PlotSC.py [result_dir] [--fps N] [--dpi N]
                             [--obs X0 X1 Y0 Y1]

result_dir: directory name under ioRes/, or an absolute path.
            Defaults to the most recent exSC* run.
"""

import os
import sys
import time
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.patches as mpatches
from pathlib import Path


# ── parsing ─────────────────────────────────────────────────────────────────

def parse_probe(fpath):
    """
    Read a Map or Field CSV probe.

    CSV layout: rows = time frames; columns = Time, (x0 y0), (x0 y1), …
    The C++ probe writes outer loop x, inner loop y — so y varies fastest
    (ny consecutive columns share the same x value).

    Returns
    -------
    frames : ndarray (n_frames, nx, ny)
    times  : ndarray (n_frames,)
    xcoords: ndarray (nx,)   — from header (Faces positions)
    ycoords: ndarray (ny,)   — from header (Faces positions)
    """
    t0 = time.time()
    name = Path(fpath).name
    print(f"  Reading {name} ...", end=" ", flush=True)

    with open(fpath) as f:
        header = f.readline().rstrip("\n")

    cols = header.split(",")
    pairs = [c.split(" ") for c in cols[1:]]   # each entry = ["x", "y"]
    all_x = [float(p[0]) for p in pairs]
    all_y = [float(p[1]) for p in pairs]

    # y is the inner (fast) loop: count how many columns share the same x
    x0 = all_x[0]
    ny = 1
    while ny < len(all_x) and all_x[ny] == x0:
        ny += 1
    nx = len(all_x) // ny

    xcoords = np.array([all_x[i * ny] for i in range(nx)])
    ycoords = np.array(all_y[:ny])

    raw = np.loadtxt(fpath, delimiter=",", skiprows=1)
    times  = raw[:, 0]
    frames = raw[:, 1:].reshape(len(times), nx, ny)

    print(f"{len(times)} frames, grid {nx}×{ny}  ({time.time()-t0:.1f}s)")
    return frames, times, xcoords, ycoords


# ── velocity interpolation ───────────────────────────────────────────────────

def speed_frames(u_frames, v_frames, u_bc_east=1.0, v_bc_north=0.0):
    """
    Interpolate face-centred u, v to pressure-cell centres and return |V|.

    The probe writes nx values in x for u (missing the East boundary face).
    Padding with the outlet BC and averaging adjacent faces gives the p-cell
    centre value.  Same logic applies to v in y.

    u_frames : (n, nx, ny)
    v_frames : (n, nx, ny)
    """
    n, nx, ny = u_frames.shape
    east_pad = np.full((n, 1, ny), u_bc_east)
    u_pad    = np.concatenate([u_frames, east_pad], axis=1)  # (n, nx+1, ny)
    u_p      = 0.5 * (u_pad[:, :-1, :] + u_pad[:, 1:, :])   # (n, nx, ny)

    north_pad = np.full((n, nx, 1), v_bc_north)
    v_pad     = np.concatenate([v_frames, north_pad], axis=2)  # (n, nx, ny+1)
    v_p       = 0.5 * (v_pad[:, :, :-1] + v_pad[:, :, 1:])    # (n, nx, ny)

    return np.sqrt(u_p**2 + v_p**2)


# ── animation ────────────────────────────────────────────────────────────────

def animate(result_path, fps=30, dpi=150,
            obs_x=(0.4, 0.6), obs_y=(0.65, 0.85)):

    result_path = Path(result_path)
    print(f"\nResult directory: {result_path.name}")

    p_frames, times, xP, yP = parse_probe(result_path / "Probe_1_Map.csv")
    u_frames, _,     xU, yU = parse_probe(result_path / "Probe_2u_Field.csv")
    v_frames, _,     xV, yV = parse_probe(result_path / "Probe_2v_Field.csv")

    spd = speed_frames(u_frames, v_frames)

    n_frames = len(times)
    print(f"  Frames: {n_frames}  |  t = [{times[0]:.3f}, {times[-1]:.3f}] s")

    # Fixed colour limits (consistent across all frames)
    p_min, p_max = float(p_frames.min()), float(p_frames.max())
    p_mid  = 0.5 * (p_min + p_max)
    p_half = max(abs(p_min - p_mid), abs(p_max - p_mid), 1e-6)
    s_min, s_max = 0.0, float(spd.max())

    # Spatial extent for imshow (use pressure grid positions)
    ext = [xP[0], xP[-1], yP[0], yP[-1]]

    # ── figure ────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.2))
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.12, top=0.88, wspace=0.35)

    # Panel 1: pressure
    ax1 = axes[0]
    im1 = ax1.imshow(p_frames[0].T, cmap="RdBu_r", origin="lower",
                     aspect="auto", extent=ext,
                     vmin=p_mid - p_half, vmax=p_mid + p_half)
    cb1 = fig.colorbar(im1, ax=ax1, fraction=0.03, pad=0.02)
    cb1.set_label("p (Pa)", fontsize=9)
    ax1.set_title("Pressure", fontsize=10)
    ax1.set_xlabel("x (m)"); ax1.set_ylabel("y (m)")
    _add_obstacle(ax1, obs_x, obs_y)

    # Panel 2: speed
    ax2 = axes[1]
    im2 = ax2.imshow(spd[0].T, cmap="viridis", origin="lower",
                     aspect="auto", extent=ext,
                     vmin=s_min, vmax=s_max)
    cb2 = fig.colorbar(im2, ax=ax2, fraction=0.03, pad=0.02)
    cb2.set_label("|V| (m/s)", fontsize=9)
    ax2.set_title("Velocity Magnitude", fontsize=10)
    ax2.set_xlabel("x (m)"); ax2.set_ylabel("y (m)")
    _add_obstacle(ax2, obs_x, obs_y)

    title = fig.suptitle(_title(times[0]), fontsize=11, y=0.97)

    def update(k):
        im1.set_data(p_frames[k].T)
        im2.set_data(spd[k].T)
        title.set_text(_title(times[k]))
        return im1, im2, title

    print(f"\nRendering {n_frames} frames at {fps} fps, dpi={dpi} ...")
    t0 = time.time()
    ani = anim.FuncAnimation(fig, update, frames=n_frames,
                             interval=1000 / fps, blit=True, repeat=False)

    out = result_path / "Animation_SC.mp4"
    ani.save(str(out), writer="ffmpeg", fps=fps, dpi=dpi,
             extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    print(f"Saved → {out}  ({time.time()-t0:.0f}s)")
    plt.close(fig)


def _add_obstacle(ax, obs_x, obs_y):
    rect = mpatches.Rectangle(
        (obs_x[0], obs_y[0]),
        obs_x[1] - obs_x[0],
        obs_y[1] - obs_y[0],
        linewidth=0.8, edgecolor="k", facecolor="dimgray", zorder=5,
    )
    ax.add_patch(rect)


def _title(t):
    return f"Square Cylinder — Re = 50   |   t = {t:.3f} s"


# ── static snapshot ──────────────────────────────────────────────────────────

def snapshot(result_path, t_target=-1, obs_x=(0.4, 0.6), obs_y=(0.65, 0.85), dpi=150):
    """Save a single 4-panel snapshot (p, u, v, |V|) at the closest time step."""

    result_path = Path(result_path)
    print(f"\nSnapshot from: {result_path.name}")

    p_frames, times, xP, yP = parse_probe(result_path / "Probe_1_Map.csv")
    u_frames, _,     xU, yU = parse_probe(result_path / "Probe_2u_Field.csv")
    v_frames, _,     xV, yV = parse_probe(result_path / "Probe_2v_Field.csv")
    spd = speed_frames(u_frames, v_frames)

    k = -1 if t_target < 0 else int(np.argmin(np.abs(times - t_target)))
    t_actual = times[k]
    print(f"  Using frame {k}  (t = {t_actual:.4f} s)")

    ext_p = [xP[0], xP[-1], yP[0], yP[-1]]
    ext_u = [xU[0], xU[-1], yU[0], yU[-1]]
    ext_v = [xV[0], xV[-1], yV[0], yV[-1]]

    fig, axes = plt.subplots(2, 2, figsize=(14, 7))
    fig.suptitle(f"Square Cylinder — Re = 50   |   t = {t_actual:.4f} s", fontsize=12)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.08, top=0.92,
                        hspace=0.38, wspace=0.35)

    panels = [
        (axes[0, 0], p_frames[k],  ext_p, "RdBu_r", "Pressure  p (Pa)"),
        (axes[0, 1], spd[k],       ext_p, "viridis", "|V|  (m/s)"),
        (axes[1, 0], u_frames[k],  ext_u, "RdBu_r", "u  (m/s)"),
        (axes[1, 1], v_frames[k],  ext_v, "RdBu_r", "v  (m/s)"),
    ]
    for ax, data, ext, cmap, label in panels:
        vabs = max(abs(data.min()), abs(data.max()), 1e-6)
        sym  = cmap == "RdBu_r"
        vmin = -vabs if sym else 0
        vmax = vabs  if sym else data.max()
        im = ax.imshow(data.T, cmap=cmap, origin="lower", aspect="auto",
                       extent=ext, vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(label, fontsize=8)
        ax.set_xlabel("x (m)", fontsize=8); ax.set_ylabel("y (m)", fontsize=8)
        ax.set_title(label.split("(")[0].strip(), fontsize=9)
        _add_obstacle(ax, obs_x, obs_y)

    out = result_path / f"Snapshot_SC_t{t_actual:.3f}.png"
    fig.savefig(str(out), dpi=dpi)
    print(f"Saved → {out}")
    plt.close(fig)


# ── CLI ──────────────────────────────────────────────────────────────────────

def _latest_sc(ioRes):
    matches = sorted(d for d in os.listdir(ioRes) if "exSC" in d)
    if not matches:
        print("No exSC* directory found in ioRes/"); sys.exit(1)
    return ioRes / matches[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result_dir", nargs="?", default=None,
                    help="Result dir name (under ioRes/) or absolute path")
    ap.add_argument("--fps",  type=int,   default=30,  help="Frames per second (default 30)")
    ap.add_argument("--dpi",  type=int,   default=150, help="Output DPI (default 150)")
    ap.add_argument("--obs",  type=float, nargs=4,
                    metavar=("X0", "X1", "Y0", "Y1"),
                    default=[0.4, 0.6, 0.65, 0.85],
                    help="Obstacle bounds (default 0.4 0.6 0.65 0.85)")
    ap.add_argument("--snapshot", type=float, default=None,
                    metavar="T",
                    help="Save a 4-panel snapshot at time T instead of animating "
                         "(use -1 for last frame)")
    args = ap.parse_args()

    ioRes = Path.cwd() / "ioRes"
    if args.result_dir is None:
        result = _latest_sc(ioRes)
    elif Path(args.result_dir).is_absolute():
        result = Path(args.result_dir)
    else:
        result = ioRes / args.result_dir

    obs_x = (args.obs[0], args.obs[1])
    obs_y = (args.obs[2], args.obs[3])

    if args.snapshot is not None:
        snapshot(result, t_target=args.snapshot, obs_x=obs_x, obs_y=obs_y, dpi=args.dpi)
    else:
        animate(result, fps=args.fps, dpi=args.dpi, obs_x=obs_x, obs_y=obs_y)


if __name__ == "__main__":
    main()
