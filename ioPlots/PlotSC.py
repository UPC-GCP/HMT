#!/usr/bin/env python3
"""
PlotSC.py — Animate Square Cylinder (SC) solver output.

Reads Probe_1_Map.csv (pressure), Probe_2u_Field.csv, Probe_2v_Field.csv
from a result directory and writes Animation_SC.mp4.

Usage:
    cd /path/to/HMT
    source .venv/bin/activate
    python ioPlots/PlotSC.py [result_dir] [--fps N] [--dpi N]
                             [--obs X0 X1 Y0 Y1] [--re N]
                             [--snapshot T] [--pressure]

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
    The C++ probe writes outer loop x, inner loop y — so y varies fastest.

    Returns
    -------
    frames : ndarray (n_frames, nx, ny)
    times  : ndarray (n_frames,)
    xcoords: ndarray (nx,)
    ycoords: ndarray (ny,)
    """
    t0 = time.time()
    name = Path(fpath).name
    print(f"  Reading {name} ...", end=" ", flush=True)

    with open(fpath) as f:
        header = f.readline().rstrip("\n")

    cols = header.split(",")
    pairs = [c.split(" ") for c in cols[1:]]
    all_x = [float(p[0]) for p in pairs]
    all_y = [float(p[1]) for p in pairs]

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

def speed_frames(u_frames, v_frames):
    """
    Interpolate face-centred u, v to cell centres (Neumann outflow at East/North)
    and return |V|.

    u_frames : (n, nx, ny)   — u-face values written by Field probe
    v_frames : (n, nx, ny)   — v-face values written by Field probe
    """
    n, nx, ny = u_frames.shape
    # Neumann outflow: duplicate last face value for the boundary
    u_pad = np.concatenate([u_frames, u_frames[:, -1:, :]], axis=1)
    u_p   = 0.5 * (u_pad[:, :-1, :] + u_pad[:, 1:, :])

    v_pad = np.concatenate([v_frames, v_frames[:, :, -1:]], axis=2)
    v_p   = 0.5 * (v_pad[:, :, :-1] + v_pad[:, :, 1:])

    return np.sqrt(u_p**2 + v_p**2)


# ── helpers ──────────────────────────────────────────────────────────────────

def _add_obstacle(ax, obs_x, obs_y):
    rect = mpatches.Rectangle(
        (obs_x[0], obs_y[0]),
        obs_x[1] - obs_x[0],
        obs_y[1] - obs_y[0],
        linewidth=0.8, edgecolor="k", facecolor="dimgray", zorder=5,
    )
    ax.add_patch(rect)


def _pcol(fig, ax, xc, yc, data, cmap, label, vmin, vmax, obs_x, obs_y):
    """pcolormesh panel — correctly handles non-uniform meshes."""
    im = ax.pcolormesh(xc, yc, data.T, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading='nearest')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(label, fontsize=9)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.set_ylabel("y (m)", fontsize=8)
    ax.set_xlim(xc[0], xc[-1])
    ax.set_ylim(yc[0], yc[-1])
    _add_obstacle(ax, obs_x, obs_y)
    return im


def _title(t, re):
    return f"Square Cylinder — Re = {re}   |   t = {t:.3f} s"


# ── animation ────────────────────────────────────────────────────────────────

def animate(result_path, plot_dir, fps=30, dpi=150,
            obs_x=(0.4, 0.6), obs_y=(0.65, 0.85),
            show_pressure=False, show_u=False, re=50):

    result_path = Path(result_path)
    print(f"\nResult directory: {result_path.name}")

    u_frames, times, xU, yU = parse_probe(result_path / "Probe_2u_Field.csv")
    v_frames, _,     _,  _  = parse_probe(result_path / "Probe_2v_Field.csv")
    if show_pressure:
        p_frames, _, xP, yP = parse_probe(result_path / "Probe_1_Map.csv")

    spd = speed_frames(u_frames, v_frames)

    n_frames = len(times)
    print(f"  Frames: {n_frames}  |  t = [{times[0]:.3f}, {times[-1]:.3f}] s")

    s_max = float(spd.max())
    if show_u:
        u_abs = max(abs(u_frames.min()), abs(u_frames.max()), 1e-6)

    ncols = 1 + int(show_pressure) + int(show_u)
    fig, axes = plt.subplots(1, ncols, figsize=(8 * ncols, 4.5), squeeze=False)
    fig.subplots_adjust(left=0.10, right=0.92, bottom=0.12, top=0.87, wspace=0.45)

    col = 0
    if show_pressure:
        p_abs = max(abs(p_frames.min()), abs(p_frames.max()), 1e-6)
        ax_p  = axes[0, col]
        im_p  = _pcol(fig, ax_p, xP, yP, p_frames[0], "RdBu_r", "p (Pa)",
                      -p_abs, p_abs, obs_x, obs_y)
        ax_p.set_title("Pressure", fontsize=10)
        col += 1

    if show_u:
        ax_u = axes[0, col]
        im_u = _pcol(fig, ax_u, xU, yU, u_frames[0], "RdBu_r", "u (m/s)",
                     -u_abs, u_abs, obs_x, obs_y)
        ax_u.set_title("u velocity", fontsize=10)
        col += 1

    ax_s = axes[0, col]
    im_s = _pcol(fig, ax_s, xU, yU, spd[0], "viridis", "|V| (m/s)",
                 0.0, s_max, obs_x, obs_y)
    ax_s.set_title("Velocity Magnitude", fontsize=10)

    title = fig.suptitle(_title(times[0], re), fontsize=11, y=0.97)

    def update(k):
        artists = []
        if show_pressure:
            im_p.set_array(p_frames[k].T.ravel())
            artists.append(im_p)
        if show_u:
            im_u.set_array(u_frames[k].T.ravel())
            artists.append(im_u)
        im_s.set_array(spd[k].T.ravel())
        title.set_text(_title(times[k], re))
        return artists + [im_s, title]

    print(f"\nRendering {n_frames} frames at {fps} fps, dpi={dpi} ...")
    t0 = time.time()
    ani_obj = anim.FuncAnimation(fig, update, frames=n_frames,
                                 interval=1000 / fps, blit=True, repeat=False)

    out = plot_dir / "Animation_SC.mp4"
    ani_obj.save(str(out), writer="ffmpeg", fps=fps, dpi=dpi,
                 extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    print(f"Saved → {out}  ({time.time()-t0:.0f}s)")
    plt.close(fig)


# ── static snapshot ──────────────────────────────────────────────────────────

def snapshot(result_path, plot_dir, t_target=-1,
             obs_x=(0.4, 0.6), obs_y=(0.65, 0.85), dpi=150, re=50):
    """Save a 4-panel snapshot (p, |V|, u, v) at the closest time step."""
    result_path = Path(result_path)
    print(f"\nSnapshot from: {result_path.name}")

    p_frames, times, xP, yP = parse_probe(result_path / "Probe_1_Map.csv")
    u_frames, _,     xU, yU = parse_probe(result_path / "Probe_2u_Field.csv")
    v_frames, _,     xV, yV = parse_probe(result_path / "Probe_2v_Field.csv")
    spd = speed_frames(u_frames, v_frames)

    k = -1 if t_target < 0 else int(np.argmin(np.abs(times - t_target)))
    t_actual = times[k]
    print(f"  Using frame {k}  (t = {t_actual:.4f} s)")

    fig, axes = plt.subplots(2, 2, figsize=(14, 7))
    fig.suptitle(f"Square Cylinder — Re = {re}   |   t = {t_actual:.4f} s", fontsize=12)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.08, top=0.92,
                        hspace=0.38, wspace=0.35)

    panels = [
        (axes[0, 0], xP, yP, p_frames[k],  "RdBu_r",  "Pressure  p (Pa)", True),
        (axes[0, 1], xU, yU, spd[k],        "viridis",  "|V|  (m/s)",       False),
        (axes[1, 0], xU, yU, u_frames[k],   "RdBu_r",  "u  (m/s)",         True),
        (axes[1, 1], xV, yV, v_frames[k],   "RdBu_r",  "v  (m/s)",         True),
    ]
    titles = ["Pressure", "|V|", "u", "v"]
    for (ax, xc, yc, data, cmap, label, sym), ttl in zip(panels, titles):
        vabs = max(abs(data.min()), abs(data.max()), 1e-6)
        vmin = -vabs if sym else 0.0
        vmax =  vabs if sym else data.max()
        _pcol(fig, ax, xc, yc, data, cmap, label, vmin, vmax, obs_x, obs_y)
        ax.set_title(ttl, fontsize=9)

    out = plot_dir / f"Snapshot_SC_t{t_actual:.3f}.png"
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
    ap.add_argument("--re",   type=int,   default=50,  help="Reynolds number for title (default 50)")
    ap.add_argument("--obs",  type=float, nargs=4,
                    metavar=("X0", "X1", "Y0", "Y1"),
                    default=[0.4, 0.6, 0.65, 0.85],
                    help="Obstacle bounds (default 0.4 0.6 0.65 0.85)")
    ap.add_argument("--snapshot", type=float, default=None,
                    metavar="T",
                    help="Save a 4-panel snapshot at time T instead of animating "
                         "(use -1 for last frame)")
    ap.add_argument("--pressure", action="store_true",
                    help="Add a pressure panel to the animation")
    ap.add_argument("--u", action="store_true",
                    help="Add a u-velocity panel to the animation")
    args = ap.parse_args()

    ioRes   = Path.cwd() / "ioRes"
    ioPlots = Path.cwd() / "ioPlots"
    if args.result_dir is None:
        result = _latest_sc(ioRes)
    elif Path(args.result_dir).is_absolute():
        result = Path(args.result_dir)
    else:
        result = ioRes / args.result_dir

    plot_dir = ioPlots / result.name
    plot_dir.mkdir(parents=True, exist_ok=True)

    obs_x = (args.obs[0], args.obs[1])
    obs_y = (args.obs[2], args.obs[3])

    if args.snapshot is not None:
        snapshot(result, plot_dir, t_target=args.snapshot,
                 obs_x=obs_x, obs_y=obs_y, dpi=args.dpi, re=args.re)
    else:
        animate(result, plot_dir, fps=args.fps, dpi=args.dpi,
                obs_x=obs_x, obs_y=obs_y,
                show_pressure=args.pressure, show_u=args.u, re=args.re)


if __name__ == "__main__":
    main()
