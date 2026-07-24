#!/usr/bin/env python3
"""
PlotSC_MeshStudy.py — Mesh refinement comparison for the SC solver.

Reads Probe_1_Map.csv (pressure) from the three mesh result directories and
produces a side-by-side snapshot and a pressure-profile comparison along the
domain centreline (y=0.75).

Usage:
    cd /path/to/HMT
    source .venv/bin/activate
    python ioPlots/PlotSC_MeshStudy.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path


MESHES = [
    ("Coarse (60×20)",  "20260724045633_exSC_coarse_implicit_Hybrid"),
    ("Medium (120×40)", "20260724045633_exSC_medium_implicit_Hybrid"),
    ("Fine (240×80)",   "20260724045633_exSC_fine_implicit_Hybrid"),
]

OBS_X = (0.4, 0.6)
OBS_Y = (0.65, 0.85)


def parse_map(fpath):
    with open(fpath) as f:
        header = f.readline()
    cols = header.strip().split(",")
    pairs = [c.split() for c in cols[1:]]
    all_x = [float(p[0]) for p in pairs]
    all_y = [float(p[1]) for p in pairs]
    x0 = all_x[0]
    ny = 1
    while ny < len(all_x) and all_x[ny] == x0:
        ny += 1
    nx = len(all_x) // ny
    xc = np.array([all_x[i * ny] for i in range(nx)])
    yc = np.array(all_y[:ny])
    raw = np.loadtxt(fpath, delimiter=",", skiprows=1)
    times = raw[:, 0]
    frames = raw[:, 1:].reshape(len(times), nx, ny)
    return frames, times, xc, yc


def add_obstacle(ax):
    rect = mpatches.Rectangle(
        (OBS_X[0], OBS_Y[0]), OBS_X[1] - OBS_X[0], OBS_Y[1] - OBS_Y[0],
        linewidth=0.8, edgecolor="k", facecolor="dimgray", zorder=5,
    )
    ax.add_patch(rect)


def main():
    ioRes   = Path.cwd() / "ioRes"
    ioPlots = Path.cwd() / "ioPlots"
    out_dir = ioPlots / "MeshStudy_SC"
    out_dir.mkdir(parents=True, exist_ok=True)

    data = []
    for label, dname in MESHES:
        fpath = ioRes / dname / "Probe_1_Map.csv"
        frames, times, xc, yc = parse_map(fpath)
        data.append((label, frames, times, xc, yc))
        print(f"  {label}: {len(times)} frames, t_end={times[-1]:.3f} s, grid {len(xc)}×{len(yc)}")

    # ── 1. Side-by-side pressure field at final frame ─────────────────────────
    p_all = np.concatenate([d[1][-1].ravel() for d in data])
    vabs  = max(abs(p_all.min()), abs(p_all.max()))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
    fig.suptitle("Pressure field — mesh refinement study   (Re = 50)", fontsize=11)
    fig.subplots_adjust(left=0.05, right=0.97, bottom=0.12, top=0.88, wspace=0.35)

    for ax, (label, frames, times, xc, yc) in zip(axes, data):
        ext = [xc[0], xc[-1], yc[0], yc[-1]]
        im = ax.imshow(frames[-1].T, cmap="RdBu_r", origin="lower",
                       aspect="auto", extent=ext, vmin=-vabs, vmax=vabs)
        fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("p (Pa)", fontsize=8)
        ax.set_title(f"{label}\nt = {times[-1]:.2f} s", fontsize=9)
        ax.set_xlabel("x (m)", fontsize=8)
        ax.set_ylabel("y (m)", fontsize=8)
        add_obstacle(ax)

    out = out_dir / "MeshComparison_Pressure.png"
    fig.savefig(str(out), dpi=150)
    print(f"\nSaved → {out}")
    plt.close(fig)

    # ── 2. Centreline pressure profile p(x) at y = 0.75 ─────────────────────
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    colors = ["#2166ac", "#f46d43", "#1a9641"]
    ls     = ["--", "-.", "-"]

    for (label, frames, times, xc, yc), col, l in zip(data, colors, ls):
        j_mid = int(np.argmin(np.abs(yc - 0.75)))
        p_cl  = frames[-1][:, j_mid]
        ax2.plot(xc, p_cl, color=col, linestyle=l, linewidth=1.4, label=label)

    ax2.axvline(OBS_X[0], color="k", linewidth=0.7, linestyle=":")
    ax2.axvline(OBS_X[1], color="k", linewidth=0.7, linestyle=":", label="Cylinder faces")
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("p (Pa)")
    ax2.set_title("Centreline pressure profile at y = 0.75 m   (Re = 50)", fontsize=10)
    ax2.legend(fontsize=8)
    ax2.grid(True, linewidth=0.4, alpha=0.5)

    out2 = out_dir / "MeshComparison_CentrelinePressure.png"
    fig2.tight_layout()
    fig2.savefig(str(out2), dpi=150)
    print(f"Saved → {out2}")
    plt.close(fig2)


if __name__ == "__main__":
    main()
