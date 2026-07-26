#!/usr/bin/env python3
"""
PlotSC_MeshStudy.py — Mesh refinement comparison for the SC solver.

Generates per-mesh 4-panel snapshots (p, |V|, u, v) and side-by-side
comparison figures (pressure field, velocity magnitude, centreline profiles).

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
    ("Coarse (60×20)",  "20260724062905_exSC_coarse_implicit_Hybrid"),
    ("Medium (120×40)", "20260724062905_exSC_medium_implicit_Hybrid"),
    ("Fine (240×80)",   "20260724062905_exSC_fine_implicit_Hybrid"),
]

OBS_X  = (0.4, 0.6)
OBS_Y  = (0.65, 0.85)
COLORS = ["#2166ac", "#f46d43", "#1a9641"]


# ── I/O ─────────────────────────────────────────────────────────────────────

def parse_map(fpath):
    with open(fpath) as f:
        header = f.readline()
    cols = header.strip().split(",")
    pairs = [c.split() for c in cols[1:]]
    all_x = [float(p[0]) for p in pairs]
    all_y = [float(p[1]) for p in pairs]
    x0 = all_x[0]; ny = 1
    while ny < len(all_x) and all_x[ny] == x0:
        ny += 1
    nx = len(all_x) // ny
    xc = np.array([all_x[i * ny] for i in range(nx)])
    yc = np.array(all_y[:ny])
    raw = np.loadtxt(fpath, delimiter=",", skiprows=1)
    times  = raw[:, 0]
    frames = raw[:, 1:].reshape(len(times), nx, ny)
    return frames, times, xc, yc


def speed_frames(u_frames, v_frames, u_bc_east=1.0, v_bc_north=0.0):
    n, nx, ny = u_frames.shape
    east_pad  = np.full((n, 1, ny), u_bc_east)
    u_pad     = np.concatenate([u_frames, east_pad], axis=1)
    u_p       = 0.5 * (u_pad[:, :-1, :] + u_pad[:, 1:, :])

    north_pad = np.full((n, nx, 1), v_bc_north)
    v_pad     = np.concatenate([v_frames, north_pad], axis=2)
    v_p       = 0.5 * (v_pad[:, :, :-1] + v_pad[:, :, 1:])

    return np.sqrt(u_p**2 + v_p**2)


def load_mesh(result_dir):
    d = Path(result_dir)
    p_frames, times, xp, yp = parse_map(d / "Probe_1_Map.csv")
    u_frames, _,     xu, yu = parse_map(d / "Probe_2u_Field.csv")
    v_frames, _,     xv, yv = parse_map(d / "Probe_2v_Field.csv")
    spd = speed_frames(u_frames, v_frames)
    return dict(p=p_frames, u=u_frames, v=v_frames, spd=spd,
                times=times, xp=xp, yp=yp, xu=xu, yu=yu)


# ── helpers ──────────────────────────────────────────────────────────────────

def add_obstacle(ax):
    rect = mpatches.Rectangle(
        (OBS_X[0], OBS_Y[0]), OBS_X[1] - OBS_X[0], OBS_Y[1] - OBS_Y[0],
        linewidth=0.8, edgecolor="k", facecolor="dimgray", zorder=5,
    )
    ax.add_patch(rect)


def imshow_panel(fig, ax, data, ext, cmap, label, sym=True):
    vabs = max(abs(data.min()), abs(data.max()), 1e-6)
    vmin = -vabs if sym else 0.0
    vmax =  vabs if sym else data.max()
    im = ax.imshow(data.T, cmap=cmap, origin="lower", aspect="auto",
                   extent=ext, vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(label, fontsize=8)
    ax.set_xlabel("x (m)", fontsize=8)
    ax.set_ylabel("y (m)", fontsize=8)
    add_obstacle(ax)
    return im


# ── 1. Individual 4-panel snapshots ──────────────────────────────────────────

def plot_snapshots(datasets, out_dir, ioRes):
    for (label, dname), ds in zip(MESHES, datasets):
        t_end = ds["times"][-1]
        ext_p = [ds["xp"][0], ds["xp"][-1], ds["yp"][0], ds["yp"][-1]]
        ext_u = [ds["xu"][0], ds["xu"][-1], ds["yu"][0], ds["yu"][-1]]

        fig, axes = plt.subplots(2, 2, figsize=(14, 7))
        fig.suptitle(
            f"Square Cylinder — Re = 50   |   {label}   |   t = {t_end:.3f} s",
            fontsize=11)
        fig.subplots_adjust(left=0.07, right=0.97, bottom=0.08, top=0.92,
                            hspace=0.38, wspace=0.35)

        panels = [
            (axes[0, 0], ds["p"][-1],   ext_p, "RdBu_r", "p (Pa)",    True),
            (axes[0, 1], ds["spd"][-1], ext_p, "viridis", "|V| (m/s)", False),
            (axes[1, 0], ds["u"][-1],   ext_u, "RdBu_r", "u (m/s)",   True),
            (axes[1, 1], ds["v"][-1],   ext_u, "RdBu_r", "v (m/s)",   True),
        ]
        titles = ["Pressure", "Velocity Magnitude", "u", "v"]
        for (ax, data, ext, cmap, lbl, sym), title in zip(panels, titles):
            imshow_panel(fig, ax, data, ext, cmap, lbl, sym)
            ax.set_title(title, fontsize=9)

        tag = dname.split("_")[2]  # coarse / medium / fine
        out = out_dir / f"Snapshot_{tag}.png"
        fig.savefig(str(out), dpi=150)
        print(f"  Saved → {out}")
        plt.close(fig)


# ── 2. Side-by-side comparison: pressure + |V| ───────────────────────────────

def plot_comparison(datasets, out_dir):
    p_all   = np.concatenate([ds["p"][-1].ravel()   for ds in datasets])
    spd_all = np.concatenate([ds["spd"][-1].ravel() for ds in datasets])
    p_vabs  = max(abs(p_all.min()), abs(p_all.max()))
    s_max   = spd_all.max()

    for row_idx, (field_key, cmap, sym, vmin, vmax, cb_label, title_base) in enumerate([
        ("p",   "RdBu_r",  True,  -p_vabs, p_vabs, "p (Pa)",    "Pressure"),
        ("spd", "viridis", False, 0.0,     s_max,  "|V| (m/s)", "Velocity magnitude"),
    ]):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
        fig.suptitle(
            f"{title_base} — mesh refinement study   (Re = 50)", fontsize=11)
        fig.subplots_adjust(left=0.05, right=0.97, bottom=0.12,
                            top=0.88, wspace=0.35)

        for ax, (label, _), ds in zip(axes, MESHES, datasets):
            data = ds[field_key][-1]
            ext  = [ds["xp"][0], ds["xp"][-1], ds["yp"][0], ds["yp"][-1]]
            im = ax.imshow(data.T, cmap=cmap, origin="lower", aspect="auto",
                           extent=ext, vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label(cb_label, fontsize=8)
            ax.set_title(f"{label}\nt = {ds['times'][-1]:.2f} s", fontsize=9)
            ax.set_xlabel("x (m)", fontsize=8)
            ax.set_ylabel("y (m)", fontsize=8)
            add_obstacle(ax)

        fname = "MeshComparison_Pressure.png" if field_key == "p" \
                else "MeshComparison_Speed.png"
        out = out_dir / fname
        fig.savefig(str(out), dpi=150)
        print(f"  Saved → {out}")
        plt.close(fig)


# ── 3. Centreline profiles ────────────────────────────────────────────────────

def plot_centreline(datasets, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
    fig.suptitle(
        "Centreline profiles at y = 0.75 m — mesh refinement study   (Re = 50)",
        fontsize=11)
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.12, top=0.88, wspace=0.38)

    quantities = [
        ("p",   "p (Pa)",    True),
        ("u",   "u (m/s)",   True),
        ("spd", "|V| (m/s)", False),
    ]

    for ax, (key, ylabel, _) in zip(axes, quantities):
        for (label, _), ds, col in zip(MESHES, datasets, COLORS):
            xc = ds["xp"] if key in ("p", "spd") else ds["xu"]
            yc = ds["yp"] if key in ("p", "spd") else ds["yu"]
            j  = int(np.argmin(np.abs(yc - 0.75)))
            profile = ds[key][-1][:, j]
            ax.plot(xc, profile, color=col, linewidth=1.4, label=label)

        ax.axvline(OBS_X[0], color="k", linewidth=0.7, linestyle=":")
        ax.axvline(OBS_X[1], color="k", linewidth=0.7, linestyle=":")
        ax.set_xlabel("x (m)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(ylabel, fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, linewidth=0.4, alpha=0.5)

    out = out_dir / "MeshComparison_Centreline.png"
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(str(out), dpi=150)
    print(f"  Saved → {out}")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ioRes   = Path.cwd() / "ioRes"
    ioPlots = Path.cwd() / "ioPlots"
    out_dir = ioPlots / "MeshStudy_SC"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading data ...")
    datasets = []
    for label, dname in MESHES:
        print(f"  {label} ...")
        ds = load_mesh(ioRes / dname)
        print(f"    {len(ds['times'])} frames, t_end = {ds['times'][-1]:.3f} s")
        datasets.append(ds)

    print("\n[1/3] Individual snapshots ...")
    plot_snapshots(datasets, out_dir, ioRes)

    print("\n[2/3] Side-by-side field comparisons ...")
    plot_comparison(datasets, out_dir)

    print("\n[3/3] Centreline profiles ...")
    plot_centreline(datasets, out_dir)

    print("\nDone. All plots saved to", out_dir)


if __name__ == "__main__":
    main()
