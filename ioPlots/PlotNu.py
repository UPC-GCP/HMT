"""
Nusselt number post-processing and mesh refinement analysis for DHC simulations.

Usage:
    python PlotNu.py <dir1> [dir2 dir3 ...]

Each dir is a subdirectory of ioRes/ produced by DHCSolver.
Multiple dirs are used for the mesh refinement study (ordered coarse → fine).

The T-map CSV (Probe_1_Map.csv) has:
  header:  Time, x0 y0, x0 y1, ..., x0 yN-1, x1 y0, ...   (outer=x, inner=y)
  rows:    t, T[0][0], T[0][1], ..., T[Nx-1][Ny-1]
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob


# ── De Vahl Davis (1983) reference values ─────────────────────────────────────
REFERENCE = {
    1e3: {"Nu_avg": 1.118, "Nu_max": 1.505, "v_max": 3.649,  "v_max_y": 0.813},
    1e4: {"Nu_avg": 2.243, "Nu_max": 3.528, "v_max": 19.617, "v_max_y": 0.823},
    1e5: {"Nu_avg": 4.519, "Nu_max": 7.717, "v_max": 68.590, "v_max_y": 0.855},
    1e6: {"Nu_avg": 8.800, "Nu_max": 17.925,"v_max": 220.560,"v_max_y": 0.850},
}


def load_T_map(res_dir):
    """Return (x_nodes, y_nodes, T_final) from the last row of Probe_1_Map.csv."""
    import subprocess
    path = glob.glob(os.path.join("ioRes", res_dir, "Probe_1_Map.csv"))
    if not path:
        raise FileNotFoundError(f"No Probe_1_Map.csv in ioRes/{res_dir}")

    # Read header only to get column names (avoids loading large file into RAM)
    with open(path[0]) as f:
        header_line = f.readline().rstrip("\n")
    columns = header_line.split(",")[1:]   # drop the index column name

    # Read only the last data row via tail (memory-efficient for large N=128 files)
    last_line = subprocess.check_output(["tail", "-1", path[0]], text=True).rstrip("\n")
    values = [float(v) for v in last_line.split(",")[1:]]   # drop timestamp

    # Parse node positions from column headers "x y"
    coords  = [tuple(map(float, c.split())) for c in columns]
    x_nodes = sorted(set(c[0] for c in coords))
    y_nodes = sorted(set(c[1] for c in coords))
    Nx, Ny  = len(x_nodes), len(y_nodes)

    T_flat = np.array(values)
    T = T_flat.reshape(Nx, Ny)           # T[i][j]

    return np.array(x_nodes), np.array(y_nodes), T


def compute_Nu(x_nodes, y_nodes, T, T_h=1.0, T_c=0.0, gamma=1.0):
    """
    Average Nusselt number on the hot (left) wall.

      Nu_local(j) = gamma * (T_h - T[0,j]) / dX_w
      Nu_avg      = integral Nu_local dy  (trapezoid)

    dX_w = x_nodes[0] - 0.0  (half-cell width at hot wall, since Faces[0]=0)
    """
    dX_w = x_nodes[0]                          # distance from wall to first node
    Nu_local = gamma * (T_h - T[0, :]) / dX_w  # shape (Ny,)
    Nu_avg = np.trapz(Nu_local, y_nodes) / (y_nodes[-1] - y_nodes[0])
    Nu_max = Nu_local.max()
    return Nu_avg, Nu_max, Nu_local


def richardson_analysis(Nu_list, r=2.0):
    """
    Given Nu values [coarse, medium, fine] and refinement ratio r,
    return (order p, Richardson extrapolation, GCI_fine).
    """
    f1, f2, f3 = Nu_list[-3], Nu_list[-2], Nu_list[-1]
    if abs(f1 - f2) < 1e-14 or abs(f2 - f3) < 1e-14:
        return None, f3, 0.0
    p   = np.log(abs(f1 - f2) / abs(f2 - f3)) / np.log(r)
    f_h = f3 + (f3 - f2) / (r**p - 1)          # Richardson extrapolation
    GCI = 1.25 * abs(f2 - f3) / abs(f3) / (r**p - 1)
    return p, f_h, GCI


# ── Main (refinement mode) ────────────────────────────────────────────────────
def main(dirs, out_dir="ioRes/DataNew"):
    os.makedirs(out_dir, exist_ok=True)
    results = []
    for d in dirs:
        x, y, T = load_T_map(d)
        Nx = len(x)
        Nu_avg, Nu_max, Nu_local = compute_Nu(x, y, T)
        results.append({"dir": d, "N": Nx, "Nu_avg": Nu_avg, "Nu_max": Nu_max,
                         "x": x, "y": y, "T": T, "Nu_local": Nu_local})
        print(f"  {d}  N={Nx:4d}  Nu_avg={Nu_avg:.4f}  Nu_max={Nu_max:.4f}")

    # ── Refinement convergence table ──────────────────────────────────────────
    if len(results) >= 3:
        print("\n─── Richardson Extrapolation ───")
        Nu_vals = [r["Nu_avg"] for r in results]
        p, f_h, GCI = richardson_analysis(Nu_vals)
        ref = REFERENCE.get(1e5, {}).get("Nu_avg")  # adjust Ra if needed
        print(f"  Apparent order p  = {p:.2f}  (expected ~1 for Hybrid/UDS)")
        print(f"  Richardson extrap = {f_h:.4f}")
        print(f"  GCI (fine)        = {GCI*100:.2f}%")
        if ref:
            print(f"  De Vahl Davis ref = {ref:.3f}  error = {abs(f_h-ref)/ref*100:.2f}%")

    # ── Plot 1: Nu_avg vs mesh size ───────────────────────────────────────────
    if len(results) > 1:
        Ns    = [r["N"] for r in results]
        Navgs = [r["Nu_avg"] for r in results]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(Ns, Navgs, "o-", label="Computed")
        ref = REFERENCE.get(1e5, {}).get("Nu_avg")
        if ref:
            ax.axhline(ref, ls="--", color="gray", label=f"De Vahl Davis {ref}")
        ax.set_xlabel("N (nodes per side)"); ax.set_ylabel("Nu_avg")
        ax.set_title("Mesh refinement — Ra=1e5")
        ax.legend(); ax.grid(True, ls=":")
        out = os.path.join(out_dir, "Nu_refinement.png")
        plt.tight_layout(); plt.savefig(out, dpi=150)
        print(f"\nSaved {out}")

    # ── Plot 2: local Nu profile on hot wall (finest mesh) ───────────────────
    finest = results[-1]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(finest["Nu_local"], finest["y"], label=f"N={finest['N']}")
    ax.set_xlabel("Nu_local(y)"); ax.set_ylabel("y")
    ax.set_title("Local Nusselt — hot wall (x=0)")
    ax.grid(True, ls=":"); ax.legend()
    out = os.path.join(out_dir, "Nu_local.png")
    plt.tight_layout(); plt.savefig(out, dpi=150)
    print(f"Saved {out}")

    # ── Plot 3: T field (finest mesh) ────────────────────────────────────────
    T  = finest["T"]
    xn = finest["x"]; yn = finest["y"]
    fig, ax = plt.subplots(figsize=(5, 4))
    cf = ax.contourf(xn, yn, T.T, levels=20, cmap="RdBu_r")
    plt.colorbar(cf, ax=ax, label="T")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"Temperature field — N={finest['N']}")
    out = os.path.join(out_dir, "T_field.png")
    plt.tight_layout(); plt.savefig(out, dpi=150)
    print(f"Saved {out}")

    plt.show()


def mode_benchmark(dirs_by_ra, out_dir="ioRes/DataNew"):
    """
    Benchmark comparison: one dir per Ra level (Ra = 1e3, 1e4, 1e5, 1e6).
    Produces a bar chart of Nu_avg with De Vahl Davis reference and a 2x2 T-field grid.
    dirs_by_ra: list of (Ra, dir) tuples in increasing Ra order.
    """
    os.makedirs(out_dir, exist_ok=True)
    Ra_vals, computed_avg, computed_max = [], [], []
    ref_avg, ref_max = [], []
    fields = []

    for Ra, d in dirs_by_ra:
        x, y, T = load_T_map(d)
        Nu_avg, Nu_max, _ = compute_Nu(x, y, T)
        Ra_vals.append(Ra)
        computed_avg.append(Nu_avg)
        computed_max.append(Nu_max)
        ref_avg.append(REFERENCE[Ra]["Nu_avg"])
        ref_max.append(REFERENCE[Ra]["Nu_max"])
        fields.append((x, y, T))
        print(f"  Ra={Ra:.0e}  N={len(x):3d}  Nu_avg={Nu_avg:.4f} (ref {REFERENCE[Ra]['Nu_avg']:.3f})"
              f"  err={abs(Nu_avg-REFERENCE[Ra]['Nu_avg'])/REFERENCE[Ra]['Nu_avg']*100:.2f}%")

    # Bar chart Nu_avg vs Ra
    Ra_labels = [f"$10^{int(np.log10(R))}$" for R in Ra_vals]
    x_pos = np.arange(len(Ra_vals))
    w = 0.35
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(x_pos - w/2, computed_avg, w, label="Computed",    color="tab:blue",  alpha=0.85)
    ax.bar(x_pos + w/2, ref_avg,      w, label="De Vahl Davis", color="tab:orange", alpha=0.85)
    ax.set_xticks(x_pos); ax.set_xticklabels(Ra_labels)
    ax.set_xlabel("Ra"); ax.set_ylabel("$\\overline{Nu}$")
    ax.set_title("Average Nusselt number — DHC benchmark")
    ax.legend(); ax.grid(True, ls=":", axis="y")
    out1 = os.path.join(out_dir, "Nu_benchmark.png")
    plt.tight_layout(); plt.savefig(out1, dpi=150)
    print(f"Saved {out1}")

    # 2×2 T-field grid
    fig2, axs = plt.subplots(2, 2, figsize=(9, 8))
    for idx, ((x, y, T), Ra) in enumerate(zip(fields, Ra_vals)):
        ax = axs[idx // 2][idx % 2]
        cf = ax.contourf(x, y, T.T, levels=20, cmap="RdBu_r", vmin=0, vmax=1)
        plt.colorbar(cf, ax=ax, label="T")
        ax.set_title(f"Ra = $10^{{{int(np.log10(Ra))}}}$")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.suptitle("Temperature field — DHC benchmark")
    out2 = os.path.join(out_dir, "T_benchmark.png")
    plt.tight_layout(); plt.savefig(out2, dpi=150)
    print(f"Saved {out2}")
    plt.show()


def _load_field_ns(res_dir, component="u"):
    """Load Probe_Nu_Field.csv or Probe_Nv_Field.csv (last row, memory-safe)."""
    import subprocess
    pattern = os.path.join("ioRes", res_dir, f"Probe_*{component}_Field.csv")
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No Field CSV for {component} in ioRes/{res_dir}")
    with open(paths[0]) as f:
        header = f.readline().rstrip("\n")
    col_names = header.split(",")[1:]
    last = subprocess.check_output(["tail", "-1", paths[0]], text=True).rstrip("\n")
    vals = np.array([float(v) for v in last.split(",")[1:]])
    coords = np.array([[float(v) for v in c.split()] for c in col_names])
    xs = np.unique(coords[:, 0])
    ys = np.unique(coords[:, 1])
    phi = vals.reshape(len(xs), len(ys))
    return xs, ys, phi


def mode_velocity(dirs_by_ra, out_dir="ioRes/DataNew"):
    """
    Velocity field maps for DHC: u, v, and streamlines for each Ra.
    dirs_by_ra: list of (Ra, dir) tuples (dirs must have Field probes).
    Produces DHC_ufield.png, DHC_vfield.png, DHC_streamlines.png.
    """
    from scipy.interpolate import RegularGridInterpolator
    os.makedirs(out_dir, exist_ok=True)

    n = len(dirs_by_ra)
    ncols = 2
    nrows = (n + 1) // 2

    fig_u,  axs_u  = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows))
    fig_v,  axs_v  = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows))
    fig_st, axs_st = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows))
    axs_u  = np.array(axs_u).reshape(nrows, ncols)
    axs_v  = np.array(axs_v).reshape(nrows, ncols)
    axs_st = np.array(axs_st).reshape(nrows, ncols)

    for idx, (Ra, d) in enumerate(dirs_by_ra):
        row, col = idx // ncols, idx % ncols
        Ra_label = f"Ra = $10^{{{int(np.log10(Ra))}}}$"

        x_u, y_u, U = _load_field_ns(d, "u")
        x_v, y_v, V = _load_field_ns(d, "v")

        # Project both onto a common collocated grid
        Ng = max(len(x_u) * 2, 64)
        xg = np.linspace(0.0, 1.0, Ng)
        yg = np.linspace(0.0, 1.0, Ng)
        XG, YG = np.meshgrid(xg, yg, indexing="ij")
        pts = np.column_stack([XG.ravel(), YG.ravel()])
        iU = RegularGridInterpolator((x_u, y_u), U, method="linear",
                                     bounds_error=False, fill_value=0.0)
        iV = RegularGridInterpolator((x_v, y_v), V, method="linear",
                                     bounds_error=False, fill_value=0.0)
        Ug = iU(pts).reshape(Ng, Ng)
        Vg = iV(pts).reshape(Ng, Ng)
        speed = np.sqrt(Ug**2 + Vg**2)

        # u panel
        ax = axs_u[row, col]
        cf = ax.contourf(xg, yg, Ug.T, levels=20, cmap="RdBu_r")
        fig_u.colorbar(cf, ax=ax, label="u")
        ax.set_title(Ra_label); ax.set_xlabel("x"); ax.set_ylabel("y")

        # v panel
        ax = axs_v[row, col]
        cf = ax.contourf(xg, yg, Vg.T, levels=20, cmap="RdBu_r")
        fig_v.colorbar(cf, ax=ax, label="v")
        ax.set_title(Ra_label); ax.set_xlabel("x"); ax.set_ylabel("y")

        # streamlines panel
        ax = axs_st[row, col]
        cf = ax.contourf(xg, yg, speed.T, levels=20, cmap="viridis")
        fig_st.colorbar(cf, ax=ax, label="|V|")
        ax.streamplot(xg, yg, Ug.T, Vg.T, color="white", linewidth=0.6,
                      density=1.2, arrowsize=0.8)
        ax.set_title(Ra_label); ax.set_xlabel("x"); ax.set_ylabel("y")
        print(f"  Ra={Ra:.0e}  u_max={Ug.max():.3f}  v_max={Vg.max():.3f}")

    # Hide unused panels
    for idx in range(n, nrows * ncols):
        axs_u.flat[idx].set_visible(False)
        axs_v.flat[idx].set_visible(False)
        axs_st.flat[idx].set_visible(False)

    for fig, fname, title in [
        (fig_u,  "DHC_ufield.png",      "Horizontal velocity u — DHC"),
        (fig_v,  "DHC_vfield.png",      "Vertical velocity v — DHC"),
        (fig_st, "DHC_streamlines.png", "Velocity streamlines — DHC"),
    ]:
        fig.suptitle(title, fontsize=12)
        fig.tight_layout()
        out = os.path.join(out_dir, fname)
        fig.savefig(out, dpi=150)
        print(f"Saved {out}")
    plt.show()


def mode_centreline(dirs_by_ra, out_dir="ioRes/DataNew"):
    """
    Centreline velocity profiles for DHC vs De Vahl Davis (1983) benchmark.
    Left panel : u(y) at x=0.5  — horizontal velocity, vertical centreline.
    Right panel: v(x) at y=0.5  — vertical velocity, horizontal centreline.
    dirs_by_ra : list of (Ra, dir) tuples with Field probes.
    """
    from scipy.interpolate import RegularGridInterpolator
    os.makedirs(out_dir, exist_ok=True)

    # De Vahl Davis (1983) Table 1 benchmark values:
    #   u_cl  = max |u| along vertical centreline x=0.5,  at y-position u_y
    #   v_cl  = max |v| along horizontal centreline y=0.5, at x-position v_x
    DVS = {
        1e3: {"u_cl": 3.544,   "u_y": 0.813, "v_cl": 3.649,   "v_x": 0.178},
        1e4: {"u_cl": 16.178,  "u_y": 0.823, "v_cl": 19.617,  "v_x": 0.119},
        1e5: {"u_cl": 34.730,  "u_y": 0.855, "v_cl": 68.590,  "v_x": 0.066},
        1e6: {"u_cl": 64.630,  "u_y": 0.850, "v_cl": 220.560, "v_x": 0.038},
    }

    Ng = 300
    xg = np.linspace(0.0, 1.0, Ng)
    yg = np.linspace(0.0, 1.0, Ng)
    XG, YG = np.meshgrid(xg, yg, indexing="ij")
    pts = np.column_stack([XG.ravel(), YG.ravel()])
    ix_half = np.argmin(np.abs(xg - 0.5))
    iy_half = np.argmin(np.abs(yg - 0.5))

    prop_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = {Ra: prop_cycle[i % len(prop_cycle)] for i, (Ra, _) in enumerate(dirs_by_ra)}

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    ref_handle = None

    for Ra, d in dirs_by_ra:
        c = colors[Ra]
        Ra_label = f"Ra = $10^{{{int(np.log10(Ra))}}}$"

        x_u, y_u, U = _load_field_ns(d, "u")
        x_v, y_v, V = _load_field_ns(d, "v")

        iU = RegularGridInterpolator((x_u, y_u), U, method="linear",
                                     bounds_error=False, fill_value=0.0)
        iV = RegularGridInterpolator((x_v, y_v), V, method="linear",
                                     bounds_error=False, fill_value=0.0)
        Ug = iU(pts).reshape(Ng, Ng)
        Vg = iV(pts).reshape(Ng, Ng)

        u_vcl = Ug[ix_half, :]   # u(y) at x=0.5
        v_hcl = Vg[:, iy_half]   # v(x) at y=0.5

        # Numerical profiles
        axes[0].plot(u_vcl, yg, color=c, linewidth=1.5, label=Ra_label)
        axes[1].plot(xg, v_hcl, color=c, linewidth=1.5, label=Ra_label)

        # De Vahl Davis reference markers (antisymmetric pair each)
        ref = DVS[Ra]
        h, = axes[0].plot(
            [ ref["u_cl"], -ref["u_cl"]],
            [ ref["u_y"],  1.0 - ref["u_y"]],
            "x", color=c, markersize=7, markeredgewidth=1.8)
        axes[1].plot(
            [ ref["v_x"],  1.0 - ref["v_x"]],
            [ ref["v_cl"], -ref["v_cl"]],
            "x", color=c, markersize=7, markeredgewidth=1.8)
        if ref_handle is None:
            ref_handle = h

    # Left panel
    ax = axes[0]
    ax.axvline(0, color="gray", linewidth=0.6, linestyle="--", zorder=0)
    ax.set_xlabel("u  [m/s]")
    ax.set_ylabel("y  [—]")
    ax.set_title("Horizontal velocity at $x = 0.5$")
    ax.set_ylim(0, 1)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(ref_handle)
    labels.append("De Vahl Davis (1983)")
    ax.legend(handles, labels, fontsize=8, loc="center left")

    # Right panel
    ax = axes[1]
    ax.axhline(0, color="gray", linewidth=0.6, linestyle="--", zorder=0)
    ax.set_xlabel("x  [—]")
    ax.set_ylabel("v  [m/s]")
    ax.set_title("Vertical velocity at $y = 0.5$")
    ax.set_xlim(0, 1)
    ax.legend(fontsize=8, loc="center right")

    fig.suptitle("DHC — Centreline velocity profiles vs De Vahl Davis (1983)", fontsize=11)
    fig.tight_layout()
    out = os.path.join(out_dir, "DHC_centreline.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dirs", nargs="+")
    parser.add_argument("--benchmark", action="store_true",
                        help="4-Ra benchmark mode: pass 4 dirs in Ra=1e3,1e4,1e5,1e6 order")
    parser.add_argument("--velocity", action="store_true",
                        help="velocity field mode: pass dirs with Field probes in Ra order")
    parser.add_argument("--centreline", action="store_true",
                        help="centreline profile mode: pass dirs with Field probes in Ra order")
    parser.add_argument("--ra", nargs="+", type=float,
                        help="Ra values matching --benchmark/--velocity/--centreline dirs")
    parser.add_argument("--out-dir", default="ioRes/DataNew",
                        help="output directory for PNG files (default: ioRes/DataNew)")
    args = parser.parse_args()

    if args.benchmark:
        Ra_list = args.ra if args.ra else [1e3, 1e4, 1e5, 1e6]
        mode_benchmark(list(zip(Ra_list, args.dirs)), out_dir=args.out_dir)
    elif args.velocity:
        Ra_list = args.ra if args.ra else [1e3, 1e4, 1e5, 1e6]
        mode_velocity(list(zip(Ra_list, args.dirs)), out_dir=args.out_dir)
    elif args.centreline:
        Ra_list = args.ra if args.ra else [1e3, 1e4, 1e5, 1e6]
        mode_centreline(list(zip(Ra_list, args.dirs)), out_dir=args.out_dir)
    else:
        main(args.dirs, out_dir=args.out_dir)
