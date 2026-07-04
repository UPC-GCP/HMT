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
def main(dirs):
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
        plt.tight_layout(); plt.savefig("Nu_refinement.png", dpi=150)
        print("\nSaved Nu_refinement.png")

    # ── Plot 2: local Nu profile on hot wall (finest mesh) ───────────────────
    finest = results[-1]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(finest["Nu_local"], finest["y"], label=f"N={finest['N']}")
    ax.set_xlabel("Nu_local(y)"); ax.set_ylabel("y")
    ax.set_title("Local Nusselt — hot wall (x=0)")
    ax.grid(True, ls=":"); ax.legend()
    plt.tight_layout(); plt.savefig("Nu_local.png", dpi=150)
    print("Saved Nu_local.png")

    # ── Plot 3: T field (finest mesh) ────────────────────────────────────────
    T  = finest["T"]
    xn = finest["x"]; yn = finest["y"]
    fig, ax = plt.subplots(figsize=(5, 4))
    cf = ax.contourf(xn, yn, T.T, levels=20, cmap="RdBu_r")
    plt.colorbar(cf, ax=ax, label="T")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"Temperature field — N={finest['N']}")
    plt.tight_layout(); plt.savefig("T_field.png", dpi=150)
    print("Saved T_field.png")

    plt.show()


def mode_benchmark(dirs_by_ra):
    """
    Benchmark comparison: one dir per Ra level (Ra = 1e3, 1e4, 1e5, 1e6).
    Produces a bar chart of Nu_avg with De Vahl Davis reference and a 2x2 T-field grid.
    dirs_by_ra: list of (Ra, dir) tuples in increasing Ra order.
    """
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
    plt.tight_layout(); plt.savefig("Nu_benchmark.png", dpi=150)
    print("Saved Nu_benchmark.png")

    # 2×2 T-field grid
    fig2, axs = plt.subplots(2, 2, figsize=(9, 8))
    for idx, ((x, y, T), Ra) in enumerate(zip(fields, Ra_vals)):
        ax = axs[idx // 2][idx % 2]
        cf = ax.contourf(x, y, T.T, levels=20, cmap="RdBu_r", vmin=0, vmax=1)
        plt.colorbar(cf, ax=ax, label="T")
        ax.set_title(f"Ra = $10^{{{int(np.log10(Ra))}}}$")
        ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.suptitle("Temperature field — DHC benchmark")
    plt.tight_layout(); plt.savefig("T_benchmark.png", dpi=150)
    print("Saved T_benchmark.png")
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dirs", nargs="+")
    parser.add_argument("--benchmark", action="store_true",
                        help="4-Ra benchmark mode: pass 4 dirs in Ra=1e3,1e4,1e5,1e6 order")
    parser.add_argument("--ra", nargs="+", type=float,
                        help="Ra values matching --benchmark dirs (default: 1e3 1e4 1e5 1e6)")
    args = parser.parse_args()

    if args.benchmark:
        Ra_list = args.ra if args.ra else [1e3, 1e4, 1e5, 1e6]
        mode_benchmark(list(zip(Ra_list, args.dirs)))
    else:
        main(args.dirs)
