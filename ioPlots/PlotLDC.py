"""
LDC post-processing: velocity profiles, Ghia validation, mesh refinement, delta sensitivity.

Usage:
  Refinement study  (pass dirs coarse → fine, all same delta):
    python PlotLDC.py --refine <dir1> <dir2> <dir3> <dir4>

  Delta sensitivity (pass dirs ordered by increasing delta, all same N):
    python PlotLDC.py --delta  <dir_d01> <dir_d1> <dir_d2> <dir_d3>

  Single-case validation:
    python PlotLDC.py --single <dir>

Each <dir> is a subdirectory of ioRes/ produced by LDCSolver.

Field probe coordinate system (from MeshNS.cpp lines 224–233)
─────────────────────────────────────────────────────────────
  u-field x: header = u.Faces[0][j] = p.Nodes[0][j-1]  (shifted ~Δx from true node)
             true   = u.Nodes[0][j] = p.Faces[0][j]
  u-field y: header = u.Faces[1][k] = p.Faces[1][k]    (face positions)
             true   = u.Nodes[1][k] = p.Nodes[1][k]    (≈ midpoint of k-th and k+1-th face)
  v-field x: header = v.Faces[0][j] = p.Faces[0][j]    (face positions)
             true   = v.Nodes[0][j] = p.Nodes[0][j]    (≈ midpoint)
  v-field y: header = v.Faces[1][k] = p.Nodes[1][k-1]  (shifted ~Δy from true node)
             true   = v.Nodes[1][k] = p.Faces[1][k]

We reconstruct true positions from header faces using midpoints where needed.
"""

import sys
import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ── Ghia et al. (1982) reference data ────────────────────────────────────────
# Table 1: u-velocity along x=0.5 (vertical centreline)
GHIA_U = {
    100:  {"y": [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,
                 0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0],
           "u": [1.0000,0.84123,0.78871,0.73722,0.68717,0.23151,0.00332,-0.13641,
                 -0.20581,-0.21090,-0.15662,-0.10150,-0.06434,-0.04775,-0.04192,-0.03717,0.0]},
    400:  {"y": [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,
                 0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0],
           "u": [1.0000,0.75837,0.68439,0.61756,0.55892,0.29093,0.16256,0.02135,
                 -0.11477,-0.17119,-0.32726,-0.24299,-0.14612,-0.10338,-0.09266,-0.08186,0.0]},
    1000: {"y": [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,
                 0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0],
           "u": [1.0000,0.65928,0.57492,0.51117,0.46604,0.33304,0.18719,0.05702,
                 -0.06080,-0.10648,-0.27805,-0.38289,-0.29730,-0.22220,-0.20196,-0.18109,0.0]},
    3200: {"y": [1.0000,0.9766,0.9688,0.9609,0.9531,0.8516,0.7344,0.6172,
                 0.5000,0.4531,0.2813,0.1719,0.1016,0.0703,0.0625,0.0547,0.0],
           "u": [1.0000,0.39560,0.33117,0.27612,0.23151,0.19188,0.12508,0.08344,
                 0.03717,0.00239,-0.12372,-0.21568,-0.31966,-0.42665,-0.51550,-0.39188,0.0]},
}
# Table 2: v-velocity along y=0.5 (horizontal centreline)
# Convention: v > 0 = upward (left side of clockwise primary vortex)
GHIA_V = {
    100:  {"x": [0.0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,
                 0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0],
           "v": [0.0,0.05906,0.07391,0.08864,0.10313,0.16914,0.22445,0.23186,
                 -0.05454,-0.17527,-0.17507,-0.16077,-0.12031,-0.10719,-0.09233,-0.07640,0.0]},
    400:  {"x": [0.0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,
                 0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0],
           "v": [0.0,0.12146,0.15663,0.19254,0.27669,0.29730,0.17117,0.10003,
                 -0.06198,-0.34779,-0.27305,-0.15991,-0.05765,-0.04294,-0.02826,-0.01385,0.0]},
    1000: {"x": [0.0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,
                 0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0],
           "v": [0.0,0.27485,0.29012,0.30353,0.32627,0.37095,0.33075,0.32235,
                 0.02526,-0.31966,-0.42665,-0.51550,-0.39188,-0.33714,-0.27669,-0.21388,0.0]},
    3200: {"x": [0.0,0.0625,0.0703,0.0781,0.0938,0.1563,0.2266,0.2344,
                 0.5,0.8047,0.8594,0.9063,0.9453,0.9531,0.9609,0.9688,1.0],
           "v": [0.0,0.20918,0.23827,0.27224,0.30054,0.30203,0.28124,0.27454,
                 0.00831,-0.36724,-0.47244,-0.39570,-0.28093,-0.26148,-0.24376,-0.21568,0.0]},
}


def _faces_to_nodes(faces):
    """Convert face positions to node (cell-centre) positions via midpoints."""
    f = np.asarray(faces)
    nodes = 0.5 * (f[:-1] + f[1:])
    return nodes


def load_field(res_dir, field="u"):
    """
    Load u or v Field probe CSV from ioRes/<res_dir>.

    Returns (x_nodes, y_nodes, phi) where phi[i, j] is the last saved row
    (steady state), and x_nodes / y_nodes are the TRUE node positions
    (corrected from the probe's face-based header).

    Uses tail -1 to avoid loading the full CSV into RAM (files can be 160 MB+).
    """
    import subprocess
    pattern = os.path.join("ioRes", res_dir, f"Probe_*{field}_Field.csv")
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No {field}-Field CSV in ioRes/{res_dir}")

    path = paths[0]
    with open(path) as f:
        header_line = f.readline().rstrip("\n")

    # Column names (drop the first "Time" index column)
    col_names = header_line.split(",")[1:]
    if not col_names:
        raise ValueError(f"Field probe is empty in {res_dir}")

    # Read only the last data row (memory-safe for large N=128 files)
    last_line = subprocess.check_output(["tail", "-1", path], text=True).rstrip("\n")
    phi_flat = np.array([float(v) for v in last_line.split(",")[1:]])

    # Parse header: every column is "x_header y_header"
    coords = np.array([[float(v) for v in c.split()] for c in col_names])
    x_hdr_uniq = np.unique(coords[:, 0])
    y_hdr_uniq = np.unique(coords[:, 1])
    Nx_hdr = len(x_hdr_uniq)
    Ny_hdr = len(y_hdr_uniq)

    phi = phi_flat.reshape(Nx_hdr, Ny_hdr)  # phi[i_x, i_y]

    if field == "u":
        # x-header: u.Faces[0][j] = p.Nodes[0][j-1]  (one slot left of true u.Nodes)
        #   → true u x-nodes ≈ one step to the right; extract from the POSITION of x=0.5
        #   We keep x_hdr as-is for column selection (argmin still finds correct column).
        x_hdr = x_hdr_uniq          # used only for selecting the centreline column
        # y-header: u.Faces[1][k] = p.Faces[1][k]  (face positions)
        #   → true u y-nodes = p.Nodes[1][k] = midpoints of consecutive faces
        y_nodes = _faces_to_nodes(np.append(y_hdr_uniq, 1.0))  # append domain top face
    else:  # field == "v"
        # x-header: v.Faces[0][j] = p.Faces[0][j]  (face positions)
        #   → true v x-nodes = p.Nodes[0][j] = midpoints
        x_hdr = x_hdr_uniq
        x_nodes_v = _faces_to_nodes(np.append(x_hdr_uniq, 1.0))
        # y-header: v.Faces[1][k] = p.Nodes[1][k-1]  (one slot left of true v.Nodes)
        #   → true v y-nodes ≈ face positions p.Faces[1][k]; keep y_hdr for row selection.
        y_nodes = y_hdr_uniq  # y_hdr used for centreline selection (argmin finds correct row)
        return x_nodes_v, y_hdr_uniq, phi

    return x_hdr, y_nodes, phi


def u_centreline(res_dir):
    """
    u-velocity along x=0.5 (vertical centreline).
    Returns (y_true, u_values).
    """
    x_hdr, y_nodes, phi = load_field(res_dir, "u")
    ic = int(np.argmin(np.abs(x_hdr - 0.5)))
    return y_nodes, phi[ic, :]


def v_centreline(res_dir):
    """
    v-velocity along y=0.5 (horizontal centreline).
    Returns (x_true, v_values).
    """
    x_nodes, y_hdr, phi = load_field(res_dir, "v")
    jc = int(np.argmin(np.abs(y_hdr - 0.5)))
    return x_nodes, phi[:, jc]


def scalar_QoI(res_dir):
    """u_min on vertical centreline and v_max/v_min on horizontal centreline."""
    y, u = u_centreline(res_dir)
    x, v = v_centreline(res_dir)
    return {"u_min": float(u.min()), "u_min_y": float(y[np.argmin(u)]),
            "v_max": float(v.max()), "v_min": float(v.min())}


def N_from_dir(res_dir):
    """Infer mesh size from the u-field: Nx = number of u-nodes in x."""
    x_hdr, _, phi = load_field(res_dir, "u")
    return phi.shape[0]


def richardson(vals, r=2.0):
    """Richardson extrapolation from last 3 values [coarse, medium, fine]."""
    f1, f2, f3 = vals[-3], vals[-2], vals[-1]
    if abs(f1 - f2) < 1e-14 or abs(f2 - f3) < 1e-14:
        return None, f3, 0.0
    p   = np.log(abs(f1 - f2) / abs(f2 - f3)) / np.log(r)
    f_h = f3 + (f3 - f2) / (r**p - 1)
    GCI = 1.25 * abs(f2 - f3) / (abs(f3) * (r**p - 1))
    return p, f_h, GCI


# ── Plotting helpers ──────────────────────────────────────────────────────────

def plot_centrelines(res_dirs, labels=None, Re=1000, title=""):
    """Centreline profiles with Ghia reference."""
    g_u = GHIA_U.get(Re); g_v = GHIA_V.get(Re)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    if g_u:
        ax1.plot(g_u["u"], g_u["y"], "kx", ms=5, label="Ghia et al. (1982)", zorder=5)
    if g_v:
        ax2.plot(g_v["x"], g_v["v"], "kx", ms=5, label="Ghia et al. (1982)", zorder=5)

    for i, d in enumerate(res_dirs):
        lbl = labels[i] if labels else d
        y, u = u_centreline(d)
        x, v = v_centreline(d)
        ax1.plot(u, y, label=lbl)
        ax2.plot(x, v, label=lbl)

    ax1.set_xlabel("u"); ax1.set_ylabel("y"); ax1.set_title(f"u @ x=0.5  {title}")
    ax1.legend(fontsize=7); ax1.grid(True, ls=":")
    ax2.set_xlabel("x"); ax2.set_ylabel("v"); ax2.set_title(f"v @ y=0.5  {title}")
    ax2.legend(fontsize=7); ax2.grid(True, ls=":")
    plt.tight_layout()
    return fig


# ── Mode: single validation ───────────────────────────────────────────────────

def mode_single(dirs, Re=1000, out_dir="ioRes/DataNew"):
    os.makedirs(out_dir, exist_ok=True)
    d = dirs[0]
    qoi = scalar_QoI(d)
    g = GHIA_U.get(Re, {})
    ref_umin = min(g.get("u", [0])) if g else None
    print(f"\n  {d}")
    print(f"  u_min = {qoi['u_min']:.5f}  at y = {qoi['u_min_y']:.4f}")
    if ref_umin:
        print(f"  Ghia  = {ref_umin:.5f}   error = {abs(qoi['u_min']-ref_umin)/abs(ref_umin)*100:.2f}%")
    print(f"  v_max = {qoi['v_max']:.5f}")
    print(f"  v_min = {qoi['v_min']:.5f}")
    fig = plot_centrelines([d], Re=Re, title=f"Re={Re}")
    out = os.path.join(out_dir, "LDC_single.png")
    fig.savefig(out, dpi=150); print(f"\nSaved {out}")
    plt.show()


# ── Mode: mesh refinement ─────────────────────────────────────────────────────

def mode_refine(dirs, Re=1000, out_dir="ioRes/DataNew", prefix="LDC_refine"):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{'N':>6}  {'u_min':>9}  {'v_max':>9}  {'v_min':>9}")
    Ns, u_mins = [], []
    for d in dirs:
        N   = N_from_dir(d)
        qoi = scalar_QoI(d)
        Ns.append(N); u_mins.append(qoi["u_min"])
        print(f"  {N:4d}   {qoi['u_min']:9.5f}   {qoi['v_max']:9.5f}   {qoi['v_min']:9.5f}")

    ref_umin = min(GHIA_U.get(Re, {}).get("u", [0]))

    if len(dirs) >= 3:
        print("\n─── Richardson Extrapolation (u_min) ───")
        p, f_h, GCI = richardson(u_mins)
        print(f"  Apparent order p  = {p:.2f}  (CDS → expected ~2)")
        print(f"  Richardson extrap = {f_h:.5f}")
        print(f"  GCI (fine)        = {GCI*100:.3f}%")
        if ref_umin:
            print(f"  Ghia reference    = {ref_umin:.5f}")
            print(f"  Error (extrap)    = {abs(f_h - ref_umin)/abs(ref_umin)*100:.2f}%")

    # Plot 1: centreline profiles
    labels = [f"N={N_from_dir(d)}" for d in dirs]
    fig1 = plot_centrelines(dirs, labels=labels, Re=Re, title=f"Re={Re} refinement")
    out1 = os.path.join(out_dir, f"{prefix}_profiles.png")
    fig1.savefig(out1, dpi=150)

    # Plot 2: QoI convergence
    fig2, ax = plt.subplots(figsize=(5, 4))
    ax.plot(Ns, u_mins, "o-", label="|u_min| computed")
    if ref_umin:
        ax.axhline(ref_umin, ls="-", color="gray", label=f"Ghia {ref_umin:.4f}")
        all_vals = u_mins + [ref_umin]
        margin = 0.02 * (max(all_vals) - min(all_vals))
        ax.set_ylim(min(all_vals) - margin, max(all_vals) + margin)
    ax.set_xlabel("N"); ax.set_ylabel("u_min"); ax.set_title(f"Mesh convergence — Re={Re}")
    ax.legend(); ax.grid(True, ls=":")
    out2 = os.path.join(out_dir, f"{prefix}_convergence.png")
    plt.tight_layout(); fig2.savefig(out2, dpi=150)

    print(f"\nSaved {out1}, {out2}")
    plt.show()


# ── Mode: delta sensitivity ───────────────────────────────────────────────────

def mode_delta(dirs, deltas, Re=1000, out_dir="ioRes/DataNew"):
    """deltas is a list of floats matching dirs, e.g. [0.1, 1, 2, 3]."""
    os.makedirs(out_dir, exist_ok=True)
    ref_umin = min(GHIA_U.get(Re, {}).get("u", [0]))
    print(f"\n{'delta':>7}  {'u_min':>9}  {'err %':>7}")
    u_mins = []
    for d, delta in zip(dirs, deltas):
        qoi = scalar_QoI(d)
        u_mins.append(qoi["u_min"])
        err = abs(qoi["u_min"] - ref_umin) / abs(ref_umin) * 100 if ref_umin else float("nan")
        print(f"  {delta:5.2f}   {qoi['u_min']:9.5f}   {err:7.2f}%")

    # Plot 1: centreline profiles
    labels = [f"δ={d}" for d in deltas]
    fig1 = plot_centrelines(dirs, labels=labels, Re=Re, title=f"Re={Re} delta study")
    out1 = os.path.join(out_dir, "LDC_delta_profiles.png")
    fig1.savefig(out1, dpi=150)

    # Plot 2: u_min vs delta
    fig2, ax = plt.subplots(figsize=(5, 4))
    ax.plot(deltas, u_mins, "s-")
    if ref_umin:
        ax.axhline(ref_umin, ls="-", color="gray", label=f"Ghia {ref_umin:.4f}")
        all_vals = u_mins + [ref_umin]
        margin = 0.02 * (max(all_vals) - min(all_vals))
        ax.set_ylim(min(all_vals) - margin, max(all_vals) + margin)
    ax.set_xlabel("delta (clustering parameter)"); ax.set_ylabel("u_min")
    ax.set_title(f"Non-uniform mesh sensitivity — Re={Re}, N=64")
    ax.legend(); ax.grid(True, ls=":")
    out2 = os.path.join(out_dir, "LDC_delta_sensitivity.png")
    plt.tight_layout(); fig2.savefig(out2, dpi=150)

    print(f"\nSaved {out1}, {out2}")
    plt.show()


# ── CLI ───────────────────────────────────────────────────────────────────────

# ── Mode: scheme comparison (CDS vs Hybrid at same N) ────────────────────────

def mode_scheme(dirs, labels, Re=1000, out_dir="ioRes/DataNew"):
    """
    Compare spatial schemes (e.g., CDS and Hybrid) at the same mesh.
    dirs:   list of result directories, one per scheme
    labels: list of scheme names matching dirs
    """
    os.makedirs(out_dir, exist_ok=True)
    ref_umin = min(GHIA_U.get(Re, {}).get("u", [0]))
    print(f"\n{'Scheme':>10}  {'u_min':>9}  {'err vs Ghia':>12}")
    for d, lbl in zip(dirs, labels):
        qoi = scalar_QoI(d)
        err = abs(qoi["u_min"] - ref_umin) / abs(ref_umin) * 100 if ref_umin else float("nan")
        print(f"  {lbl:8s}   {qoi['u_min']:9.5f}   {err:10.2f}%")

    # Centreline profiles for all schemes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    g_u = GHIA_U.get(Re); g_v = GHIA_V.get(Re)
    if g_u:
        ax1.plot(g_u["u"], g_u["y"], "kx", ms=5, label="Ghia et al. (1982)", zorder=5)
    if g_v:
        ax2.plot(g_v["x"], g_v["v"], "kx", ms=5, label="Ghia et al. (1982)", zorder=5)

    for d, lbl in zip(dirs, labels):
        y, u = u_centreline(d)
        x, v = v_centreline(d)
        ax1.plot(u, y, label=lbl)
        ax2.plot(x, v, label=lbl)

    ax1.set_xlabel("u"); ax1.set_ylabel("y")
    ax1.set_title(f"u @ x=0.5  Re={Re} scheme comparison")
    ax1.legend(fontsize=8); ax1.grid(True, ls=":")
    ax2.set_xlabel("x"); ax2.set_ylabel("v")
    ax2.set_title(f"v @ y=0.5  Re={Re} scheme comparison")
    ax2.legend(fontsize=8); ax2.grid(True, ls=":")
    plt.tight_layout()
    out = os.path.join(out_dir, "LDC_scheme_comparison.png")
    fig.savefig(out, dpi=150)
    print(f"\nSaved {out}")
    plt.show()


# ── Mode: 2D field plots (contours + streamlines) ────────────────────────────

def mode_field(dirs, Re=1000, out_dir="ioRes/DataNew"):
    """
    2D colour-map plots of u, v, and streamlines for one or more cases.
    Produces one PNG per directory: LDC_field_N<NNN>.png
    """
    from scipy.interpolate import RegularGridInterpolator
    os.makedirs(out_dir, exist_ok=True)

    for d in dirs:
        N = N_from_dir(d)
        x_u, y_u, U = load_field(d, "u")   # U[ix, iy], shape (Nx_u, Ny_u)
        x_v, y_v, V = load_field(d, "v")   # V[ix, iy], shape (Nx_v, Ny_v)

        # Common collocated grid for streamplot interpolation
        Ng = max(N * 2, 64)
        xg = np.linspace(0.0, 1.0, Ng)
        yg = np.linspace(0.0, 1.0, Ng)
        XG, YG = np.meshgrid(xg, yg, indexing="ij")   # shape (Ng, Ng)
        pts = np.column_stack([XG.ravel(), YG.ravel()])

        interp_u = RegularGridInterpolator((x_u, y_u), U,
                                           method="linear", bounds_error=False, fill_value=0.0)
        interp_v = RegularGridInterpolator((x_v, y_v), V,
                                           method="linear", bounds_error=False, fill_value=0.0)
        U_g = interp_u(pts).reshape(Ng, Ng)   # U_g[ix, iy]
        V_g = interp_v(pts).reshape(Ng, Ng)
        speed = np.sqrt(U_g**2 + V_g**2)

        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

        # ── Panel 1: velocity magnitude ───────────────────────────────────────
        ax = axes[0]
        cf = ax.contourf(xg, yg, speed.T, levels=20, cmap="viridis",
                         vmin=0, vmax=speed.max())
        plt.colorbar(cf, ax=ax, label="|u|  [m/s]")
        ax.set_title("Velocity magnitude"); ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

        # ── Panel 2: streamlines on velocity-magnitude background ─────────────
        ax = axes[1]
        cf = ax.contourf(xg, yg, speed.T, levels=20, cmap="viridis",
                         vmin=0, vmax=speed.max())
        plt.colorbar(cf, ax=ax, label="|u|  [m/s]")
        # streamplot expects u[j,i] (y-row, x-col)
        ax.streamplot(xg, yg, U_g.T, V_g.T,
                      color="white", linewidth=0.8, density=1.8, arrowsize=0.9)
        ax.set_title("Streamlines"); ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

        plt.suptitle(f"LDC velocity field — Re={Re}, N={N}", fontsize=11)
        plt.tight_layout()
        out = os.path.join(out_dir, f"LDC_field_N{N:03d}.png")
        fig.savefig(out, dpi=150)
        print(f"Saved {out}")
        plt.show()


def mode_multi_re(dirs_by_re, out_dir="ioRes/DataNew"):
    """
    Multi-Re validation: one dir per Re level. 2×N grid of (u, v) centreline plots,
    each panel overlaying Ghia reference. dirs_by_re: list of (Re, dir) tuples.
    """
    os.makedirs(out_dir, exist_ok=True)
    n = len(dirs_by_re)
    fig, axes = plt.subplots(n, 2, figsize=(10, 3.5 * n))
    if n == 1:
        axes = axes[np.newaxis, :]

    for row, (Re, d) in enumerate(dirs_by_re):
        g_u = GHIA_U.get(Re); g_v = GHIA_V.get(Re)
        y, u = u_centreline(d)
        x, v = v_centreline(d)

        ax1 = axes[row, 0]
        ax2 = axes[row, 1]
        if g_u:
            ax1.plot(g_u["u"], g_u["y"], "kx", ms=5, label="Ghia et al.", zorder=5)
        if g_v:
            ax2.plot(g_v["x"], g_v["v"], "kx", ms=5, label="Ghia et al.", zorder=5)

        ax1.plot(u, y, label=f"Computed (Re={Re})")
        ax2.plot(x, v, label=f"Computed (Re={Re})")
        ax1.set_xlabel("u"); ax1.set_ylabel("y")
        ax1.set_title(f"Re = {Re} — u @ x=0.5")
        ax1.legend(fontsize=7); ax1.grid(True, ls=":")
        ax2.set_xlabel("x"); ax2.set_ylabel("v")
        ax2.set_title(f"Re = {Re} — v @ y=0.5")
        ax2.legend(fontsize=7); ax2.grid(True, ls=":")

        qoi = scalar_QoI(d)
        ref_umin = min(g_u["u"]) if g_u else None
        err = abs(qoi["u_min"] - ref_umin) / abs(ref_umin) * 100 if ref_umin else float("nan")
        print(f"  Re={Re:5d}  u_min={qoi['u_min']:.5f}  err={err:.2f}%")

    plt.suptitle("LDC centreline profiles vs Ghia et al. (1982)", fontsize=12)
    plt.tight_layout()
    out = os.path.join(out_dir, "LDC_multi_re.png")
    fig.savefig(out, dpi=150)
    print(f"\nSaved {out}")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--refine", nargs="+", metavar="DIR")
    grp.add_argument("--delta",  nargs="+", metavar="DIR")
    grp.add_argument("--single", nargs=1,   metavar="DIR")
    grp.add_argument("--scheme", nargs="+", metavar="DIR",
                     help="dirs for scheme comparison (one per scheme, ordered CDS Hybrid ...)")
    grp.add_argument("--field", nargs="+", metavar="DIR",
                     help="2D field plots (u contour, v contour, streamlines) for each dir")
    grp.add_argument("--multi-re", nargs="+", metavar="DIR",
                     help="multi-Re validation: pass dirs in matching --re-vals order")
    parser.add_argument("--Re", type=int, default=1000, help="Reynolds number (default 1000)")
    parser.add_argument("--re-vals", nargs="+", type=int, default=[100, 400, 1000, 3200],
                        help="Re values matching --multi-re dirs (default: 100 400 1000 3200)")
    parser.add_argument("--deltas", nargs="+", type=float, default=[0.1, 1.0, 2.0, 3.0],
                        help="delta values matching --delta dirs (default: 0.1 1 2 3)")
    parser.add_argument("--labels", nargs="+", default=["CDS", "Hybrid"],
                        help="scheme labels matching --scheme dirs (default: CDS Hybrid)")
    parser.add_argument("--out-dir", default="ioRes/DataNew",
                        help="output directory for PNG files (default: ioRes/DataNew)")
    parser.add_argument("--prefix", default=None,
                        help="filename prefix for --refine plots (default: LDC_CDS_refine or LDC_refine)")
    args = parser.parse_args()

    if args.single:
        mode_single(args.single, Re=args.Re, out_dir=args.out_dir)
    elif args.refine:
        prefix = args.prefix if args.prefix else "LDC_refine"
        mode_refine(args.refine, Re=args.Re, out_dir=args.out_dir, prefix=prefix)
    elif args.delta:
        mode_delta(args.delta, args.deltas, Re=args.Re, out_dir=args.out_dir)
    elif args.scheme:
        mode_scheme(args.scheme, args.labels, Re=args.Re, out_dir=args.out_dir)
    elif args.field:
        mode_field(args.field, Re=args.Re, out_dir=args.out_dir)
    elif getattr(args, "multi_re"):
        dirs = getattr(args, "multi_re")
        mode_multi_re(list(zip(args.re_vals, dirs)), out_dir=args.out_dir)
