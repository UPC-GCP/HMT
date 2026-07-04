"""
Smith-Hutton post-processing: 3-Pe scheme comparison and mesh refinement study.

Smith-Hutton benchmark (Smith & Hutton, 1982):
  Domain   : x ∈ [-1,1], y ∈ [0,1]
  Velocity : u = 2y(1-x²),  v = -2x(1-y²)   [divergence-free, rotating]
  Inlet    : φ = 1 + tanh(α(2x+1))  at y=0, x<0     (α=10)
  Outlet   : ∂φ/∂n = 0              at y=0, x>0
  Walls    : φ = 1 − tanh(10)       at x=±1, y=1
  Materials: ρ/γ = 10 (Pe1), 10³ (Pe2), 10⁶ (Pe3)

Analytical outlet reference for Pe→∞ (streamline symmetry, −x_in → x_out):
  φ_ref(x) = 1 + tanh(10·(1 − 2x))     for x ∈ [0, 1]

Usage:
  Scheme comparison (3 Pe levels, 3 schemes each):
    python PlotSmithHutton.py --compare \\
        <Pe1_UDS> <Pe1_Hybrid> <Pe1_CDS> \\
        <Pe2_UDS> <Pe2_Hybrid> <Pe2_CDS> \\
        <Pe3_UDS> <Pe3_Hybrid> <Pe3_CDS>

  Mesh refinement study (Hybrid, coarse → fine):
    python PlotSmithHutton.py --refine <N016_dir> <N032_dir> <N064_dir> <N128_dir>

Each <dir> is a subdirectory of ioRes/ produced by MainSolver.
"""

import sys
import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ── Reference values ─────────────────────────────────────────────────────────
def inlet_phi(x, alpha=10):
    """Dirichlet inlet at y=0, x∈[-1,0]."""
    return 1.0 + np.tanh(alpha * (2 * x + 1))

def outlet_ref_inf(x, alpha=10):
    """
    Analytical outlet profile for Pe→∞ (ρ/γ=10⁶).
    Streamline symmetry: a particle entering at x_in exits at x_out = -x_in.
    φ_inlet(x_in) = 1 + tanh(α(2x_in+1))
    → φ_outlet(x) = 1 + tanh(α(1−2x))     for x ∈ [0,1].
    """
    return 1.0 + np.tanh(alpha * (1.0 - 2.0 * x))

PHI_WALL = 1.0 - np.tanh(10)   # value on all non-inlet boundaries ≈ -0.9999
PE_LABELS = {0: r"$\rho/\gamma = 10$", 1: r"$\rho/\gamma = 10^3$", 2: r"$\rho/\gamma = 10^6$"}
SCHEME_COLORS = {"UDS": "tab:blue", "Hybrid": "tab:orange", "CDS": "tab:red"}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_boundary(res_dir):
    """
    Load Probe_*_Boundary.csv — boundary node values at y=0.
    Returns (x_all, phi_all) for the full bottom edge.
    Boundary probe saves bcPhi (Msh.boundaryConditions[k].Phi[node_index]).
    """
    pattern = os.path.join("ioRes", res_dir, "Probe_*_Boundary.csv")
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No Boundary CSV in ioRes/{res_dir}")

    df = pd.read_csv(paths[0], header=0, index_col=0)
    if df.empty:
        raise ValueError(f"Boundary probe empty in {res_dir} — check probe time range")

    coords = np.array([[float(v) for v in c.split()] for c in df.columns])
    x_all  = coords[:, 0]
    phi    = df.iloc[-1].values.astype(float)   # last row = steady state
    sort   = np.argsort(x_all)
    return x_all[sort], phi[sort]


def outlet_profile(res_dir):
    """phi(x, y=0) for x in [0, 1] from the Boundary probe."""
    x, phi = load_boundary(res_dir)
    mask = x >= -1e-9           # outlet: x >= 0 (include x=0 boundary)
    return x[mask], phi[mask]


def load_map(res_dir):
    """
    Load Probe_1_Map.csv — full interior phi field.
    Header uses true node positions (Msh.Nodes[0][j], Msh.Nodes[1][k]).
    Returns (x_nodes, y_nodes, phi[Nx, Ny]).
    """
    pattern = os.path.join("ioRes", res_dir, "Probe_1_Map.csv")
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No Map CSV in ioRes/{res_dir}")

    df = pd.read_csv(paths[0], header=0, index_col=0)
    coords   = np.array([[float(v) for v in c.split()] for c in df.columns])
    x_unique = np.unique(coords[:, 0])
    y_unique = np.unique(coords[:, 1])
    Nx, Ny   = len(x_unique), len(y_unique)
    phi_flat = df.iloc[-1].values.astype(float)
    phi      = phi_flat.reshape(Nx, Ny)
    return x_unique, y_unique, phi


def N_from_dir(res_dir):
    """Infer Nx from map probe."""
    x, _, phi = load_map(res_dir)
    return phi.shape[0]


# ── Richardson extrapolation ──────────────────────────────────────────────────

def l2_error(phi_coarse, x_coarse, phi_fine, x_fine):
    """L2 norm of phi_coarse interpolated to x_fine positions vs phi_fine."""
    phi_interp = np.interp(x_fine, x_coarse, phi_coarse)
    return float(np.sqrt(np.mean((phi_interp - phi_fine) ** 2)))


def richardson(vals, r=2.0):
    f1, f2, f3 = vals[-3], vals[-2], vals[-1]
    if abs(f1 - f2) < 1e-15 or abs(f2 - f3) < 1e-15:
        return None, f3, 0.0
    p   = np.log(abs(f1 - f2) / abs(f2 - f3)) / np.log(r)
    f_h = f3 + (f3 - f2) / (r**p - 1)
    GCI = 1.25 * abs(f2 - f3) / (abs(f3) * (r**p - 1))
    return p, f_h, GCI


# ── Mode: 3-Pe scheme comparison ─────────────────────────────────────────────

def _infer_scheme(dirname):
    for s in ("UDS", "Hybrid", "CDS"):
        if s.lower() in dirname.lower():
            return s
    return dirname


def mode_compare(dirs_flat):
    """
    dirs_flat: 9 directories in order Pe1_UDS Pe1_Hybrid Pe1_CDS
                                       Pe2_UDS Pe2_Hybrid Pe2_CDS
                                       Pe3_UDS Pe3_Hybrid Pe3_CDS
    Produces one figure with 3 subplots (one per Pe), each showing 3 schemes.
    """
    if len(dirs_flat) != 9:
        raise ValueError(f"--compare needs exactly 9 dirs (3 Pe × 3 schemes), got {len(dirs_flat)}")

    schemes  = ["UDS", "Hybrid", "CDS"]
    pe_names = [r"$\rho/\gamma = 10$", r"$\rho/\gamma = 10^3$", r"$\rho/\gamma = 10^6$"]
    pe_keys  = ["Pe1", "Pe2", "Pe3"]   # cosmetic only

    # Analytical reference (valid only for Pe3 → ∞)
    x_ref = np.linspace(0, 1, 300)
    phi_ref_inf = outlet_ref_inf(x_ref)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=False)

    for pe_idx in range(3):
        ax = axes[pe_idx]
        for sch_idx, scheme in enumerate(schemes):
            d    = dirs_flat[pe_idx * 3 + sch_idx]
            x, phi = outlet_profile(d)
            ax.plot(x, phi, color=list(SCHEME_COLORS.values())[sch_idx],
                    lw=1.8, label=scheme)

        # Analytical reference for Pe3 only
        if pe_idx == 2:
            ax.plot(x_ref, phi_ref_inf, "k--", lw=1.2, label="Ref. (Pe→∞)")

        ax.axhline(PHI_WALL, ls=":", color="gray", lw=0.8)
        ax.set_xlabel("x");  ax.set_title(pe_names[pe_idx])
        ax.set_xlim(0, 1);   ax.set_ylim(PHI_WALL - 0.1, 2.1)
        ax.legend(fontsize=8); ax.grid(True, ls=":", alpha=0.5)

    axes[0].set_ylabel(r"$\varphi(x,\, y=0)$")
    plt.suptitle("Smith-Hutton outlet profile — scheme comparison", fontsize=12)
    plt.tight_layout()
    fig.savefig("SH_scheme_comparison.png", dpi=150)
    print("Saved SH_scheme_comparison.png")

    # Separate 3×3 grid of 2-D φ fields
    fig2, axs = plt.subplots(3, 3, figsize=(13, 10))
    for pe_idx in range(3):
        for sch_idx, scheme in enumerate(schemes):
            d = dirs_flat[pe_idx * 3 + sch_idx]
            ax = axs[pe_idx][sch_idx]
            try:
                x_n, y_n, phi = load_map(d)
                cf = ax.contourf(x_n, y_n, phi.T, levels=20,
                                 cmap="RdYlBu_r", vmin=PHI_WALL, vmax=2.0)
                plt.colorbar(cf, ax=ax)
            except Exception as e:
                ax.text(0.5, 0.5, str(e), ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"{scheme} | {pe_names[pe_idx]}", fontsize=8)
            ax.set_xlabel("x", fontsize=7); ax.set_ylabel("y", fontsize=7)
    plt.suptitle("Smith-Hutton φ fields", fontsize=11)
    plt.tight_layout()
    fig2.savefig("SH_scheme_fields.png", dpi=150)
    print("Saved SH_scheme_fields.png")
    plt.show()


# ── Mode: mesh refinement ─────────────────────────────────────────────────────

def mode_refine(dirs):
    """Convergence study: outlet L2 error vs N (Hybrid)."""
    Ns        = []
    outlets   = []   # (x, phi) per case
    for d in dirs:
        Ns.append(N_from_dir(d))
        outlets.append(outlet_profile(d))

    # L2 error relative to finest mesh
    x_fine, phi_fine = outlets[-1]
    errors = []
    for x_c, phi_c in outlets[:-1]:
        errors.append(l2_error(phi_c, x_c, phi_fine, x_fine))

    print(f"\n{'Nx':>6}  {'L2 err vs finest':>18}")
    for N, e in zip(Ns[:-1], errors):
        print(f"  {N:4d}   {e:.6f}")
    print(f"  {Ns[-1]:4d}   (reference)")

    if len(dirs) >= 3:
        print("\n─── Richardson (L2 error) ───")
        p, f_h, GCI = richardson(errors + [0.0])
        # Use last 3 non-reference errors (errors has N-1 entries for N dirs)
        if len(errors) >= 2:
            r = 2.0
            p_simple = np.log(errors[-2] / errors[-1]) / np.log(r) if errors[-1] > 0 else None
            if p_simple:
                print(f"  Apparent order p  = {p_simple:.2f}  (Hybrid → expected ~1 at high Pe, ~2 at low Pe)")

    # Plot 1: outlet profiles for all N
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    for (x_c, phi_c), N in zip(outlets, Ns):
        ax1.plot(x_c, phi_c, label=f"Nx={N}")
    ax1.axhline(PHI_WALL, ls=":", color="gray", lw=1, label="φ_wall")
    ax1.set_xlabel("x"); ax1.set_ylabel("φ(x, y=0)")
    ax1.set_title("Outlet profile vs mesh size (Hybrid)")
    ax1.legend(fontsize=8); ax1.grid(True, ls=":")

    # Plot 2: convergence (L2 error vs h)
    hs = [2.0 / N for N in Ns[:-1]]
    ax2.loglog(hs, errors, "o-", label="L2 error vs finest")
    # Reference slopes
    h_ref = np.array([hs[0], hs[-1]])
    ax2.loglog(h_ref, errors[0] * (h_ref / h_ref[0]) ** 1, "k--", lw=1, label="slope 1")
    ax2.loglog(h_ref, errors[0] * (h_ref / h_ref[0]) ** 2, "k:",  lw=1, label="slope 2")
    ax2.set_xlabel("h = 2/Nx"); ax2.set_ylabel("L2 error (outlet)")
    ax2.set_title("Mesh convergence — Smith-Hutton Hybrid")
    ax2.legend(); ax2.grid(True, which="both", ls=":")
    ax2.invert_xaxis()

    plt.suptitle("Smith-Hutton mesh refinement", fontsize=11)
    plt.tight_layout()
    fig.savefig("SH_refine.png", dpi=150)
    print("\nSaved SH_refine.png")
    plt.show()


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--compare", nargs=9, metavar="DIR",
                     help="9 dirs: Pe1_UDS Pe1_Hybrid Pe1_CDS Pe2_UDS ... Pe3_CDS")
    grp.add_argument("--refine", nargs="+", metavar="DIR",
                     help="dirs for mesh refinement study (coarse → fine, Hybrid)")
    args = parser.parse_args()

    if args.compare:
        mode_compare(args.compare)
    elif args.refine:
        mode_refine(args.refine)
