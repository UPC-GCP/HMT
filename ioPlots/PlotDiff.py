"""
PlotDiff.py — Post-processing for Diffusion chapter (1D Rod + 4Materials)

Usage examples:
  python3 PlotDiff.py --solver-cmp CG_dir GS_dir --out-dir ioRes/DataNew
  python3 PlotDiff.py --N-sweep dir1 dir2 ... --out-dir ioRes/DataNew
  python3 PlotDiff.py --dt-sweep-impl d1 d2 ... --dt-sweep-cn d1 d2 ... --dt-sweep-expl d --out-dir ioRes/DataNew
  python3 PlotDiff.py --bc-val ND_dir CD_dir DN_dir DC_dir --out-dir ioRes/DataNew
  python3 PlotDiff.py --mat4-refine dir1 dir2 dir3 dir4 --out-dir ioRes/DataNew
  python3 PlotDiff.py --mat4-dt-sweep impl10 impl50 impl500 cn10 cn50 cn500 expl --out-dir ioRes/DataNew
  python3 PlotDiff.py --mat4-ref dir --out-dir ioRes/DataNew
  python3 PlotDiff.py --mat4-colormap dir --out-dir ioRes/DataNew
"""

import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from glob import glob

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "lines.linewidth": 1.5,
    "figure.dpi": 150,
})

COLORS   = plt.rcParams["axes.prop_cycle"].by_key()["color"]
RES_BASE = Path("ioRes")

# ── I/O helpers ───────────────────────────────────────────────────────────────
def res_path(d):
    """Return absolute path to a result directory."""
    p = Path(d) if Path(d).is_absolute() else RES_BASE / d
    if not p.exists():
        raise FileNotFoundError(f"Result dir not found: {p}")
    return p

def load_debug(res_dir):
    """Return (t, iters, res) from Probe_*_Bug.csv"""
    p = res_path(res_dir)
    f = next(iter(sorted(p.glob("Probe_*_Bug.csv"))), None)
    if f is None:
        raise FileNotFoundError(f"No Bug probe in {p}")
    d = np.genfromtxt(f, delimiter=",", skip_header=1)
    if d.ndim == 1:
        d = d[np.newaxis, :]
    return d[:, 0], d[:, 1].astype(int), d[:, 2]

def load_map(res_dir, probe_idx=None):
    """
    Load a Map probe CSV.  probe_idx=None → first Map probe found.
    Returns (t_arr, x_arr, T_2d) where T_2d[i] = T profile at time i.
    Columns interleaved: x0y0, x0y1, x1y0, x1y1, ...
    Returns the y-averaged 1D profile.
    """
    p = res_path(res_dir)
    if probe_idx is not None:
        f = p / f"Probe_{probe_idx}_Map.csv"
    else:
        candidates = sorted(p.glob("Probe_*_Map.csv"))
        if not candidates:
            raise FileNotFoundError(f"No Map probe in {p}")
        f = candidates[0]

    with open(f) as fh:
        header = fh.readline().strip().split(",")
    # header[0] = "Time", header[1:] = "x y" strings
    coords = [(float(h.split()[0]), float(h.split()[1])) for h in header[1:]]
    x_all  = np.array([c[0] for c in coords])
    y_all  = np.array([c[1] for c in coords])

    data = np.genfromtxt(f, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    t_arr = data[:, 0]

    # Average over y to get 1D profile
    x_unique = np.unique(x_all)
    T_1d = np.zeros((len(t_arr), len(x_unique)))
    for k, xv in enumerate(x_unique):
        mask = x_all == xv
        T_1d[:, k] = data[:, 1:][:, mask].mean(axis=1)

    return t_arr, x_unique, T_1d

def load_all_maps(res_dir):
    """Load all Map probes in a directory, return list of (t, x, T_1d)."""
    p = res_path(res_dir)
    maps = []
    for f in sorted(p.glob("Probe_*_Map.csv")):
        with open(f) as fh:
            header = fh.readline().strip().split(",")
        coords = [(float(h.split()[0]), float(h.split()[1])) for h in header[1:]]
        x_all  = np.array([c[0] for c in coords])
        y_all  = np.array([c[1] for c in coords])
        data = np.genfromtxt(f, delimiter=",", skip_header=1)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        t_arr = data[:, 0]
        x_unique = np.unique(x_all)
        T_1d = np.zeros((len(t_arr), len(x_unique)))
        for k, xv in enumerate(x_unique):
            mask = x_all == xv
            T_1d[:, k] = data[:, 1:][:, mask].mean(axis=1)
        maps.append((t_arr, x_unique, T_1d))
    return maps

def load_point(res_dir):
    """Return (t, x_labels, T_matrix) from Probe_0_Point.csv."""
    f = res_path(res_dir) / "Probe_0_Point.csv"
    with open(f) as fh:
        header = fh.readline().strip().split(",")
    labels = header[1:]
    data = np.genfromtxt(f, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data[:, 0], labels, data[:, 1:]

def config_name(res_dir):
    """Extract short config tag from directory name."""
    return Path(str(res_dir)).name.split("_", 1)[1] if "_" in Path(str(res_dir)).name else str(res_dir)

def save(fig, out_dir, fname):
    os.makedirs(out_dir, exist_ok=True)
    path = Path(out_dir) / fname
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")

# ── Analytical solutions ───────────────────────────────────────────────────────
# Source term used in all 1D-rod configs
ROD_S = 200.0   # [W/m³]  internal heat generation

def T_rod_ss(x, T_L=100.0, T_R=0.0, L=1.0, S=ROD_S, gamma=1.0):
    """Dirichlet–Dirichlet steady state with uniform source S.
    Solves γ T'' + S = 0, T(0)=T_L, T(L)=T_R.
    T_ss = T_L + (T_R-T_L)*x/L + S/(2γ) * (x/L)*(1 - x/L)
    """
    xi = x / L
    return T_L + (T_R - T_L) * xi + (S / (2.0 * gamma)) * xi * (1.0 - xi)

def T_rod_transient_fast(x, t, T_L=100.0, T_R=0.0, T0=0.0, L=1.0, alpha=1.0,
                         S=ROD_S, gamma=1.0, n_terms=200):
    """Vectorised transient solution for 1D rod (Dirichlet–Dirichlet) with source S.

    T(x,t) = T_ss(x) + Σ_n B_n sin(nπx/L) exp(-α(nπ/L)²t)

    T_ss = T_L + (T_R-T_L)x/L + S/(2γ) x/L (1-x/L)   [parabolic]

    B_n = (2/L) ∫₀ᴸ [T₀ - T_ss(x)] sin(nπx/L) dx
        = B_n_linear + B_n_source
    B_n_linear  = closed-form from T_L, T_R, T0
    B_n_source  = -S/γ · 2[1-(-1)ⁿ] / (nπ)³   (from ∫x(L-x)sin(nπx/L)dx)
    """
    T_ss = T_rod_ss(x, T_L, T_R, L, S, gamma)
    ns   = np.arange(1, n_terms + 1, dtype=float)[:, None]  # (n,1)
    xv   = x[None, :]                                         # (1, Nx)
    lam  = ns * np.pi / L
    cos_n = (-1.0) ** ns

    # B_n for the linear part of T_ss (T_L + (T_R-T_L)x - T0)
    B_n = (
        -2.0 * T_L * (1.0 - cos_n) / (ns * np.pi)
        - 2.0 * (T_R - T_L) * (-cos_n) / (ns * np.pi)
        + 2.0 * T0 * (1.0 - cos_n) / (ns * np.pi)
    )
    # Source correction: ∫₀¹ x(1-x) sin(nπx) dx = 2[1-(-1)ⁿ]/(nπ)³
    B_n += -(S / gamma) * 2.0 * (1.0 - cos_n) / (ns * np.pi) ** 3

    mode = B_n * np.sin(lam * xv) * np.exp(-alpha * lam**2 * t)
    return T_ss + mode.sum(axis=0)

# ── 4Materials analytical steady state ───────────────────────────────────────
MATS_4 = [
    {"gamma": 170, "x0": 0.00, "x1": 0.25},
    {"gamma": 140, "x0": 0.25, "x1": 0.50},
    {"gamma": 200, "x0": 0.50, "x1": 0.75},
    {"gamma": 140, "x0": 0.75, "x1": 1.00},
]

def T_mat4_ss(x, T_L=100.0, T_R=0.0):
    """Piecewise-linear steady state for 4-material rod."""
    # Compute heat flux: q = ΔT / R_total, R_i = L_i/γ_i
    R_total = sum((m["x1"] - m["x0"]) / m["gamma"] for m in MATS_4)
    q = (T_L - T_R) / R_total  # heat flux [W/m²]

    T_out = np.empty_like(x, dtype=float)
    T_prev = T_L
    for m in MATS_4:
        mask = (x >= m["x0"]) & (x <= m["x1"])
        slope = q / m["gamma"]
        T_out[mask] = T_prev - slope * (x[mask] - m["x0"])
        T_prev = T_prev - slope * (m["x1"] - m["x0"])  # always advance
    return T_out

def mat4_interfaces():
    return [m["x0"] for m in MATS_4[1:]]


# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — CG vs GS solver performance
# ══════════════════════════════════════════════════════════════════════════════
def mode_solver_cmp(dirs, out_dir):
    labels = ["CG", "GS"]
    styles = ["-", "--"]
    colors = [COLORS[0], COLORS[1]]

    fig, axes = plt.subplots(2, 1, figsize=(7, 5.5), sharex=True)
    ax_iter, ax_res = axes

    for i, d in enumerate(dirs):
        t, iters, res = load_debug(d)
        mask = iters > 0
        ax_iter.plot(t[mask], iters[mask], styles[i], color=colors[i], label=labels[i])
        ax_res.semilogy(t[mask], res[mask],  styles[i], color=colors[i], label=labels[i])

    ax_iter.set_ylabel("Iterations per step")
    ax_iter.legend()
    ax_iter.grid(True, ls=":", alpha=0.5)

    ax_res.set_ylabel("Residual")
    ax_res.set_xlabel("Time [s]")
    ax_res.legend()
    ax_res.grid(True, ls=":", which="both", alpha=0.5)

    fig.suptitle("1D Rod — CG vs GS solver comparison", fontsize=12)
    fig.tight_layout()
    save(fig, out_dir, "Diff_solver_cmp.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — Numerical parameter study (N sweep + dt sweep + scheme comparison)
# ══════════════════════════════════════════════════════════════════════════════
def mode_N_sweep(dirs, N_vals, out_dir):
    """Spatial convergence: L2 error vs N.
    We compare against the analytical transient T(x, t_sim) so the spatial
    discretisation error (O(Δx²)) is visible even for a near-linear profile.
    """
    errors = []
    for d in dirs:
        t_arr, x, T_1d = load_map(d)
        T_num  = T_1d[-1]
        t_sim  = t_arr[-1]
        T_anal = T_rod_transient_fast(x, t_sim)
        l2 = np.sqrt(np.mean((T_num - T_anal)**2))
        errors.append(l2)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.loglog(N_vals, errors, "o-", color=COLORS[0], label="L₂ error (implicit)")

    # Reference slopes
    N_ref = np.array([N_vals[0], N_vals[-1]], dtype=float)
    for p, ls, lab in [(1, "--", "O(h)"), (2, ":", "O(h²)")]:
        scale = errors[0] * (N_vals[0] / N_ref[0]) ** p
        ax.loglog(N_ref, scale * (N_vals[0] / N_ref) ** p, ls, color="gray", label=lab)

    ax.set_xlabel("N (mesh nodes)")
    ax.set_ylabel("L₂ error vs analytical T(x,t)")
    ax.set_title("Spatial convergence — 1D Rod (implicit, t=0.1 s)")
    ax.legend()
    ax.grid(True, ls=":", which="both", alpha=0.5)
    fig.tight_layout()
    save(fig, out_dir, "Diff_N_sweep.png")

    # Print table
    print("\n  Spatial convergence:")
    print(f"  {'N':>6}  {'L2 error':>12}  {'Rate':>6}")
    for i, (n, e) in enumerate(zip(N_vals, errors)):
        rate = np.log2(errors[i-1]/e) / np.log2(N_vals[i]/N_vals[i-1]) if i > 0 else float("nan")
        print(f"  {n:>6}  {e:>12.4e}  {rate:>6.2f}")


def mode_dt_sweep(impl_dirs, impl_dts, cn_dirs, cn_dts, expl_dir, out_dir):
    """Temporal convergence at t*=0.1, N=128.
    Implicit: O(Δt) — error grows with dt.
    CN: flat at spatial error floor (CN temporal error << spatial error for r<0.5).
    """
    alpha, T_L, T_R, T0 = 1.0, 100.0, 0.0, 0.0

    def get_error(d):
        t_arr, x, T_1d = load_map(d)
        T_num  = T_1d[-1]
        t_sim  = t_arr[-1]
        T_anal = T_rod_transient_fast(x, t_sim, T_L, T_R, T0, alpha=alpha)
        return np.sqrt(np.mean((T_num - T_anal)**2))

    err_impl = [get_error(d) for d in impl_dirs]
    err_cn   = [get_error(d) for d in cn_dirs]

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(impl_dts, err_impl, "o-",  color=COLORS[0], label="Implicit")
    ax.loglog(cn_dts,   err_cn,   "s--", color=COLORS[1], label="Crank-Nicolson (r<0.5, spatial floor)")

    # O(Δt) reference line through smallest implicit point
    dt_ref = np.array([min(impl_dts), max(impl_dts)])
    ax.loglog(dt_ref, err_impl[0]*(dt_ref/impl_dts[0])**1, ":", color=COLORS[0], alpha=0.6, label="O(Δt)")

    # CN spatial error floor (dashed horizontal)
    e_floor = np.mean(err_cn)
    ax.axhline(e_floor, color=COLORS[1], ls=":", alpha=0.5)

    ax.set_xlabel("Δt [s]")
    ax.set_ylabel("L₂ error vs T_analytical(t_sim)")
    ax.set_title("Temporal convergence — 1D Rod (N=128, t*=0.1 s)")
    ax.legend()
    ax.grid(True, ls=":", which="both", alpha=0.5)
    fig.tight_layout()
    save(fig, out_dir, "Diff_dt_sweep.png")

    print("\n  Temporal convergence (Implicit, N=128, t*=0.1):")
    print(f"  {'dt':>10}  {'L2 error':>12}  {'Rate':>6}")
    for i, (dt, e) in enumerate(zip(impl_dts, err_impl)):
        rate = np.log(err_impl[i]/err_impl[i-1])/np.log(impl_dts[i]/impl_dts[i-1]) if i > 0 else float("nan")
        print(f"  {dt:>10.1e}  {e:>12.4e}  {rate:>6.2f}")
    print(f"\n  CN spatial error floor: {e_floor:.4e} (temporal error negligible for r<0.5)")


def mode_scheme_residue(impl_dir, cn_dir, expl_dir, out_dir):
    """Iterations and residue per timestep for 3 temporal schemes."""
    dirs   = [impl_dir, cn_dir, expl_dir]
    labels = ["Implicit (1st order)", "Crank-Nicolson (2nd order)", "Explicit"]
    styles = ["-", "--", ":"]
    colors = [COLORS[0], COLORS[1], COLORS[2]]

    fig, axes = plt.subplots(2, 1, figsize=(7, 5.5), sharex=False)
    ax_iter, ax_res = axes

    for d, lab, sty, col in zip(dirs, labels, styles, colors):
        t, iters, res = load_debug(d)
        mask = iters > 0
        ax_iter.plot(t[mask], iters[mask], sty, color=col, label=lab)
        ax_res.semilogy(t[mask], res[mask],  sty, color=col, label=lab)

    ax_iter.set_ylabel("Iterations per step")
    ax_iter.set_xlabel("Time [s]")
    ax_iter.legend()
    ax_iter.grid(True, ls=":", alpha=0.5)

    ax_res.set_ylabel("Residual")
    ax_res.set_xlabel("Time [s]")
    ax_res.legend()
    ax_res.grid(True, ls=":", which="both", alpha=0.5)

    fig.suptitle("1D Rod — Temporal scheme comparison (N=64)", fontsize=12)
    fig.tight_layout()
    save(fig, out_dir, "Diff_scheme_residue.png")


def mode_param_study(impl_dirs, impl_dts, cn_dirs, cn_dts, expl_dir, N_dirs, N_vals, out_dir):
    """Combined 3-panel numerical parameter study figure."""
    alpha, T_L, T_R, T0 = 1.0, 100.0, 0.0, 0.0
    t_star = 0.5

    def err_at(d):
        t_arr, x, T_1d = load_map(d)
        T_num  = T_1d[-1]
        t_sim  = t_arr[-1]
        T_anal = T_rod_transient_fast(x, t_sim, T_L, T_R, T0, alpha=alpha)
        return np.sqrt(np.mean((T_num - T_anal)**2))

    err_N    = [err_at(d)  for d in N_dirs]
    err_impl = [err_at(d)  for d in impl_dirs]
    err_cn   = [err_at(d)  for d in cn_dirs]
    err_expl =  err_at(expl_dir)
    dt_expl  = (1.0/64)**2 / (2.0*alpha)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # ── Panel A: N sweep ──────────────────────────────────────────────────────
    ax = axes[0]
    ax.loglog(N_vals, err_N, "o-", color=COLORS[0], label="Implicit")
    N_ref = np.array([N_vals[0], N_vals[-1]], dtype=float)
    ax.loglog(N_ref, err_N[0]*(N_vals[0]/N_ref)**2, "--", color="gray", alpha=0.7, label="O(Δx²)")
    ax.loglog(N_ref, err_N[0]*(N_vals[0]/N_ref)**1, ":",  color="gray", alpha=0.7, label="O(Δx)")
    ax.set_xlabel("N"); ax.set_ylabel("L₂ error"); ax.set_title("Spatial convergence")
    ax.legend(); ax.grid(True, ls=":", which="both", alpha=0.5)

    # ── Panel B: dt sweep ─────────────────────────────────────────────────────
    ax = axes[1]
    ax.loglog(impl_dts, err_impl, "o-",  color=COLORS[0], label="Implicit")
    ax.loglog(cn_dts,   err_cn,   "s--", color=COLORS[1], label="Crank-Nicolson")
    ax.loglog([dt_expl], [err_expl], "^", color=COLORS[2], label="Explicit", ms=8)
    dt_ref = np.array([min(impl_dts), max(impl_dts)])
    ax.loglog(dt_ref, err_impl[0]*(impl_dts[0]/dt_ref)**1, ":", color=COLORS[0], alpha=0.5, label="O(Δt)")
    ax.loglog(dt_ref, err_cn[0]  *(cn_dts[0]  /dt_ref)**2, ":", color=COLORS[1], alpha=0.5, label="O(Δt²)")
    ax.set_xlabel("Δt [s]"); ax.set_ylabel("L₂ error"); ax.set_title(f"Temporal convergence (t*={t_star} s)")
    ax.legend(fontsize=8); ax.grid(True, ls=":", which="both", alpha=0.5)

    # ── Panel C: iterations + residue ─────────────────────────────────────────
    ax = axes[2]
    # Use medium-dt dirs for scheme comparison
    dirs   = [impl_dirs[2], cn_dirs[2], expl_dir]   # dt≈0.01
    labels = ["Implicit", "Crank-Nicolson", "Explicit"]
    colors = [COLORS[0], COLORS[1], COLORS[2]]
    styles = ["-", "--", ":"]
    ax2    = ax.twinx()
    lines  = []
    for d, lab, col, sty in zip(dirs, labels, colors, styles):
        t, iters, res = load_debug(d)
        mask = iters > 0
        l1, = ax.plot(t[mask], iters[mask], sty, color=col, label=lab+" iter")
        l2, = ax2.semilogy(t[mask], res[mask], sty, color=col, alpha=0.5, label=lab+" res")
        lines.extend([l1, l2])
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Iterations")
    ax2.set_ylabel("Residual")
    ax.set_title("Scheme comparison (N=64, Δt≈0.01 s)")
    ax.grid(True, ls=":", alpha=0.3)
    labs = [l.get_label() for l in lines]
    ax.legend(lines, labs, fontsize=7, ncol=2)

    fig.suptitle("1D Rod — Numerical parameter study", fontsize=13)
    fig.tight_layout()
    save(fig, out_dir, "Diff_param_study.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 — BC combination validation
# ══════════════════════════════════════════════════════════════════════════════
# Analytical steady states with Q=200 W/m³ source term (γ=1 for all BC cases).
# Each BVP: γ T'' + Q = 0  →  T(x) = C1 + C2·x − Q/(2γ)·x²
#
# ND: q_in=100 at x=0 (→ T'(0)=−100), T(1)=0
#     C2=−100, C1=200  →  T_ss = 200 − 100x − 100x²
#
# CD: Robin h=10, T∞=100 at x=0 (→ −T'(0)=10(100−T(0))), T(1)=0
#     With Q=200: C2=0, C1=100  →  T_ss = 100(1−x²)
#
# DN: T(0)=100, q_exit=50 at x=1 (→ −T'(1)=50 → T'(1)=−50)
#     C1=100, C2=150  →  T_ss = 100 + 150x − 100x²  (dome, peak at x=0.75)
#
# DC: T(0)=100, Robin h=10, T∞=0 at x=1 (→ −T'(1)=10·T(1))
#     C1=100, C2=200/11  →  T_ss = 100 + (200/11)x − 100x²
ANALYTICAL_SS = {
    "ND": lambda x: 200.0 - 100.0*x - 100.0*x**2,
    "CD": lambda x: 100.0 * (1.0 - x**2),
    "DN": lambda x: 100.0 + 150.0*x - 100.0*x**2,
    "DC": lambda x: 100.0 + (200.0/11.0)*x - 100.0*x**2,
}
BC_LABELS = {
    "ND": "Neumann–Dirichlet\n(q_in=100, T_R=0)",
    "CD": "Robin–Dirichlet\n(h=10, T∞=100; T_R=0)",
    "DN": "Dirichlet–Neumann\n(T_L=100, q_out=50)",
    "DC": "Dirichlet–Robin\n(T_L=100; h=10, T∞=0)",
}

def mode_bc_val(dirs, out_dir):
    """4 BC combinations — numerical vs analytical steady state."""
    keys   = ["ND", "CD", "DN", "DC"]
    titles = [BC_LABELS[k] for k in keys]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    axes = axes.ravel()

    for i, (d, key) in enumerate(zip(dirs, keys)):
        t_arr, x, T_1d = load_map(d)
        T_num  = T_1d[-1]
        T_anal = ANALYTICAL_SS[key](x)
        l2 = np.sqrt(np.mean((T_num - T_anal)**2))

        ax = axes[i]
        ax.plot(x, T_num,  "o", color=COLORS[i], ms=4, label=f"Numerical (L₂={l2:.2e})")
        ax.plot(x, T_anal, "-", color="k", lw=1.2, label="Analytical (steady state)")
        ax.set_xlabel("x [m]"); ax.set_ylabel("T [°C]")
        ax.set_title(titles[i])
        ax.legend(fontsize=8)
        ax.grid(True, ls=":", alpha=0.5)

    fig.suptitle("1D Rod — Boundary condition validation", fontsize=13)
    fig.tight_layout()
    save(fig, out_dir, "Diff_BC_validation.png")

    print("\n  BC validation summary:")
    for d, key in zip(dirs, keys):
        t_arr, x, T_1d = load_map(d)
        T_num  = T_1d[-1]
        T_anal = ANALYTICAL_SS[key](x)
        l2 = np.sqrt(np.mean((T_num - T_anal)**2))
        print(f"  {key}: L2={l2:.4e}")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — 4Materials numerical parameter study
# ══════════════════════════════════════════════════════════════════════════════
def mode_mat4_param_study(N_dirs, N_vals, dt_impl_dirs, dt_impl_vals, dt_cn_dirs, dt_cn_vals, expl_dir, out_dir):
    """Iterations/residue for N sweep and dt/scheme sweep (4Materials)."""

    def get_debug_stats(d):
        t, iters, res = load_debug(d)
        mask = iters > 0
        return t[mask], iters[mask], res[mask]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # ── Panel A: N sweep (iterations and residue vs N) ────────────────────────
    ax = axes[0]
    avg_iters = []
    avg_res   = []
    for d in N_dirs:
        t, iters, res = get_debug_stats(d)
        avg_iters.append(np.mean(iters))
        avg_res.append(np.mean(res))

    ax2 = ax.twinx()
    ax.plot(N_vals, avg_iters, "o-", color=COLORS[0], label="Avg iterations")
    ax2.semilogy(N_vals, avg_res, "s--", color=COLORS[1], label="Avg residual")
    ax.set_xlabel("N (mesh nodes)"); ax.set_ylabel("Average iterations/step", color=COLORS[0])
    ax2.set_ylabel("Average residual", color=COLORS[1])
    ax.set_title("4-Material rod — N sweep (implicit, Δt=50 s)")
    ax.tick_params(axis="y", labelcolor=COLORS[0])
    ax2.tick_params(axis="y", labelcolor=COLORS[1])
    ax.grid(True, ls=":", alpha=0.5)
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labs1+labs2, fontsize=8)

    # ── Panel B: dt/scheme sweep (iterations vs time for 3 representative cases)
    ax = axes[1]
    scheme_dirs  = [dt_impl_dirs[1], dt_cn_dirs[1], expl_dir]   # dt=50 for impl/CN, auto for expl
    scheme_labs  = ["Implicit (Δt=50 s)", "Crank-Nicolson (Δt=50 s)", "Explicit (auto-Δt)"]
    scheme_cols  = [COLORS[0], COLORS[1], COLORS[2]]
    scheme_sty   = ["-", "--", ":"]
    ax2 = ax.twinx()
    lines = []
    for d, lab, col, sty in zip(scheme_dirs, scheme_labs, scheme_cols, scheme_sty):
        t, iters, res = get_debug_stats(d)
        # Downsample for readability
        step = max(1, len(t)//200)
        l1, = ax.plot(t[::step], iters[::step], sty, color=col, label=lab+" (iter)")
        l2, = ax2.semilogy(t[::step], res[::step], sty, color=col, alpha=0.5, label=lab+" (res)")
        lines.extend([l1, l2])
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Iterations/step")
    ax2.set_ylabel("Residual")
    ax.set_title("4-Material rod — Scheme comparison (N=64)")
    ax.grid(True, ls=":", alpha=0.3)
    ax.legend(lines, [l.get_label() for l in lines], fontsize=7, ncol=2)

    fig.suptitle("4-Material rod — Numerical parameter study", fontsize=13)
    fig.tight_layout()
    save(fig, out_dir, "Diff_mat4_param_study.png")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 5 — 4Materials CN vs analytical benchmark
# ══════════════════════════════════════════════════════════════════════════════
def mode_mat4_bench(ref_dir, out_dir):
    maps = load_all_maps(ref_dir)
    # Last map = final steady state
    t_arr, x, T_1d = maps[-1]
    T_num  = T_1d[-1]
    T_anal = T_mat4_ss(x)
    l2 = np.sqrt(np.mean((T_num - T_anal)**2))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(x, T_num,  "o", color=COLORS[0], ms=4, label=f"Crank-Nicolson (L₂={l2:.3e})")
    ax.plot(x, T_anal, "-k", lw=1.5,              label="Analytical (steady state)")
    for xi in mat4_interfaces():
        ax.axvline(xi, color="gray", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlabel("x [m]"); ax.set_ylabel("T [°C]")
    ax.set_title("4-Material rod — Steady-state CN vs analytical benchmark")
    ax.legend()
    ax.grid(True, ls=":", alpha=0.5)

    # Add material labels
    mid_x = [0.125, 0.375, 0.625, 0.875]
    mat_labs = ["Mat 1\nγ=170", "Mat 2\nγ=140", "Mat 3\nγ=200", "Mat 4\nγ=140"]
    ymin, ymax = ax.get_ylim()
    for mx, ml in zip(mid_x, mat_labs):
        ax.text(mx, ymin + 0.05*(ymax-ymin), ml, ha="center", va="bottom",
                fontsize=8, color="gray")

    # Annotate interface temperatures (analytical)
    T_interfaces = [T_mat4_ss(np.array([xi]))[0] for xi in [0.25, 0.5, 0.75]]
    for xi, Ti in zip([0.25, 0.5, 0.75], T_interfaces):
        ax.annotate(f"{Ti:.1f}°C", xy=(xi, Ti), xytext=(xi+0.04, Ti+3),
                    fontsize=8, color="k",
                    arrowprops=dict(arrowstyle="->", color="k", lw=0.8))

    fig.tight_layout()
    save(fig, out_dir, "Diff_mat4_benchmark.png")

    print(f"\n  4Mat benchmark L2 error = {l2:.4e}")
    print("  Interface temperatures (analytical):")
    for xi, Ti in zip([0.25, 0.5, 0.75], T_interfaces):
        Ti_num = np.interp(xi, x, T_num)
        print(f"    x={xi:.2f}: T_anal={Ti:.2f}, T_num={Ti_num:.2f}, err={abs(Ti_num-Ti):.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Plot 6 — 4Materials temperature colormap at multiple times
# ══════════════════════════════════════════════════════════════════════════════
def mode_mat4_colormap(ref_dir, out_dir):
    maps = load_all_maps(ref_dir)
    n_snap = min(len(maps), 6)
    maps   = maps[:n_snap]

    fig, axes = plt.subplots(1, n_snap, figsize=(3.2*n_snap, 3.0), squeeze=False)
    axes = axes[0]

    vmin = min(T_1d.min() for _, _, T_1d in maps)
    vmax = max(T_1d.max() for _, _, T_1d in maps)

    for i, (t_arr, x, T_1d) in enumerate(maps):
        ax = axes[i]
        T_snap = T_1d[-1]           # last row in this probe window
        t_snap = t_arr[-1]

        # Pseudo-2D colormap: duplicate in y for visual effect
        T_2d = np.vstack([T_snap, T_snap])
        im = ax.imshow(T_2d, aspect="auto",
                       extent=[x[0], x[-1], 0, 1],
                       cmap="hot", vmin=vmin, vmax=vmax,
                       origin="lower")
        # Overlay the 1D profile
        ax.plot(x, 0.5 + 0.4 * (T_snap - vmin)/(vmax - vmin + 1e-12), "w-", lw=1)
        for xi in mat4_interfaces():
            ax.axvline(xi, color="cyan", ls="--", lw=0.8, alpha=0.7)
        ax.set_xlabel("x [m]")
        ax.set_title(f"t = {t_snap:.0f} s")
        if i == 0:
            ax.set_ylabel("T profile")
        ax.set_yticks([])

    fig.colorbar(im, ax=axes[-1], label="T [°C]", shrink=0.85)
    fig.suptitle("4-Material rod — Temperature evolution (CN, N=64)", fontsize=12)
    fig.tight_layout()
    save(fig, out_dir, "Diff_mat4_colormap.png")


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--solver-cmp",    nargs=2,  metavar="DIR", help="CG dir, GS dir")
    p.add_argument("--N-sweep",       nargs="+",metavar="DIR")
    p.add_argument("--N-vals",        nargs="+",type=int, default=[8,16,32,64,128])
    p.add_argument("--dt-sweep-impl", nargs="+",metavar="DIR")
    p.add_argument("--dt-vals-impl",  nargs="+",type=float, default=[0.001,0.005,0.01,0.05,0.1])
    p.add_argument("--dt-sweep-cn",   nargs="+",metavar="DIR")
    p.add_argument("--dt-vals-cn",    nargs="+",type=float, default=[0.001,0.005,0.01,0.05,0.1])
    p.add_argument("--dt-sweep-expl", metavar="DIR")
    p.add_argument("--scheme-dirs",   nargs=3,  metavar="DIR", help="impl, CN, expl dirs")
    p.add_argument("--param-study",   action="store_true", help="Generate combined 3-panel param study")
    p.add_argument("--bc-val",        nargs=4,  metavar="DIR", help="ND CD DN DC dirs")
    p.add_argument("--mat4-refine",   nargs="+",metavar="DIR")
    p.add_argument("--mat4-N-vals",   nargs="+",type=int, default=[16,32,64,128])
    p.add_argument("--mat4-dt-impl",  nargs="+",metavar="DIR")
    p.add_argument("--mat4-dt-impl-vals", nargs="+", type=float, default=[10,50,500])
    p.add_argument("--mat4-dt-cn",    nargs="+",metavar="DIR")
    p.add_argument("--mat4-dt-cn-vals",   nargs="+", type=float, default=[10,50,500])
    p.add_argument("--mat4-expl",     metavar="DIR")
    p.add_argument("--mat4-ref",      metavar="DIR")
    p.add_argument("--mat4-colormap", metavar="DIR")
    p.add_argument("--out-dir",       default="ioRes/DataNew")
    args = p.parse_args()

    out = args.out_dir

    if args.solver_cmp:
        mode_solver_cmp(args.solver_cmp, out)

    if args.N_sweep:
        mode_N_sweep(args.N_sweep, args.N_vals, out)

    if args.dt_sweep_impl and args.dt_sweep_cn:
        mode_dt_sweep(args.dt_sweep_impl, args.dt_vals_impl,
                      args.dt_sweep_cn,   args.dt_vals_cn,
                      None, out)

    if args.scheme_dirs:
        mode_scheme_residue(args.scheme_dirs[0], args.scheme_dirs[1],
                            args.scheme_dirs[2], out)

    if args.param_study and args.N_sweep and args.dt_sweep_impl and args.dt_sweep_cn and args.dt_sweep_expl:
        mode_param_study(
            args.dt_sweep_impl, args.dt_vals_impl,
            args.dt_sweep_cn,   args.dt_vals_cn,
            args.dt_sweep_expl,
            args.N_sweep, args.N_vals, out)

    if args.bc_val:
        mode_bc_val(args.bc_val, out)

    if args.mat4_refine and args.mat4_dt_impl and args.mat4_dt_cn and args.mat4_expl:
        mode_mat4_param_study(
            args.mat4_refine,   args.mat4_N_vals,
            args.mat4_dt_impl,  args.mat4_dt_impl_vals,
            args.mat4_dt_cn,    args.mat4_dt_cn_vals,
            args.mat4_expl, out)

    if args.mat4_ref:
        mode_mat4_bench(args.mat4_ref, out)

    if args.mat4_colormap:
        mode_mat4_colormap(args.mat4_colormap, out)


if __name__ == "__main__":
    main()
