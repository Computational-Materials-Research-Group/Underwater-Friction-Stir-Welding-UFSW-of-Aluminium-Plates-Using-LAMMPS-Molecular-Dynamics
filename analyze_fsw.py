#!/usr/bin/env python3
"""
Post-process the LAMMPS underwater-FSW runs (fsw_underwater.in).

Run in the folder with the log + dump files:
    python3 analyze_fsw.py                       # log file auto-detected (fsw.log / log.lammps)
    python3 analyze_fsw.py --cu                  # if plate B was Cu
Compare runs (e.g. wet vs dry, each in its own folder):
    python3 analyze_fsw.py wet=/mnt/d/extras/wet dry=/mnt/d/extras/dry

Writes (in the current folder):
    fsw_timeseries.png   temperature, forces, torque, mixing vs time
    fsw_Tmaps.png        surface temperature maps (from dumps, flow removed)
    fsw_peakT.png        peak surface temperature + crystal structure in the stir zone vs time
    fsw_flow.png         stirring flow (averaged velocity arrows) of water and metal around the tool
    fsw_summary.csv      the numbers behind the plots
Needs: numpy, pandas, matplotlib
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

KB = 8.617333e-5          # eV/K
MVV2E = 1.0364269e-4      # amu*(A/ps)^2 -> eV
EVA2NN = 1.602177         # eV/A -> nN
EV2NNNM = 0.1602177       # eV   -> nN*nm

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="*", default=["run=."], help="label=folder (default: current folder)")
ap.add_argument("--log", default=None, help="log file name inside each folder")
ap.add_argument("--cu", action="store_true", help="type 4 is Cu")
ap.add_argument("--bin", type=float, default=8.0, help="temperature-map bin size (A)")
ap.add_argument("--depth", type=float, default=15.0, help="surface layer depth for maps (A)")
ap.add_argument("--fbin", type=float, default=4.0, help="flow-field bin size (A)")
ap.add_argument("--fmax", type=int, default=25, help="max dumps averaged per flow field")
ap.add_argument("--every", type=int, default=1, help="use every N-th dump for peak-T curve")
a = ap.parse_args()
MASS = {1: 26.982, 4: 63.546 if a.cu else 26.982}
runs = [r.split("=", 1) if "=" in r else (os.path.basename(os.path.abspath(r)) or r, r) for r in a.runs]

# ---------------------------------------------------------------- log
def read_weld_thermo(folder):
    cands = [a.log] if a.log else ["fsw.log", "log.lammps"] + sorted(glob.glob(os.path.join(folder, "*.log")))
    for c in cands:
        f = c if os.path.isabs(c) else os.path.join(folder, c)
        if f and os.path.isfile(f):
            break
    else:
        sys.exit(f"no log file found in {folder}")
    lines = open(f, errors="ignore").read().splitlines()
    blocks, i = [], 0
    while i < len(lines):
        if lines[i].split()[:1] == ["Step"]:
            hdr = lines[i].split(); rows = []; i += 1
            while i < len(lines) and not lines[i].startswith("Loop time"):
                p = lines[i].split()
                if len(p) == len(hdr):
                    try: rows.append([float(v) for v in p])
                    except ValueError: pass
                i += 1
            blocks.append(pd.DataFrame(rows, columns=hdr))
        i += 1
    weld = [b for b in blocks if "v_t" in b.columns]
    if not weld:
        sys.exit(f"{f}: no welding thermo block (column v_t) found")
    print(f"[{folder}] log: {f}  ({len(weld[-1])} thermo rows)")
    return weld[-1]

# ---------------------------------------------------------------- dumps
def read_dump(fn):
    with open(fn) as fh:
        head = [next(fh) for _ in range(9)]
    step = int(head[1]); cols = head[8].split()[2:]
    box = [list(map(float, head[i].split()[:2])) for i in (5, 6, 7)]
    df = pd.read_csv(fn, sep=r"\s+", skiprows=9, header=None, names=cols, engine="c")
    return step, box, df

def bin_temperature(df, x0, x1, y0, y1, zmin, b):
    m = df[df.type.isin([1, 4]) & (df.z > zmin)]
    ix = ((m.x.values - x0) // b).astype(int); iy = ((m.y.values - y0) // b).astype(int)
    nx, ny = int(np.ceil((x1 - x0) / b)), int(np.ceil((y1 - y0) / b))
    ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    key = ix[ok] * ny + iy[ok]
    mass = m.type.map(MASS).values[ok]
    v = m[["vx", "vy", "vz"]].values[ok]
    n = np.bincount(key, minlength=nx * ny)
    M = np.bincount(key, mass, minlength=nx * ny)
    P = np.stack([np.bincount(key, mass * v[:, k], minlength=nx * ny) for k in range(3)], 1)
    vcom = P / np.maximum(M, 1e-12)[:, None]
    dv = v - vcom[key]
    ke = 0.5 * np.bincount(key, mass * (dv ** 2).sum(1), minlength=nx * ny) * MVV2E
    T = np.where(n >= 6, 2 * ke / (3 * np.maximum(n - 1, 1) * KB), np.nan)   # flow (COM) removed
    return T.reshape(nx, ny)

def tool_geom(d):
    tl = d[d.type == 2]
    xc, yc = tl.x.mean(), tl.y.mean()
    r = np.hypot(tl.x - xc, tl.y - yc).values
    z = tl.z.values; zmin = z.min()
    Rp = r[z < zmin + 4].max(); Rs = r.max()
    zsb = z[r > Rp + 1.5].min() if (r > Rp + 1.5).any() else z.max()
    return dict(xc=xc, yc=yc, Rp=Rp, Rs=Rs, ztip=zmin, zsb=zsb, ztop=z.max())

def flow_field(fnames, kind, box, b, half):
    """velocity averaged in x-y bins around the (moving) tool; random thermal motion averages out"""
    Lx, Ly = box[0][1] - box[0][0], box[1][1] - box[1][0]
    n = int(np.ceil(2 * half / b)); N = np.zeros(n * n); U = np.zeros(n * n); V = np.zeros(n * n)
    for fn in fnames:
        _, _, d = read_dump(fn); g = tool_geom(d)
        if kind == "water":   # water layer beside the shoulder
            m = d[(d.type == 3) & (d.z > g["zsb"] - 2) & (d.z < g["ztop"] + 4)]
        else:                 # metal in the pin zone, below the shoulder
            m = d[d.type.isin([1, 4]) & (d.z > g["ztip"] + 1) & (d.z < g["zsb"])]
        dx = m.x.values - g["xc"]; dx -= Lx * np.round(dx / Lx)
        dy = m.y.values - g["yc"]; dy -= Ly * np.round(dy / Ly)
        ix = ((dx + half) // b).astype(int); iy = ((dy + half) // b).astype(int)
        ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        k = ix[ok] * n + iy[ok]
        N += np.bincount(k, minlength=n * n)
        U += np.bincount(k, m.vx.values[ok], minlength=n * n)
        V += np.bincount(k, m.vy.values[ok], minlength=n * n)
    good = N >= 8 * max(1, len(fnames) // 4)
    U = np.where(good, U / np.maximum(N, 1), np.nan).reshape(n, n)
    V = np.where(good, V / np.maximum(N, 1), np.nan).reshape(n, n)
    return U, V, g

# ---------------------------------------------------------------- per run
res = {}
for label, folder in runs:
    th = read_weld_thermo(folder)
    t = th["v_t"].values
    zsh = th["v_zsh"].values if "v_zsh" in th else None
    xsh = th["v_xsh"].values if "v_xsh" in th else None
    tp = t[np.argmax(zsh <= zsh.min() + 1e-9)] if zsh is not None else None
    t1 = t[np.argmax(xsh > 1e-9)] if xsh is not None and (xsh > 0).any() else None

    files = sorted(glob.glob(os.path.join(folder, "dump.fsw.*")), key=lambda s: int(re.findall(r"\d+", s)[-1]))
    if not files:
        sys.exit(f"no dump.fsw.* files in {folder}")
    _, box, d0 = read_dump(files[0])
    metal0 = d0[d0.type.isin([1, 4])]
    zAl = np.percentile(metal0.z, 99.5)
    tool0 = d0[d0.type == 2]
    xt0, yt = tool0.x.mean(), tool0.y.mean()
    Rs = np.hypot(tool0.x - xt0, tool0.y - yt).max()
    (x0, x1), (y0, y1) = box[0], box[1]
    dt_step = t[-1] / th["Step"].values[-1] if th["Step"].values[-1] > 0 else 0.001
    print(f"[{label}] {len(files)} dumps, Al top z={zAl:.1f} A, tool at ({xt0:.0f},{yt:.0f}), shoulder R={Rs:.0f} A")

    # time series from dumps: peak surface T, crystal structure in stir zone
    rows = []
    for fn in files[::a.every]:
        step, _, d = read_dump(fn)
        tt = step * dt_step
        T = bin_temperature(d, x0, x1, y0, y1, zAl - a.depth, a.bin)
        row = dict(t=tt, Tpeak=np.nanmax(T), T95=np.nanpercentile(T, 95))
        if "c_cna" in d:
            sz = d[d.type.isin([1, 4]) & (np.abs(d.y - yt) < Rs) & (d.z > zAl - a.depth)]
            cna = sz.c_cna.values
            for k, nm in ((1, "fcc"), (2, "hcp"), (5, "other")):
                row[nm] = (cna == k).mean() if len(cna) else np.nan
        rows.append(row)
        print(f"\r[{label}] dumps processed: {len(rows)}", end="")
    print()
    dumpts = pd.DataFrame(rows)

    # maps at 3 moments: end of plunge, mid-traverse, final
    def nearest(target):
        steps = np.array([int(re.findall(r"\d+", s)[-1]) for s in files]) * dt_step
        return files[int(np.argmin(np.abs(steps - target)))]
    moments = [("end of plunge", tp if tp else t[-1] / 3),
               ("mid traverse", (t1 + t[-1]) / 2 if t1 else 2 * t[-1] / 3),
               ("final", t[-1])]
    maps = []
    for name, tm in moments:
        step, _, d = read_dump(nearest(tm))
        maps.append((f"{name} (t = {step*dt_step:.0f} ps)", bin_temperature(d, x0, x1, y0, y1, zAl - a.depth, a.bin)))
    # flow fields: average many dumps in the tool's frame (dwell = tool spinning in place)
    ftimes = np.array([int(re.findall(r"\d+", s)[-1]) for s in files]) * dt_step
    windows = []
    if tp is not None and t1 is not None:
        dw = [f for f, tt in zip(files, ftimes) if tp + 1 < tt <= t1]
        if dw: windows.append(("dwell", dw[-a.fmax:]))
        tr = [f for f, tt in zip(files, ftimes) if tt > t1 + 2]
        if tr: windows.append(("traverse", tr[-a.fmax:]))
    if not windows:
        windows.append(("last frames", files[-min(a.fmax, len(files)):]))
    half = Rs + 22
    flows = []
    for wname, fl in windows:
        for kind in ("water", "metal"):
            U, V, g = flow_field(fl, kind, box, a.fbin, half)
            if np.isfinite(U).any():
                flows.append((f"{kind} – {wname} ({len(fl)} dumps)", U, V, g, half))
        print(f"[{label}] flow field: {wname}, {len(fl)} dumps averaged")
    res[label] = dict(th=th, tp=tp, t1=t1, dumpts=dumpts, maps=maps, ext=(x0, x1, y0, y1), yt=yt, flows=flows)

# ---------------------------------------------------------------- plots
def phases(ax, r):
    for v in (r["tp"], r["t1"]):
        if v: ax.axvline(v, color="0.6", ls=":", lw=1)

panels = [("c_tM", "Metal temperature (K)", 1), ("c_Ft[3]", "Plunge force Fz (nN)", EVA2NN),
          ("c_Ft[1]", "Traverse force Fx (nN)", EVA2NN), ("c_Tq", "Torque (nN·nm)", EV2NNNM),
          ("c_mix", "Plate-B atoms across seam", 1)]
panels = [p for p in panels if any(p[0] in r["th"] for r in res.values())]
fig, axs = plt.subplots(len(panels), 1, figsize=(9, 2.3 * len(panels)), sharex=True)
for ax, (col, lab, f) in zip(np.atleast_1d(axs), panels):
    for label, r in res.items():
        if col in r["th"]:
            y = r["th"][col].values * f
            ax.plot(r["th"]["v_t"], y, lw=0.6, alpha=0.35)
            w = max(1, len(y) // 60)
            ax.plot(r["th"]["v_t"], pd.Series(y).rolling(w, center=True, min_periods=1).mean(),
                    lw=1.8, color=ax.lines[-1].get_color(), label=label)
        phases(ax, r)
    ax.set_ylabel(lab, fontsize=9); ax.grid(alpha=0.3)
np.atleast_1d(axs)[0].legend(fontsize=8)
np.atleast_1d(axs)[0].set_title("plunge | dwell | traverse  (dotted lines)", fontsize=9)
np.atleast_1d(axs)[-1].set_xlabel("time (ps)")
fig.tight_layout(); fig.savefig("fsw_timeseries.png", dpi=150); plt.close(fig)

allT = np.concatenate([m[1][np.isfinite(m[1])] for r in res.values() for m in r["maps"]])
vmax = np.percentile(allT, 99)
fig, axs = plt.subplots(len(res), 3, figsize=(15, 3.2 * len(res)), squeeze=False)
for i, (label, r) in enumerate(res.items()):
    x0, x1, y0, y1 = r["ext"]
    for j, (title, T) in enumerate(r["maps"]):
        im = axs[i, j].imshow(T.T, origin="lower", extent=(x0, x1, y0, y1), cmap="inferno", vmin=250, vmax=vmax)
        axs[i, j].axhline(r["yt"], color="w", ls="--", lw=0.6)
        axs[i, j].set_title(f"{label}: {title}", fontsize=9); axs[i, j].set_xlabel("x (Å)"); axs[i, j].set_ylabel("y (Å)")
fig.colorbar(im, ax=axs, label="T (K), flow removed", shrink=0.8)
fig.savefig("fsw_Tmaps.png", dpi=150); plt.close(fig)

has_cna = any("fcc" in r["dumpts"] for r in res.values())
fig, axs = plt.subplots(2 if has_cna else 1, 1, figsize=(9, 6 if has_cna else 3.2), sharex=True, squeeze=False)
for label, r in res.items():
    d = r["dumpts"]
    axs[0, 0].plot(d.t, d.T95, lw=1.8, label=f"{label} (95th pct)")
    axs[0, 0].plot(d.t, d.Tpeak, lw=0.8, alpha=0.5, color=axs[0, 0].lines[-1].get_color(), label=f"{label} (max bin)")
    phases(axs[0, 0], r)
    if has_cna and "fcc" in d:
        for nm, ls in (("fcc", "-"), ("hcp", "--"), ("other", ":")):
            axs[1, 0].plot(d.t, 100 * d[nm], ls=ls, label=f"{label} {nm}")
        phases(axs[1, 0], r)
axs[0, 0].set_ylabel("Surface temperature (K)"); axs[0, 0].legend(fontsize=8); axs[0, 0].grid(alpha=0.3)
if has_cna:
    axs[1, 0].set_ylabel("Stir-zone atoms (%)"); axs[1, 0].legend(fontsize=8, ncol=3); axs[1, 0].grid(alpha=0.3)
    axs[1, 0].set_title("hcp = stacking faults / dislocations, other = disorder/hot/surface", fontsize=8)
axs[-1, 0].set_xlabel("time (ps)")
fig.tight_layout(); fig.savefig("fsw_peakT.png", dpi=150); plt.close(fig)

from matplotlib.patches import Circle
ncol = max(len(r["flows"]) for r in res.values())
if ncol:
    fig, axs = plt.subplots(len(res), ncol, figsize=(4.6 * ncol, 4.4 * len(res)), squeeze=False)
    for i, (label, r) in enumerate(res.items()):
        for j in range(ncol):
            ax = axs[i, j]
            if j >= len(r["flows"]): ax.axis("off"); continue
            title, U, V, g, half = r["flows"][j]
            sp = np.hypot(U, V); n = U.shape[0]
            c = (np.arange(n) + 0.5) * a.fbin - half
            vmax = np.nanpercentile(sp, 98)
            im = ax.imshow(sp.T, origin="lower", extent=(-half, half, -half, half), cmap="viridis", vmin=0, vmax=vmax)
            X, Y = np.meshgrid(c, c, indexing="ij")
            ax.quiver(X, Y, U, V, color="w", scale=vmax * n * 0.9, width=0.004)
            ax.add_patch(Circle((0, 0), g["Rp"], fc="0.15", ec="k"))
            ax.add_patch(Circle((0, 0), g["Rs"], fc="none", ec="k" if "water" in title else "w", ls="--"))
            ax.set_title(f"{label}: {title}", fontsize=9)
            ax.set_xlabel("x – tool x (Å)"); ax.set_ylabel("y – tool y (Å)")
            ax.annotate("", xy=(half * 0.95, -half * 0.88), xytext=(half * 0.55, -half * 0.88),
                        arrowprops=dict(arrowstyle="->", color="r", lw=2))
            ax.text(half * 0.55, -half * 0.8, "weld dir.", color="r", fontsize=7)
            fig.colorbar(im, ax=ax, shrink=0.75, label="mean flow speed (Å/ps)")
    fig.suptitle("Stirring flow: velocities averaged in %.0f Å bins over many dumps (thermal noise removed)\n"
                 "dark disc = pin, dashed = shoulder edge; water = layer beside shoulder, metal = pin zone"
                 % a.fbin, fontsize=10)
    fig.tight_layout(); fig.savefig("fsw_flow.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- summary
out = []
for label, r in res.items():
    th, d = r["th"], r["dumpts"]
    trav = th[th.v_t > r["t1"]] if r["t1"] else th.iloc[len(th) // 2:]
    s = dict(run=label, peak_surface_T_K=d.T95.max(), final_surface_T_K=d.T95.iloc[-1])
    for col, nm, f in (("c_Ft[3]", "mean_Fz_traverse_nN", EVA2NN), ("c_Ft[1]", "mean_Fx_traverse_nN", EVA2NN),
                       ("c_Tq", "mean_torque_traverse_nNnm", EV2NNNM), ("c_mix", "final_mixed_atoms", 1)):
        if col in th:
            s[nm] = trav[col].mean() * f if nm.startswith("mean") else th[col].iloc[-1]
    if "fcc" in d: s["final_stirzone_fcc_pct"] = 100 * d.fcc.iloc[-1]
    out.append(s)
    d.assign(run=label).to_csv(f"fsw_dumpseries_{label}.csv", index=False)
pd.DataFrame(out).to_csv("fsw_summary.csv", index=False)
print(pd.DataFrame(out).round(2).to_string(index=False))
print("wrote fsw_timeseries.png, fsw_Tmaps.png, fsw_peakT.png, fsw_flow.png, fsw_summary.csv")
