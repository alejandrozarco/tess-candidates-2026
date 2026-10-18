#!/usr/bin/env python3
"""Regenerate the figures in figures/ from the tables in data/ (light and dark SVG variants).
usage: python3 figures/make_figures.py        (needs numpy and matplotlib only)
Figures: fold_<tic>  (TTV-corrected fold of the detrended epochs with 20-min bins and the transit model),
         epochs_<tic> (one panel per transit epoch), oc_<tic> (observed minus calculated mid-transit times),
         diffimage_<tic> (per-sector difference images with the catalogue positions of nearby sources).
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); OUT = os.path.join(ROOT, "figures")
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9", axis="#c3c2b7",
                  pts="#b9b7ae", mark="#2a78d6", mark2="#d4661f", mark3="#2f9e6a", cmap="RdBu_r"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a", axis="#383835",
                 pts="#5a5953", mark="#3987e5", mark2="#e8813a", mark3="#43b27e", cmap="RdBu_r"),
}
SAVE = dict(metadata={"Date": None})
TICS = [402710280, 324429952, 89247735]
NAMES = {402710280: "TIC 402710280", 324429952: "TIC 324429952", 89247735: "TIC 89247735"}
YLIM_FOLD = {89247735: (0.983, 1.006)}       # default: (min(0.3rd percentile, 0.9965), 1.0025)
YLIM_EPOCHS = {89247735: (0.981, 1.008)}     # default: (0.9955, 1.0035)


def read(name):
    with open(os.path.join(DATA, name)) as f:
        rows = [r for r in csv.reader(l for l in f if not l.startswith("#"))]
    head = rows[0]; body = rows[1:]
    out = {}
    for i, h in enumerate(head):
        col = [r[i] for r in body]
        try:
            out[h] = np.array([float(x) if x != "" else np.nan for x in col])
        except ValueError:
            out[h] = np.array(col)
    return out


def style(ax, t):
    ax.set_facecolor(t["surface"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(t["axis"])
    ax.tick_params(colors=t["muted"], labelsize=9, length=0, which="both")
    ax.xaxis.label.set_color(t["ink2"]); ax.yaxis.label.set_color(t["ink2"])
    ax.grid(color=t["grid"], linewidth=0.6); ax.set_axisbelow(True)


def title(ax, t, s):
    ax.set_title(s, color=t["ink"], fontsize=10, loc="left")


def save(fig, t, name):
    fig.savefig(os.path.join(OUT, name), facecolor=t["surface"], **SAVE); plt.close(fig)


def binned(x, y, w):
    b = np.round(x / w); ub = np.unique(b)
    return np.array([x[b == u].mean() for u in ub]), np.array([y[b == u].mean() for u in ub]), np.array([y[b == u].std() / np.sqrt((b == u).sum()) for u in ub])


def fig_fold(t, name, tic):
    F = read(f"folded_{tic}.csv"); M = read(f"model_{tic}.csv"); m = F["in_ephemeris"] == 1
    x = F["dt_days"][m] * 24; y = F["flux"][m]
    fig, ax = plt.subplots(figsize=(8, 4.2)); fig.patch.set_facecolor(t["surface"])
    ax.plot(x, y, ".", color=t["pts"], ms=2, rasterized=False, label="detrended cadences")
    xb, yb, eb = binned(x, y, 20 / 60)
    ax.errorbar(xb, yb, eb, fmt="o", color=t["ink"], ms=3.5, lw=0.8, label="20-min bins")
    ax.plot(M["dt_days"] * 24, M["flux"], "-", color=t["mark2"], lw=1.6, label="transit model")
    ax.set_xlim(-12, 12); lo = np.percentile(y, 0.3); ax.set_ylim(*YLIM_FOLD.get(tic, (min(lo, 0.9965), 1.0025)))
    ax.set_xlabel("hours from fitted mid-transit"); ax.set_ylabel("relative flux")
    n = len(np.unique(F["epoch"][m])); style(ax, t)
    title(ax, t, f"{NAMES[tic]}: {n} transits folded on the per-epoch mid-transit times")
    ax.legend(loc="lower right", fontsize=8, frameon=False, labelcolor=t["ink2"])
    fig.tight_layout(); save(fig, t, name)


def fig_epochs(t, name, tic):
    F = read(f"folded_{tic}.csv"); M = read(f"model_{tic}.csv"); T = read(f"transit_times_{tic}.csv")
    eps = [int(e) for e in np.unique(F["epoch"])]; eps.sort(key=lambda e: F["time_btjd"][F["epoch"] == e].min())
    n = len(eps); nc = 4 if n > 4 else n; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.3 * nc, 2.9 * nr), squeeze=False, sharey=True); fig.patch.set_facecolor(t["surface"])
    for k, ep in enumerate(eps):
        ax = axs[k // nc, k % nc]; m = F["epoch"] == ep; x = F["dt_days"][m] * 24; y = F["flux"][m]
        ax.plot(x, y, ".", color=t["pts"], ms=2); xb, yb, eb = binned(x, y, 30 / 60); ax.plot(xb, yb, "o", color=t["ink"], ms=3)
        ax.plot(M["dt_days"] * 24, M["flux"], "-", color=t["mark2"], lw=1.3); ax.set_xlim(-12, 12)
        j = np.where(T["epoch"] == ep)[0][0]; note = "" if T["note"][j].startswith("used") else " (scattered light)"
        style(ax, t); title(ax, t, f"S{int(T['sector'][j])}, epoch {ep}{note}\n$t_c$ = BTJD {T['tc_btjd'][j]:.4f}, depth {T['depth_ppm'][j]:.0f} ppm")
        if k // nc == nr - 1: ax.set_xlabel("hours from mid-transit")
        if k % nc == 0: ax.set_ylabel("relative flux")
    for k in range(n, nr * nc):
        axs[k // nc, k % nc].axis("off")
    axs[0, 0].set_ylim(*YLIM_EPOCHS.get(tic, (0.9955, 1.0035)))
    prods = sorted(set(F["product"])); lab = "TESScut r = 2.0 px aperture photometry" + ("; sector 13: SPOC 2-min SAP" if any(p.startswith("SPOC") for p in prods) else "")
    fig.suptitle(f"{NAMES[tic]}: individual transits ({lab})", color=t["ink"], fontsize=10, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


def fig_oc(t, name, tic):
    T = read(f"transit_times_{tic}.csv"); used = np.array([s.startswith("used") for s in T["note"]])
    fig, ax = plt.subplots(figsize=(6.5, 3.6)); fig.patch.set_facecolor(t["surface"])
    ax.axhline(0, color=t["axis"], lw=0.8)
    ax.errorbar(T["epoch"][used], T["oc_min"][used], T["tc_err_min"][used], fmt="o", color=t["mark"], ms=5, lw=1, capsize=2, label="used in the ephemeris")
    if (~used).any():
        ax.errorbar(T["epoch"][~used], T["oc_min"][~used], T["tc_err_min"][~used], fmt="s", color=t["muted"], ms=5, lw=1, capsize=2, label="scattered light, not used")
    for e, o, s in zip(T["epoch"], T["oc_min"], T["sector"]):
        ax.annotate(f"S{int(s)}", (e, o), textcoords="offset points", xytext=(6, 4), fontsize=8, color=t["ink2"])
    ax.set_xlabel("epoch (cycles from $T_0$)"); ax.set_ylabel("O $-$ C (min)")
    chi = np.sum((T["oc_min"][used] / 1440 / T["tc_err_d"][used]) ** 2); dof = used.sum() - 2
    style(ax, t); title(ax, t, f"{NAMES[tic]}: mid-transit times minus the linear ephemeris")
    ax.text(0.99, 0.03, f"$\\chi^2$ = {chi:.1f}, {dof} degrees of freedom", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=t["ink2"])
    ax.legend(loc="best", fontsize=8, frameon=False, labelcolor=t["ink2"]); fig.tight_layout(); save(fig, t, name)


def fig_diff(t, name, tic):
    D = read(f"diffimage_{tic}.csv"); S = read(f"diffimage_sources_{tic}.csv")
    key = "epoch" if "epoch" in D else "sector"       # one image per epoch (two transits in one sector for TIC 89247735)
    keys = sorted(int(s) for s in np.unique(D[key])); n = len(keys); nc = 4 if n > 4 else n; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.2 * nr), squeeze=False); fig.patch.set_facecolor(t["surface"])
    for k, s in enumerate(keys):
        ax = axs[k // nc, k % nc]; m = D[key] == s; img = np.zeros((21, 21)); img[D["y"][m].astype(int), D["x"][m].astype(int)] = D["diff"][m]
        yy, xx = np.indices(img.shape); tx, ty = [S[c][(S[key] == s) & (S["source"] == "target")][0] for c in ("x_px", "y_px")]
        img = img - np.median(img[np.hypot(xx - tx, yy - ty) > 6])      # display only: remove the uniform sky offset of the difference image
        v = np.nanpercentile(np.abs(img), 99.5); ax.imshow(img, origin="lower", cmap=t["cmap"], vmin=-v, vmax=v, interpolation="nearest")
        ms = S[key] == s
        for src, x, y in zip(S["source"][ms], S["x_px"][ms], S["y_px"][ms]):
            if src == "target":
                ax.plot(x, y, "+", color=t["mark3"], ms=14, mew=1.8)
            elif 0 <= x <= 20 and 0 <= y <= 20:
                ax.plot(x, y, "x", color=t["mark2"], ms=7, mew=1.2)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor(t["surface"])
        for sp in ax.spines.values():
            sp.set_color(t["axis"])
        title(ax, t, f"S{s}" if key == "sector" else f"S{int(D['sector'][m][0])}, epoch {s}")
    for k in range(n, nr * nc):
        axs[k // nc, k % nc].axis("off")
    fig.suptitle(f"{NAMES[tic]}: difference images (out-of-transit minus in-transit, uniform offset removed for display; 21 arcsec pixels); + target, x other Gaia DR3 sources of diffimage_sources_{tic}.csv", color=t["ink"], fontsize=9, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


if __name__ == "__main__":
    for th, t in THEMES.items():
        plt.rcParams.update({"text.color": t["ink"], "axes.edgecolor": t["axis"], "axes.titlecolor": t["ink"]})
        for tic in TICS:
            fig_fold(t, f"fold_{tic}_{th}.svg", tic)
            fig_epochs(t, f"epochs_{tic}_{th}.svg", tic)
            fig_oc(t, f"oc_{tic}_{th}.svg", tic)
            fig_diff(t, f"diffimage_{tic}_{th}.svg", tic)
    print("figures written to", OUT)
