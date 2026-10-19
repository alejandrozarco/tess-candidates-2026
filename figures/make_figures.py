#!/usr/bin/env python3
"""Regenerate the figures in figures/ from the tables in data/ (light and dark SVG variants).
usage: python3 figures/make_figures.py [tic ...] [other]   (needs numpy and matplotlib only; default: all)
Figures: fold_<tic>  (TTV-corrected fold of the detrended epochs with 20-min bins and the transit model),
         epochs_<tic> (one panel per transit epoch), oc_<tic> (observed minus calculated mid-transit times),
         diffimage_<tic> (per-sector difference images with the catalogue positions of nearby sources).
'other' (README section 6): eb_history_11613676 (eclipse depth, duration and primary O-C against time, with the
         model of eclipse_model_11613676.csv), eb_profiles_11613676 (per-sector eclipse profiles),
         toi540_sectors (per-sector folded transits of TOI-540 b, sector 98 with and without an injected transit),
         toi540_depths (depth per sector, and sector-98 depths of the control planets).
'c311179599': epochs_311179599c (the three events of the candidate second transit signal with the joint model).
'eb124751008': eb_depths_124751008 (eclipse depth of TIC 124751008 against time, with the inclination models).
"""
import csv, os, sys
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
TICS = [402710280, 324429952, 89247735, 455017950, 311179599]
NAMES = {402710280: "TIC 402710280", 324429952: "TIC 324429952", 89247735: "TIC 89247735", 455017950: "TIC 455017950 (HD 134147)", 311179599: "TIC 311179599"}
FOLD_ALL = {311179599}                       # fold every fitted epoch (the sector-39 transit is in the shape fit, not in the ephemeris)
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
    F = read(f"folded_{tic}.csv"); M = read(f"model_{tic}.csv"); m = (F["in_ephemeris"] == 1) | (tic in FOLD_ALL)
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
        j = np.where(T["epoch"] == ep)[0][0]; note = "" if T["note"][j].startswith("used") else (" (scattered light)" if "scattered" in T["note"][j] else " (not in the ephemeris)")
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
        lab = "scattered light, not used" if any("scattered" in n for n in T["note"][~used]) else "not used in the ephemeris"
        ax.errorbar(T["epoch"][~used], T["oc_min"][~used], T["tc_err_min"][~used], fmt="s", color=t["muted"], ms=5, lw=1, capsize=2, label=lab)
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


SECCOL = {3: "muted", 30: "mark3", 97: "mark", 104: "mark2", 106: "ink"}


def fig_eb_history(t, name):
    D = read("depths_11613676.csv"); E = read("eclipse_times_11613676.csv"); M = read("model_history_11613676.csv")
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.9)); fig.patch.set_facecolor(t["surface"])
    mods = [m for m in ("N4", "MA") if m in set(M["model"])]
    for kind, col, mk, lab in (("pri", t["mark"], "o", "primary"), ("sec", t["mark2"], "s", "secondary")):
        m = (D["kind"] == kind) & np.isnan(D["limit3_pct"]); l = (D["kind"] == kind) & ~np.isnan(D["limit3_pct"])
        axs[0].errorbar(D["year"][m], D["depth_pct"][m], D["depth_err_pct"][m], fmt=mk, color=col, ms=5, label=lab)
        axs[0].errorbar(D["year"][l], D["limit3_pct"][l], yerr=0.25, uplims=True, fmt=mk, mfc="none", color=col, ms=5)
        axs[1].errorbar(D["year"][m], D["T14_h"][m], D["T14_err_h"][m], fmt=mk, color=col, ms=5, label=lab)
        for mo, ls in zip(mods, ("-", ":")):
            mm = M["model"] == mo
            axs[0].plot(M["year"][mm], M[f"depth_{kind}_pct"][mm], ls, color=col, lw=1.1, alpha=0.8)
            axs[1].plot(M["year"][mm], np.where(M[f"T14_{kind}_h"][mm] > 0, M[f"T14_{kind}_h"][mm], np.nan), ls, color=col, lw=1.1, alpha=0.8)
    u = (E["timing_set"] == 1) & (E["kind"] == "pri"); yr = 2000 + (E["t0_btjd"][u] + 2457000 - 2451544.5) / 365.25
    axs[2].axhline(0, color=t["axis"], lw=0.8); axs[2].errorbar(yr, E["oc_min"][u], E["t0_err_min"][u], fmt="o", color=t["mark"], ms=5, capsize=2)
    for ax, yl, ti in ((axs[0], "eclipse depth (%)", "depth (open: 3-sigma limits, sector 3)"), (axs[1], "duration $T_{14}$ (h)", "total duration"), (axs[2], "O $-$ C (min)", "primary mid-eclipse times minus linear ephemeris")):
        ax.set_xlabel("year"); ax.set_ylabel(yl); style(ax, t); title(ax, t, ti); ax.set_xlim(2018, 2029 if ax is not axs[2] else 2027)
    axs[0].legend(loc="upper left", fontsize=8, frameon=False, labelcolor=t["ink2"])
    if mods:
        axs[1].text(0.02, 0.03, "lines: " + ", ".join(f"{m} ({'solid' if i == 0 else 'dotted'})" for i, m in enumerate(mods)) + " models,\nextrapolated outside 2018.7-2026.6", transform=axs[1].transAxes, fontsize=7.5, color=t["ink2"])
    fig.suptitle("TIC 11613676 (CD-30 448): eclipse history, TESS sectors 3, 30, 97, 104, 106", color=t["ink"], fontsize=10, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


def fig_eb_profiles(t, name):
    F = read("profiles_11613676.csv"); D = read("depths_11613676.csv")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True); fig.patch.set_facecolor(t["surface"])
    for ax, kind, lab in ((axs[0], "pri", "primary eclipse"), (axs[1], "sec", "secondary eclipse")):
        for s in sorted(set(int(x) for x in F["sector"])):
            m = (F["sector"] == s) & (F["kind"] == kind)
            if not m.any(): continue
            yr = D["year"][(D["sector"] == s) & (D["kind"] == kind)][0]
            ax.errorbar(F["dt_h"][m], F["flux"][m], F["err"][m], fmt="o-", ms=2.5, lw=0.8, color=t[SECCOL.get(s, "ink")], label=f"S{s} ({yr:.1f})")
        ax.set_xlim(-4, 4); ax.set_xlabel("hours from mid-eclipse"); style(ax, t); title(ax, t, lab)
    axs[0].set_ylabel("relative flux"); axs[1].legend(loc="lower right", fontsize=8, frameon=False, labelcolor=t["ink2"])
    fig.suptitle("TIC 11613676 (CD-30 448): eclipse profiles per sector (10-min bins; SPOC 2-min SAP in S3, S30; TESScut r = 2.5 px otherwise; S3 at the predicted times)", color=t["ink"], fontsize=9.5, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


def fig_toi540_sectors(t, name):
    F = read("folded_200322593.csv"); T = read("transit_depths_200322593.csv")
    sets = [x for x in ("4", "5", "6", "31", "32", "98", "98_injected") if (F["set"] == x).any()]
    fig, axs = plt.subplots(2, 4, figsize=(13, 5.6), sharex=True, sharey=True, squeeze=False); fig.patch.set_facecolor(t["surface"])
    for k, sname in enumerate(sets):
        ax = axs[k // 4, k % 4]; m = F["set"] == sname
        ax.errorbar(F["dt_min"][m], F["flux"][m], F["err"][m], fmt="o", ms=2.5, lw=0.6, color=t["ink"] if not sname.startswith("98") else t["mark2"])
        ax.plot(F["dt_min"][m], F["model"][m], "-", color=t["mark"], lw=1.2)
        s = int(sname.split("_")[0]); j = np.where(T["sector"] == s)[0][0]
        lab = f"S{s} ({T['year'][j]:.1f}), {int(T['n_transits'][j])} transits\ndepth {T['depth_own_ppm'][j]:.0f} $\\pm$ {T['depth_own_err_ppm'][j]:.0f} ppm"
        if sname == "98_injected": lab = "S98 with the reference transit injected\nat the predicted times"
        style(ax, t); title(ax, t, lab); ax.set_xlim(-90, 90); ax.set_ylim(0.9965, 1.0015)
        if k // 4 == 1 or k >= 4: ax.set_xlabel("minutes from predicted mid-transit")
        if k % 4 == 0: ax.set_ylabel("relative flux")
    for k in range(len(sets), 8): axs[k // 4, k % 4].axis("off")
    fig.suptitle("TOI-540 b (TIC 200322593): transits folded on the sectors 4-32 ephemeris (own TPF photometry, r = 2.0 px, 2-min bins; line: reference transit)", color=t["ink"], fontsize=10, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


def fig_toi540_depths(t, name):
    T = read("transit_depths_200322593.csv"); Cc = read("control_planets_s98.csv")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw=dict(width_ratios=[1, 1.25])); fig.patch.set_facecolor(t["surface"]); ax = axs[0]
    for col, key, dx, lab in ((t["mark"], "own", -0.06, "own TPF photometry"), (t["mark2"], "pdc", 0.06, "SPOC PDCSAP")):
        ax.errorbar(T["year"] + dx, T[f"depth_{key}_ppm"], T[f"depth_{key}_err_ppm"], fmt="o", color=col, ms=5, capsize=2, label=lab)
    ax.axhline(0, color=t["axis"], lw=0.8); ax.set_xlabel("year"); ax.set_ylabel("transit depth at the predicted times (ppm)")
    style(ax, t); title(ax, t, "TOI-540 b: depth per TESS sector"); ax.legend(loc="lower left", fontsize=8, frameon=False, labelcolor=t["ink2"])
    for x, y, s in zip(T["year"], T["depth_own_ppm"], T["sector"]): ax.annotate(f"S{int(s)}", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8, color=t["ink2"])
    ax = axs[1]; tois = []; skipped = []
    for toi in [x for x in dict.fromkeys(Cc["toi"]) if abs(x - 540.01) > 1e-6]:
        m = (Cc["toi"] == toi) & ~np.isnan(Cc["snr"]) & (Cc["centroid_frac"] > 0)      # defective SPOC light curves (no centroids) are not shown
        if not m.any(): continue
        pre = m & (Cc["sector"] < 98); s98 = m & (Cc["sector"] == 98)
        if not s98.any():
            skipped.append(toi); continue
        w = 1 / Cc["depth_err_ppm"][pre]**2; i = len(tois); tois.append(toi)
        ref = np.sum(w * Cc["best_depth_ppm"][pre]) / w.sum() if pre.any() else Cc["toi_depth_ppm"][s98][0]   # no earlier 2-min sector: TOI-table depth
        ax.errorbar(i, Cc["depth_ppm"][s98][0] / ref, Cc["depth_err_ppm"][s98][0] / ref, fmt="o", mfc="none", color=t["muted"], ms=5)
        ax.errorbar(i, Cc["best_depth_ppm"][s98][0] / ref, Cc["depth_err_ppm"][s98][0] / ref, fmt="o", color=t["mark3"], ms=5)
    o = np.where(T["sector"] == 98)[0][0]; pre = T["sector"] < 98; w = 1 / T["depth_own_err_ppm"][pre]**2; ref = np.sum(w * T["depth_own_ppm"][pre]) / w.sum()
    ax.errorbar(len(tois), T["depth_own_ppm"][o] / ref, T["depth_own_err_ppm"][o] / ref, fmt="D", color=t["mark2"], ms=6); tois.append("540.01")
    ax.axhline(1, color=t["axis"], lw=0.8); ax.axhline(0, color=t["axis"], lw=0.8, ls=":")
    ax.set_xticks(range(len(tois))); ax.set_xticklabels([f"TOI-{x:.2f}" if not isinstance(x, str) else f"TOI-{x}" for x in tois], rotation=45, ha="right", fontsize=8)
    if skipped: ax.text(0.01, 0.03, "not shown (sector-98 SPOC light curve without centroids): " + ", ".join(f"TOI-{x:.2f}" for x in skipped) + "; TOI-4981.01: relative to the TOI-table depth", transform=ax.transAxes, fontsize=7.5, color=t["ink2"])
    ax.set_ylabel("sector-98 depth / earlier depth"); style(ax, t)
    title(ax, t, "sector-98 depth relative to earlier sectors: control planets (filled: best\nwithin the ephemeris window; open: at the TOI-table ephemeris) and TOI-540 b")
    fig.tight_layout(); save(fig, t, name)


def fig_epochs_c(t, name):
    F = read("folded_311179599c.csv"); T = read("transit_times_311179599c.csv")
    evs = list(T["event"]); fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True); fig.patch.set_facecolor(t["surface"])
    for ax, ev in zip(axs, evs):
        m = F["event"] == ev; x = F["dt_days"][m] * 24; y = F["flux"][m]; o = np.argsort(x); j = list(T["event"]).index(ev)
        ax.plot(x, y, ".", color=t["pts"], ms=2); xb, yb, eb = binned(x, y, 30 / 60); ax.errorbar(xb, yb, eb, fmt="o", color=t["ink"], ms=3, lw=0.7)
        ax.plot(x[o], F["model"][m][o], "-", color=t["mark2"], lw=1.4); ax.set_xlim(-14, 14); ax.set_xlabel("hours from fitted mid-transit")
        style(ax, t); title(ax, t, f"event {ev}: S{int(T['sector'][j])}, {T['cadence_min'][j]:.1f}-min cadence\n$t_c$ = BTJD {T['tc_btjd'][j]:.4f}, depth {T['depth_ppm'][j]:.0f} $\\pm$ {T['depth_err_ppm'][j]:.0f} ppm")
    axs[0].set_ylabel("relative flux"); axs[0].set_ylim(0.9965, 1.0025)
    fig.suptitle("TIC 311179599: the three events of the candidate second transit signal (TESScut r = 2.0 px; 30-min bins; line: joint model)", color=t["ink"], fontsize=10, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


def fig_eb_depths_124(t, name):
    E = read("eclipse_times_124751008.csv"); D = read("depths_124751008.csv"); M = read("model_history_124751008.csv")
    fig, ax = plt.subplots(figsize=(9, 4.4)); fig.patch.set_facecolor(t["surface"])
    ok = E["depth_err_pct"] < 0.2; ax.errorbar(E["year"][ok], E["depth_pct"][ok], E["depth_err_pct"][ok], fmt="o", color=t["mark"], ms=4.5, capsize=2, label="primary eclipses")
    if (~ok).any(): ax.plot(E["year"][~ok], E["depth_pct"][~ok], "o", mfc="none", color=t["mark"], ms=4.5, label="primary (error > 0.2 % or not determined)")
    s = (D["kind"] == "sec") & np.isnan(D["limit3_pct"]); l = (D["kind"] == "sec") & ~np.isnan(D["limit3_pct"])
    ax.errorbar(D["year"][s], D["depth_pct"][s], D["depth_err_pct"][s], fmt="s", color=t["mark2"], ms=5, capsize=2, label="secondary eclipse")
    ax.errorbar(D["year"][l], D["limit3_pct"][l], yerr=0.08, uplims=True, fmt="s", mfc="none", color=t["mark2"], ms=5, label="secondary, 3-sigma limit")
    for (l3, mo), ls, lab in ((("0.00", "sinusoid_free"), "-", "sinusoid, third light 0"), (("0.70", "sinusoid_free"), "--", "sinusoid, third light 0.70"), (("0.00", "quadratic"), ":", "quadratic, third light 0"), (("0.00", "linear"), "-.", "linear (rejected), third light 0")):
        mm = (M["third_light"] == float(l3)) & (M["model"] == mo)
        ax.plot(M["year"][mm], M["pri_depth_pct"][mm], ls, color=t["mark"], lw=1.0, alpha=0.75, label=f"model primary: {lab}")
        if mo != "linear": ax.plot(M["year"][mm], M["sec_depth_pct"][mm], ls, color=t["mark2"], lw=1.0, alpha=0.75)
    ax.set_xlim(2018.5, 2030); ax.set_ylim(-0.05, 3.4); ax.set_xlabel("year"); ax.set_ylabel("eclipse depth (%)"); style(ax, t)
    title(ax, t, "TIC 124751008: eclipse depth against time (TESS sectors 6-109; lines: inclination models,\nextrapolated after 2026.8; blue: primary, orange: secondary)")
    ax.legend(loc="upper center", fontsize=7.5, frameon=False, labelcolor=t["ink2"], ncol=2)
    fig.tight_layout(); save(fig, t, name)


if __name__ == "__main__":
    for th, t in THEMES.items():
        plt.rcParams.update({"text.color": t["ink"], "axes.edgecolor": t["axis"], "axes.titlecolor": t["ink"]})
        sel = sys.argv[1:] or [str(x) for x in TICS] + ["other", "c311179599", "eb124751008"]
        if "other" in sel:
            fig_eb_history(t, f"eb_history_11613676_{th}.svg"); fig_eb_profiles(t, f"eb_profiles_11613676_{th}.svg")
            fig_toi540_sectors(t, f"toi540_sectors_{th}.svg"); fig_toi540_depths(t, f"toi540_depths_{th}.svg")
        if "c311179599" in sel:
            fig_epochs_c(t, f"epochs_311179599c_{th}.svg")
        if "eb124751008" in sel:
            fig_eb_depths_124(t, f"eb_depths_124751008_{th}.svg")
        for tic in [int(x) for x in sel if x not in ("other", "c311179599", "eb124751008")]:
            fig_fold(t, f"fold_{tic}_{th}.svg", tic)
            fig_epochs(t, f"epochs_{tic}_{th}.svg", tic)
            fig_oc(t, f"oc_{tic}_{th}.svg", tic)
            fig_diff(t, f"diffimage_{tic}_{th}.svg", tic)
    print("figures written to", OUT)
