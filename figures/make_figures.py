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
'tdv': tdv_307210830 (per-sector transit duration of L 98-59 d with the JWST fits and the linear fit; controls b and c),
       tdv_296739893 (per-sector transit duration and width of TOI-620 b in three photometries).
'toi6482': transits_189300053 (the six TESS QLP transits of TOI-6482.01 with the common model, and the O-C diagram).
'jwstchrom': jwst_chromatic_307210830 (JWST transit duration of L 98-59 d per wavelength bin and per detector, the change
       between the two NIRSpec visits per bin, and the per-detector and detector-averaged durations per visit).
'toi2638': depth_t14_458342769 (transit depth of TOI-2638.01 per TESS sector in every product, and the per-sector
       trapezoid duration, with model M2 and the extrapolation to 2027.0).
'toi1117': fold_295541511 (the 4.5787-d signal of TOI-1117 c folded on the selected model sample of the adopted fit: SPOC 2-min sectors 39, 100-103 and
       TESS-SPOC sector 13, with the adopted and the free-density transit models).
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


def _lin(x, y, e, x0=2022.0):
    A = np.vstack([np.ones_like(x), x - x0]).T; W = 1 / e**2
    cov = np.linalg.inv(A.T @ (A * W[:, None])); c = cov @ (A.T @ (W * y)); return c


def fig_tdv_l9859(t, name):
    D = read("durations_307210830.csv"); S = read("duration_slopes_307210830.csv"); J = read("jwst_fits_307210830.csv")
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1.35, 1])); fig.patch.set_facecolor(t["surface"]); ax = axs[0]
    m = (D["planet"] == "d") & (D["variant"] == "2m_pdc_ttv") & (D["in_regression"] == 1)
    x, y, e = D["year"][m], D["T14_min"][m], D["T14_err_min"][m]
    ax.errorbar(x, y, e, fmt="o", color=t["mark"], ms=4.5, capsize=2, label="TESS, one point per sector (2-min PDCSAP; template a/R* = 34.8)")
    c = _lin(x, y, e); xx = np.array([2018.4, 2026.8]); j = np.where((S["planet"] == "d") & (S["variant"] == "2m_pdc_ttv") & (S["quantity"] == "T14"))[0][0]
    ax.plot(xx, c[0] + c[1] * (xx - 2022), "-", color=t["mark"], lw=1.2, label=f"linear fit, {S['slope_per_yr'][j]:+.2f} $\\pm$ {S['slope_err_per_yr'][j]:.2f} min/yr")
    for arm, fmt, col, lab in (("free", "D", t["mark2"], "JWST, a/R* free"), ("33.6", "s", t["mark3"], "JWST, a/R* = 33.6")):
        k = (J["planet"] == "d") & (J["mode"] == "single") & (J["aR_mode"] == arm)
        ax.errorbar(J["year"][k], J["T14_min"][k], [J["T14_minus"][k], J["T14_plus"][k]], fmt=fmt, color=col, ms=6, capsize=2, label=lab + " (GTO 1201, GTO 1224, GO 4098)")
    ax.set_xlim(2018.4, 2026.8); ax.set_xlabel("year"); ax.set_ylabel("transit duration T14 (min)"); style(ax, t)
    title(ax, t, "L 98-59 d: transit duration (first to fourth contact) per TESS sector and per JWST visit\n(absolute values depend on the a/R* of each fit; the slopes are the compared quantity)")
    ax.legend(loc="upper right", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    ax = axs[1]
    for p, col, off in (("d", t["mark"], 0.0), ("c", t["mark3"], 0.05), ("b", t["mark2"], -0.05)):
        v = "2m_pdc_ttv" if p != "b" else "2m_pdc_lin"
        m = (D["planet"] == p) & (D["variant"] == v) & (D["in_regression"] == 1); x, y, e = D["year"][m], D["T14_min"][m], D["T14_err_min"][m]
        mu = np.sum(y / e**2) / np.sum(1 / e**2); j = np.where((S["planet"] == p) & (S["variant"] == v) & (S["quantity"] == "T14"))[0][0]
        ax.errorbar(x + off, y - mu, e, fmt="o", color=col, ms=3.5, capsize=1.5, alpha=0.9, label=f"{p}: {S['slope_per_yr'][j]:+.2f} $\\pm$ {S['slope_err_per_yr'][j]:.2f} min/yr")
        c = _lin(x, y - mu, e); ax.plot(xx, c[0] + c[1] * (xx - 2022), "-", color=col, lw=1.0)
    ax.axhline(0, color=t["axis"], lw=0.8); ax.set_xlim(2018.4, 2026.8); ax.set_xlabel("year"); ax.set_ylabel("T14 minus weighted mean (min)"); style(ax, t)
    title(ax, t, "same data and code: planets b, c and d\n(per-sector T14 relative to the mean; lines: linear fits)")
    ax.legend(loc="lower left", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    fig.tight_layout(); save(fig, t, name)


def fig_tdv_toi620(t, name):
    D = read("durations_296739893.csv")
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2)); fig.patch.set_facecolor(t["surface"])
    for ax, q, lab in ((axs[0], "T14", "transit duration T14 (min), b-only fit"), (axs[1], "w", "width w relative to the template")):
        for v, col, fmt, off, nm in (("2m_pdc_sec", t["mark"], "o", 0.0, "SPOC PDCSAP"), ("2m_sap_pt", t["mark2"], "s", 0.06, "SPOC SAP"),
                                     ("own_r1.5_sec", t["mark3"], "D", -0.06, "own TPF, r = 1.5 px"), ("own_r3.0_sec", t["muted"], "^", 0.12, "own TPF, r = 3.0 px")):
            if q == "T14" and v == "own_r3.0_sec": continue
            m = D["variant"] == v; y = D[q if q == "w" else "T14_min"][m]
            lo = D["w_lo" if q == "w" else "T14_lo_min"][m]; hi = D["w_hi" if q == "w" else "T14_hi_min"][m]
            ax.errorbar(D["year"][m] + off, y, [y - lo, hi - y], fmt=fmt, color=col, ms=5, capsize=2, label=nm)
        ax.set_xlabel("year"); ax.set_ylabel(lab); style(ax, t); ax.legend(loc="lower left", fontsize=8, frameon=False, labelcolor=t["ink2"])
    title(axs[0], t, "TOI-620 b: transit duration per TESS sector (8, 35, 62, 89, 99)"); title(axs[1], t, "TOI-620 b: time-stretch factor w of the transit template (free depth)")
    fig.tight_layout(); save(fig, t, name)


def fig_jwst_chrom(t, name):
    B = read("jwst_chromatic_307210830.csv"); Q = read("jwst_chromatic_summary_307210830.csv")
    fig = plt.figure(figsize=(12.5, 9)); fig.patch.set_facecolor(t["surface"]); gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
    ax = fig.add_subplot(gs[0, :]); bx = fig.add_subplot(gs[1, 0]); cx = fig.add_subplot(gs[1, 1])
    vis = (("d1224_prior", t["mark"], "o", -0.012, "GTO 1224 obs. 1 (NIRSpec, 2023.48)"), ("d4098_prior", t["mark2"], "s", 0.012, "GO 4098 obs. 2 (NIRSpec, 2024.20)"),
           ("d1201_prior", t["mark3"], "D", 0.0, "GTO 1201 obs. 311 (NIRISS SOSS order 1, 2022.99)"))
    for f, col, fmt, off, lab in vis:
        m = (B["fit"] == f) & (B["row"] == "bin") & (B["used"] == 1); lam = 0.5 * (B["lam_lo_um"][m] + B["lam_hi_um"][m])
        ax.errorbar(lam + off, B["T14_min"][m], B["T14_err_min"][m], fmt=fmt, color=col, ms=4.5, capsize=2, label=lab)
        for ch in ("NRS1", "NRS2", "SOSS1"):
            k = (B["fit"] == f) & (B["row"] == "channel, shared b") & (B["channel"] == ch)
            if k.any():
                y, e = B["T14_min"][k][0], B["T14_err_min"][k][0]; x0, x1 = B["lam_lo_um"][k][0], B["lam_hi_um"][k][0]
                ax.fill_between([x0, x1], y - e, y + e, color=col, alpha=0.25, lw=0); ax.plot([x0, x1], [y, y], "-", color=col, lw=1.6)
    for a in (ax, bx):
        a.axvspan(3.71, 3.83, color=t["grid"], alpha=0.9, lw=0)
    ax.text(3.77, 47.35, "detector\ngap", ha="center", va="top", fontsize=7.5, color=t["muted"])
    ax.set_xlim(0.8, 5.25); ax.set_ylim(42.6, 47.4); ax.set_xlabel("wavelength (µm)"); ax.set_ylabel("T14 per bin (min)"); style(ax, t)
    title(ax, t, "L 98-59 d: JWST transit duration per wavelength bin (Kurucz-prior limb darkening; P, mid-time and a/R* = 36.25 fixed per visit)\n"
              "(bands: one b per detector, T14 with its approximate conditional error (see section 6.6.11); NRS1 2.87-3.71 µm, NRS2 3.83-5.15 µm)")
    ax.legend(loc="lower left", fontsize=8, frameon=False, labelcolor=t["ink2"])
    a1 = (B["fit"] == "d1224_prior") & (B["row"] == "bin"); a2 = (B["fit"] == "d4098_prior") & (B["row"] == "bin")
    lam = 0.5 * (B["lam_lo_um"][a1] + B["lam_hi_um"][a1]); d = B["T14_min"][a2] - B["T14_min"][a1]; e = np.hypot(B["T14_err_min"][a1], B["T14_err_min"][a2])
    bx.errorbar(lam, d, e, fmt="o", color=t["ink"], ms=4, capsize=2, label="per bin")
    for ch, x0, x1 in (("NRS1", 2.87, 3.71), ("NRS2", 3.83, 5.15)):
        k = (Q["ld"] == "prior") & (Q["quantity"] == f"{ch}, weighted mean of the per-bin delta T14 (min)")
        y, s = Q["value"][k][0], Q["error"][k][0]
        bx.fill_between([x0, x1], y - s, y + s, color=t["mark"], alpha=0.25, lw=0); bx.plot([x0, x1], [y, y], "-", color=t["mark"], lw=1.6)
        bx.text(0.5 * (x0 + x1), 1.75, f"{ch} mean: {y:+.2f} $\\pm$ {s:.2f} min", ha="center", va="top", fontsize=8, color=t["ink2"])
    bx.axhline(0, color=t["axis"], lw=0.8); bx.axhline(-0.85, color=t["mark3"], lw=1.2, ls="--", label="TESS trend, -0.85 $\\pm$ 0.08 min")
    bx.set_xlim(2.8, 5.22); bx.set_ylim(-2.9, 1.9); bx.set_xlabel("wavelength (µm)"); bx.set_ylabel("T14(2024.20) - T14(2023.48) (min)"); style(bx, t)
    title(bx, t, "change between the two NIRSpec visits, per bin\n(bands: weighted mean per detector, formal error)")
    bx.legend(loc="lower right", fontsize=8, frameon=False, labelcolor=t["ink2"])
    yr = {"d1224_prior": 2023.481, "d4098_prior": 2024.195}
    for ch, fmt, col, off in (("NRS1", "o", t["mark"], -0.012), ("NRS2", "s", t["mark2"], 0.012)):
        xs, ys, es = [], [], []
        for f in ("d1224_prior", "d4098_prior"):
            k = (B["fit"] == f) & (B["row"] == "channel, shared b") & (B["channel"] == ch); xs.append(yr[f] + off); ys.append(B["T14_min"][k][0]); es.append(B["T14_err_min"][k][0])
        cx.errorbar(xs, ys, es, fmt=fmt, color=col, ms=5, capsize=2, lw=1, ls=":", label=f"{ch} alone (approximate conditional error)")
    xs, ys, es = [], [], []
    for y0 in (2023.481, 2024.195):
        k = (Q["ld"] == "prior") & (Q["quantity"] == f"detector-averaged T14, {y0} (min)"); xs.append(y0); ys.append(Q["value"][k][0]); es.append(Q["error"][k][0])
    k = (B["fit"] == "d1201_prior") & (B["row"] == "channel, shared b") & (B["channel"] == "SOSS1"); fl = Q["value"][(Q["ld"] == "prior") & (Q["quantity"] == "empirical systematic floor per detector and visit (min)")][0]
    cx.errorbar(xs, ys, es, fmt="D", color=t["ink"], ms=6, capsize=3, label="NIRSpec detector average (error incl. floor / $\\sqrt{2}$)")
    cx.errorbar([2022.992], B["T14_min"][k], [np.hypot(B["T14_err_min"][k][0], fl)], fmt="D", color=t["mark3"], ms=6, capsize=3, label="NIRISS SOSS order 1 (error incl. floor)")
    xa = np.array(xs + [2022.992]); ya = np.array(ys + [B["T14_min"][k][0]]); wa = 1 / np.array(es + [np.hypot(B["T14_err_min"][k][0], fl)]) ** 2
    xx = np.array([2022.85, 2024.35]); y0 = np.sum(wa * (ya + 1.185 * (xa - 2023.5))) / np.sum(wa)
    cx.plot(xx, y0 - 1.185 * (xx - 2023.5), "--", color=t["mark3"], lw=1.0, label="slope of the TESS trend (offset fitted)")
    cx.set_xlim(2022.85, 2024.35); cx.set_xlabel("year"); cx.set_ylabel("T14 (min)"); style(cx, t)
    title(cx, t, "per-detector and detector-averaged T14 per visit\n(averages and floor defined after the detector difference was seen)")
    cx.legend(loc="lower left", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    fig.tight_layout(); save(fig, t, name)


def fig_toi6482(t, name):
    L = read("transits_189300053.csv"); T = read("transit_times_189300053.csv"); E = read("ephemeris_189300053.csv")
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(8, 10.5), gridspec_kw=dict(height_ratios=[3.3, 1])); fig.patch.set_facecolor(t["surface"])
    tess = T["in_fit"] == 1
    for k, (s, e) in enumerate(zip(T["sector"][tess], T["epoch"][tess])):
        m = (L["sector"] == s) & (L["epoch"] == e); x = L["dt_lin_h"][m]; o = np.argsort(x); off = -0.022 * k
        ax.plot(x, L["flux"][m] + off, ".", color=t["pts"], ms=1.6)
        xb, yb, _ = binned(x, L["flux"][m], 0.5); ax.plot(xb, yb + off, "o", color=t["mark"], ms=3.2)
        ax.plot(x[o], L["model"][m][o] + off, "-", color=t["ink"], lw=1.1)
        ax.text(-15.6, 1.0075 + off, f"sector {int(s)}, epoch {int(e)}", color=t["ink2"], fontsize=8)
    ax.set_xlim(-16, 16); ax.set_ylim(0.87, 1.014); ax.set_xlabel("hours from mid-transit (linear ephemeris)"); ax.set_ylabel("relative flux (baseline divided out) + offset")
    style(ax, t); title(ax, t, f"TOI-6482.01 (TIC 189300053): six TESS QLP transits, 2019-2026, P = {E['P_d'][0]:.6f} d\n(grey: cadences; blue: 30-min bins; line: common exposure-averaged trapezoid at the fitted mid-time)")
    bx.axhline(0, color=t["axis"], lw=0.8)
    bx.errorbar(T["epoch"][tess], T["o_c_min"][tess], T["tc_err_min"][tess], fmt="o", color=t["ink"], ms=5, capsize=2, label="TESS QLP (fitted)")
    g = ~tess
    if g.any():
        bx.errorbar(T["epoch"][g], T["o_c_min"][g], T["o_c_err_min"][g], fmt="s", color=t["mark2"], ms=5, capsize=2, label="TRAPPIST-South 2025-04-21, reduction v2 (not in the fit; error includes the ephemeris uncertainty)")
    bx.set_xlabel("epoch"); bx.set_ylabel("O - C (min)"); style(bx, t)
    bx.legend(loc="upper left", fontsize=8, frameon=False, labelcolor=t["ink2"])
    fig.tight_layout(); save(fig, t, name)


def fig_toi2638(t, name):
    D = read("depths_458342769.csv"); S = read("shape_458342769.csv"); X = read("pixels_458342769.csv"); R = read("predictions_458342769.csv"); TR = read("trends_458342769.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 4.8)); fig.patch.set_facecolor(t["surface"])
    for prod, fmt, col, off, lab in (("SPOC PDCSAP", "o", t["mark"], 0.0, "SPOC 2-min PDCSAP"), ("SPOC SAP crowding-corrected", "s", t["ink2"], 0.06, "SPOC 2-min SAP / CROWDSAP"),
                                     ("TESS-SPOC PDCSAP", "D", t["mark3"], -0.06, "TESS-SPOC FFI"), ("QLP KSPSAP", "^", t["mark2"], 0.12, "QLP FFI (not deblended)"),
                                     ("TGLC cal_psf_flux", "v", t["muted"], -0.12, "TGLC PSF (sector 10)")):
        m = D["product"] == prod
        ax.errorbar(D["year"][m] + off, D["depth_ppm"][m] / 1e3, D["err_ppm"][m] / 1e3, fmt=fmt, color=col, ms=5, capsize=2, label=lab)
    yx = np.array([D["year"][(D["product"] == "SPOC PDCSAP") & (D["sector"] == s)][0] for s in X["sector"]])
    ax.errorbar(yx - 0.18, X["ap1.5_undiluted_ppm"] / 1e3, X["ap1.5_err_ppm"] / X["ap1.5_target_fraction"] / 1e3, fmt="P", color=t["mark"], mfc="none", ms=6, capsize=2,
                label="own TPF photometry, r = 1.5 px, scene-model dilution removed")
    m = np.char.startswith(R["source"].astype(str), "M2 (SPOC PDCSAP + QLP sector 10)") & np.isfinite(R["depth_ppm"])
    ax.plot(R["year"][m], R["depth_ppm"][m] / 1e3, "-", color=t["ink"], lw=1.0, label="model M2 (impact parameter linear in time), peak depth")
    j = np.where((TR["quantity"] == "depth slope") & (TR["variant"] == "SPOC PDCSAP + QLP sector 10"))[0][0]
    ax.set_xlim(2018.7, 2027.6); ax.set_xlabel("year"); ax.set_ylabel("transit depth (ppt)"); style(ax, t)
    title(ax, t, f"TOI-2638.01: transit depth per TESS sector and product\n(common transit shape; linear trend {TR['value'][j]:+.0f} $\\pm$ {TR['err'][j]:.0f} ppm/yr, SPOC PDCSAP + QLP sector 10)")
    ax.legend(loc="upper left", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    ok = np.isfinite(S["trapezoid_T14_min"])
    bx.errorbar(S["year"][ok], S["trapezoid_T14_min"][ok], S["trapezoid_T14_err_min"][ok], fmt="o", color=t["mark"], ms=5, capsize=2,
                label="trapezoid T14 per sector, depth free (SPOC PDCSAP; bootstrap errors)")
    m2 = m & np.isfinite(R["T14_min"])
    bx.plot(R["year"][m2], R["T14_min"][m2], "-", color=t["ink"], lw=1.0, label="model M2, first to fourth contact")
    k = np.char.startswith(R["source"].astype(str), "linear trapezoid") & (np.abs(R["year"] - 2027.0) < 1e-6)
    bx.errorbar(R["year"][k] + 0.05, R["T14_min"][k], R["T14_err_min"][k], fmt="s", mfc="none", color=t["mark2"], ms=6, capsize=2,
                label="linear-trend extrapolation to 2027.0 (TESS sectors 111-113)")
    bx.set_xlim(2018.7, 2027.6); bx.set_xlabel("year"); bx.set_ylabel("transit duration T14 (min)"); style(bx, t)
    title(bx, t, "TOI-2638.01: transit duration per SPOC sector\n(trapezoid with free depth: independent of dilution)")
    bx.legend(loc="lower right", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    fig.tight_layout(); save(fig, t, name)


def fig_toi1117(t, name):
    F = read("folded_295541511.csv"); M = read("model_295541511.csv"); E = read("ephemeris_295541511.csv"); R = read("fit_295541511.csv")
    bval = lambda cfg: R["value"][(R["config"] == cfg) & (R["quantity"] == "b")][0]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.6, 1])); fig.patch.set_facecolor(t["surface"])
    for ax, sel, bw, mcol, ttl in ((axs[0], F["sector"] != 13, 15 / 60, "rho_2min", "SPOC 2-min PDCSAP, sectors 39, 100-103"),
                                   (axs[1], F["sector"] == 13, 30 / 60, "rho_30min", "TESS-SPOC FFI (1800 s), sector 13: not in the adopted fit")):
        x = F["dt_days"][sel] * 24; y = (F["flux"][sel] - 1) * 1e6
        ax.plot(x, y, ".", color=t["pts"], ms=1.2 if bw < 0.3 else 3, rasterized=bw < 0.3, label="detrended cadences")
        xb, yb, eb = binned(x, y, bw)
        ax.errorbar(xb, yb, eb, fmt="o", color=t["ink"], ms=3.2, lw=0.8, label=f"{int(round(bw * 60))}-min bins")
        ax.plot(M["dt_days"] * 24, (M[mcol] - 1) * 1e6, "-", color=t["mark2"], lw=1.6, label=f"adopted fit (stellar-density prior, b = {bval('rho_2min_pdc'):.2f})")
        if bw < 0.3:
            ax.plot(M["dt_days"] * 24, (M["free_2min"] - 1) * 1e6, "--", color=t["mark"], lw=1.3, label="free-density model (uniform q1, q2)")
        ax.set_xlim(-10.8, 10.8); ax.set_ylim(-1500, 900); ax.set_xlabel("hours from mid-transit (selected model sample)"); ax.set_ylabel("relative flux - 1 (ppm)")
        style(ax, t); title(ax, t, ttl); ax.legend(loc="lower left", fontsize=7.5, frameon=False, labelcolor=t["ink2"])
    i = list(E["role"]).index("adopted")
    fig.suptitle(f"TOI-1117 (TIC 295541511): transit-like signal at the period of planet c, P = {E['P_d'][i]:.7f} d (planet b masked)", color=t["ink"], fontsize=11, x=0.01, ha="left")
    fig.tight_layout(); save(fig, t, name)


if __name__ == "__main__":
    for th, t in THEMES.items():
        plt.rcParams.update({"text.color": t["ink"], "axes.edgecolor": t["axis"], "axes.titlecolor": t["ink"]})
        sel = sys.argv[1:] or [str(x) for x in TICS] + ["other", "c311179599", "eb124751008", "tdv", "toi6482", "jwstchrom", "toi2638", "toi1117"]
        if "other" in sel:
            fig_eb_history(t, f"eb_history_11613676_{th}.svg"); fig_eb_profiles(t, f"eb_profiles_11613676_{th}.svg")
            fig_toi540_sectors(t, f"toi540_sectors_{th}.svg"); fig_toi540_depths(t, f"toi540_depths_{th}.svg")
        if "c311179599" in sel:
            fig_epochs_c(t, f"epochs_311179599c_{th}.svg")
        if "eb124751008" in sel:
            fig_eb_depths_124(t, f"eb_depths_124751008_{th}.svg")
        if "tdv" in sel:
            fig_tdv_l9859(t, f"tdv_307210830_{th}.svg"); fig_tdv_toi620(t, f"tdv_296739893_{th}.svg")
        if "toi6482" in sel:
            fig_toi6482(t, f"transits_189300053_{th}.svg")
        if "jwstchrom" in sel:
            fig_jwst_chrom(t, f"jwst_chromatic_307210830_{th}.svg")
        if "toi2638" in sel:
            fig_toi2638(t, f"depth_t14_458342769_{th}.svg")
        if "toi1117" in sel:
            fig_toi1117(t, f"fold_295541511_{th}.svg")
        for tic in [int(x) for x in sel if x not in ("other", "c311179599", "eb124751008", "tdv", "toi6482", "jwstchrom", "toi2638", "toi1117")]:
            fig_fold(t, f"fold_{tic}_{th}.svg", tic)
            fig_epochs(t, f"epochs_{tic}_{th}.svg", tic)
            fig_oc(t, f"oc_{tic}_{th}.svg", tic)
            fig_diff(t, f"diffimage_{tic}_{th}.svg", tic)
    print("figures written to", OUT)
