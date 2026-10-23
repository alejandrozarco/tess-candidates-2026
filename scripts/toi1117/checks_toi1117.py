"""TOI-1117 c (TIC 295541511), README section 6.10: false-alarm tests, systematics per transit, odd/even and secondary-eclipse
tests, the data of sectors 13 and 39 alone, transit timing with an exposure-integrated trapezoid (c and b), sector-13 timing
injections, and the TCEs of the data-validation reports.

usage: python scripts/toi1117/checks_toi1117.py [--work ./work] [--n 1000]
input:  <work>/295541511/dl/ (download_toi1117.py), data/ephemeris_295541511.csv (ephemeris_toi1117.py)
output: data/falsealarm_295541511.csv  S/N maxima of the false-alarm tests (one row per trial)
        data/epochs_295541511.csv      per predicted transit: coverage, quality flags, momentum dumps, gaps, background, centroids
        data/checks_295541511.csv      summary of all tests
        data/timing_295541511.csv      per-product trapezoid mid-times of c and b, and the sector-13 injection test
Search light curve of the false-alarm tests: TESS-SPOC sector 13 (1800 s) and SPOC 2-min sectors 39, 100-103 binned to 10 min,
b masked, biweight 1.0 d with c NOT masked, 5-sigma clip; astropy BoxLeastSquares (objective 'snr'), six durations
1.0-3.2 h scaled as P^(1/3); the RV window is P = 4.579 +- 0.0114 d (3 sigma of the RV period of Lockley et al. 2025), all phases.
Pseudo windows: random centres log-uniform in 2.5-12 d with the same relative width, avoiding n/m aliases (n, m <= 4) of b, c and
d, searched on the data with c's transits removed (+-D) and on the inverted flux (2 - f). Shuffle: each sector shifted by a
random time in [0, 4.579) d (c kept). Seeded (numpy default_rng(42) and (7)).
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from astropy.timeseries import BoxLeastSquares
from astropy.io import fits
from wotan import flatten

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--n', type=int, default=1000); args = p.parse_args()
W = args.work; DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
LABS6 = ['TSPOC_13_pdc'] + [f'SPOC_{s}_pdc' for s in SPOC2]
summ = []


def add(check, variant, value, err='', note=''):
    summ.append((check, variant, value, err, note))


def binned(t, f, dt):
    k = np.floor((t - t[0]) / dt).astype(np.int64); u, inv, c = np.unique(k, return_inverse=True, return_counts=True)
    return np.bincount(inv, t) / c, np.bincount(inv, f) / c, c


def build(invert=False):
    T, R, S = [], [], []
    for lab in LABS6:
        L = load(W, lab); t, f = L['t'], L['f']
        if invert:
            f = 2 - f
        k = mask_b(t) & trim_gaps(t); t, f = t[k], f[k]
        r = detrend(t, f, 1.0) - 1; ok = np.isfinite(r); t, r = t[ok], r[ok]
        if L['cad'] < 0.005:
            t, r, c = binned(t, r, 10 / 1440); t, r = t[c >= 3], r[c >= 3]
        sd = robust_sd(r); cl = np.abs(r) < 5 * sd
        T.append(t[cl]); R.append(r[cl]); S.append(np.full(cl.sum(), sd))
    return np.concatenate(T), np.concatenate(R), np.concatenate(S), [len(x) for x in T]


DUR0 = np.array([1.0, 1.3, 1.6, 2.0, 2.5, 3.2]) / 24


def win_search(t, r, dy, Pc, rel=0.0114 / 4.579):
    T = t.max() - t.min(); durs = DUR0 * (Pc / 4.579) ** (1 / 3)
    P = np.arange(Pc * (1 - rel), Pc * (1 + rel), 0.3 * durs.min() * Pc / T)
    res = BoxLeastSquares(t, r, dy).power(P, durs, objective='snr'); i = int(np.argmax(res.depth_snr))
    return float(res.depth_snr[i]), float(res.period[i]), float(res.transit_time[i]), float(res.duration[i]), float(res.depth[i])


def c_removed(t, D=DBOX):
    return np.abs(phase(t, PC, TC)) > D


# 1. false-alarm tests
t, r, dy, ns = build(); ti, ri, dyi, _ = build(invert=True)
real = win_search(t, r, dy, 4.579); inv = win_search(ti, ri, dyi, 4.579)
add('RV-window search (P = 4.579 +- 0.0114 d), real data', 'BLS S/N; period; depth', f'{real[0]:.2f}', '', f'P = {real[1]:.6f} d, duration {real[3] * 24:.2f} h, depth {real[4] * 1e6:.0f} ppm')
add('RV-window search, inverted flux (2 - f)', 'BLS S/N', f'{inv[0]:.2f}', '', f'P = {inv[1]:.6f} d')
print('real', real, 'inverted', inv, flush=True)
rng = np.random.default_rng(42); bad = [k * n / m for k in [PB, 4.579, PD] for n in range(1, 5) for m in range(1, 5)]
kc = c_removed(t); kci = c_removed(ti); fa = []
while len(fa) < args.n:
    Pc = float(np.exp(rng.uniform(np.log(2.5), np.log(12))))
    if min(abs(Pc / x - 1) for x in bad) < 0.006:
        continue
    fa.append((Pc, win_search(t[kc], r[kc], dy[kc], Pc)[0], win_search(ti[kci], ri[kci], dyi[kci], Pc)[0]))
    if len(fa) % 100 == 0:
        print(len(fa), flush=True)
fa = np.array(fa)
for j, lab in ((1, 'pseudo windows, real data with c removed'), (2, 'pseudo windows, inverted flux with c removed')):
    x = fa[:, j]
    add(f'false alarms: {lab}', f'{args.n} windows: max; 99th percentile; median', f'{x.max():.2f}', '', f'p99 {np.percentile(x, 99):.2f}, median {np.median(x):.2f}; '
        f'{int((x >= real[0]).sum())} of {args.n} reach the observed {real[0]:.2f}; {int((x >= 8).sum())} reach 8')
edges = np.cumsum([0] + ns); sh = []
for k in range(args.n // 2):
    tt = t.copy()
    for i in range(len(ns)):
        tt[edges[i]:edges[i + 1]] += rng.uniform(0, 4.579)
    sh.append(win_search(tt, r, dy, 4.579)[0])
sh = np.array(sh)
add('false alarms: per-sector random time shifts (c kept), RV window', f'{len(sh)} shuffles: max; 99th percentile; median', f'{sh.max():.2f}', '',
    f'p99 {np.percentile(sh, 99):.2f}, median {np.median(sh):.2f}; {int((sh >= real[0]).sum())} of {len(sh)} reach {real[0]:.2f}')
write_csv(os.path.join(DATA, 'falsealarm_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): false-alarm tests of the RV-window search; scripts/toi1117/checks_toi1117.py, {DATE}',
           f'observed: BLS S/N {real[0]:.2f} at P = {real[1]:.6f} d in the RV window (4.579 +- 0.0114 d); search light curve and windows as in the script docstring',
           'test pseudo_real / pseudo_inverted: window centre P_centre_d, maximum S/N in the window (c removed); test shuffle: per-sector random time shifts, RV window (P_centre_d = 4.579)'],
          ['test', 'P_centre_d', 'max_snr'],
          [('pseudo_real', f'{a:.5f}', f'{b:.3f}') for a, b, c in fa] + [('pseudo_inverted', f'{a:.5f}', f'{c:.3f}') for a, b, c in fa] + [('shuffle', '4.579', f'{x:.3f}') for x in sh])

# 2. data of sectors 13 and 39 only (the TESS data available to Lockley et al. 2025)
for invert in (False, True):
    tt, rr, dd, _ = build(invert=invert); k = tt < 3000; tt, rr, dd = tt[k], rr[k], dd[k]
    w = win_search(tt, rr, dd, 4.579)
    T = tt.max() - tt.min(); fq = np.arange(1 / 15, 1 / 1.0, 0.3 * (1.0 / 24) / T / 1.5); P = np.sort(1 / fq)
    res = BoxLeastSquares(tt, rr, dd).power(P, DUR0, objective='snr'); o = np.argsort(res.depth_snr)[::-1]; top = []
    for i in o:
        if all(abs(res.period[i] / x[0] - 1) > 0.002 for x in top):
            top.append((float(res.period[i]), float(res.depth_snr[i])))
        if len(top) == 5:
            break
    lab = 'inverted flux' if invert else 'real data'
    add(f'sectors 13 + 39 only, {lab}: RV-window search', 'BLS S/N', f'{w[0]:.2f}', '', f'P = {w[1]:.5f} d')
    add(f'sectors 13 + 39 only, {lab}: blind BLS 1-15 d ({len(P)} periods)', 'highest S/N; top five peaks (P d: S/N)', f'{top[0][1]:.2f}', '', '; '.join(f'{a:.4f}: {b:.2f}' for a, b in top))
X = [prep(W, l) for l in LABS6]; tt = np.concatenate([x['t'] for x in X]); rr = np.concatenate([x['r'] for x in X])
for lab, sel in (('sectors 13 + 39', tt < 3000), ('sectors 100-103', tt > 3000)):
    d, e, k = box(tt[sel], rr[sel], PC, TBOX, DBOX); add('stacked box depth', lab, f'{d * 1e6:.0f}', f'{e * 1e6:.0f}', f'S/N {d / e:.1f}')

# 3. systematics per predicted transit (all cadences of the SPOC / TESS-SPOC files)
files = {13: glob.glob(os.path.join(dldir(W, 'lc'), 'hlsp_tess-spoc*s0013_tess_v1_lc.fits'))[0]}
for s in SPOC2:
    files[s] = glob.glob(os.path.join(dldir(W, 'lc'), f'tess*-s{s:04d}-0000000295541511-*-s_lc.fits'))[0]
ep, cen, dumps_all = [], [], {}
for s, fp in files.items():
    with fits.open(fp) as h:
        d = h[1].data; tq = np.array(d['TIME'], float); q = np.array(d['QUALITY']); fl = np.array(d['PDCSAP_FLUX'], float)
        bk = np.array(d['SAP_BKG'], float); c1 = np.array(d['MOM_CENTR1'], float); c2 = np.array(d['MOM_CENTR2'], float)
    ok = np.isfinite(tq); tq, q, fl, bk, c1, c2 = tq[ok], q[ok], fl[ok], bk[ok], c1[ok], c2[ok]
    dumps = tq[(q & 32) > 0]; dumps_all[s] = dumps
    good = (q == 0) & np.isfinite(fl); tg = tq[good]; gaps = np.where(np.diff(tg) > 0.15)[0]; edg = np.r_[tg[0], tg[gaps], tg[gaps + 1], tg[-1]]
    for n in np.unique(np.round((tq - TBOX) / PC)):
        tc = TBOX + n * PC
        if tc < tg[0] - DBOX or tc > tg[-1] + DBOX:
            continue
        win = np.abs(tq - tc) < DBOX / 2; wg = win & good; loc = (np.abs(tq - tc) < 0.6) & (np.abs(tq - tc) > 0.75 * DBOX) & good
        flags = int(np.bitwise_or.reduce(q[win])) if win.any() else -1
        bz = (np.nanmean(bk[wg]) - np.nanmedian(bk[loc])) / (robust_sd(bk[loc]) / np.sqrt(max(wg.sum(), 1))) if wg.sum() > 2 and loc.sum() > 10 else np.nan
        cz = []
        for c in (c1, c2):
            if wg.sum() > 2 and loc.sum() > 10:
                A = np.vstack([np.ones(loc.sum()), tq[loc] - tc]).T; cf = np.linalg.lstsq(A, c[loc], rcond=None)[0]
                cz.append((float(np.nanmean(c[wg] - (cf[0] + cf[1] * (tq[wg] - tc)))), float(robust_sd(c[loc] - A @ cf) / np.sqrt(wg.sum()))))
            else:
                cz.append((np.nan, np.nan))
        cen.append([cz[0][0], cz[0][1], cz[1][0], cz[1][1]])
        orb = float(tc - edg[edg <= tc + 1e-9].max()) if (edg <= tc).any() else np.nan
        ep.append([s, int(n), tc, wg.sum() / max(win.sum(), 1) if win.any() else 0.0, flags, float(np.min(np.abs(dumps - tc))) * 24 if len(dumps) else np.nan,
                   float(np.min(np.abs(edg - tc))), bz, cz[0][0] * 1e3, cz[0][1] * 1e3, cz[1][0] * 1e3, cz[1][1] * 1e3, orb])
cen = np.array(cen, float); okc = np.isfinite(cen).all(1)
def wm(v, e):
    w = 1 / e ** 2; return np.sum(w * v) / np.sum(w), 1 / np.sqrt(np.sum(w))
mc, ec = wm(cen[okc, 0], cen[okc, 1]); mr, er = wm(cen[okc, 2], cen[okc, 3])
add('flux-weighted centroid in transit minus local baseline, stacked over transits (SPOC / TESS-SPOC MOM_CENTR)', 'column (mpx)', f'{mc * 1e3:.2f}', f'{ec * 1e3:.2f}', f'{okc.sum()} transits')
add('flux-weighted centroid in transit minus local baseline, stacked over transits (SPOC / TESS-SPOC MOM_CENTR)', 'row (mpx)', f'{mr * 1e3:.2f}', f'{er * 1e3:.2f}', f'{okc.sum()} transits')
for s, v in dumps_all.items():
    u_ = np.unique(np.round(v, 2)); add('momentum dumps (QUALITY bit 32)', f'sector {s}', f'{np.median(np.diff(u_)):.2f}' if len(u_) > 2 else '', '', f'median spacing (d); {len(u_)} dumps')
used = [e_ for e_ in ep if e_[3] > 0.5]
add('momentum dumps', 'smallest time between a dump and a covered transit (h)', f'{min(e_[5] for e_ in used):.1f}', '', f'median {np.median([e_[5] for e_ in used]):.1f} h; covered = more than half of the box with QUALITY = 0')
flagged = sorted({e_[1] for e_ in ep if e_[3] > 0.5 and ((abs(e_[7]) > 3 if np.isfinite(e_[7]) else False) or e_[6] < 0.35 or e_[5] < 8 or
                  (np.isfinite(e_[8]) and (abs(e_[8] / e_[9]) > 3 or abs(e_[10] / e_[11]) > 3)))})
write_csv(os.path.join(DATA, 'epochs_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): systematics at every predicted transit (P = {PC} d, T0 = BTJD {TBOX:.4f}, box D = {DBOX * 24:.3f} h) in the SPOC 2-min and TESS-SPOC sector-13 files; scripts/toi1117/checks_toi1117.py, {DATE}',
           'coverage = fraction of the box cadences with QUALITY = 0 and finite PDCSAP; flags_in_box = bitwise OR of QUALITY in the box; dump_dt_h = time to the nearest momentum dump (bit 32)',
           'gap_dt_d = time to the nearest edge of a gap > 0.15 d; orbit_phase_d = time since the preceding gap end or sector start; bkg_z = (mean SAP_BKG in the box - local median) / (robust scatter / sqrt(n))',
           'dcol, drow = MOM_CENTR1, MOM_CENTR2 in the box minus a local linear baseline (0.75 D to 0.6 d), in milli-pixels, with errors',
           f'flagged (background |z| > 3, gap edge < 0.35 d, dump < 8 h or centroid > 3 sigma) among covered transits: epochs {", ".join(str(x) for x in flagged)}'],
          ['sector', 'epoch', 'tc_btjd', 'coverage', 'flags_in_box', 'dump_dt_h', 'gap_dt_d', 'bkg_z', 'dcol_mpx', 'dcol_err_mpx', 'drow_mpx', 'drow_err_mpx', 'orbit_phase_d'],
          [(a, b, f'{c:.4f}', f'{d:.2f}', f_, f'{g:.1f}', f'{h:.2f}', f'{i:.1f}', f'{j:.3f}', f'{k:.3f}', f'{l:.3f}', f'{m:.3f}', f'{n:.2f}') for a, b, c, d, f_, g, h, i, j, k, l, m, n in ep])

# 4. odd/even, secondary eclipse, dropping the flagged transits
n = np.round((tt - TBOX) / PC)
for nm, sel in (('odd', n % 2 == 1), ('even', n % 2 == 0)):
    d, e, k = box(tt[sel], rr[sel], 2 * PC, TBOX + (PC if nm == 'odd' else 0), DBOX); add('odd and even transits (box, sectors 13-103)', nm, f'{d * 1e6:.0f}', f'{e * 1e6:.0f}')
    if nm == 'odd':
        do, eo = d, e
    else:
        de, ee_ = d, e
add('odd and even transits (box, sectors 13-103)', 'difference / error', f'{abs(do - de) / np.hypot(eo, ee_):.2f}')
for lab, sec_ in (('sector 13', [13]), ('sector 39', [39]), ('sectors 100-101', [100, 101]), ('sectors 102-103', [102, 103])):
    m = np.isin(np.concatenate([np.full(len(x['t']), x['sec']) for x in X]), sec_); d_ = []
    for par in (1, 0):
        sel = m & (n % 2 == par); d, e, k = box(tt[sel], rr[sel], 2 * PC, TBOX + (PC if par else 0), DBOX); d_.append((d, e))
    add('odd and even transits (box), per group', lab, f'{abs(d_[0][0] - d_[1][0]) / np.hypot(d_[0][1], d_[1][1]):.2f}', '', f'odd {d_[0][0] * 1e6:.0f} +- {d_[0][1] * 1e6:.0f}, even {d_[1][0] * 1e6:.0f} +- {d_[1][1] * 1e6:.0f} ppm (sigma of the difference in value)')
d, e, k = box(tt, rr, PC, TBOX + 0.5 * PC, DBOX); add('secondary eclipse (box D at phase 0.5)', 'depth (ppm)', f'{d * 1e6:.0f}', f'{e * 1e6:.0f}')
sc = np.array([(f_, box(tt, rr, PC, TBOX + f_ * PC, DBOX)) for f_ in np.arange(0.30, 0.701, 0.005)], dtype=object)
snr = np.array([x[0] / x[1] for x in sc[:, 1]]); add('secondary eclipse, phase scan 0.30-0.70', 'largest S/N; phase', f'{snr.max():.2f}', '', f'phase {sc[np.argmax(snr), 0]:.3f}')
keep = ~np.isin(n, flagged) | (np.abs(phase(tt, PC, TBOX)) > DBOX)
d, e, k = box(tt[keep], rr[keep], PC, TBOX, DBOX); ee2, _ = emp_err(tt[keep], rr[keep], PC, TBOX, DBOX)
add('stacked box depth without the flagged transits', f'{len(flagged)} epochs dropped', f'{d * 1e6:.0f}', f'{e * 1e6:.0f}', f'S/N_emp {d / ee2:.1f}')

# 5. trapezoid timing (exposure-integrated; c and b) per product, and sector-13 injections
E = read_csv(os.path.join(DATA, 'ephemeris_295541511.csv')); ia = list(E['role']).index('adopted'); Tc_ref, Pc_ = float(E['Tref_btjd'][ia]), float(E['P_d'][ia])
T14c, T12c = 1.85 / 24, 0.25 / 24


def trap_int(dt, T14, T12, cad):
    tg, nn = subgrid(dt, cad); y = np.clip((T14 / 2 - np.abs(tg)) / T12, 0, 1); return y.reshape(-1, nn).mean(1)


def timing(t, r, cad, P, T0, T14, T12):
    n = np.round((t - T0) / P); k = np.abs(t - (T0 + n * P)) < 0.4; t, r, n = t[k], r[k], n[k]; sd = robust_sd(r)
    dts = np.arange(-0.06, 0.06001, 0.0005); chi = []
    for d in dts:
        s = trap_int(t - (T0 + n * P + d), T14, T12, cad); A = s @ s; dep = -(s @ r) / A if A > 0 else 0; chi.append(np.sum((r + dep * s) ** 2) / sd ** 2)
    chi = np.array(chi); i = np.argmin(chi); ok = chi - chi[i] <= 1
    return float(dts[i] * 1440), float((dts[ok].max() - dts[ok].min()) / 2 * 1440), int(np.round(np.median(n)))


def det_excl(t, f, P, T0, T14, other):
    keep = mask_b(t) if other == 'b' else np.abs(phase(t, Pc_, Tc_ref)) > T14c / 2 + 1 / 24
    t, f = t[keep], f[keep]; inn = np.abs(phase(t, P, T0)) < 0.75 * T14 + 0.5 / 24
    tr = flatten(t[~inn], f[~inn], method='biweight', window_length=1.0, edge_cutoff=0, break_tolerance=0.4, return_trend=True)[1]
    ok = np.isfinite(tr); return t, f / np.interp(t, t[~inn][ok], tr[ok]) - 1


tim = []
for lab in ['TSPOC_13_pdc', 'TSPOC_13_sap', 'QLP_13_ksp', 'ELE_13_corr', 'SPOC_39_pdc', 'FAST_39_pdc', 'TSPOC_39_pdc', 'QLP_39_ksp'] + \
           [f'SPOC_{s}_pdc' for s in (100, 101, 102, 103)] + [f'QLP_{s}_ksp' for s in (100, 101, 102, 103)]:
    L = load(W, lab); k = trim_gaps(L['t']); t_, f_ = L['t'][k], L['f'][k]
    for pl, P, T0, T14, T12, other in (('c', Pc_, Tc_ref, T14c, T12c, 'b'), ('b', PB, TB, DB, 0.06 / 24, 'c')):
        a, b_ = det_excl(t_, f_, P, T0, T14, other); oc, e, nm = timing(a, b_, L['cad'], P, T0, T14, T12)
        tim.append((lab, pl, nm, f'{oc:+.1f}', f'{e:.1f}', '', '', '')); print(lab, pl, oc, e, flush=True)
for lab in ['TSPOC_13_pdc', 'QLP_13_ksp', 'SPOC_39_pdc', 'SPOC_101_pdc']:
    L = load(W, lab); k = trim_gaps(L['t']) & mask_b(L['t']); t_, f_ = L['t'][k], L['f'][k]
    a, b_ = det_excl(t_, f_, Pc_, Tc_ref, T14c, 'b'); oc, e, nm = timing(a, b_, L['cad'], Pc_, Tc_ref, T14c, T12c)
    away = np.abs(phase(t_, Pc_, Tc_ref)) > T14c; rng = np.random.default_rng(7); rec = []
    for i in range(150):
        T0i = Tc_ref + rng.uniform(0.2, 0.8) * Pc_; fi = f_.copy(); fi[~away] = np.nan
        fi = fi * (1 - 550e-6 * trap_int(phase(t_, Pc_, T0i), T14c, T12c, L['cad'])); g = np.isfinite(fi)
        inn = np.abs(phase(t_[g], Pc_, T0i)) < 0.75 * T14c + 0.5 / 24
        tr = flatten(t_[g][~inn], fi[g][~inn], method='biweight', window_length=1.0, edge_cutoff=0, break_tolerance=0.4, return_trend=True)[1]
        ok = np.isfinite(tr); rr_ = fi[g] / np.interp(t_[g], t_[g][~inn][ok], tr[ok]) - 1
        rec.append(timing(t_[g], rr_, L['cad'], Pc_, T0i, T14c, T12c)[0])
    rec = np.array(rec)
    tim.append((lab, 'c (injection test)', nm, f'{oc:+.1f}', f'{e:.1f}', f'{np.mean(rec):+.1f}', f'{np.std(rec):.1f}', f'{oc / np.std(rec):+.2f}'))
    print(tim[-1], flush=True)
write_csv(os.path.join(DATA, 'timing_295541511.csv'),
          [f'TOI-1117 (TIC 295541511): per-product mid-transit offsets of c and b from their linear ephemerides with an exposure-integrated trapezoid; scripts/toi1117/checks_toi1117.py, {DATE}',
           f'c: adopted ephemeris (ephemeris_295541511.csv, Tref = BTJD {Tc_ref:.5f}, P = {Pc_:.7f} d), trapezoid T14 = 1.85 h, ingress 0.25 h; b: P = {PB:.7f} d, T0 = BTJD {TB:.5f}, T14 = {DB * 24:.3f} h, ingress 0.06 h',
           'each light curve: 0.25 d trimmed at gaps, the other planet masked, biweight 1.0 d with the planet\'s own transits masked; one offset for all transits of the light curve (0.72-min grid, chi2 + 1 half-width = err_min, white errors)',
           'injection rows: 150 trapezoids of 550 ppm injected at random phases (0.2-0.8 P from c) into the same light curve with c removed; inj_mean_min, inj_sd_min = mean and standard deviation of the recovered offsets; z = o_c_min / inj_sd_min',
           'epoch = median epoch number of the transits in the light curve (c: from Tref; b: from its T0); the trapezoid is averaged over each exposure (midpoint sampling)'],
          ['product', 'planet', 'epoch', 'o_c_min', 'err_min', 'inj_mean_min', 'inj_sd_min', 'z'], tim)

# 6. TCEs in the data-validation reports
import xml.etree.ElementTree as ET
for fp in sorted(glob.glob(os.path.join(dldir(W, 'dv'), '*dvr.xml'))):
    root = ET.parse(fp).getroot(); per = []
    for el in root.iter():
        if el.tag.endswith('}planetResults'):
            for sub in el.iter():
                if sub.tag.endswith('}modelParameter') and sub.get('name') == 'orbitalPeriodDays':
                    per.append(float(sub.get('value'))); break
    add('data-validation report TCEs', os.path.basename(fp), str(len(per)), '', 'periods (d): ' + ', '.join(f'{x:.4f}' for x in per))

# 7. sectors in which tess-point places the star on silicon
import tess_stars2px as tp
r = tp.tess_stars2px_function_entry(TIC, RA, DEC)
for s, cam, ccd, x, y in zip(r[3], r[4], r[5], r[6], r[7]):
    add('tess-point: sectors with the star on silicon', f'sector {s}', f'camera {cam}, CCD {ccd}', '', f'column {x:.1f}, row {y:.1f}')

write_csv(os.path.join(DATA, 'checks_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): summary of the checks; scripts/toi1117/checks_toi1117.py, {DATE}',
           f'box measurements at P = {PC} d, T0 = BTJD {TBOX:.4f}, D = {DBOX * 24:.3f} h on TESS-SPOC sector 13 + SPOC 2-min sectors 39, 100-103 PDCSAP (c masked during detrending), as in depths_295541511.csv',
           'value and err in the units of the variant (ppm for depths); see falsealarm_295541511.csv, epochs_295541511.csv and timing_295541511.csv for the details'],
          ['check', 'variant', 'value', 'err', 'note'], summ)
