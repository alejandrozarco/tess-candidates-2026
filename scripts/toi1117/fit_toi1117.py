"""TOI-1117 c (TIC 295541511), README section 6.10: box depths per product and sector, per-transit depths, detrending tests,
and the detrended light curves for the transit fits.

usage: python scripts/toi1117/fit_toi1117.py [--work ./work]
input:  <work>/295541511/dl/lc/ (download_toi1117.py)
output: data/depths_295541511.csv      box depth and S/N per product and sector, c masked and not masked during detrending
        data/pertransit_295541511.csv  depth of every transit (local linear baseline), SPOC 2-min PDCSAP and TESS-SPOC sector 13
        data/depth_tests_295541511.csv duration scan, per-sector and leave-one-sector-out S/N, detrending-window and odd/even tests
        <work>/295541511/lc_fit_{pdc,sap}.npz  detrended windows (+-0.45 d) for mcmc_toi1117.py and triceratops_toi1117.py
Box measurements use the starting ephemeris of the search (common.PC, common.TBOX) and a box of D = 1.625 h (best of the
duration scan); the white error is the robust out-of-box scatter / sqrt(n_in); the empirical error is the standard deviation of
the box depths at 300 random phase shifts of the same light curve (|shift| > 2 D).
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); args = p.parse_args()
W = args.work; DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
NAME = {'SPOC': 'SPOC 2-min', 'FAST': 'SPOC 20-s', 'TSPOC': 'TESS-SPOC FFI', 'QLP': 'QLP FFI', 'ELE': 'eleanor-lite FFI', 'TARS': 'TARS FFI'}

# 1. duration and mid-time scan on the SPOC 2-min PDCSAP light curves
S2 = [prep(W, f'SPOC_{s}_pdc', T0=TC, D=1.6 / 24) for s in SPOC2]
t = np.concatenate([x['t'] for x in S2]); r = np.concatenate([x['r'] for x in S2])
grid = []
for Dh in np.arange(0.75, 4.01, 0.125):
    for dT in np.arange(-0.03, 0.0301, 0.0025):
        d, e, k = box(t, r, PC, TC + dT, Dh / 24); grid.append((Dh, dT, d, e, d / e))
grid = np.array(grid); best = grid[np.nanargmax(grid[:, 4])]
tests = [('duration scan (SPOC 2-min PDCSAP)', f'best D = {best[0]:.3f} h, T0 offset {best[1] * 1440:+.1f} min', best[2] * 1e6, best[3] * 1e6, np.nan, best[4], np.nan)]
for Dh in (1.0, 1.25, 1.5, 1.625, 1.75, 2.0, 2.5, 3.25):
    g = grid[np.isclose(grid[:, 0], Dh)]; i = np.nanargmax(g[:, 4])
    tests.append(('duration scan (SPOC 2-min PDCSAP)', f'D = {g[i, 0]:.3f} h (best T0 offset)', g[i, 2] * 1e6, g[i, 3] * 1e6, np.nan, g[i, 4], np.nan))
print('duration scan best', best)

# 2. every product and sector, c masked and not masked during detrending
labs = ['TSPOC_13_pdc', 'TSPOC_13_sap', 'QLP_13_ksp', 'QLP_13_sap', 'ELE_13_corr', 'ELE_13_pca', 'TARS_13',
        'SPOC_39_pdc', 'SPOC_39_sap', 'FAST_39_pdc', 'FAST_39_sap', 'TSPOC_39_pdc', 'TSPOC_39_sap', 'QLP_39_ksp', 'QLP_39_sap', 'TARS_39']
for s in (100, 101, 102, 103):
    labs += [f'SPOC_{s}_pdc', f'SPOC_{s}_sap', f'QLP_{s}_ksp', f'QLP_{s}_sap']
rows = []
for lab in labs:
    prov, sec, *kind = lab.split('_'); kind = kind[0] if kind else ''
    for maskc in (1, 0):
        X = prep(W, lab, maskc=bool(maskc))
        if X is None:
            continue
        d, e, k = box(X['t'], X['r'], PC, TBOX, DBOX); ee, _ = emp_err(X['t'], X['r'], PC, TBOX, DBOX)
        rows.append((int(sec), NAME[prov], X['col'], maskc, f'{X["cad"] * 1440:.2f}', k, f'{d * 1e6:.0f}', f'{e * 1e6:.0f}', f'{ee * 1e6:.0f}',
                     f'{d / e:.2f}', f'{d / ee:.2f}', f'{robust_sd(X["r"]) * 1e6:.0f}'))
        print(lab, maskc, rows[-1][6:], flush=True)
write_csv(os.path.join(DATA, 'depths_295541511.csv'),
          [f'TOI-1117 = TIC 295541511 = Gaia DR3 {GAIA}: box depth of the 4.5787-d signal (TOI-1117 c) per TESS product and sector; scripts/toi1117/fit_toi1117.py, {DATE}',
           f'fixed ephemeris P = {PC} d, T0 = BTJD {TBOX:.4f} (starting ephemeris of the search), box duration D = {DBOX * 24:.3f} h; planet b masked (+-(T14_b/2 + 1 h)); 0.25 d trimmed at gaps > 0.15 d',
           'detrending: wotan biweight, 1.0-d window; c_masked = 1: windows of c (+-(0.75 D + 0.5 h)) excluded from the trend fit; two-sided 5-sigma clip of the residuals',
           'depth_ppm = -mean residual in the box; err_white_ppm = robust out-of-box scatter / sqrt(n_in); err_emp_ppm = std of the box depth at 300 random phase shifts (|shift| > 2 D) of the same light curve',
           'not dilution corrected beyond the corrections of each pipeline (SPOC and TESS-SPOC PDCSAP: CROWDSAP; SAP and QLP SAP: none); TARS has no quality column (all cadences used)'],
          ['sector', 'product', 'flux', 'c_masked', 'cadence_min', 'n_in', 'depth_ppm', 'err_white_ppm', 'err_emp_ppm', 'snr_white', 'snr_emp', 'scatter_ppm'], rows)

# 3. per-transit depths (local linear baseline on the detrended residuals)
S13 = prep(W, 'TSPOC_13_pdc'); S2 = [prep(W, f'SPOC_{s}_pdc') for s in SPOC2]
ep = []
for X in [S13] + S2:
    n = np.round((X['t'] - TBOX) / PC)
    for e_ in np.unique(n):
        tc = TBOX + e_ * PC; w = np.abs(X['t'] - tc) < 0.5; m = w & (np.abs(X['t'] - tc) < DBOX / 2); o = w & (np.abs(X['t'] - tc) > 0.75 * DBOX)
        if m.sum() < max(3, 0.5 * DBOX / X['cad']) or o.sum() < 10:
            continue
        A = np.vstack([np.ones(o.sum()), X['t'][o] - tc]).T; cf = np.linalg.lstsq(A, X['r'][o], rcond=None)[0]
        d = float(-(X['r'][m] - (cf[0] + cf[1] * (X['t'][m] - tc))).mean()); sd = robust_sd(X['r'][o] - A @ cf)
        ep.append((X['sec'], 'SPOC 2-min PDCSAP' if X['sec'] != 13 else 'TESS-SPOC FFI PDCSAP', int(e_), tc, d * 1e6, sd / np.sqrt(m.sum()) * 1e6, int(m.sum()), m.sum() * X['cad'] / DBOX))
dd = np.array([x[4] for x in ep]); ee = np.array([x[5] for x in ep]); w = 1 / ee ** 2; mu = np.sum(w * dd) / np.sum(w); emu = 1 / np.sqrt(w.sum())
chi2 = float(np.sum((dd - mu) ** 2 * w)); snr2 = (dd / ee) ** 2 * (dd > 0)
write_csv(os.path.join(DATA, 'pertransit_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): depth of every transit at the fixed ephemeris P = {PC} d, T0 = BTJD {TBOX:.4f}, box D = {DBOX * 24:.3f} h; scripts/toi1117/fit_toi1117.py, {DATE}',
           'detrended residuals as in depths_295541511.csv (c masked); a linear baseline fitted to the points 0.75 D to 0.5 d from mid-transit is subtracted; depth = -mean residual in the box',
           'err_ppm = robust scatter of the baseline residuals / sqrt(n_in); coverage = in-box cadences x cadence / D',
           f'weighted mean {mu:.0f} +- {emu:.0f} ppm over {len(ep)} transits, chi2 = {chi2:.1f} for {len(ep) - 1} degrees of freedom about it; negative depths: {int((dd < 0).sum())}; largest single-transit share of the summed positive S/N^2: {snr2.max() / snr2.sum():.3f}'],
          ['sector', 'product', 'epoch', 'tc_btjd', 'depth_ppm', 'err_ppm', 'n_in', 'coverage'],
          [(s, pr, e_, f'{tc:.5f}', f'{d:.0f}', f'{er:.0f}', k, f'{c:.2f}') for s, pr, e_, tc, d, er, k, c in ep])
print('per-transit', len(ep), mu, emu, chi2)
sec = np.array([x[0] for x in ep])
for s in SECTORS:
    m = sec == s; ws = 1 / ee[m] ** 2; ms = np.sum(ws * dd[m]) / ws.sum(); es = 1 / np.sqrt(ws.sum())
    tests.append(('per-sector S/N (per-transit depths, weighted mean)', f'sector {s}, {m.sum()} transits', ms, es, np.nan, ms / es, np.nan))
    m = ~m; ws = 1 / ee[m] ** 2; ms = np.sum(ws * dd[m]) / ws.sum(); es = 1 / np.sqrt(ws.sum())
    tests.append(('leave one sector out (per-transit depths, weighted mean)', f'without sector {s}', ms, es, np.nan, ms / es, np.nan))
tests.append(('all transits (per-transit depths, weighted mean)', f'{len(ep)} transits', mu, emu, np.nan, mu / emu, np.nan))

# 4. stacked box depth: detrending window, odd and even transits, first and second half of the data
labs6 = ['TSPOC_13_pdc'] + [f'SPOC_{s}_pdc' for s in SPOC2]
for win in (0.5, 0.75, 1.0, 1.5, 2.5):
    X = [prep(W, l, win=win) for l in labs6]; tt = np.concatenate([x['t'] for x in X]); rr = np.concatenate([x['r'] for x in X])
    d, e, k = box(tt, rr, PC, TBOX, DBOX); ee_, _ = emp_err(tt, rr, PC, TBOX, DBOX)
    tests.append(('stacked box depth, sectors 13 (TESS-SPOC) + 39, 100-103 (SPOC 2-min) PDCSAP', f'biweight window {win} d', d * 1e6, e * 1e6, ee_ * 1e6, d / e, d / ee_))

# 5. detrended windows for the transit fits
for kind in ('pdc', 'sap'):
    T, R, E, S, C = [], [], [], [], []
    for lab in [f'TSPOC_13_{kind}'] + [f'SPOC_{s}_{kind}' for s in SPOC2]:
        X = prep(W, lab); ph = phase(X['t'], PC, TC); w = np.abs(ph) < 0.45; sd = robust_sd(X['r'][np.abs(ph) > 0.1])
        T.append(X['t'][w]); R.append(X['r'][w]); E.append(np.full(w.sum(), sd)); S.append(np.full(w.sum(), X['sec'])); C.append(np.full(w.sum(), X['cad']))
    np.savez(os.path.join(workdir(W), f'lc_fit_{kind}.npz'), t=np.concatenate(T), f=1 + np.concatenate(R), e=np.concatenate(E), sec=np.concatenate(S), cad=np.concatenate(C))
    print(kind, 'fit windows:', sum(len(x) for x in T), 'points')

write_csv(os.path.join(DATA, 'depth_tests_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): depth and S/N tests of the 4.5787-d signal; scripts/toi1117/fit_toi1117.py, {DATE}',
           f'box measurements at P = {PC} d, T0 = BTJD {TBOX:.4f} (duration scan: T0 = BTJD {TC} plus the offset given); see depths_295541511.csv and pertransit_295541511.csv for the definitions',
           'err_emp_ppm, snr_emp: empirical error from 300 random phase shifts (empty where not computed)'],
          ['test', 'variant', 'depth_ppm', 'err_ppm', 'err_emp_ppm', 'snr', 'snr_emp'],
          [(a, b, f'{c:.0f}', f'{d:.0f}', '' if not np.isfinite(e) else f'{e:.0f}', f'{f:.2f}', '' if not np.isfinite(g) else f'{g:.2f}') for a, b, c, d, e, f, g in tests])
