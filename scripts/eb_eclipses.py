"""TIC 11613676 (CD-30 448): photometry, per-eclipse fits, ephemeris, depth history, difference images, predictions.

usage: python scripts/eb_eclipses.py [--work ./work] [--out data] [--start 2026-10-19] [--months 12]
input:  <work>/11613676/cut/s<sector>.npz (other_download.py) and the SPOC 2-min light curves of sectors 3 and 30 in
        <work>/11613676/dl/mastDownload/TESS/
output: data/eclipse_times_11613676.csv, depths_11613676.csv, sectors_11613676.csv, diffimage_11613676.csv,
        checks_11613676.csv, profiles_11613676.csv, predicted_eclipses_11613676.csv (formats in the file headers);
        <work>/11613676/lc_s<sector>.npz (aperture light curves) and events.json.
Photometry: circular apertures r = 1.5 and 2.5 px around the proper-motion-corrected Gaia DR3 position, background =
a plane fitted per cadence to the faint pixels (below the 40th percentile of the median image) more than 5 px from the
target; quality-0 cadences. Contamination: Gaussian PSF with the sigma fitted to the median image of each sector,
the six nearest Gaia DR3 sources of targets.py (30-97"). Eclipse fits: trapezoid (mid-time, depth, T14, T12) times a
linear baseline over +-0.4 d, exposure-time supersampling; errors = max(residual-scaled covariance, cyclic-shift
bootstrap of the residuals, 200 draws). Eclipses with incomplete ingress or egress coverage are skipped. Timing and depth
set: SPOC 2-min SAP flux for sectors 3 and 30, cutout r = 2.5 px otherwise; r = 1.5 px and PDCSAP are checks.
Sector-3 limits: the sector-30 eclipse shape (SPOC SAP) with the mid-time fixed at the ephemeris prediction.
"""
import argparse, csv, glob, json, os
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.io.fits import Header
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from targets import OTHER

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
p.add_argument('--start', default='2026-10-19'); p.add_argument('--months', type=float, default=12)
args = p.parse_args(); TIC = 11613676; c = OTHER[TIC]; wd = os.path.join(args.work, str(TIC)); os.makedirs(args.out, exist_ok=True)
T0G, PG = c['eph0']; PH2 = c['ph2']; T14G = 3.0 / 24; N = 21; yy, xx = np.indices((N, N)); rng = np.random.default_rng(1)
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
hdr = lambda lines: ''.join('# ' + l + '\n' for l in lines)
year = lambda t: 2000 + (t + 2457000 - 2451544.5) / 365.25
tmag = lambda g, bprp: g - 0.00522555 * bprp**3 + 0.0891337 * bprp**2 - 0.633923 * bprp + 0.0324473   # Stassun et al. 2019


# ---------------------------------------------------------------- photometry
def photometry(sec):
    z = np.load(os.path.join(wd, 'cut', f's{sec}.npz')); tm = z['tm']; f = z['f'].astype(float); q = z['q']
    ok = np.isfinite(tm) & np.all(np.isfinite(f), axis=(1, 2)); tm, f, q = tm[ok], f[ok], q[ok]
    tx, ty = float(z['tx']), float(z['ty']); w = WCS(Header.fromstring(str(z['wcs'])))
    O0 = np.median(f[q == 0], 0); far = (np.hypot(xx - tx, yy - ty) > 5) & (O0 < np.percentile(O0, 40))
    A = np.c_[np.ones(far.sum()), xx[far], yy[far]]; Af = np.c_[np.ones(N * N), xx.ravel(), yy.ravel()]
    g = f - (Af @ np.linalg.lstsq(A, f[:, far].T, rcond=None)[0]).T.reshape(-1, N, N)
    O = np.median(g[q == 0], 0); m5 = (abs(xx - tx) <= 3) & (abs(yy - ty) <= 3)
    gm = lambda p_: p_[3] * np.exp(-((xx - p_[0])**2 + (yy - p_[1])**2) / (2 * p_[2]**2)) + p_[4]
    x0, y0, sig = least_squares(lambda p_: (gm(p_) - O)[m5], [tx, ty, 0.7, O.max(), 0]).x[:3]
    nb = [(*[float(v) for v in w.world_to_pixel(SkyCoord(ra, de, unit='deg'))], tmag(G, br)) for ra, de, G, br in c['neighbours']]
    tt = tmag(11.041, 1.025); row = dict(sector=sec, camera=int(z['cam']), ccd=int(z['ccd']), cadence_min=float(np.median(np.diff(tm)) * 1440),
                                         start=float(tm.min()), end=float(tm.max()), n_q0=int((q == 0).sum()), psf_sigma=float(sig), psf_offset=float(np.hypot(x0 - tx, y0 - ty)))
    lc = dict(t=tm, q=q)
    for rad in (1.5, 2.5):
        msk = np.hypot(xx - tx, yy - ty) <= rad; lc[f'r{rad}'] = g[:, msk].sum(1)
        fin = lambda px, py: np.exp(-((xx[msk] - px)**2 + (yy[msk] - py)**2) / (2 * sig**2)).sum() / (2 * np.pi * sig**2)
        ft = fin(tx, ty); row[f'npix_r{rad}'] = int(msk.sum()); row[f'target_frac_r{rad}'] = float(ft)
        row[f'contamination_r{rad}'] = float(sum(fin(px, py) * 10**(-0.4 * (tn - tt)) for px, py, tn in nb) / ft)
    np.savez(os.path.join(wd, f'lc_s{sec}.npz'), **lc)
    return row, dict(t=tm, q=q, g=g, O=O, tx=tx, ty=ty, sig=sig, nb=nb)


def spoc(sec):
    fn = glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{sec:04d}-{TIC:016d}*', '*-s_lc.fits'))[0]
    d = fits.open(fn)[1].data
    return dict(t=np.array(d['TIME'], float), q=np.array(d['QUALITY']), sap=np.array(d['SAP_FLUX'], float), pdcsap=np.array(d['PDCSAP_FLUX'], float))


# ---------------------------------------------------------------- per-eclipse trapezoid fits
def trap(p_, x):
    t0, dep, T14, T12, a, b = p_; a_ = np.abs(x - t0); h = T14 / 2; r = max(T12, 1e-4)
    return (a + b * (x - t0)) * (1 - np.where(a_ < h - r, dep, np.where(a_ < h, dep * (h - a_) / r, 0)))


def fit_event(t, f, tp, texp):
    nsub = max(1, int(round(texp / 0.002))); off = ((np.arange(nsub) + 0.5) / nsub - 0.5) * texp
    model = lambda p_, x: trap(p_, (x[:, None] + off[None, :]).ravel()).reshape(-1, nsub).mean(1)
    lo = [tp - 0.8 * T14G, -0.02, 0.3 * T14G, 0.004, 0.5, -1]; hi = [tp + 0.8 * T14G, 0.5, 2.5 * T14G, 1.2 * T14G, 1.5, 1]
    best = None
    for d0 in (0.005, 0.03):
        for T in (0.7 * T14G, T14G, 1.3 * T14G):
            r = least_squares(lambda p_: model(p_, t) - f, [tp, d0, T, 0.25 * T, 1, 0], bounds=(lo, hi))
            if best is None or r.cost < best.cost: best = r
    res = best.fun; s2 = (res**2).sum() / max(len(t) - 6, 1)
    try: e = np.sqrt(np.diag(np.linalg.inv(best.jac.T @ best.jac) * s2))
    except np.linalg.LinAlgError: e = np.full(6, np.nan)
    mfit = model(best.x, t)
    bs = []
    for _ in range(200):
        fb = mfit + np.roll(res, rng.integers(1, len(t))); bs.append(least_squares(lambda p_: model(p_, t) - fb, best.x, bounds=(lo, hi)).x)
    bs = np.array(bs)
    return best.x, np.fmax(e, bs.std(0)), np.sqrt(s2)


def events(t, f, q, src):
    g = (q == 0) & np.isfinite(f) & np.isfinite(t); tt, ff = t[g], f[g]; texp = float(np.median(np.diff(tt))); out = []
    for kind, ph in (('pri', 0.0), ('sec', PH2)):
        for n in range(int(np.floor((tt.min() - T0G) / PG)) - 1, int(np.ceil((tt.max() - T0G) / PG)) + 2):
            tp = T0G + (n + ph) * PG; m = np.abs(tt - tp) < 0.4
            if m.sum() < 15 or np.sum(np.abs(tt - tp) < T14G / 2) < 0.5 * T14G / texp: continue
            need = 0.3 * (0.75 * T14G / texp)
            if np.sum((tt > tp - T14G) & (tt < tp - T14G / 4)) <= need or np.sum((tt < tp + T14G) & (tt > tp + T14G / 4)) <= need: continue
            partial = not (tt[m].min() < tp - 0.25 and tt[m].max() > tp + 0.25)
            x, e, rms = fit_event(tt[m], ff[m] / np.median(ff[m][np.abs(tt[m] - tp) > T14G]), tp, texp)
            out.append(dict(src=src, kind=kind, n=n, t_pred=tp, t0=x[0], t0_err=e[0], depth=x[1], depth_err=e[1], T14_h=x[2] * 24, T14_err_h=e[2] * 24,
                            T12_h=x[3] * 24, T12_err_h=e[3] * 24, rms=rms, partial=partial))
            print(f"{src:14s} {kind} n={n:5d} t0={x[0]:.5f}+-{e[0]:.5f} depth={x[1]*100:.3f}+-{e[1]*100:.3f}% T14={x[2]*24:.2f}h{' partial' if partial else ''}", flush=True)
    return out


def wmean(v, e):
    w = 1 / np.asarray(e)**2; return float(np.sum(w * v) / w.sum()), float(1 / np.sqrt(w.sum()))


def wlin(x, y, e):
    w = 1 / e**2; A = np.c_[np.ones_like(x), x]; cov = np.linalg.inv(A.T @ (A * w[:, None])); return cov @ (A.T @ (w * y)), np.sqrt(np.diag(cov))


# ---------------------------------------------------------------- run
SECS = c['sectors']; rows_sec = []; cubes = {}; EV = []; LCS = {}
for s in SECS:
    r, cube = photometry(s); rows_sec.append(r); cubes[s] = cube; L = np.load(os.path.join(wd, f'lc_s{s}.npz'))
    LCS[(s, 'r2.5')] = (L['t'], L['r2.5'], L['q']); LCS[(s, 'r1.5')] = (L['t'], L['r1.5'], L['q'])
for s in c['spoc_2min']:
    d = spoc(s); LCS[(s, 'sap')] = (d['t'], d['sap'], d['q']); LCS[(s, 'pdcsap')] = (d['t'], d['pdcsap'], d['q'])
for (s, col), (t, f, q) in LCS.items():
    for e in events(t, f, q, f'{col}'): e['sector'] = s; EV.append(e)
json.dump(EV, open(os.path.join(wd, 'events.json'), 'w'), indent=1, default=float)
TIMING = lambda e: (e['src'] == 'sap') if e['sector'] in c['spoc_2min'] else (e['src'] == 'r2.5')
S = [e for e in EV if TIMING(e) and not e['partial']]
det = [e for e in S if e['depth'] / e['depth_err'] > 5 and e['depth'] > 0.001]
pri = [e for e in det if e['kind'] == 'pri']; sec = [e for e in det if e['kind'] == 'sec']
cyc = np.array([round((e['t0'] - T0G) / PG) for e in pri]); tp = np.array([e['t0'] for e in pri]); ep = np.array([e['t0_err'] for e in pri])
(T0, P), (eT0, eP) = wlin(cyc, tp, ep); oc = (tp - T0 - P * cyc) * 1440; chi2 = float(np.sum((oc / (ep * 1440))**2))
print(f'primary ephemeris: T0 = {T0:.5f} +- {eT0:.5f}, P = {P:.6f} +- {eP:.6f} d, chi2 = {chi2:.0f} / {len(pri) - 2}')
for e in EV:
    e['cycle'] = int(round((e['t0'] - T0 - (0 if e['kind'] == 'pri' else PH2 * P)) / P)); e['phase'] = ((e['t0'] - T0) / P) % 1
    e['oc_min'] = (e['t0'] - T0 - P * e['cycle']) * 1440 if e['kind'] == 'pri' else np.nan
cs = np.array([e['cycle'] for e in sec]); ts = np.array([e['t0'] for e in sec]); es = np.array([e['t0_err'] for e in sec])
(T0s, Ps), (eT0s, ePs) = wlin(cs, ts, es); ys = year(ts) - 2023.0
(ph0, dph), (eph0, edph) = wlin(ys, np.array([e['phase'] for e in sec]), es / P)
print(f'secondary: P_sec - P = {Ps - P:+.6f} +- {np.hypot(ePs, eP):.6f} d; d(phase)/dt = {dph:+.5f} +- {edph:.5f} /yr')

# ---- tables: eclipse times
with open(os.path.join(args.out, 'eclipse_times_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448 = Gaia DR3 {c["gaia_dr3"]}: per-eclipse trapezoid fits (scripts/eb_eclipses.py), {DATE}',
                  f'Primary linear ephemeris (complete primaries, depth > 5 sigma, timing set): T0 = BTJD {T0:.5f} +- {eT0:.5f}, P = {P:.6f} +- {eP:.6f} d, chi2 = {chi2:.0f} for {len(pri) - 2} dof',
                  'kind = pri / sec; cycle = cycle number relative to T0; t0_btjd = BJD_TDB - 2457000; errors = max(covariance, cyclic-shift bootstrap); oc_min = primary O-C (secondaries: empty)',
                  'phase = (t0 - T0)/P mod 1; depth = trapezoid depth (fraction of the flux in the aperture); T14_h, T12_h = total and ingress duration (h)',
                  'source: sap / pdcsap = SPOC 2-min SAP_FLUX / PDCSAP_FLUX (sectors 3, 30); r2.5 / r1.5 = TESScut circular aperture (px); timing_set = 1 for the set used for the ephemeris and the depth history',
                  'partial = 1 if the +-0.25 d window around the eclipse is not fully covered (not used)']))
    w = csv.writer(fo); w.writerow(['sector', 'kind', 'cycle', 'source', 't0_btjd', 't0_err_d', 't0_err_min', 'oc_min', 'phase', 'depth', 'depth_err', 'T14_h', 'T14_err_h', 'T12_h', 'T12_err_h', 'rms', 'partial', 'timing_set'])
    for e in sorted(EV, key=lambda e: (e['sector'], e['src'], e['t0'])):
        w.writerow([e['sector'], e['kind'], e['cycle'], e['src'], f"{e['t0']:.6f}", f"{e['t0_err']:.6f}", f"{e['t0_err'] * 1440:.2f}", '' if np.isnan(e['oc_min']) else f"{e['oc_min']:+.2f}",
                    f"{e['phase']:.5f}", f"{e['depth']:.5f}", f"{e['depth_err']:.5f}", f"{e['T14_h']:.3f}", f"{e['T14_err_h']:.3f}", f"{e['T12_h']:.3f}", f"{e['T12_err_h']:.3f}", f"{e['rms']:.5f}",
                    int(e['partial']), int(TIMING(e) and not e['partial'] and e['depth'] / e['depth_err'] > 5 and e['depth'] > 0.001)])

# ---- sector-3 limits with the sector-30 shape
def trapz(x, t0, T14, T12):
    a_ = np.abs(x - t0); h = T14 / 2; return np.where(a_ < h - T12, 1.0, np.where(a_ < h, (h - a_) / T12, 0))
lim = {}; PH30 = np.mean([e['phase'] for e in sec if e['sector'] == 30])
for kind, ph in (('pri', 0.0), ('sec', PH30)):
    g30 = [e for e in S if e['sector'] == 30 and e['kind'] == kind]
    T14 = wmean([e['T14_h'] for e in g30], [e['T14_err_h'] for e in g30])[0] / 24; T12 = wmean([e['T12_h'] for e in g30], [e['T12_err_h'] for e in g30])[0] / 24
    for col in ('sap', 'pdcsap', 'r2.5', 'r1.5'):
        t, f, q = LCS[(3, col)]; gm = (q == 0) & np.isfinite(f); t, f = t[gm], f[gm]; vals = []
        for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2):
            tq = T0 + (n + ph) * P; m = np.abs(t - tq) < 0.4
            if np.sum(np.abs(t[m] - tq) < T14 / 2) < 0.5 * T14 / np.median(np.diff(t)): continue
            x = t[m]; y = f[m] / np.median(f[m]); A = np.c_[np.ones_like(x), x - tq, -trapz(x, tq, T14, T12)]
            cc = np.linalg.lstsq(A, y, rcond=None)[0]; res = y - A @ cc; cov = np.linalg.inv(A.T @ A) * (res**2).sum() / (len(x) - 3)
            vals.append((cc[2] / cc[0], np.sqrt(cov[2, 2]) / cc[0]))
        v = np.array(vals); lim[(kind, col)] = (*wmean(v[:, 0], v[:, 1]), len(v), T14 * 24, T12 * 24)
        print(f'S3 {kind} {col}: depth {lim[(kind, col)][0]*100:+.3f} +- {lim[(kind, col)][1]*100:.3f} % ({len(v)} eclipses), 3-sigma limit {3*lim[(kind, col)][1]*100:.3f} %')

# ---- depth history
with open(os.path.join(args.out, 'depths_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: eclipse depth and duration per sector and eclipse (scripts/eb_eclipses.py), {DATE}',
                  'depth_pct = weighted mean trapezoid depth (% of the aperture flux) of the complete eclipses; for sector 3 the depth with the sector-30 shape at the predicted time (no eclipse detected) and limit3_pct = 3-sigma upper limit',
                  'T14_h, T12_h = weighted mean total and ingress durations; T23_h = T14 - 2 T12 (mean); phase_mean = mean phase (secondary) relative to the primary ephemeris of eclipse_times_11613676.csv',
                  'source = timing set (SPOC 2-min SAP in sectors 3, 30; TESScut r = 2.5 px otherwise); depth_r1.5_pct, depth_pdcsap_pct = the same with the TESScut r = 1.5 px aperture / SPOC PDCSAP flux (checks)']))
    w = csv.writer(fo); w.writerow(['sector', 'year', 'kind', 'n', 'source', 'depth_pct', 'depth_err_pct', 'limit3_pct', 'T14_h', 'T14_err_h', 'T12_h', 'T23_h', 'phase_mean', 'depth_r1.5_pct', 'depth_pdcsap_pct'])
    for s in SECS:
        for kind in ('pri', 'sec'):
            g = [e for e in S if e['sector'] == s and e['kind'] == kind]
            alt = lambda col: [e for e in EV if e['sector'] == s and e['kind'] == kind and e['src'] == col and not e['partial']]
            fmt = lambda gg: f"{wmean([e['depth'] for e in gg], [e['depth_err'] for e in gg])[0] * 100:.3f}" if gg else ''
            src = 'SPOC 2-min SAP' if s in c['spoc_2min'] else 'TESScut r=2.5 px'
            if s == 3:
                d, e_, n, T14, T12 = lim[(kind, 'sap')]
                w.writerow([s, f'{year(np.mean(LCS[(3, "sap")][0][np.isfinite(LCS[(3, "sap")][0])])):.2f}', kind, n, src + ' (sector-30 shape fixed)', f'{d*100:.3f}', f'{e_*100:.3f}', f'{3*e_*100:.3f}',
                            f'{T14:.2f}', '', f'{T12:.2f}', '', '', f"{lim[(kind, 'r1.5')][0]*100:.3f}", f"{lim[(kind, 'pdcsap')][0]*100:.3f}"])
                continue
            dm, de = wmean([e['depth'] for e in g], [e['depth_err'] for e in g]); T14, eT14 = wmean([e['T14_h'] for e in g], [e['T14_err_h'] for e in g])
            w.writerow([s, f"{year(np.mean([e['t0'] for e in g])):.2f}", kind, len(g), src, f'{dm*100:.3f}', f'{de*100:.3f}', '', f'{T14:.2f}', f'{eT14:.2f}',
                        f"{np.mean([e['T12_h'] for e in g]):.2f}", f"{np.mean([e['T14_h'] - 2 * e['T12_h'] for e in g]):.2f}",
                        f"{np.mean([e['phase'] for e in g]):.4f}" if kind == 'sec' else '', fmt(alt('r1.5')), fmt(alt('pdcsap'))])

# ---- sectors
with open(os.path.join(args.out, 'sectors_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: TESS sectors (21x21 TESScut cutouts) and aperture properties (scripts/eb_eclipses.py), {DATE}',
                  'span = first/last cadence (BTJD); psf_sigma = Gaussian sigma (px) fitted to the median image; npix, target_frac = pixels in and fraction of the target PSF inside the r = 1.5 / 2.5 px apertures',
                  'contamination = flux of the six nearest Gaia DR3 sources (30-97") in the aperture relative to the target flux in the aperture (Gaussian PSF, T from G and BP-RP)']))
    w = csv.writer(fo); keys = ['sector', 'camera', 'ccd', 'cadence_min', 'start', 'end', 'n_q0', 'psf_sigma', 'psf_offset', 'npix_r1.5', 'target_frac_r1.5', 'contamination_r1.5', 'npix_r2.5', 'target_frac_r2.5', 'contamination_r2.5']
    w.writerow([k.replace('start', 'span_start_btjd').replace('end', 'span_end_btjd').replace('psf_offset', 'psf_offset_px') for k in keys])
    for r in rows_sec:
        w.writerow([f'{r[k]:.4f}' if isinstance(r[k], float) and 'contamination' in k else (f'{r[k]:.3f}' if isinstance(r[k], float) else r[k]) for k in keys])

# ---- difference images (stacked per sector at primary and secondary phase)
rows_d = []
for s in SECS:
    C = cubes[s]; t, q, g = C['t'], C['q'], C['g']; sig = C['sig']
    for kind in ('pri', 'sec'):
        ph = 0.0 if kind == 'pri' else (np.mean([e['phase'] for e in sec if e['sector'] == s]) if any(e['sector'] == s for e in sec) else PH2)
        D = np.zeros((N, N)); E = np.zeros((N, N)); nev = 0
        for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2):
            dt = t - (T0 + (n + ph) * P); win = (np.abs(dt) < 0.09) & (q == 0); oot = (np.abs(dt) > 0.15) & (np.abs(dt) < 0.6) & (q == 0)
            if win.sum() < 5 or oot.sum() < 20: continue
            nev += 1; X = np.c_[np.ones(oot.sum()), dt[oot], dt[oot]**2]; cf = np.linalg.lstsq(X, g[oot].reshape(oot.sum(), -1), rcond=None)[0]
            res = g.reshape(len(t), -1) - np.c_[np.ones(len(t)), dt, dt**2] @ cf
            D += -res[win].mean(0).reshape(N, N); E += (res[oot].std(0) / np.sqrt(win.sum())).reshape(N, N)**2
        if nev == 0: continue
        D /= nev; E = np.sqrt(E) / nev
        dres = lambda p_, pos=None: ((p_[0] * np.exp(-((xx - (pos or p_[1:])[0])**2 + (yy - (pos or p_[1:])[1])**2) / (2 * sig**2)) - D) / E).ravel()
        rf = least_squares(lambda p_: dres(p_), [D.max(), C['tx'], C['ty']]); rt = least_squares(lambda p_: dres(p_, (C['tx'], C['ty'])), [D.max()])
        nbc = [((least_squares(lambda p_: dres(p_, (px, py)), [D.max()]).fun**2).sum(), tn) for px, py, tn in C['nb'] if 0 <= px < N and 0 <= py < N]
        rows_d.append([s, kind, nev, f'{(D / E).max():.1f}', f'{np.hypot(rf.x[1] - C["tx"], rf.x[2] - C["ty"]):.2f}', f'{((D / E)**2).sum():.1f}', f'{(rt.fun**2).sum():.1f}',
                       f'{min(nbc)[0]:.1f}' if nbc else '', f'{min(nbc)[1]:.2f}' if nbc else ''])
        print('diff', rows_d[-1], flush=True)
with open(os.path.join(args.out, 'diffimage_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: stacked difference images (out-of-eclipse quadratic model per pixel minus in-eclipse mean, +-0.09 d) per sector and eclipse (scripts/eb_eclipses.py), {DATE}',
                  'peak_snr = largest pixel S/N; free_offset_px = position of a free single-PSF fit minus the target position; chi2_null / chi2_target / chi2_best_neighbour = chi2 of no signal / a PSF at the target / a PSF at the best Gaia neighbour on the cutout (its T magnitude)']))
    w = csv.writer(fo); w.writerow(['sector', 'kind', 'n_eclipses', 'peak_snr', 'free_offset_px', 'chi2_null', 'chi2_target', 'chi2_best_neighbour', 'best_neighbour_tmag']); w.writerows(rows_d)

# ---- checks: odd/even, phase 0.5, secondary phase
rows_c = []
for kind in ('pri', 'sec'):
    for s in SECS:
        g = [e for e in det if e['sector'] == s and e['kind'] == kind]; o = [e for e in g if e['cycle'] % 2]; ev = [e for e in g if not e['cycle'] % 2]
        if o and ev:
            (do, eo), (de, ee) = wmean([e['depth'] for e in o], [e['depth_err'] for e in o]), wmean([e['depth'] for e in ev], [e['depth_err'] for e in ev])
            rows_c.append([s, f'odd/even {kind}', f'{do*100:.3f} / {de*100:.3f}', f'{(do - de) / np.hypot(eo, ee):+.1f} sigma', f'{len(o)} + {len(ev)} eclipses'])
for s in SECS:
    col = 'sap' if s in c['spoc_2min'] else 'r2.5'; t, f, q = LCS[(s, col)]; gm = (q == 0) & np.isfinite(f); t, f = t[gm], f[gm]; v = []
    for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2):
        tq = T0 + (n + 0.5) * P; m = np.abs(t - tq) < 0.3
        if m.sum() < 20: continue
        x, y = t[m], f[m]; oot = np.abs(x - tq) > 0.1; win = np.abs(x - tq) < 0.0625
        if win.sum() < 3: continue
        r = y / np.polyval(np.polyfit(x[oot] - tq, y[oot], 1), x - tq); v.append((1 - r[win].mean(), r[oot].std() / np.sqrt(win.sum())))
    if v:
        v = np.array(v); d, e_ = wmean(v[:, 0], v[:, 1]); rows_c.append([s, 'depth at phase 0.5 (3-h box)', f'{d*100:+.3f} +- {e_*100:.3f} %', '', f'{len(v)} windows'])
for s in SECS:
    g = [e for e in sec if e['sector'] == s]
    if g: rows_c.append([s, 'secondary phase', f"{np.mean([e['phase'] for e in g]):.4f}", f"scatter {np.std([e['phase'] for e in g]):.4f}", f'{len(g)} eclipses'])
rows_c.append(['all', 'secondary period minus primary period (d)', f'{Ps - P:+.6f}', f'+- {np.hypot(ePs, eP):.6f}', f'{len(sec)} secondaries'])
rows_c.append(['all', 'd(secondary phase)/dt (per yr)', f'{dph:+.5f}', f'+- {edph:.5f}', ''])
with open(os.path.join(args.out, 'checks_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: period and eclipse checks from the timing set (scripts/eb_eclipses.py), {DATE}']))
    w = csv.writer(fo); w.writerow(['sector', 'check', 'value', 'significance_or_error', 'note']); w.writerows(rows_c)

# ---- per-sector eclipse profiles (for the figures): each eclipse divided by its fitted linear baseline, 10-min bins
with open(os.path.join(args.out, 'profiles_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: per-sector eclipse profiles, timing set (scripts/eb_eclipses.py), {DATE}',
                  'dt_h = hours from the fitted mid-eclipse time of each eclipse (sector 3: from the predicted time); flux = mean of the eclipses of the sector, each divided by a linear out-of-eclipse baseline; 10-min bins; err = standard error of the bin']))
    w = csv.writer(fo); w.writerow(['sector', 'kind', 'dt_h', 'flux', 'err', 'n'])
    for s in SECS:
        col = 'sap' if s in c['spoc_2min'] else 'r2.5'; t, f, q = LCS[(s, col)]; gm = (q == 0) & np.isfinite(f); t, f = t[gm], f[gm]
        for kind in ('pri', 'sec'):
            cent = [e['t0'] for e in S if e['sector'] == s and e['kind'] == kind] if s != 3 else \
                   [T0 + (n + (0 if kind == 'pri' else PH30)) * P for n in range(int((t.min() - T0) / P) - 1, int((t.max() - T0) / P) + 2)]
            X, Y = [], []
            for t0 in cent:
                m = np.abs(t - t0) < 0.3
                if m.sum() < 20: continue
                x, y = t[m] - t0, f[m]; oot = np.abs(x) > 0.11
                if oot.sum() < 10: continue
                X.append(x); Y.append(y / np.polyval(np.polyfit(x[oot], y[oot], 1), x))
            if not X: continue
            X, Y = np.concatenate(X) * 24, np.concatenate(Y); b = np.round(X * 6)
            for k in np.unique(b):
                m = b == k
                if m.sum() >= 2: w.writerow([s, kind, f'{X[m].mean():.3f}', f'{Y[m].mean():.5f}', f'{Y[m].std() / np.sqrt(m.sum()):.5f}', int(m.sum())])

# ---- predicted eclipses (geocentric UTC as in predict_transits.py)
COORD = SkyCoord(c['ra'] * u.deg, c['dec'] * u.deg); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
def utc(bjd):
    t = Time(bjd, format='jd', scale='tdb', location=GEO)
    for _ in range(2): t = Time(bjd, format='jd', scale='tdb', location=GEO) - Time(t, location=GEO).light_travel_time(COORD, kind='barycentric')
    return t.utc.iso[:16]
phs = np.mean([e['phase'] for e in sec if e['sector'] == max(e['sector'] for e in sec)])
t1 = Time(args.start).jd - 2457000; t2 = t1 + args.months * 30.5; n0 = int(np.ceil((t1 - T0) / P))
with open(os.path.join(args.out, 'predicted_eclipses_11613676.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 11613676 = CD-30 448: predicted eclipses for the {args.months:g} months from {args.start} (scripts/eb_eclipses.py), {DATE}',
                  f'primary: T0 = BTJD {T0:.5f} + {P:.6f} E (linear ephemeris of eclipse_times_11613676.csv); the observed O-C spans {oc.min():+.0f} to {oc.max():+.0f} min, so a window of at least +-30 min should be used',
                  f'secondary: phase {phs:.4f} (sector {max(e["sector"] for e in sec)}; the phase drifts by about {dph:+.4f} per yr); *_utc = geocentric UTC (barycentric light-travel time removed)',
                  'the eclipse depths are changing (README); last measured T14 about 4.0 h (primary) and 3.8 h (secondary)']))
    w = csv.writer(fo); w.writerow(['cycle', 'kind', 't_btjd', 't_bjd_tdb', 't_utc'])
    for n in range(n0, n0 + int((t2 - t1) / P) + 2):
        for kind, ph in (('pri', 0.0), ('sec', phs)):
            tq = T0 + (n + ph) * P
            if t1 <= tq <= t2: w.writerow([n, kind, f'{tq:.4f}', f'{tq + 2457000:.4f}', utc(tq + 2457000)])
print('tables written to', args.out)
