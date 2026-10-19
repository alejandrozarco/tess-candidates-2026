"""TIC 124751008 (2MASS J06162084-3509125): photometry, per-eclipse fits, ephemeris, depth history, secondary-eclipse
limits and predictions (README section 6.4).

usage: python scripts/eb_124751008.py [--work ./work] [--out data] [--start 2026-10-19] [--months 12]
input:  <work>/124751008/cut/s<sector>.npz (other_download.py 124751008)
output: data/sectors_124751008.csv, eclipse_times_124751008.csv, depths_124751008.csv, predicted_eclipses_124751008.csv
        (formats in the file headers); <work>/124751008/lc_s<sector>.npz (light curves) and eclipses.json.
Photometry: the saturated G = 4.1 star kappa Col (158") bleeds along columns 1.5-5.7 px from the target; pixels of the
median image above `bleed` (targets.py) and their horizontal neighbours are masked from the aperture and the background.
Background = a plane fitted per cadence to the unmasked faint pixels (below the 40th percentile) more than 5 px from the
proper-motion-corrected Gaia DR3 position. Light curves: r = 2.0 px (adopted) and r = 2.5 px circular apertures without
the masked pixels, and an unmasked r = 2.0 px aperture (check); quality-0 cadences; in 200-s data, cadences more than
5 robust sigma from a 7-cadence running median are removed. Contamination: Gaussian PSF of the sigma (0.5-1.5 px) fitted to the
unmasked pixels of the median image, Gaia DR3 sources within 60" (targets.py). Primary eclipses with ingress and egress coverage are fitted
with a trapezoid (mid-time, depth, T14, T12) times a linear baseline over +-0.35 d, exposure-time supersampling; errors =
max(residual-scaled covariance, cyclic-shift bootstrap of the residuals, 60 draws). The linear ephemeris uses the
primaries whose window scatter is below twice the sector scatter and whose depth exceeds 5 sigma. Secondary eclipse:
the same fit at phase `ph2` in the 2019 sectors; in the other sectors the depth of a fixed trapezoid (T14 = 1.9 h,
T12 = 0.75 h) at phase `ph2` of the sector's own primaries, inverse-variance weighted over the covered epochs.
"""
import argparse, csv, json, os
import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from astropy.io.fits import Header
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from targets import OTHER

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
p.add_argument('--start', default='2026-10-19'); p.add_argument('--months', type=float, default=12)
args = p.parse_args(); TIC = 124751008; c = OTHER[TIC]; wd = os.path.join(args.work, str(TIC)); os.makedirs(args.out, exist_ok=True)
T0G, PG = c['eph0']; PH2 = c['ph2']; N = 21; yy, xx = np.indices((N, N)); rng = np.random.default_rng(2); rng_check = np.random.default_rng(3)
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
hdr = lambda lines: ''.join('# ' + l + '\n' for l in lines)
year = lambda t: Time(np.asarray(t) + 2457000, format='jd').jyear
tmag = lambda g, bprp: g - 0.00522555 * bprp**3 + 0.0891337 * bprp**2 - 0.633923 * bprp + 0.0324473   # Stassun et al. 2019


# ---------------------------------------------------------------- photometry
def photometry(sec):
    z = np.load(os.path.join(wd, 'cut', f's{sec}.npz')); tm = z['tm']; f = z['f'].astype(float); q = z['q']
    ok = np.isfinite(tm) & np.all(np.isfinite(f), axis=(1, 2)); tm, f, q = tm[ok], f[ok], q[ok]
    tx, ty = float(z['tx']), float(z['ty']); w = WCS(Header.fromstring(str(z['wcs'])))
    O0 = np.median(f[q == 0], 0); bleed = O0 > c['bleed']; bad = bleed.copy(); bad[:, 1:] |= bleed[:, :-1]; bad[:, :-1] |= bleed[:, 1:]
    r = np.hypot(xx - tx, yy - ty); far = (r > 5) & (O0 < np.percentile(O0[~bad], 40)) & ~bad
    A = np.c_[np.ones(far.sum()), xx[far], yy[far]]; Af = np.c_[np.ones(N * N), xx.ravel(), yy.ravel()]
    g = f - (Af @ np.linalg.lstsq(A, f[:, far].T, rcond=None)[0]).T.reshape(-1, N, N)
    O = np.median(g[q == 0], 0); m5 = (abs(xx - tx) <= 3) & (abs(yy - ty) <= 3) & ~bad
    gm = lambda p_: p_[3] * np.exp(-((xx - p_[0])**2 + (yy - p_[1])**2) / (2 * p_[2]**2)) + p_[4]
    x0, y0, sig = least_squares(lambda p_: (gm(p_) - O)[m5], [tx, ty, 0.7, O.max(), 0], bounds=([tx - 1, ty - 1, 0.5, 0, -np.inf], [tx + 1, ty + 1, 1.5, np.inf, np.inf])).x[:3]
    nb = [(*[float(v) for v in w.world_to_pixel(SkyCoord(ra, de, unit='deg'))], tmag(G, br)) for ra, de, G, br in c['neighbours']]
    tt = tmag(c['G'], c['bprp'])
    row = dict(sector=sec, camera=int(z['cam']), ccd=int(z['ccd']), cadence_min=float(np.median(np.diff(tm)) * 1440), start=float(tm.min()), end=float(tm.max()),
               n_q0=int((q == 0).sum()), psf_sigma=float(sig), nearest_bleed_px=float(r[bleed].min()) if bleed.any() else np.nan)
    lc = dict(t=tm, q=q)
    for lab, msk in (('m2.0', (r <= 2.0) & ~bad), ('m2.5', (r <= 2.5) & ~bad), ('u2.0', r <= 2.0)):
        lc[lab] = g[:, msk].sum(1)
        fin = lambda px, py: np.exp(-((xx[msk] - px)**2 + (yy[msk] - py)**2) / (2 * sig**2)).sum() / (2 * np.pi * sig**2)
        ft = fin(tx, ty); row[f'npix_{lab}'] = int(msk.sum()); row[f'target_frac_{lab}'] = float(ft)
        row[f'contamination_{lab}'] = float(sum(fin(px, py) * 10**(-0.4 * (tn - tt)) for px, py, tn in nb) / ft)
        row[f'flux_{lab}'] = float(np.median(lc[lab][q == 0]))
    row['n_masked_r2.0'] = int(((r <= 2.0) & bad).sum())
    np.savez(os.path.join(wd, f'lc_s{sec}.npz'), **lc)
    return row, lc


def clean(lc, key):
    t, f, q = lc['t'], lc[key].astype(float), lc['q']; g = (q == 0) & np.isfinite(f); t, f = t[g], f[g] / np.median(f[g])
    if np.median(np.diff(t)) * 1440 < 5:
        r_ = f - median_filter(f, size=7, mode='nearest'); s_ = 1.4826 * np.median(np.abs(r_)); k = np.abs(r_) < 5 * s_; t, f = t[k], f[k]
    return t, f


# ---------------------------------------------------------------- trapezoid fits
def trap(p_, x):
    t0, dep, T14, T12, a, b = p_; a_ = np.abs(x - t0); h = T14 / 2; rr = max(T12, 1e-4)
    return (a + b * (x - t0)) * (1 - np.where(a_ < h - rr, dep, np.where(a_ < h, dep * (h - a_) / rr, 0)))


def model(p_, x, texp):
    nsub = max(1, int(round(texp / 0.002))); off = ((np.arange(nsub) + 0.5) / nsub - 0.5) * texp
    return trap(p_, (x[:, None] + off[None, :]).ravel()).reshape(-1, nsub).mean(1)


def fitfree(t, f, tp, texp, rg=None):
    rg = rng if rg is None else rg
    lo = [tp - 0.1, -0.02, 0.4 / 24, 0.1 / 24, 0.9, -0.1]; hi = [tp + 0.1, 0.2, 3.0 / 24, 1.5 / 24, 1.1, 0.1]; best = None
    for dt0 in (-0.04, -0.02, 0, 0.02, 0.04):
        for T in (1.0 / 24, 1.5 / 24, 2.0 / 24):
            r_ = least_squares(lambda p_: model(p_, t, texp) - f, [tp + dt0, 0.005, T, 0.3 * T, 1, 0], bounds=(lo, hi))
            if best is None or r_.cost < best.cost: best = r_
    res = best.fun; s2 = (res**2).sum() / max(len(t) - 6, 1)
    try: e = np.sqrt(np.diag(np.linalg.inv(best.jac.T @ best.jac) * s2))
    except np.linalg.LinAlgError: e = np.full(6, np.nan)
    mf = model(best.x, t, texp); bs = []
    for _ in range(60):
        fb = mf + np.roll(res, rg.integers(1, len(t))); bs.append(least_squares(lambda p_: model(p_, t, texp) - fb, best.x, bounds=(lo, hi)).x)
    return best.x, np.fmax(e, np.nanstd(np.array(bs), 0)), np.sqrt(s2)


def fitfixed(t, f, tc, texp, T14, T12):
    shape = 1 - model([tc, 1.0, T14, T12, 1, 0], t, texp); A = np.c_[np.ones_like(t), t - tc, -shape]; cc = np.linalg.lstsq(A, f, rcond=None)[0]
    res = f - A @ cc; s = 1.4826 * np.median(np.abs(res - np.median(res))); return cc[2], np.sqrt((np.linalg.inv(A.T @ A) * s**2)[2, 2])


def fit_window(t, f, tp, texp, rg=None):
    m = np.abs(t - tp) < 0.35
    if m.sum() < 10: return None
    cin = np.sum((t > tp - 0.12) & (t < tp - 0.01)) >= max(2, 0.3 * 0.11 / texp); ceg = np.sum((t < tp + 0.12) & (t > tp + 0.01)) >= max(2, 0.3 * 0.11 / texp)
    if not (cin and ceg): return None
    tt = t[m]; ff = f[m] / np.median(f[m][np.abs(tt - tp) > 0.12]); x, e, rms = fitfree(tt, ff, tp, texp, rg)
    return x, e, rms, 1.4826 * np.median(np.abs(np.diff(ff))) / np.sqrt(2), int(m.sum())


SEC, LC = [], {}
for s in c['sectors']:
    row, lc = photometry(s); SEC.append(row); LC[s] = lc
    print(s, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
CL = {s: {k: clean(LC[s], k) for k in ('m2.0', 'm2.5', 'u2.0')} for s in LC}
PRI = []
for s in c['sectors']:
    t, f = CL[s]['m2.0']; texp = float(np.median(np.diff(t))); srms = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
    for n in range(int(np.floor((t.min() - T0G) / PG)) - 1, int(np.ceil((t.max() - T0G) / PG)) + 2):
        tp = T0G + n * PG; R = fit_window(t, f, tp, texp)
        if R is None: continue
        x, e, rms, wrms, npts = R; row = dict(sector=s, cycle=n, t0=float(x[0]), t0_err=float(e[0]), depth=float(x[1]), depth_err=float(e[1]), T14_h=float(x[2] * 24), T14_err_h=float(e[2] * 24),
                                              T12_h=float(x[3] * 24), rms_ppm=float(rms * 1e6), noisy=bool(wrms > 2 * srms), npts=npts, cadence_min=texp * 1440)
        for k in ('m2.5', 'u2.0'):
            t2, f2 = CL[s][k]; R2 = fit_window(t2, f2, x[0], texp, rng_check); row[f'depth_{k}'] = float(R2[0][1]) if R2 else np.nan
        PRI.append(row); print('primary', {k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
good = [r for r in PRI if not r['noisy'] and r['depth'] / r['depth_err'] > 5]
ok_err = lambda r: r['depth_err'] < r['depth']      # False where the bootstrap error diverged (shown as empty in the tables)
n = np.array([r['cycle'] for r in good]); tc = np.array([r['t0'] for r in good]); te = np.array([r['t0_err'] for r in good])
A = np.c_[np.ones_like(n), n]; W = 1 / te**2; cov = np.linalg.inv(A.T @ (A * W[:, None])); cc = cov @ (A.T @ (W * tc))
T0, P = float(cc[0]), float(cc[1]); eT0, eP = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
for r in PRI:
    r['oc_min'] = (r['t0'] - (T0 + r['cycle'] * P)) * 1440; r['in_ephemeris'] = int(r in good)
chi2 = float(np.sum(((tc - A @ cc) / te)**2)); print('ephemeris', T0, eT0, P, eP, 'chi2', chi2, len(n) - 2)
# secondary eclipses
SECD = []
for s in c['sectors']:
    t, f = CL[s]['m2.0']; texp = float(np.median(np.diff(t))); pr = [r for r in good if r['sector'] == s] or [r for r in PRI if r['sector'] == s]
    ref = pr[0]['t0'] - pr[0]['cycle'] * PG if pr else T0
    ks = [k for k in range(int(np.floor((t.min() - ref) / PG)) - 1, int(np.ceil((t.max() - ref) / PG)) + 1) if np.sum(np.abs(t - (ref + (k + PH2) * PG)) < 0.3) >= 8]
    free = []
    if s in (6, 7):
        for k in ks:
            R = fit_window(t, f, ref + (k + PH2) * PG, texp, rng_check)
            if R: free.append(dict(t0=float(R[0][0]), depth=float(R[0][1]), depth_err=float(R[1][1]), T14_h=float(R[0][2] * 24), phase=float(((R[0][0] - ref) / PG) % 1)))
    ds, es = [], []
    for k in ks:
        tcand = ref + (k + PH2) * PG; m = np.abs(t - tcand) < 0.3
        if np.sum(np.abs(t[m] - tcand) < 0.95 / 24) < max(2, 0.5 * 1.9 / 24 / texp): continue
        d_, e_ = fitfixed(t[m], f[m] / np.median(f[m]), tcand, texp, 1.9 / 24, 0.75 / 24); ds.append(d_); es.append(e_)
    if ds:
        w_ = 1 / np.array(es)**2; SECD.append(dict(sector=s, n_epochs=len(ds), depth=float(np.sum(w_ * np.array(ds)) / w_.sum()), err=float(1 / np.sqrt(w_.sum())), free=free))
    print('secondary', SECD[-1] if ds else (s, 'no coverage'), flush=True)
json.dump(dict(sectors=SEC, primaries=PRI, ephemeris=dict(T0=T0, eT0=eT0, P=P, eP=eP, chi2=chi2, dof=len(n) - 2), secondaries=SECD), open(os.path.join(wd, 'eclipses.json'), 'w'), indent=1)

# ---------------------------------------------------------------- tables
with open(os.path.join(args.out, 'sectors_124751008.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {TIC} = Gaia DR3 {c["gaia_dr3"]} ({c["name"]}): TESS sectors used (TESScut 21x21 cutouts at the Gaia DR3 position), {DATE}',
                  f'nearest_bleed_px = distance from the target to the nearest pixel of the kappa Col saturation bleed column (median image > {c["bleed"]:.0f} e-/s); n_masked_r2.0 = bleed and neighbouring pixels removed from the r = 2.0 px aperture',
                  'flux_* = median aperture flux (e-/s) of the masked r = 2.0 and 2.5 px apertures and the unmasked r = 2.0 px aperture; contamination_* = Gaia DR3 neighbour flux / target flux in the aperture (Gaussian PSF of sigma psf_sigma px; kappa Col not included)']))
    keys = ['sector', 'camera', 'ccd', 'cadence_min', 'start', 'end', 'n_q0', 'psf_sigma', 'nearest_bleed_px', 'n_masked_r2.0', 'npix_m2.0', 'flux_m2.0', 'flux_m2.5', 'flux_u2.0', 'contamination_m2.0', 'contamination_m2.5', 'contamination_u2.0']
    w = csv.writer(fo); w.writerow(['sector', 'camera', 'ccd', 'cadence_min', 'start_utc', 'end_utc', 'n_q0', 'psf_sigma_px', 'nearest_bleed_px', 'n_masked_r2.0', 'npix_m2.0', 'flux_m2.0', 'flux_m2.5', 'flux_u2.0', 'contamination_m2.0', 'contamination_m2.5', 'contamination_u2.0'])
    for r in SEC:
        w.writerow([r['sector'], r['camera'], r['ccd'], f'{r["cadence_min"]:.2f}', Time(r['start'] + 2457000, format='jd').iso[:10], Time(r['end'] + 2457000, format='jd').iso[:10], r['n_q0'], f'{r["psf_sigma"]:.2f}',
                    f'{r["nearest_bleed_px"]:.2f}', r['n_masked_r2.0'], r['npix_m2.0'], f'{r["flux_m2.0"]:.0f}', f'{r["flux_m2.5"]:.0f}', f'{r["flux_u2.0"]:.0f}', f'{r["contamination_m2.0"]:.4f}', f'{r["contamination_m2.5"]:.4f}', f'{r["contamination_u2.0"]:.4f}'])
with open(os.path.join(args.out, 'eclipse_times_124751008.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {TIC}: primary eclipses fitted with a trapezoid (scripts/eb_124751008.py; masked r = 2.0 px aperture), {DATE}',
                  f'linear ephemeris of the primaries with in_ephemeris = 1: T_pri = BTJD {T0:.5f} + {P:.7f} E (formal errors {eT0:.5f} d and {eP:.7f} d; chi2 = {chi2:.0f} for {len(n) - 2} degrees of freedom); cycle = E',
                  't0_err_min, depth_err_pct = max(formal, cyclic-shift bootstrap; empty where the bootstrap diverged); noisy = window scatter above twice the sector scatter; depth_m2.5_pct, depth_u2.0_pct = the same fit on the masked r = 2.5 px and the unmasked r = 2.0 px apertures (checks)']))
    w = csv.writer(fo); w.writerow(['sector', 'cycle', 'year', 't0_btjd', 't0_err_min', 'oc_min', 'depth_pct', 'depth_err_pct', 'T14_h', 'T14_err_h', 'T12_h', 'depth_m2.5_pct', 'depth_u2.0_pct', 'noisy', 'in_ephemeris'])
    for r in PRI:
        w.writerow([r['sector'], r['cycle'], f'{float(year(r["t0"])):.3f}', f'{r["t0"]:.5f}', f'{r["t0_err"] * 1440:.2f}', f'{r["oc_min"]:+.1f}', f'{r["depth"] * 100:.3f}', f'{r["depth_err"] * 100:.3f}' if ok_err(r) else '',
                    f'{r["T14_h"]:.2f}', f'{r["T14_err_h"]:.2f}', f'{r["T12_h"]:.2f}', f'{r["depth_m2.5"] * 100:.3f}', f'{r["depth_u2.0"] * 100:.3f}', int(r['noisy']), r['in_ephemeris']])
with open(os.path.join(args.out, 'depths_124751008.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {TIC}: eclipse depth and duration per sector (scripts/eb_124751008.py; masked r = 2.0 px aperture), {DATE}',
                  'kind = pri: weighted mean of the primaries of the sector (noisy windows and diverged errors excluded); T14_err_h = scatter of T14 between the eclipses; kind = sec: weighted mean of the free trapezoid fits (sectors 6 and 7; T14_err_h = scatter) or, elsewhere, the fixed-trapezoid depth (T14 = 1.9 h, T12 = 0.75 h) at the secondary phase of the sector, weighted over n_events',
                  'depth_fixed_pct, depth_fixed_err_pct = the fixed-trapezoid depth at the secondary phase in every sector (in sectors 6 and 7 lower than the free fits because the fixed shape differs)',
                  f'limit3_pct = 3-sigma upper limit where the secondary is not detected (depth < 3 err); phase = secondary phase relative to the sector primaries (nominal {PH2})']))
    w = csv.writer(fo); w.writerow(['sector', 'year', 'kind', 'n_events', 'depth_pct', 'depth_err_pct', 'limit3_pct', 'T14_h', 'T14_err_h', 'phase', 'depth_fixed_pct', 'depth_fixed_err_pct'])
    for s in c['sectors']:
        pr = [r for r in PRI if r['sector'] == s and not r['noisy'] and ok_err(r)] or [r for r in PRI if r['sector'] == s and ok_err(r)]
        if pr:
            wt = 1 / np.array([r['depth_err'] for r in pr])**2; dm = np.sum(wt * [r['depth'] for r in pr]) / wt.sum(); de = 1 / np.sqrt(wt.sum())
            w.writerow([s, f'{float(np.mean(year([r["t0"] for r in pr]))):.3f}', 'pri', len(pr), f'{dm * 100:.3f}', f'{de * 100:.3f}', '', f'{np.mean([r["T14_h"] for r in pr]):.2f}', f'{np.std([r["T14_h"] for r in pr]):.2f}', '', '', ''])
        sd = [x for x in SECD if x['sector'] == s]
        if sd:
            x = sd[0]; yr_ = float(np.mean(year([r['t0'] for r in pr]))) if pr else float(year(np.mean([LC[s]['t'].min(), LC[s]['t'].max()])))
            fr = [y for y in x['free'] if y['depth_err'] < y['depth']]
            if fr:
                wt = 1 / np.array([y['depth_err'] for y in fr])**2; dm = np.sum(wt * [y['depth'] for y in fr]) / wt.sum(); de = 1 / np.sqrt(wt.sum())
                w.writerow([s, f'{yr_:.3f}', 'sec', len(fr), f'{dm * 100:.3f}', f'{de * 100:.3f}', '', f'{np.mean([y["T14_h"] for y in fr]):.2f}', f'{np.std([y["T14_h"] for y in fr]):.2f}', f'{np.mean([y["phase"] for y in fr]):.4f}', f'{x["depth"] * 100:.3f}', f'{x["err"] * 100:.3f}'])
            else:
                det = x['depth'] > 3 * x['err']
                w.writerow([s, f'{yr_:.3f}', 'sec', x['n_epochs'], f'{x["depth"] * 100:.3f}', f'{x["err"] * 100:.3f}', '' if det else f'{max(x["depth"], 0) * 100 + 3 * x["err"] * 100:.3f}', '1.90 (fixed)', '', f'{PH2}', f'{x["depth"] * 100:.3f}', f'{x["err"] * 100:.3f}'])
COORD = SkyCoord(c['ra'] * u.deg, c['dec'] * u.deg); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
def utc(bjd):
    t = Time(bjd, format='jd', scale='tdb', location=GEO)
    for _ in range(2): t = Time(bjd, format='jd', scale='tdb', location=GEO) - Time(t, location=GEO).light_travel_time(COORD, kind='barycentric')
    return t.utc.iso[:16]
t_start = Time(args.start).jd - 2457000; t_end = t_start + args.months * 30.5
with open(os.path.join(args.out, 'predicted_eclipses_124751008.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {TIC}: predicted primary and secondary eclipses for the {args.months:g} months from {args.start} (linear ephemeris of eclipse_times_124751008.csv; secondary at phase {PH2}), {DATE}',
                  '*_utc = geocentric UTC (barycentric light-travel time toward the target removed); the observed O-C spans -22 to +68 min with changes of up to 43 min between consecutive cycles, so a window of at least +-70 min should be used',
                  'the secondary eclipse has not been detected since 2019 (README section 6.4); its predicted times are given for the test of its return']))
    w = csv.writer(fo); w.writerow(['kind', 'cycle', 'tc_btjd', 'tc_bjd_tdb', 'tc_utc', 'ephemeris_err_min'])
    rows = []
    for kind, ph in (('pri', 0.0), ('sec', PH2)):
        for k in range(int(np.ceil((t_start - T0) / P - ph)), int(np.floor((t_end - T0) / P - ph)) + 1):
            tp = T0 + (k + ph) * P; rows.append((tp, kind, k))
    for tp, kind, k in sorted(rows):
        w.writerow([kind, k, f'{tp:.4f}', f'{tp + 2457000:.4f}', utc(tp + 2457000), f'{np.hypot(eT0, k * eP) * 1440:.1f}'])
print('tables written to', args.out)
