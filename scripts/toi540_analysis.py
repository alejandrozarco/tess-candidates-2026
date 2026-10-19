"""TOI-540 b (TIC 200322593): transit depths per sector, ephemeris, sector-98 tests, b(t) fits, ground-based windows.

usage: python scripts/toi540_analysis.py [--work ./work] [--out data] [--start 2026-10-19] [--months 3] [--skip-injections]
input:  SPOC 2-min target pixel files and light curves of sectors 4, 5, 6, 31, 32, 98 and the 20-s light curves of
        sectors 31, 32, 98 in <work>/200322593/dl/mastDownload/TESS/ (other_download.py 200322593)
output: data/transit_depths_200322593.csv, ephemeris_200322593.csv, s98_tests_200322593.csv, injection_200322593.csv,
        bt_fits_200322593.csv, perturber_200322593.csv, ground_windows_200322593.csv and folded_200322593.csv
        (formats in the file headers); <work>/200322593/phot_s<sector>.npz.
Photometry ("own"): sum of the target-pixel-file flux with the SPOC background added back and a per-cadence
background (median of the pixels more than 4 px from the target and below the 50th percentile of the median image)
subtracted, in circular apertures of r = 1.5-3.0 px around the proper-motion-corrected Gaia DR3 position and in the
SPOC aperture. Light curves: quality 0, divided by a biweight trend (wotan, window 0.3 d unless stated), points more
than 4 sigma above or 8 sigma below the trend removed. Transit model: batman, quadratic limb darkening u = (0.25, 0.35),
circular orbit, 2-min exposure integration; a/R* has the prior 13.87 +- 0.71 from the mass and radius of Ment et al.
(2021). Depths are peak depths of the common transit shape (template amplitude), with errors scaled by the red-noise
factor beta of the residuals at 28 min. Reference depth = weighted mean of sectors 4-32 (own, r = 2.0 px, 0.3 d).
"""
import argparse, csv, glob, json, os
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as CHI2, norm
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import astropy.units as u
from wotan import flatten
from targets import OTHER, import_batman

batman = import_batman()
p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
p.add_argument('--start', default='2026-10-19'); p.add_argument('--months', type=float, default=3); p.add_argument('--skip-injections', action='store_true')
args = p.parse_args(); TIC = 200322593; c = OTHER[TIC]; wd = os.path.join(args.work, str(TIC)); os.makedirs(args.out, exist_ok=True)
DL = os.path.join(wd, 'dl', 'mastDownload', 'TESS'); PRE = c['pre']; YR = {4: 2018.84, 5: 2018.92, 6: 2018.99, 31: 2020.84, 32: 2020.91, 98: 2025.93}
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
hdr = lambda lines: ''.join('# ' + l + '\n' for l in lines)
G = 6.674e-11; MSUN = 1.989e30; RSUN = 6.957e8; TEXP = 120 / 86400
rng = np.random.default_rng(1); Ms = rng.normal(c['star']['M'], c['star']['eM'], 20000); Rs = rng.normal(c['star']['R'], c['star']['eR'], 20000)
aRs = (G * Ms * MSUN * (c['toi'][1] * 86400)**2 / (4 * np.pi**2))**(1 / 3) / (Rs * RSUN); AR0, EAR = float(np.median(aRs)), float(np.std(aRs))
print(f'a/R* prior from Ment et al. (2021) M, R: {AR0:.2f} +- {EAR:.2f}')
BAD = 1 | 2 | 4 | 8 | 32 | 128 | 1024 | 2048 | 8192 | 16384      # 'default' quality mask (scattered-light and impulsive-outlier bits kept)


# ---------------------------------------------------------------- photometry and light curves
def phot(sec):
    out = os.path.join(wd, f'phot_s{sec}.npz')
    if os.path.exists(out): return
    fn = glob.glob(os.path.join(DL, f'*s{sec:04d}-{TIC:016d}-*-s', '*-s_tp.fits'))[0]
    with fits.open(fn, memmap=True) as h:
        d = h[1].data; t = np.array(d['TIME'], float); q = np.array(d['QUALITY']); F = np.array(d['FLUX'], np.float32); B = np.array(d['FLUX_BKG'], np.float32)
        ap = np.array(h[2].data, int); w = WCS(h[2].header)
    yr = 2016.0 + (np.nanmedian(t) + 2457000 - 2457388.5) / 365.25
    px, py = w.all_world2pix([[c['ra'] + c['pmra'] * (yr - 2016) / 3.6e6 / np.cos(np.radians(c['dec'])), c['dec'] + c['pmdec'] * (yr - 2016) / 3.6e6]], 0)[0]
    yy, xx = np.mgrid[0:F.shape[1], 0:F.shape[2]]; rr = np.hypot(xx - px, yy - py); raw = F + B; med = np.nanmedian(raw, 0)
    bm = (rr > 4) & (med < np.nanpercentile(med[rr > 4], 50)); own = raw - np.nanmedian(raw[:, bm], 1)[:, None, None]
    res = dict(t=t, q=q, px=px, py=py)
    for k, m in {'spoc': (ap & 2) > 0, 'r1.5': rr <= 1.5, 'r2.0': rr <= 2.0, 'r2.5': rr <= 2.5, 'r3.0': rr <= 3.0}.items():
        res['f_' + k] = np.nansum(own[:, m], 1); res['npix_' + k] = int(m.sum())
    np.savez(out, **res); print(f'S{sec}: target pixel ({px:.2f}, {py:.2f}), SPOC aperture {res["npix_spoc"]} px', flush=True)


def detrend(t, f, win):
    tr = flatten(t, f, method='biweight', window_length=win, return_trend=True, break_tolerance=0.3, edge_cutoff=0.)[1]
    o = f / tr; ok = np.isfinite(o); s = 1.4826 * np.nanmedian(np.abs(o[ok] - 1)); ok &= (o - 1 < 4 * s) & (o - 1 > -8 * s)
    return t[ok], o[ok], s


def raw_lc(sec, src='own', ap='r2.0', qual='zero'):
    if src == 'own':
        z = np.load(os.path.join(wd, f'phot_s{sec}.npz')); t, q, f = z['t'], z['q'], z['f_' + ap]
    else:
        pat = '*-a_fast-lc.fits' if src.endswith('fast') else '*-s_lc.fits'
        d = fits.open(glob.glob(os.path.join(DL, f'*s{sec:04d}-{TIC:016d}-*', pat))[0])[1].data
        t, q, f = np.array(d['TIME'], float), np.array(d['QUALITY']), np.array(d['SAP_FLUX' if src == 'sap' else 'PDCSAP_FLUX'], float)
    m = np.isfinite(t) & np.isfinite(f) & (f > 0)
    if qual == 'zero': m &= q == 0
    elif qual == 'default': m &= (q & BAD) == 0
    return t[m], f[m] / np.nanmedian(f[m])


def lc(sec, src='own', win=0.3, ap='r2.0', qual='zero'):
    t, f = raw_lc(sec, src, ap, qual); t, f, s = detrend(t, f, win); return t, f, s, float(np.nanmedian(np.diff(t)))


def bmodel(t, tc, P, k, aR, b, texp=TEXP):
    pp = batman.TransitParams(); pp.t0 = tc; pp.per = P; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(min(b / aR, 1))); pp.ecc = 0; pp.w = 90
    pp.u = c['u']; pp.limb_dark = 'quadratic'
    return batman.TransitModel(pp, t, supersample_factor=max(1, int(round(texp * 86400 / 24))), exp_time=texp).light_curve(pp)


def t14(P, k, aR, b):
    return P / np.pi * np.arcsin(min(1, np.sqrt(max((1 + k)**2 - b**2, 0)) / (aR * np.sin(np.arccos(b / aR)))))


SHAPES = {}
def grid(P, k, aR, b, texp):
    """unit-peak transit shape on a 5-s grid of +-0.06 d (cached)"""
    key = (round(P, 9), round(k, 6), round(aR, 4), round(b, 5), round(texp * 86400))
    if key not in SHAPES:
        g = np.arange(-0.06, 0.06, 5 / 86400.); s = 1 - bmodel(g, 0.0, P, k, aR, b, texp); SHAPES[key] = (g, s / s.max() if s.max() > 0 else s)
    return SHAPES[key]


def template(t, tc, P, k, aR, b, texp):
    g, s = grid(P, k, aR, b, texp); return np.interp((((t - tc) / P + 0.5) % 1 - 0.5) * P, g, s, left=0, right=0)


def beta_factor(t, r, cad, tbin=28 / 1440):
    n = max(int(round(tbin / cad)), 1); rr = r[np.argsort(t)]; nb = len(rr) // n
    if nb < 10: return 1.0
    return max(1.0, np.std(rr[:nb * n].reshape(nb, n).mean(1)) / (np.std(r) / np.sqrt(n) * np.sqrt(nb / (nb - 1))))


def fit_depth(t, f, tc, P, shape, cad, sig, beta=1.0):
    s = template(t, tc, P, *shape, cad); ss = np.sum(s * s)
    if ss == 0: return np.nan, np.nan
    return -np.sum((f - 1) * s) / ss, sig / np.sqrt(ss) * beta


def phase_scan(t, f, tc, P, shape, cad, sig, beta, offs_min):
    """template amplitude and S/N at each offset, on 0.5-min phase bins of the residuals"""
    dt = (((t - tc) / P + 0.5) % 1 - 0.5) * P * 1440; edges = np.arange(-P * 720, P * 720 + 0.5, 0.5); cen = 0.5 * (edges[1:] + edges[:-1])
    R = np.histogram(dt, edges, weights=f - 1)[0]; Nn = np.histogram(dt, edges)[0]
    g, s = grid(P, *shape, cad)
    A, S = [], []
    for o in offs_min:
        x = ((cen - o + P * 720) % (P * 1440)) - P * 720; sb = np.interp(x / 1440, g, s, left=0, right=0); den = np.sum(sb * sb * Nn)
        a = -np.sum(sb * R) / den; A.append(a); S.append(a / (sig / np.sqrt(den) * beta))
    return np.array(A), np.array(S)


for s in c['sectors']:
    phot(s)

# ---------------------------------------------------------------- shape, per-sector times and ephemeris (sectors 4-32)
def ephemeris(src='own', win=0.3, ap='r2.0'):
    TC0, P0 = c['toi'][0] - 2457000, c['toi'][1]; data = {s: lc(s, src, win, ap) for s in PRE}
    sel = lambda t, tc, P, w: np.abs(((t - tc) / P + 0.5) % 1 - 0.5) * P < w
    tt = np.concatenate([data[s][0] for s in PRE]); ff = np.concatenate([data[s][1] for s in PRE]); ee = np.concatenate([np.full(len(data[s][0]), data[s][2]) for s in PRE])
    m = sel(tt, TC0, P0, 0.2); tt, ff, ee = tt[m], ff[m], ee[m]
    res = lambda x: np.append((ff - bmodel(tt, *x)) / ee, (x[3] - AR0) / EAR)
    best = None
    for b0 in (0.3, 0.6, 0.75, 0.85):
        for dt in (-0.006, 0, 0.006):
            r = least_squares(res, [TC0 + dt, P0, 0.044, AR0, b0], x_scale=[1e-3, 1e-6, 1e-3, 0.5, 0.05], bounds=([TC0 - 0.05, P0 - 1e-4, 0.01, 3, 0], [TC0 + 0.05, P0 + 1e-4, 0.2, 40, 1.2]))
            if best is None or r.cost < best.cost: best = r
    chi2 = 2 * best.cost; cov = np.linalg.inv(best.jac.T @ best.jac) * max(1, chi2 / (len(best.fun) - 5)); x = best.x; ex = np.sqrt(np.diag(cov))
    tc_, P_, k, aR, b = x; shape = (k, aR, b); rows = []
    for s in PRE:
        t, f, sg, cad = data[s]; E = np.round((np.median(t) - tc_) / P_); tref = tc_ + E * P_; mm = sel(t, tref, P_, 0.2); t, f = t[mm], f[mm]
        dts = np.arange(-30, 30.01, 0.25); chi = []
        for d in dts:
            stp = template(t, tref + d / 1440, P_, *shape, cad); A = -np.sum((f - 1) * stp) / np.sum(stp**2); chi.append(np.sum(((f - 1 + A * stp) / sg)**2))
        chi = np.array(chi) - min(chi); i = np.argmin(chi); beta = beta_factor(t, f - 1, cad); inside = dts[chi <= beta**2]
        A, eA = fit_depth(t, f, tref + dts[i] / 1440, P_, shape, cad, sg, beta)
        rows.append(dict(sector=s, E=E, tc=tref + dts[i] / 1440, sig=max((inside.max() - inside.min()) / 2, 0.25) / 1440, depth=A, edepth=eA, beta=beta, ntr=len(np.unique(np.round((t - tref) / P_)))))
    E = np.array([r['E'] for r in rows]); tcs = np.array([r['tc'] for r in rows]); sig = np.array([r['sig'] for r in rows]); Eref = float(np.round(np.average(E, weights=1 / sig**2)))
    X = np.c_[np.ones_like(E), E - Eref]; W = np.diag(1 / sig**2); C = np.linalg.inv(X.T @ W @ X); be = C @ X.T @ W @ tcs; oc = tcs - X @ be; ce = float(np.sum((oc / sig)**2))
    C = C * max(1, ce / (len(E) - 2))
    for r, o in zip(rows, oc): r['oc_min'] = o * 1440
    T14 = t14(P_, k, aR, b) * 1440; pk = 1 - bmodel(np.array([tc_]), tc_, P_, k, aR, b)[0]
    print(f'[{src} {ap} w={win}] shape k={k:.4f} b={b:.3f} a/R*={aR:.2f} T14={T14:.1f} min peak={pk*1e6:.0f} ppm; ephemeris T0={be[0]:.6f} (E={Eref:.0f}) P={be[1]:.8f}+-{np.sqrt(C[1,1]):.8f}; chi2={ce:.2f}/{len(E)-2}', flush=True)
    return dict(src=src, win=win, ap=ap, shape=shape, eshape=ex[2:], T14=T14, peak=pk, chi2_fit=chi2, n_fit=len(ff), T0=be[0], P=be[1], Eref=Eref, cov=C, chi2_eph=ce, rows=rows)


EPH = {key: ephemeris(*key) for key in [('own', 0.3, 'r2.0'), ('own', 0.5, 'r2.0'), ('pdc', 0.3, 'r2.0')]}
R0 = EPH[('own', 0.3, 'r2.0')]; T0, P, Eref, C, SH = R0['T0'], R0['P'], R0['Eref'], R0['cov'], R0['shape']
pred = lambda E: (T0 + (E - Eref) * P, float(np.sqrt(np.array([1, E - Eref]) @ C @ np.array([1, E - Eref]))))
REF = float(np.average([r['depth'] for r in R0['rows']], weights=[1 / r['edepth']**2 for r in R0['rows']]))
print(f'reference peak depth S4-S32 = {REF*1e6:.0f} ppm')
with open(os.path.join(args.out, 'ephemeris_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: common transit shape and linear ephemeris from TESS sectors 4, 5, 6, 31, 32 (scripts/toi540_analysis.py), {DATE}',
                  'source = own TPF photometry (r = 2.0 px) or SPOC PDCSAP; window_d = biweight detrending window; k, b, a/R* from the joint fit with the a/R* prior; T14_min, peak_ppm = model duration and peak depth',
                  'T0_btjd at cycle E_ref (cycles counted from the ExoFOP TOI epoch BJD 2459175.141647), P_d with errors scaled by max(1, chi2/dof) of the five per-sector mid-transit times; cov_T0_P = covariance (d^2)']))
    w = csv.writer(fo); w.writerow(['source', 'window_d', 'k', 'k_err', 'b', 'b_err', 'aR', 'aR_err', 'T14_min', 'peak_ppm', 'chi2_fit', 'n_points', 'E_ref', 'T0_btjd', 'T0_err_d', 'P_d', 'P_err_d', 'cov_T0_P', 'chi2_ephemeris', 'dof'])
    for (src, win, ap), R in EPH.items():
        w.writerow([src, win, f'{R["shape"][0]:.4f}', f'{R["eshape"][0]:.4f}', f'{R["shape"][2]:.3f}', f'{R["eshape"][2]:.3f}', f'{R["shape"][1]:.2f}', f'{R["eshape"][1]:.2f}', f'{R["T14"]:.1f}',
                    f'{R["peak"]*1e6:.0f}', f'{R["chi2_fit"]:.1f}', R['n_fit'], f'{R["Eref"]:.0f}', f'{R["T0"]:.6f}', f'{np.sqrt(R["cov"][0,0]):.6f}', f'{R["P"]:.8f}', f'{np.sqrt(R["cov"][1,1]):.8f}',
                    f'{R["cov"][0,1]:.3e}', f'{R["chi2_eph"]:.2f}', len(PRE) - 2])

# ---------------------------------------------------------------- per-sector depths and sector-98 tests
def analyse(t, f, sg, cad, scan=True):
    E = np.round((np.median(t) - T0) / P) + Eref; tc, stc = pred(E); beta = beta_factor(t, f - 1, cad); A, eA = fit_depth(t, f, tc, P, SH, cad, sg, beta)
    ntr = len(np.unique(np.round((t[np.abs(((t - tc) / P + 0.5) % 1 - 0.5) * P < 0.007] - tc) / P)))
    out = dict(n=len(t), rms=sg, beta=beta, E=E, tc=tc, stc=stc, ntr=ntr, depth=A, edepth=eA)
    if scan:
        offs = np.arange(-P * 720, P * 720, 2.0); a, snr = phase_scan(t, f, tc, P, SH, cad, sg, beta, offs); i = np.argmax(snr); wn = np.abs(offs) <= 5 * stc * 1440
        out.update(scan_max=snr[i], scan_off=offs[i], scan_near=snr[wn].max(), scan_rms=np.std(snr), n_half=int((a > 0.5 * REF).sum()), n_off=len(offs))
    return out


TESTS = []
def test(label, sec, src, ap, qual, win, scan=True, note=''):
    t, f, sg, cad = lc(sec, src, win, ap, qual); o = analyse(t, f, sg, cad, scan); o.update(test=label, sector=sec, src=src, ap=ap if src == 'own' else '', qual=qual, win=win, note=note, cad=cad)
    TESTS.append(o); print(f"{label} S{sec} {src} {ap} q={qual} w={win}: depth {o['depth']*1e6:.0f} +- {o['edepth']*1e6:.0f} ppm ({(REF - o['depth']) / o['edepth']:.1f} sigma below ref)" +
                           (f"; scan max S/N {o['scan_max']:.1f} at {o['scan_off']:+.0f} min" if scan else ''), flush=True)
    return o, (t, f, sg, cad)


for src in ('own', 'pdc', 'sap'):
    for win in (0.2, 0.4):
        test('S98 product and window', 98, src, 'r2.0', 'zero', win)
for ap in ('r1.5', 'r2.5', 'r3.0', 'spoc'):
    test('S98 aperture', 98, 'own', ap, 'zero', 0.3, scan=False)
for qual in ('default', 'none'):
    test('S98 quality flags', 98, 'own', 'r2.0', qual, 0.3)
test('S98 20-s cadence', 98, 'pdc_fast', '', 'zero', 0.3)
PER = {}
for s in PRE + [98]:
    PER[s] = {k: test('per sector', s, src, 'r2.0', 'zero', 0.3, scan=(k == 'own'))[0] for k, src in (('own', 'own'), ('pdc', 'pdc'), ('sap', 'sap'))}
for s in (31, 32):
    PER[s]['fast'] = test('per sector 20-s', s, 'pdc_fast', '', 'zero', 0.3, scan=False)[0]
PER[98]['fast'] = [o for o in TESTS if o['test'] == 'S98 20-s cadence'][0]
o98, (t98, f98, s98, c98) = test('S98 reference', 98, 'own', 'r2.0', 'zero', 0.3)
b98 = o98['beta']; tc98, stc98 = o98['tc'], o98['stc']
gaps = np.where(np.diff(t98) > 0.5)[0]; edges = np.r_[0, gaps + 1, len(t98)]
for a0, a1 in zip(edges[:-1], edges[1:]):
    if a1 - a0 < 500: continue
    A, eA = fit_depth(t98[a0:a1], f98[a0:a1], tc98, P, SH, c98, s98, b98)
    TESTS.append(dict(test='S98 segment', sector=98, src='own', ap='r2.0', qual='zero', win=0.3, n=a1 - a0, rms=s98, beta=b98, depth=A, edepth=eA, ntr='', note=f'BTJD {t98[a0]:.2f}-{t98[a1-1]:.2f}'))
best = (-99, 0, 0, 0)
for Pt in P + np.arange(-6e-4, 6.0001e-4, 2e-5):
    offs = np.arange(-Pt * 720, Pt * 720, 2.0); a, snr = phase_scan(t98, f98, tc98, Pt, SH, c98, s98, b98, offs); i = np.argmax(snr)
    if snr[i] > best[0]: best = (snr[i], Pt, offs[i], a[i])
TESTS.append(dict(test='S98 period x phase grid', sector=98, src='own', ap='r2.0', qual='zero', win=0.3, n=len(t98), rms=s98, beta=b98, depth=best[3], edepth=np.nan, scan_max=best[0], scan_off=best[2],
                  note=f'61 periods P +- 6e-4 d ({6e-4 / np.sqrt(C[1,1]):.0f} sigma_P), 2-min phase steps; best P = {best[1]:.6f} d'))
for b in (0.78, 0.85, 0.90, 0.95, 0.98, 1.00, 1.02):
    for src in ('own', 'pdc'):
        t, f, sg, cad = lc(98, src, 0.2); be = beta_factor(t, f - 1, cad); shb = (SH[0], SH[1], b); pk = 1 - bmodel(np.array([0.]), 0, P, SH[0], SH[1], b)[0]
        A, eA = fit_depth(t, f, tc98, P, shb, cad, sg, be)
        TESTS.append(dict(test='S98 grazing shape', sector=98, src=src, ap='r2.0' if src == 'own' else '', qual='zero', win=0.2, n=len(t), rms=sg, beta=be, depth=A, edepth=eA, expected=pk,
                          note=f'b = {b:.2f}, T14 = {t14(P, SH[0], SH[1], b) * 1440:.1f} min'))
with open(os.path.join(args.out, 's98_tests_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: transit depth at the predicted times in TESS sector 98 and checks (scripts/toi540_analysis.py), {DATE}',
                  f'ephemeris: own photometry, sectors 4-32, T0 = BTJD {T0:.6f} at cycle {Eref:.0f}, P = {P:.8f} d; predicted mid-time of the central sector-98 transit BTJD {tc98:.5f} +- {stc98*1440:.2f} min',
                  f'depth_ppm = peak depth of the common shape (k = {SH[0]:.4f}, a/R* = {SH[1]:.2f}, b = {SH[2]:.3f}, T14 = {R0["T14"]:.1f} min) fitted at the predicted times; error scaled by beta; reference depth (sectors 4-32) = {REF*1e6:.0f} ppm',
                  'sigma_below = (expected - depth)/error, expected = reference depth or, for the grazing-shape rows, the peak depth of the shape with the b given in the note',
                  'scan_max_snr, scan_offset_min = highest S/N of the same template at any phase (2-min steps, 0.5-min phase bins) and its offset from the prediction; scan_near_snr = highest S/N within +-5 sigma of the prediction',
                  'src: own = TPF photometry (aperture ap), pdc / sap = SPOC 2-min PDCSAP / SAP, pdc_fast = SPOC 20-s PDCSAP; qual: zero = quality 0, default = standard mask, none = all finite cadences']))
    w = csv.writer(fo); w.writerow(['test', 'sector', 'src', 'ap', 'qual', 'window_d', 'n_points', 'rms_ppm', 'beta', 'n_transits', 'depth_ppm', 'depth_err_ppm', 'expected_ppm', 'sigma_below', 'scan_max_snr', 'scan_offset_min', 'scan_near_snr', 'note'])
    for o in TESTS:
        if o['test'].startswith('per sector'): continue
        ex = o.get('expected', REF)
        w.writerow([o['test'], o['sector'], o['src'], o['ap'], o['qual'], o['win'], o['n'], f"{o['rms']*1e6:.0f}", f"{o['beta']:.2f}", o.get('ntr', ''), f"{o['depth']*1e6:.0f}",
                    '' if not np.isfinite(o['edepth']) else f"{o['edepth']*1e6:.0f}", f'{ex*1e6:.0f}', '' if not np.isfinite(o['edepth']) else f"{(ex - o['depth']) / o['edepth']:.1f}",
                    f"{o['scan_max']:.1f}" if 'scan_max' in o else '', f"{o['scan_off']:+.0f}" if 'scan_off' in o else '', f"{o['scan_near']:.1f}" if 'scan_near' in o else '', o.get('note', '')])
with open(os.path.join(args.out, 'transit_depths_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: transit depth and mid-transit time per TESS sector (scripts/toi540_analysis.py), {DATE}',
                  f'depth_own/pdc/sap/fast_ppm = peak depth of the common shape at the predicted times (own TPF photometry r = 2.0 px; SPOC 2-min PDCSAP; SAP; SPOC 20-s PDCSAP), quality 0, 0.3-d biweight; errors scaled by beta',
                  f'tc_btjd, tc_err_min = per-sector mid-transit time of the stacked transits (own photometry, shape fixed) and its error; oc_min = O - C against T0 = BTJD {T0:.6f} (cycle {Eref:.0f}), P = {P:.8f} d (sectors 4-32 only)',
                  'n_transits = predicted transits with data; year = mean epoch; sector 98: no mid-transit time (no transit detected)']))
    w = csv.writer(fo); w.writerow(['sector', 'year', 'cycle', 'n_transits', 'depth_own_ppm', 'depth_own_err_ppm', 'depth_pdc_ppm', 'depth_pdc_err_ppm', 'depth_sap_ppm', 'depth_sap_err_ppm', 'depth_fast_ppm', 'depth_fast_err_ppm', 'beta_own', 'tc_btjd', 'tc_err_min', 'oc_min'])
    for s in PRE + [98]:
        o = PER[s]; rr = [r for r in R0['rows'] if r['sector'] == s]; fz = o.get('fast')
        w.writerow([s, YR[s], f"{o['own']['E']:.0f}", o['own']['ntr'], f"{o['own']['depth']*1e6:.0f}", f"{o['own']['edepth']*1e6:.0f}", f"{o['pdc']['depth']*1e6:.0f}", f"{o['pdc']['edepth']*1e6:.0f}",
                    f"{o['sap']['depth']*1e6:.0f}", f"{o['sap']['edepth']*1e6:.0f}", f"{fz['depth']*1e6:.0f}" if fz else '', f"{fz['edepth']*1e6:.0f}" if fz else '', f"{o['own']['beta']:.2f}",
                    f"{rr[0]['tc']:.6f}" if rr else '', f"{rr[0]['sig']*1440:.2f}" if rr else '', f"{rr[0]['oc_min']:+.2f}" if rr else ''])

# ---------------------------------------------------------------- injections into the sector-98 photometry
tr, fr = raw_lc(98); pk0 = 1 - bmodel(np.array([0.]), 0, P, *SH)[0]; kk = SH[0] * np.sqrt(REF / pk0)
def inject(off_min):
    t2, f2, s2 = detrend(tr, fr * bmodel(tr, tc98 + off_min / 1440, P, kk, SH[1], SH[2]), 0.3); A, eA = fit_depth(t2, f2, tc98 + off_min / 1440, P, SH, c98, s2, b98)
    offs = np.arange(-P * 720, P * 720, 2.0); a, snr = phase_scan(t2, f2, tc98, P, SH, c98, s2, b98, offs); i = np.argmax(snr)
    return A, eA, offs[i], snr[i], (t2, f2)
INJ = []
if not args.skip_injections:
    A, eA, o, snr, INJ_LC = inject(0.0); INJ.append((0.0, A, eA, o, snr))
    for off in np.random.default_rng(7).uniform(-P * 720, P * 720, 20):
        A, eA, o, snr, _ = inject(off); INJ.append((off, A, eA, o, snr))
    with open(os.path.join(args.out, 'injection_200322593.csv'), 'w') as fo:
        fo.write(hdr([f'TOI-540 b = TIC 200322593: transits of the reference depth ({REF*1e6:.0f} ppm peak) injected into the raw sector-98 own photometry, then detrended and measured as the data (scripts/toi540_analysis.py), {DATE}',
                      'offset_min = injected phase offset from the prediction (first row: at the prediction; others: 20 random phases, seed 7); depth at the injected phase; blind_* = highest-S/N phase of a full-phase scan']))
        w = csv.writer(fo); w.writerow(['offset_min', 'depth_ppm', 'depth_err_ppm', 'blind_offset_min', 'blind_snr', 'recovered_within_6min'])
        for off, A, eA, o, snr in INJ:
            d = ((o - off + P * 720) % (P * 1440)) - P * 720; w.writerow([f'{off:+.1f}', f'{A*1e6:.0f}', f'{eA*1e6:.0f}', f'{o:+.0f}', f'{snr:.1f}', int(abs(d) <= 6)])

# ---------------------------------------------------------------- folded light curves for the figures (2-min bins, +-90 min)
with open(os.path.join(args.out, 'folded_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: detrended own photometry (r = 2.0 px, quality 0, 0.3 d) folded on the sectors 4-32 ephemeris, 2-min bins (scripts/toi540_analysis.py), {DATE}',
                  'set = sector, or 98_injected = sector 98 with the reference transit injected at the prediction; dt_min = minutes from the predicted mid-transit; err = standard error of the bin',
                  f'model: common shape scaled to the reference peak depth {REF*1e6:.0f} ppm (k = {SH[0]:.4f} x sqrt(ref/peak), a/R* = {SH[1]:.2f}, b = {SH[2]:.3f})']))
    w = csv.writer(fo); w.writerow(['set', 'dt_min', 'flux', 'err', 'n', 'model'])
    sets = {str(s): lc(s)[:2] for s in PRE + [98]}
    if INJ: sets['98_injected'] = INJ_LC
    for name, (t, f) in sets.items():
        dt = (((t - T0) / P + 0.5) % 1 - 0.5) * P * 1440; m = np.abs(dt) < 90; b = np.round(dt[m] / 2) * 2
        for k in np.unique(b):
            mm = b == k; w.writerow([name, f'{k:.0f}', f'{f[m][mm].mean():.6f}', f'{f[m][mm].std() / np.sqrt(mm.sum()):.6f}', int(mm.sum()), f'{bmodel(np.array([k / 1440]), 0, P, kk, SH[1], SH[2])[0]:.6f}'])

# ---------------------------------------------------------------- impact parameter per epoch group (sectors 4-6 vs 31-32)
def bt_fit(src, ap='r2.0'):
    dat = []
    for gi, ss in enumerate(([4, 5, 6], [31, 32])):
        for s in ss:
            t, f, sg, cad = lc(s, src, 0.3, ap); m = np.abs(((t - T0) / P + 0.5) % 1 - 0.5) * P < 0.15; dat.append((t[m], f[m], np.full(m.sum(), sg), np.full(m.sum(), gi), beta_factor(t, f - 1, cad)))
    t = np.concatenate([d[0] for d in dat]); f = np.concatenate([d[1] for d in dat]); e = np.concatenate([d[2] for d in dat]); gid = np.concatenate([d[3] for d in dat]); beta = np.median([d[4] for d in dat])
    out = {}
    for kind in ('const', 'bfree', 'kbfree'):
        def unpack(x):
            k = [x[3], x[3]] if kind != 'kbfree' else [x[3], x[6]]; b = [x[4], x[4]] if kind == 'const' else [x[4], x[5]]; return x[0], x[1], x[2], k, b
        def res(x):
            tc, Pp, aR, k, b = unpack(x); mo = np.empty_like(t)
            for g in (0, 1): mo[gid == g] = bmodel(t[gid == g], tc, Pp, k[g], aR, b[g])
            return np.append((f - mo) / e, (aR - AR0) / EAR)
        nx = 5 + (kind != 'const') + (kind == 'kbfree'); bestr = None
        for b0 in (0.5, 0.7, 0.8, 0.9):
            x0 = [T0, P, AR0, 0.041] + [b0] * (1 + (kind != 'const')) + ([0.041] if kind == 'kbfree' else [])
            r = least_squares(res, x0, bounds=([T0 - 0.01, P - 2e-5, 5, 0.01] + [0] * (nx - 4 - (kind == 'kbfree')) + ([0.01] if kind == 'kbfree' else []),
                                               [T0 + 0.01, P + 2e-5, 30, 0.2] + [1.2] * (nx - 4 - (kind == 'kbfree')) + ([0.2] if kind == 'kbfree' else [])))
            if bestr is None or r.cost < bestr.cost: bestr = r
        cov = np.linalg.inv(bestr.jac.T @ bestr.jac) * beta**2; out[kind] = (2 * bestr.cost, bestr.x, cov)
    chi_c, chi_b = out['const'][0], out['bfree'][0]; x, cov = out['bfree'][1], out['bfree'][2]; d1 = (chi_c - chi_b) / beta**2; d2 = (chi_c - out['kbfree'][0]) / beta**2
    rate = (x[5] - x[4]) / (2020.88 - 2018.92); erate = np.sqrt(cov[4, 4] + cov[5, 5] - 2 * cov[4, 5]) / (2020.88 - 2018.92); L = (2025.93 - 2020.88) / (2020.88 - 2018.92)
    b98x = x[5] + L * (x[5] - x[4]); eb98 = np.sqrt(((1 + L))**2 * cov[5, 5] + L**2 * cov[4, 4] - 2 * L * (1 + L) * cov[4, 5])
    bs = np.linspace(0.5, 1.1, 601); dd = np.array([1 - bmodel(np.array([0.]), 0, P, x[3], x[2], bb)[0] for bb in bs])
    r = dict(src=src if src != 'own' else f'own {ap}', beta=beta, chi2_const=chi_c, chi2_bfree=chi_b, chi2_kbfree=out['kbfree'][0], b_const=out['const'][1][4], b1=x[4], eb1=np.sqrt(cov[4, 4]), b2=x[5], eb2=np.sqrt(cov[5, 5]),
             k=x[3], aR=x[2], d1=d1, p1=CHI2.sf(d1, 1), sig1=norm.isf(CHI2.sf(d1, 1) / 2), p2=CHI2.sf(d2, 2), rate=rate, erate=erate, b98=b98x, eb98=eb98,
             pgone=1 - norm.cdf((1 + x[3] - b98x) / eb98), blim=bs[np.argmax(dd < 3 * o98['edepth'])])
    print(f"b(t) {r['src']}: b = {r['b1']:.3f} -> {r['b2']:.3f}, db/dt = {rate:+.4f} +- {erate:.4f}/yr, p = {r['p1']:.3f} ({r['sig1']:.1f} sigma); b(2025.93) = {b98x:.3f} +- {eb98:.3f}", flush=True)
    return r


BT = [bt_fit('own', 'r2.0'), bt_fit('own', 'r1.5'), bt_fit('own', 'spoc'), bt_fit('pdc'), bt_fit('sap')]
with open(os.path.join(args.out, 'bt_fits_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: joint fits of sectors 4-6 (2018.92) and 31-32 (2020.88) with a common or a per-group impact parameter (scripts/toi540_analysis.py), {DATE}',
                  'a/R* shared with the stellar-density prior; k shared (bfree) or per group (kbfree); chi2 with per-sector white-noise errors; delta chi2 divided by beta^2 (median red-noise factor) for p_b (1 dof) and p_bk (2 dof)',
                  'db_dt = (b2 - b1)/1.96 yr; b98_linear = linear extrapolation to 2025.93; p_no_transit = probability that b98 > 1 + k for a Gaussian error; b_needed = smallest b whose peak depth is below the sector-98 3-sigma limit (own photometry)']))
    w = csv.writer(fo); w.writerow(['source', 'beta', 'chi2_const', 'chi2_bfree', 'chi2_kbfree', 'b_const', 'b_2018.9', 'b_2018.9_err', 'b_2020.9', 'b_2020.9_err', 'k', 'aR', 'p_b', 'sigma_b', 'p_bk', 'db_dt_per_yr', 'db_dt_err', 'b98_linear', 'b98_err', 'p_no_transit', 'b_needed'])
    for r in BT:
        w.writerow([r['src'], f"{r['beta']:.2f}", f"{r['chi2_const']:.2f}", f"{r['chi2_bfree']:.2f}", f"{r['chi2_kbfree']:.2f}", f"{r['b_const']:.3f}", f"{r['b1']:.3f}", f"{r['eb1']:.3f}", f"{r['b2']:.3f}", f"{r['eb2']:.3f}",
                    f"{r['k']:.4f}", f"{r['aR']:.2f}", f"{r['p1']:.3f}", f"{r['sig1']:.1f}", f"{r['p2']:.3f}", f"{r['rate']:+.4f}", f"{r['erate']:.4f}", f"{r['b98']:.3f}", f"{r['eb98']:.3f}", f"{r['pgone']:.2f}", f"{r['blim']:.3f}"])

# ---------------------------------------------------------------- order-of-magnitude perturber masses for nodal precession
Gy = 3.15576e7; M = c['star']['M'] * MSUN; R = c['star']['R'] * RSUN; aR = SH[1]; ab = aR * R; nb = 2 * np.pi / (P * 86400) * Gy
bneed = BT[0]['blim']; b20 = BT[0]['b2']; rate = (bneed - b20) / aR / (2025.93 - 2020.88)
with open(os.path.join(args.out, 'perturber_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: order-of-magnitude companion masses for a nodal precession that moves b from {b20:.3f} (2020.88) to {bneed:.3f} (2025.93) (scripts/toi540_analysis.py), {DATE}',
                  f'required |di/dt| = {np.degrees(rate):.3f} deg/yr; secular rate |di/dt| <= (3/8) n_b (m_c/M*) (a_b/a_c)^3 (circular orbits, a_c >> a_b, mutual inclination 45 deg); M* = {c["star"]["M"]} Msun, a_b = {ab/1.496e11:.4f} au',
                  'K_ms = radial-velocity semi-amplitude of the star for that companion with sin i_c = 0.7; values are lower limits on m_c for the assumptions stated and are not a model of the system']))
    w = csv.writer(fo); w.writerow(['P_c_d', 'a_c_au', 'm_c_min_earth', 'm_c_min_jup', 'K_ms'])
    for Pc in (2.5, 3.5, 5, 7.5, 10, 20, 30, 60, 100, 365):
        ac = (G * M * (Pc * 86400)**2 / (4 * np.pi**2))**(1 / 3); mc = rate / (0.375 * nb * (ab / ac)**3) * M
        K = (2 * np.pi * G / (Pc * 86400))**(1 / 3) * mc / (M + mc)**(2 / 3) * 0.7
        w.writerow([Pc, f'{ac/1.496e11:.4f}', f'{mc/5.972e24:.1f}', f'{mc/1.898e27:.3f}', f'{K:.1f}'])
    Om = 2 * np.pi / (17.43 * 3600); J2 = 0.29 / 3 * Om**2 * R**3 / (G * M)
    fo.write(f'# stellar oblateness (rotation 17.43 h, Love number 0.29): J2 = {J2:.2e}, maximum |db/dt| = {aR * 1.5 * nb * J2 / aR**2 / 2:.4f} per yr (required {(bneed - b20) / 5.05:.3f})\n')

# ---------------------------------------------------------------- ground-based windows (geocentric UTC, as predict_transits.py)
COORD = SkyCoord(c['ra'] * u.deg, c['dec'] * u.deg); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
def geo(bjd):
    t = Time(bjd, format='jd', scale='tdb', location=GEO)
    for _ in range(2): t = Time(bjd, format='jd', scale='tdb', location=GEO) - Time(t, location=GEO).light_travel_time(COORD, kind='barycentric')
    return t.utc
SITES = {'CTIO / LCO-LSC (Chile)': EarthLocation(lat=-30.167 * u.deg, lon=-70.805 * u.deg, height=2200 * u.m), 'SAAO / LCO-CPT (South Africa)': EarthLocation(lat=-32.380 * u.deg, lon=20.810 * u.deg, height=1760 * u.m),
         'SSO / LCO-COJ (Australia)': EarthLocation(lat=-31.273 * u.deg, lon=149.071 * u.deg, height=1165 * u.m)}
t1 = Time(args.start).jd; t2 = t1 + args.months * 30.5; Es = np.arange(np.ceil((t1 - 2457000 - T0) / P) + Eref, np.floor((t2 - 2457000 - T0) / P) + Eref + 1)
yr = 2016.0 + (t1 - 2457388.5) / 365.25; star = SkyCoord((c['ra'] + c['pmra'] * (yr - 2016) / 3.6e6 / np.cos(np.radians(c['dec']))) * u.deg, (c['dec'] + c['pmdec'] * (yr - 2016) / 3.6e6) * u.deg)
with open(os.path.join(args.out, 'ground_windows_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: predicted transits {args.start} + {args.months:g} months that are fully observable from three southern sites (scripts/toi540_analysis.py), {DATE}',
                  f'ephemeris: own photometry, sectors 4-32 (T0 = BTJD {T0:.6f} at cycle {Eref:.0f}, P = {P:.8f} d, full covariance); tc_utc = geocentric UTC (barycentric light-travel time removed)',
                  f'window = tc +- (T14/2 + 3 sigma + 30 min), T14 = {R0["T14"]:.0f} min; observable = target altitude > 30 deg and Sun altitude < -12 deg over the whole window; min_alt_deg = lowest target altitude in the window']))
    w = csv.writer(fo); w.writerow(['site', 'cycle', 'tc_bjd_tdb', 'tc_utc', 'tc_err_min', 'window_start_utc', 'window_end_utc', 'min_alt_deg'])
    for name, loc in SITES.items():
        for E in Es:
            tc, s = pred(E); tc += 2457000; half = (R0['T14'] / 2 + 3 * s * 1440 + 30) / 1440; tu = geo(tc); tt = tu + np.linspace(-half, half, 9) * u.day
            fr_ = AltAz(obstime=tt, location=loc); alt = star.transform_to(fr_).alt.deg; sun = get_sun(tt).transform_to(fr_).alt.deg
            if alt.min() > 30 and sun.max() < -12:
                w.writerow([name, f'{E:.0f}', f'{tc:.5f}', tu.iso[:16], f'{s*1440:.1f}', (tu - half * u.day).iso[:16], (tu + half * u.day).iso[:16], f'{alt.min():.0f}'])
print('tables written to', args.out)
