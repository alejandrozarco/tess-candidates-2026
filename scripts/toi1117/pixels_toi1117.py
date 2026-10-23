"""TOI-1117 c (TIC 295541511), README section 6.10: pixel-level location of the 4.5787-d signal.

usage: python scripts/toi1117/pixels_toi1117.py [--work ./work]
input:  <work>/295541511/dl/tpf/ and dl/gaia_150as.csv, dl/tic_150as.csv (download_toi1117.py)
output: data/pixels_295541511.csv     per pixel set: transits used, difference-image S/N, target depth, free-position fit
        data/sources_295541511.csv    per Gaia DR3 source within 2' (G < 20.5): chi2 of the source hypothesis minus that of the
                                      target, summed over the six independent pixel sets, and the eclipse depth the source would need
        data/apertures_295541511.csv  contamination-corrected depth in circular apertures of radius 1-4 px
Method, per pixel set (SPOC 2-min target pixel files of sectors 39, 100-103; TESS-SPOC FFI target pixel files of sectors 13
and 39, the latter a duplicate of the sector-39 SPOC set): cadences with QUALITY = 0 and finite pixels, planet b masked; for every
transit of c (P, T0 and box D as in fit_toi1117.py) a linear baseline per pixel fitted to 0.75 D to D/2 + 0.5 d from mid-transit;
the in-transit deficit image is the mean over transits of (baseline - flux) in the box, with its error from the baseline
residuals. PSF: a pixel-integrated double Gaussian whose widths and weight, a common position offset, a background and the fluxes
of the Gaia DR3 stars brighter than G = 14 are fitted to the mean out-of-transit image, with the stars of G = 14-17 at fluxes
scaled from G (Gaia DR3 positions through the file WCS). Hypotheses: deficit image = amplitude x PSF at the position of one
source + constant; chi2 relative to the target hypothesis; implied depth = amplitude / modelled flux of that source. The
free-position fit has a jackknife error over transits. Apertures: summed deficit / summed out-of-transit flux, divided by the
modelled target fraction of the aperture flux.
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.special import erf
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); args = p.parse_args()
W = args.work; DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
P, T0, D = PC, TBOX, DBOX
gaia = pd.read_csv(os.path.join(dldir(W), 'gaia_150as.csv')); gaia = gaia[np.isfinite(gaia.Gmag)].reset_index(drop=True).drop(columns=['_r'])
from astropy.coordinates import SkyCoord
gaia['sep'] = SkyCoord(
    gaia.RA_ICRS.to_numpy(), gaia.DE_ICRS.to_numpy(), unit='deg'
).separation(SkyCoord(RA, DEC, unit='deg')).arcsec
it = int(np.where(gaia.Source == GAIA)[0][0]); Gt = gaia.Gmag[it]
tic = pd.read_csv(os.path.join(dldir(W), 'tic_150as.csv'))
TPF = os.path.join(dldir(W, 'tpf'))
SRC = {'TSPOC13': glob.glob(f'{TPF}/hlsp_tess-spoc_*-s0013_tess_v1_tp.fits')[0], 'TSPOC39': glob.glob(f'{TPF}/hlsp_tess-spoc_*-s0039_tess_v1_tp.fits')[0]}
for s in SPOC2:
    SRC[f'SPOC{s}'] = glob.glob(f'{TPF}/tess*-s{s:04d}-0000000295541511-*-s_tp.fits')[0]
INDEP = ['TSPOC13'] + [f'SPOC{s}' for s in SPOC2]


def pix_gauss(xx, yy, x0, y0, s):
    gx = 0.5 * (erf((xx + 0.5 - x0) / (np.sqrt(2) * s)) - erf((xx - 0.5 - x0) / (np.sqrt(2) * s)))
    gy = 0.5 * (erf((yy + 0.5 - y0) / (np.sqrt(2) * s)) - erf((yy - 0.5 - y0) / (np.sqrt(2) * s)))
    return gx * gy


def psf(xx, yy, x0, y0, th):
    s1, s2, eta = th; return eta * pix_gauss(xx, yy, x0, y0, s1) + (1 - eta) * pix_gauss(xx, yy, x0, y0, s2)


def load_pix(key):
    with fits.open(SRC[key]) as h:
        d = h[1].data; t = np.array(d['TIME'], float); F = np.array(d['FLUX'], float); q = np.array(d['QUALITY']); w = WCS(h[2].header)
    good = np.isfinite(t) & (q == 0) & np.isfinite(F).all(axis=(1, 2)); t, F = t[good], F[good]
    k = mask_b(t); return t[k], F[k], w


def stack(t, F):
    n = np.round((t - T0) / P); ph = t - (T0 + n * P)
    inn = np.abs(ph) < D / 2; near = (np.abs(ph) > 0.75 * D) & (np.abs(ph) < D / 2 + 0.5); cad = np.median(np.diff(t)); per = []
    for e in np.unique(n[inn]):
        mi = inn & (n == e); mo = near & (n == e)
        if mi.sum() < max(2, 0.6 * D / cad) or mo.sum() < 8:
            continue
        X = np.vstack([np.ones(mo.sum()), ph[mo]]).T; Xi = np.vstack([np.ones(mi.sum()), ph[mi]]).T
        Y = F[mo].reshape(mo.sum(), -1); cf = np.linalg.lstsq(X, Y, rcond=None)[0]
        per.append(((Xi @ cf - F[mi].reshape(mi.sum(), -1)).mean(0), (Y - X @ cf).std(0) / np.sqrt(mi.sum()), (X @ cf).mean(0)))
    return per


def fit_shape(O, xx, yy, stars):
    bright = stars.Gmag < 14; fb = stars[bright]; ff = stars[~bright & (stars.Gmag < 17)]
    def model(p):
        dx, dy, s1, s2, eta, bg = p[:6]; fl = p[6:]; th = (s1, s2, eta); M = np.full(O.shape, bg)
        for j, r in enumerate(fb.itertuples()):
            M += fl[j] * psf(xx, yy, r.px + dx, r.py + dy, th)
        for g, r in zip(fl[0] * 10 ** (-0.4 * (ff.Gmag.values - Gt)), ff.itertuples()):
            M += g * psf(xx, yy, r.px + dx, r.py + dy, th)
        return M
    p0 = [0, 0, 0.7, 2.0, 0.8, 0] + [O.sum() * 10 ** (-0.4 * (g - Gt)) for g in fb.Gmag]
    lb = [-2, -2, 0.3, 0.8, 0.05, -np.inf] + [0] * len(fb); ub = [2, 2, 2.5, 6, 1, np.inf] + [np.inf] * len(fb)
    r = least_squares(lambda q: ((model(q) - O) / np.sqrt(np.clip(O, 1, None))).ravel(), p0, bounds=(lb, ub))
    return r.x, fb


def run(key):
    t, F, w = load_pix(key); ny, nx = F.shape[1:]; yy, xx = np.mgrid[:ny, :nx]; per = stack(t, F)
    Dimg = np.mean([q[0] for q in per], 0).reshape(ny, nx); E = (np.sqrt(np.sum([q[1] ** 2 for q in per], 0)) / len(per)).reshape(ny, nx)
    O = np.mean([q[2] for q in per], 0).reshape(ny, nx)
    px, py = w.world_to_pixel_values(gaia.RA_ICRS.values, gaia.DE_ICRS.values)
    st = gaia.assign(px=px, py=py); st = st[(px > -3) & (px < nx + 2) & (py > -3) & (py < ny + 2)].copy(); st = pd.concat([st[st.index == it], st[st.index != it]])
    pp, fb = fit_shape(O, xx, yy, st); dx, dy = pp[0], pp[1]; th = tuple(pp[2:5]); Ftar = pp[6]
    st['fx'] = st.px + dx; st['fy'] = st.py + dy; st['F'] = Ftar * 10 ** (-0.4 * (st.Gmag - Gt)); st.loc[st.index.isin(fb.index), 'F'] = pp[6:6 + len(fb)]
    tx, ty = float(st.fx.iloc[0]), float(st.fy.iloc[0])
    def hyp(x0, y0, Dm, Em):
        A = np.vstack([psf(xx, yy, x0, y0, th).ravel(), np.ones(nx * ny)]).T / Em.ravel()[:, None]
        cf = np.linalg.lstsq(A, Dm.ravel() / Em.ravel(), rcond=None)[0]; cov = np.linalg.inv(A.T @ A)
        return float(np.sum((A @ cf - Dm.ravel() / Em.ravel()) ** 2)), float(cf[0]), float(np.sqrt(cov[0, 0]))
    def free(Dm, Em):
        f = lambda q: (np.vstack([psf(xx, yy, q[0], q[1], th).ravel(), np.ones(nx * ny)]).T @ q[2:4] - Dm.ravel()) / Em.ravel()
        return least_squares(f, [tx, ty, Dm.sum(), 0], bounds=([tx - 4, ty - 4, -np.inf, -np.inf], [tx + 4, ty + 4, np.inf, np.inf])).x
    c0, a0, ea0 = hyp(tx, ty, Dimg, E); H = []
    for r in st[(st.sep < 120) & (st.Gmag < 20.5)].itertuples():
        chi2, a, ea = hyp(r.fx, r.fy, Dimg, E); H.append((key, int(r.Source), float(r.Gmag), float(r.sep), chi2 - c0, a / r.F, ea / r.F))
    fr = free(Dimg, E); off = float(np.hypot(fr[0] - tx, fr[1] - ty)); jk = []
    for k in range(len(per)):
        sel = [per[j] for j in range(len(per)) if j != k]
        jk.append(free(np.mean([q[0] for q in sel], 0).reshape(ny, nx), (np.sqrt(np.sum([q[1] ** 2 for q in sel], 0)) / len(sel)).reshape(ny, nx))[:2])
    jk = np.array(jk); nj = len(jk); sx, sy = np.sqrt(np.diag((nj - 1) * np.cov(jk.T, bias=True)))
    rr = np.hypot(xx - tx, yy - ty); aps = []
    for rad in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        m = rr <= rad
        dep = Dimg[m].sum() / O[m].sum(); err = np.sqrt((E[m] ** 2).sum()) / O[m].sum()
        frac = Ftar * psf(xx, yy, tx, ty, th)[m].sum() / sum(r.F * psf(xx, yy, r.fx, r.fy, th)[m].sum() for r in st.itertuples())
        aps.append((key, rad, int(m.sum()), dep * 1e6, err * 1e6, frac, dep / frac * 1e6, err / frac * 1e6))
    summ = (key, os.path.basename(SRC[key]), len(per), Dimg.sum() / np.sqrt((E ** 2).sum()), a0 / Ftar * 1e6, ea0 / Ftar * 1e6, off, sx, sy,
            off / np.hypot(sx, sy), th[0], th[1], th[2], int(key in INDEP))
    print(key, 'n', len(per), 'offset %.2f px (%.1f sigma)' % (off, off / np.hypot(sx, sy)), 'target depth %.0f ppm' % (a0 / Ftar * 1e6), flush=True)
    return summ, H, aps


S, Hs, As = [], [], []
for key in INDEP + ['TSPOC39']:
    a, b, c = run(key); S.append(a); Hs += b; As += c
write_csv(os.path.join(DATA, 'pixels_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): pixel-level tests per pixel set; scripts/toi1117/pixels_toi1117.py, {DATE}',
           f'ephemeris P = {PC} d, T0 = BTJD {TBOX:.4f}, box D = {DBOX * 24:.3f} h; planet b masked; PSF = pixel-integrated double Gaussian (sigma1, sigma2 in px, weight eta of sigma1) fitted to the out-of-transit image (not the TESS PRF)',
           'diff_snr = summed deficit / its error; target_depth_ppm = PSF amplitude at the target position / modelled target flux; free_offset_px = distance of the free-position PSF fit from the target, jk_sx/jk_sy = jackknife errors over transits, offset_sigma = offset / hypot(jk_sx, jk_sy)',
           'independent = 1 for the six pixel sets summed in sources_295541511.csv (TSPOC39 duplicates SPOC39)'],
          ['pixel_set', 'file', 'n_transits', 'diff_snr', 'target_depth_ppm', 'target_depth_err_ppm', 'free_offset_px', 'jk_sx_px', 'jk_sy_px', 'offset_sigma', 'psf_sigma1_px', 'psf_sigma2_px', 'psf_eta', 'independent'],
          [(a, b, c, f'{d:.1f}', f'{e:.0f}', f'{f:.0f}', f'{g:.2f}', f'{h:.2f}', f'{i:.2f}', f'{j:.2f}', f'{k:.2f}', f'{l:.2f}', f'{m:.2f}', n) for a, b, c, d, e, f, g, h, i, j, k, l, m, n in S])
H = pd.DataFrame(Hs, columns=['key', 'source', 'G', 'sep', 'dchi2', 'depth', 'depth_err']); H = H[H.key.isin(INDEP)]
rows = []
for src, g in H.groupby('source'):
    w = 1 / g.depth_err ** 2; jd = float(np.sum(w * g.depth) / np.sum(w)); je = float(1 / np.sqrt(np.sum(w)))
    tm = tic[tic.GAIA == src]
    if len(tm) == 0:
        gg = gaia[gaia.Source == src].iloc[0]; dd = np.hypot((tic.ra - gg.RA_ICRS) * np.cos(np.radians(DEC)), tic.dec - gg.DE_ICRS) * 3600; tm = tic[dd < 1.5]
    rows.append((src, '' if len(tm) == 0 else int(tm.ID.iloc[0]), g.G.iloc[0], g.sep.iloc[0], len(g), g.dchi2.sum(), jd, je, ';'.join(f'{k}:{v:.1f}' for k, v in zip(g.key, g.dchi2))))
rows.sort(key=lambda r: r[5])
write_csv(os.path.join(DATA, 'sources_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): source hypotheses for the in-transit deficit images, summed over the six independent pixel sets (TSPOC13, SPOC39, SPOC100-103); scripts/toi1117/pixels_toi1117.py, {DATE}',
           'every Gaia DR3 source within 120" with G < 20.5 that falls on the pixel sets; sep_arcsec = separation from the target position (Gaia DR3, epoch 2016.0)',
           'dchi2_sum = sum over pixel sets of chi2(source hypothesis) - chi2(target hypothesis) (positive: target preferred); implied_depth = weighted mean of amplitude / modelled flux of the source (fraction; 1 = total eclipse)',
           'tic = TIC v8.2 identifier matched by Gaia identifier or position (< 1.5"); per_set = dchi2 per pixel set; the target itself has dchi2_sum = 0'],
          ['gaia_dr3', 'tic', 'G', 'sep_arcsec', 'n_sets', 'dchi2_sum', 'implied_depth', 'implied_depth_err', 'per_set'],
          [(a, b, f'{c:.2f}', f'{d:.1f}', e, f'{f:.1f}', f'{g:.5f}', f'{h:.5f}', i) for a, b, c, d, e, f, g, h, i in rows])
A = pd.DataFrame(As, columns=['key', 'r', 'npix', 'dep', 'err', 'frac', 'cdep', 'cerr']); out = []
for r in A.itertuples():
    out.append((r.key, r.r, r.npix, f'{r.dep:.0f}', f'{r.err:.0f}', f'{r.frac:.3f}', f'{r.cdep:.0f}', f'{r.cerr:.0f}'))
for rad, g in A[A.key.isin(INDEP)].groupby('r'):
    w = 1 / g.cerr ** 2; out.append(('combined (six independent sets)', rad, '', '', '', '', f'{np.sum(w * g.cdep) / np.sum(w):.0f}', f'{1 / np.sqrt(np.sum(w)):.0f}'))
write_csv(os.path.join(DATA, 'apertures_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): transit depth in circular apertures centred on the target, from the deficit and out-of-transit images; scripts/toi1117/pixels_toi1117.py, {DATE}',
           'depth_ppm = summed deficit / summed out-of-transit flux; target_fraction = modelled target flux / modelled flux of all Gaia DR3 stars (G < 17 plus fitted brighter stars) in the aperture',
           'depth_corr_ppm = depth_ppm / target_fraction (the depth on the target if the signal is on it); combined rows: weighted mean over the six independent pixel sets'],
          ['pixel_set', 'radius_px', 'n_pix', 'depth_ppm', 'err_ppm', 'target_fraction', 'depth_corr_ppm', 'err_corr_ppm'], out)
