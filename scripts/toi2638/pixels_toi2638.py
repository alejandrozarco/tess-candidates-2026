"""TOI-2638.01 (TIC 458342769): pixel-level tests on the SPOC 2-min target pixel files (README section 6.9).

usage: python scripts/toi2638/pixels_toi2638.py [--work ./work] [--out data]
input:  the target pixel files of download_toi2638.py --tpf in <work>/458342769/dl/tpf/; data/neighbours_458342769.csv
        (checks_toi2638.py), data/ephemeris_458342769.csv and data/jointfit_458342769.csv (fit_toi2638.py), read from --out
output: <out>/pixels_458342769.csv (format in the file header)
Per sector, cadences with QUALITY = 0 and all pixels finite:
 1. scene model of the out-of-transit median image (phase more than 0.06 d from mid-transit): circular Gaussian PSFs
    at the Gaia DR3 positions of the sources with G < 16.5 (propagated with their proper motions to the sector median
    time; WCS of the file), common width, free fluxes for the target and the 17" neighbour (Gaia DR3
    5353676910584550272), the other sources scaled from their G magnitude with one free factor, a constant background
    and a common position offset; weights 1 / sqrt(|image| + 100).
 2. per-pixel transit amplitude: regression of each pixel's flux on the common transit shape (exposure-averaged)
    with a quadratic baseline per transit (+-0.30 d windows); the amplitude image is fitted with the target PSF only,
    the neighbour PSF only, and both (weighted by the per-pixel amplitude errors).
 3. circular apertures of radius 1, 1.5, 2 and 3 px centred on the target: transit depth from the same regression
    on the aperture sum (background of the scene model subtracted), divided by the target's flux fraction in the
    aperture from the scene model ('undiluted' depth).
 4. location of the 7.31 c/d signal: per-pixel sinusoid at the frequency of the target's light curve in that sector
    (best within 7.311 +- 0.01 c/d; first harmonic, quadratic baseline, out of transit), projected on the phase of
    the brightest-amplitude pixel; centroid of the positive part compared with the target and neighbour positions.
Requires numpy, scipy, astropy, pandas, batman-package.
"""
import argparse, csv, glob, os, sys
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import GAIA_TARGET, GAIA_NB17, Model, peak_depth, hdr, dldir, DATE

ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--work', default='./work'); ap.add_argument('--out', default='data')
args = ap.parse_args()


def rows(fn):
    return [r for r in csv.reader(l for l in open(os.path.join(args.out, fn)) if not l.startswith('#'))]


eph = [r for r in rows('ephemeris_458342769.csv')[1:] if r[0] == 'ephemeris'][0]; T0, P = float(eph[10]), float(eph[12])
J = {}
for r in rows('jointfit_458342769.csv')[1:]:
    J.setdefault((r[0], r[1]), []).append(float(r[6]))
k0, b0, a0 = J[('SPOC PDCSAP + QLP sector 10', 'null')]; D0 = peak_depth(k0, b0, a0, P)
g = pd.read_csv(os.path.join(args.out, 'neighbours_458342769.csv'), comment='#'); g = g[g.G < 16.5].reset_index(drop=True)
OUTR = []
for fn in sorted(glob.glob(os.path.join(dldir(args.work), 'tpf', '*_tp.fits'))):
    with fits.open(fn) as h:
        d = h[1].data; t = np.array(d['TIME'], float); fl = np.array(d['FLUX'], float); q = np.array(d['QUALITY']); w = WCS(h[2].header); sec = int(h[0].header['SECTOR'])
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    ny, nx = fl.shape[1:]; yy, xx = np.mgrid[:ny, :nx]; cad = float(np.median(np.diff(t)))
    yrs = 2016.0 + (np.median(t) + 2457000 - 2457389.0) / 365.25
    dra = g.pmra.fillna(0) * (yrs - 2016) / 3.6e6 / np.cos(np.radians(g.dec)); dde = g.pmdec.fillna(0) * (yrs - 2016) / 3.6e6
    px, py = w.all_world2pix(g.ra + dra, g.dec + dde, 0)
    it = int(np.where(g.source_id == GAIA_TARGET)[0][0]); inn = int(np.where(g.source_id == GAIA_NB17)[0][0])
    inside = (px > -3) & (px < nx + 2) & (py > -3) & (py < ny + 2); others = [i for i in np.where(inside)[0] if i not in (it, inn)]
    dt = (t - T0 + 0.5 * P) % P - 0.5 * P; oot = np.abs(dt) > 0.06; med = np.median(fl[oot], 0)
    psf = lambda x0, y0, s: np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2)) / (2 * np.pi * s ** 2)
    sc_ = {i: 10 ** (-0.4 * (g.G[i] - g.G[it])) for i in others}

    def scene(p):
        s, at, an, ao, bk, dx, dy = p
        m = at * psf(px[it] + dx, py[it] + dy, s) + an * psf(px[inn] + dx, py[inn] + dy, s) + bk
        for i in others:
            m = m + ao * sc_[i] * psf(px[i] + dx, py[i] + dy, s)
        return m
    ft = med.sum()
    r = least_squares(lambda p: ((scene(p) - med) / np.sqrt(np.abs(med) + 100)).ravel(), [0.8, ft * 0.6, ft * 0.2, ft * 0.6, 0, 0, 0],
                      bounds=([0.3, 0, 0, 0, -1e4, -1, -1], [3, 1e7, 1e7, 1e7, 1e4, 1, 1]))
    s, at, an, ao, bk, dx, dy = r.x
    Pt, Pn = psf(px[it] + dx, py[it] + dy, s), psf(px[inn] + dx, py[inn] + dy, s)
    Po = sum(ao * sc_[i] * psf(px[i] + dx, py[i] + dy, s) for i in others) if others else 0 * Pt
    tm = Model(t, cad, P, T0).depth(k0, b0, a0) / D0           # unit-depth common shape
    E = np.round((t - T0) / P); cols = []
    for n in np.unique(E):
        m = (E == n) & (np.abs(dt) < 0.3)
        if m.sum() >= 200 and (np.abs(dt[m]) < 0.02).sum() >= 10:
            cols.append(m)
    sel = np.any(cols, 0); X = [tm[sel]]
    for m in cols:
        mm = m[sel] * 1.0; dd = dt[sel]; X += [mm, mm * dd, mm * dd ** 2]
    X = np.array(X).T; Y = fl[sel].reshape(sel.sum(), -1)
    c, *_ = np.linalg.lstsq(X, Y, rcond=None); rr = Y - X @ c; cov0 = np.linalg.inv(X.T @ X)[0, 0]
    dF = -c[0].reshape(ny, nx); eF = (np.std(rr, 0) * np.sqrt(cov0)).reshape(ny, nx)
    Mx = np.vstack([Pt.ravel(), Pn.ravel()]).T; Wt = 1 / eF.ravel()
    cc, *_ = np.linalg.lstsq(Mx * Wt[:, None], dF.ravel() * Wt, rcond=None); cv = np.linalg.inv((Mx * Wt[:, None]).T @ (Mx * Wt[:, None])); ecc = np.sqrt(np.diag(cv))
    ch = {}
    for lab, PP in (('target', Pt), ('nb', Pn)):
        a1 = np.sum(PP.ravel() * dF.ravel() * Wt ** 2) / np.sum(PP.ravel() ** 2 * Wt ** 2); ch[lab] = float(np.sum(((dF.ravel() - a1 * PP.ravel()) * Wt) ** 2))
    rad = np.hypot(xx - px[it] - dx, yy - py[it] - dy); apr = {}
    for R_ in (1.0, 1.5, 2.0, 3.0):
        msk = rad <= R_; lc = fl[:, msk].sum(1) - bk * msk.sum(); ya = lc[sel] / np.median(lc[oot])
        ca, *_ = np.linalg.lstsq(X, ya, rcond=None); ea = np.std(ya - X @ ca) * np.sqrt(cov0)
        tf = (at * Pt[msk]).sum() / ((at * Pt + an * Pn + Po)[msk].sum())
        apr[R_] = (-ca[0] * 1e6, ea * 1e6, tf, msk.sum())
    # 7.31 c/d location
    hdr0 = None
    lcf = glob.glob(os.path.join(dldir(args.work), 'lc', f'tess*-s{sec:04d}-0000000458342769-*-s_lc.fits'))[0]
    dl_ = fits.open(lcf)[1].data; tl = np.array(dl_['TIME'], float); fl_l = np.array(dl_['PDCSAP_FLUX'], float); okl = np.isfinite(tl) & np.isfinite(fl_l) & (dl_['QUALITY'] == 0)
    tl, fl_l = tl[okl], fl_l[okl] / np.nanmedian(fl_l[okl]) - 1; ml = np.abs((tl - T0 + 0.5 * P) % P - 0.5 * P) > 0.084; best = None
    for f0 in np.linspace(7.301, 7.321, 201):
        Xl = np.c_[np.ones(ml.sum()), np.sin(2 * np.pi * f0 * tl[ml]), np.cos(2 * np.pi * f0 * tl[ml])]; cl, *_ = np.linalg.lstsq(Xl, fl_l[ml], rcond=None)
        if best is None or np.hypot(cl[1], cl[2]) > best[0]:
            best = (np.hypot(cl[1], cl[2]), f0)
    f0 = best[1]; mo = oot & (np.abs(dt) > 0.084); to = t[mo]
    Xp = np.c_[np.ones(len(to)), to - to.mean(), (to - to.mean()) ** 2, np.sin(2 * np.pi * f0 * to), np.cos(2 * np.pi * f0 * to)]
    cp, *_ = np.linalg.lstsq(Xp, fl[mo].reshape(len(to), -1), rcond=None); amp = np.hypot(cp[3], cp[4]); j = np.argmax(amp); ph0 = np.arctan2(cp[4][j], cp[3][j])
    proj = (cp[3] * np.cos(ph0) + cp[4] * np.sin(ph0)).reshape(ny, nx).clip(0); cy, cx = np.sum(proj * yy) / proj.sum(), np.sum(proj * xx) / proj.sum()
    dtg, dnb = np.hypot(cx - px[it] - dx, cy - py[it] - dy), np.hypot(cx - px[inn] - dx, cy - py[inn] - dy)
    o = dict(sector=sec, npix=nx * ny, psf_sigma=s, flux_target=at, nb_over_target=an / at, nb_over_target_gaia=10 ** (-0.4 * (g.G[inn] - g.G[it])),
             chi2_t=ch['target'], chi2_n=ch['nb'], a_t=cc[0], ea_t=ecc[0], a_n=cc[1], ea_n=ecc[1], depth_on_target=cc[0] / at * 1e6, edepth_on_target=ecc[0] / at * 1e6,
             depth_on_nb=cc[1] / an * 1e6, f731=f0, a731=best[0] * 1e6, d731_t=dtg, d731_n=dnb, sep_tn=np.hypot(px[inn] - px[it], py[inn] - py[it]), apr=apr)
    OUTR.append(o)
    print(f"S{sec}: sigma {s:.2f} px, nb/target {an / at:.3f} (Gaia {o['nb_over_target_gaia']:.3f}); chi2 nb-only - target-only {ch['nb'] - ch['target']:.1f}; "
          f"two-source: target {cc[0] / ecc[0]:.1f} sigma, nb {cc[1] / ecc[1]:+.1f} sigma; depth if on target {o['depth_on_target']:.0f} +- {o['edepth_on_target']:.0f}; "
          f"apertures {[(k_, round(v[0]), round(v[0] / v[2])) for k_, v in apr.items()]}; 7.31 c/d centroid {dtg:.2f} px from target, {dnb:.2f} px from neighbour", flush=True)
with open(os.path.join(args.out, 'pixels_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: pixel-level tests on the SPOC 2-min target pixel files; scripts/toi2638/pixels_toi2638.py, {DATE}',
            'scene: Gaussian PSFs at the Gaia DR3 positions (G < 16.5, proper motions applied); psf_sigma_px; flux_target_e_s = fitted target flux; nb_over_target = fitted flux ratio '
            'of the 17" neighbour to the target, nb_over_target_gaia = the Gaia G-band ratio',
            'transit amplitude image (e-/s per pixel, common transit shape of jointfit_458342769.csv, quadratic baseline per transit): chi2_target_only, chi2_nb_only = chi2 of fits '
            'with one PSF; amp_target, amp_nb = two-PSF fit (e-/s) with errors; depth_if_on_target_ppm = amp_target / flux_target; depth_if_on_nb_ppm = amp_nb / fitted neighbour flux',
            'ap<r>_raw_ppm, ap<r>_err_ppm = depth in a circular aperture of r px on the target; ap<r>_target_fraction = target share of the aperture flux in the scene model; '
            'ap<r>_undiluted_ppm = raw / target_fraction',
            'f731_c_d, a731_ppm = frequency and amplitude of the 7.31 c/d signal in the PDCSAP light curve of the sector; d731_target_px, d731_nb_px = distance of the centroid of its '
            'phase-projected amplitude image from the target and from the neighbour; sep_target_nb_px = their separation',
            'sector,n_pixels,psf_sigma_px,flux_target_e_s,nb_over_target,nb_over_target_gaia,chi2_target_only,chi2_nb_only,amp_target,amp_target_err,amp_nb,amp_nb_err,'
            'depth_if_on_target_ppm,depth_if_on_target_err_ppm,depth_if_on_nb_ppm,' + ','.join(f'ap{r_:g}_raw_ppm,ap{r_:g}_err_ppm,ap{r_:g}_target_fraction,ap{r_:g}_undiluted_ppm,ap{r_:g}_npix' for r_ in (1.0, 1.5, 2.0, 3.0))
            + ',f731_c_d,a731_ppm,d731_target_px,d731_nb_px,sep_target_nb_px'])
    for o in OUTR:
        aps = ','.join(f'{v[0]:.0f},{v[1]:.0f},{v[2]:.3f},{v[0] / v[2]:.0f},{v[3]}' for v in o['apr'].values())
        f.write(f"{o['sector']},{o['npix']},{o['psf_sigma']:.3f},{o['flux_target']:.0f},{o['nb_over_target']:.3f},{o['nb_over_target_gaia']:.3f},{o['chi2_t']:.1f},{o['chi2_n']:.1f},"
                f"{o['a_t']:.3f},{o['ea_t']:.3f},{o['a_n']:.3f},{o['ea_n']:.3f},{o['depth_on_target']:.0f},{o['edepth_on_target']:.0f},{o['depth_on_nb']:.0f},{aps},"
                f"{o['f731']:.4f},{o['a731']:.0f},{o['d731_t']:.2f},{o['d731_n']:.2f},{o['sep_tn']:.2f}\n")
