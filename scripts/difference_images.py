"""Per-epoch difference images from the TESScut cutouts: centroid offsets and a joint PSF fit on the listed sources.

usage: python scripts/difference_images.py <tic> --fit <fit json> [--work ./work]
input:  <work>/<tic>/cut/s<sector>.npz (tesscut_download.py), mid-transit times from <work>/<tic>/<fit json>
output: <work>/<tic>/D_s<sector>.npy = [difference image, per-pixel error, out-of-transit mean image] (21x21 each),
        <work>/<tic>/diffimg.json with the per-sector centroid offsets (arcsec E/N) and their errors.
Difference image = (quadratic-in-time out-of-transit model of every pixel, |dt| > 0.15 d, within 0.7 d) minus the
in-transit mean (|dt| < hw). A free 2-D Gaussian is fitted to the difference image and to the out-of-transit image;
their centroid difference is the offset. Its error is the scatter of the same measurement over 80 fake windows at
random out-of-transit times (the real difference image is added so that the signal-to-noise is the same). A joint
linear fit of Gaussian PSFs placed at the catalogue positions of the sources in targets.py (sigma from the
out-of-transit fit) gives the share of the dimming per source; single-source fits give delta chi2 between sources.
"""
import argparse, json, os
import numpy as np
from scipy.optimize import least_squares
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--fit', required=True); p.add_argument('--hw', type=float, default=None, help='in-transit half-window (d); default T14/2 - 0.03')
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic)
J = json.load(open(os.path.join(wd, args.fit))); hw = args.hw or max(c['T14h'] / 48 - 0.03, 0.05)
N = 21; yy, xx = np.indices((N, N)); src = c['sources']


def g2(p):
    a, x0, y0, s, cst = p; return a * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s * s)) + cst


def cen(img, x0, y0, rad=3.5):
    w = np.hypot(xx - x0, yy - y0) < rad
    return least_squares(lambda p: (g2(p) - img)[w], [max(img[int(round(y0)), int(round(x0))], 1), x0, y0, 0.9, 0],
                         bounds=([-np.inf, x0 - 3, y0 - 3, 0.4, -np.inf], [np.inf, x0 + 3, y0 + 3, 3, np.inf])).x


def dimg(tm, f, tc):
    m = np.abs(tm - tc) < 0.7; t = tm[m] - tc; F = f[m]; win = np.abs(t) < hw; oot = np.abs(t) > 0.15
    A = np.c_[np.ones(oot.sum()), t[oot], t[oot] ** 2]; cc = np.linalg.lstsq(A, F[oot].reshape(oot.sum(), -1), rcond=None)[0]
    pred = (np.c_[np.ones(win.sum()), t[win], t[win] ** 2] @ cc).reshape(-1, N, N); r = F[oot] - (A @ cc).reshape(-1, N, N)
    return (pred - F[win]).mean(0), F[oot].mean(0), r.std(0) / np.sqrt(win.sum())


rng = np.random.default_rng(2); OUT = {}
for s_, v in sorted(J.items(), key=lambda kv: int(kv[0])):
    s = int(s_); tc = v['tc']; z = np.load(os.path.join(wd, 'cut', f's{s}.npz')); tm = z['tm']; f = z['f'].astype(float); q = z['q']
    ok = np.isfinite(tm) & np.all(np.isfinite(f), axis=(1, 2)) & (q == 0); tm, f = tm[ok], f[ok]
    w = WCS(Header.fromstring(str(z['wcs']))); yr = 2016 + (tc + 2457000 - 2457389) / 365.25; pos = {}
    for k, (ra, de, pa, pd) in src.items():
        dra = pa * (yr - 2016) / 3.6e6 / np.cos(np.radians(de)); dde = pd * (yr - 2016) / 3.6e6
        pos[k] = [float(x) for x in w.world_to_pixel(SkyCoord(ra + dra, de + dde, unit='deg'))]
    tx, ty = pos['target']; D, O, E = dimg(tm, f, tc); pD = cen(D, tx, ty); pO = cen(O, tx, ty)
    offs = []
    for _ in range(80):
        tf = rng.uniform(tm.min() + 0.8, tm.max() - 0.8)
        if abs(tf - tc) < 0.8 or (np.abs(tm - tf) < hw).sum() < 3:
            continue
        try:
            Df, _, _ = dimg(tm, f, tf); offs.append(cen(Df + D, tx, ty)[1:3])
        except Exception:
            pass
    offs = np.array(offs); sx, sy = offs[:, 0].std(), offs[:, 1].std()
    sky = lambda x, y: (lambda cc: (cc.ra.deg, cc.dec.deg))(w.pixel_to_world(x, y))
    rD, dD = sky(pD[1], pD[2]); rO, dO = sky(pO[1], pO[2]); rT, dT = sky(tx, ty)
    dE = (rD - rO) * np.cos(np.radians(dT)) * 3600; dN = (dD - dO) * 3600; Jc = np.abs(np.linalg.det(w.pixel_scale_matrix)) ** 0.5 * 3600
    print(f'S{s} tc={tc:.3f}: difference-image centroid minus out-of-transit centroid: dE={dE:+.1f}" dN={dN:+.1f}" (+-{sx * Jc:.1f}, {sy * Jc:.1f}" from {len(offs)} fake windows); OOT centroid - target = ({(rO - rT) * np.cos(np.radians(dT)) * 3600:+.1f},{(dO - dT) * 3600:+.1f})"')
    sg = pO[3]; Ev = np.where(np.isfinite(E) & (E > 0), E, np.nanmedian(E))
    psf = lambda x0, y0: (lambda G: (G / G.sum()).ravel())(np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sg * sg)))
    M = np.array([psf(*pos[k]) for k in src] + [np.ones(N * N)]).T; Mw = M / Ev.ravel()[:, None]
    b = np.linalg.lstsq(Mw, (D / Ev).ravel(), rcond=None)[0]; e = np.sqrt(np.diag(np.linalg.inv(Mw.T @ Mw))); bo = np.linalg.lstsq(M, O.ravel(), rcond=None)[0]
    for i, k in enumerate(src):
        print(f'      joint PSF fit: {k:52s} dflux {b[i]:+8.1f} +- {e[i]:5.1f} ({b[i] / e[i]:+.1f} sigma); OOT flux {bo[i]:9.0f}; implied depth {b[i] / bo[i] * 100 if bo[i] > 0 else np.nan:7.3f} %')
    ch = {}
    for i, k in enumerate(src):
        Ms = np.c_[M[:, i], M[:, -1]]; bs = np.linalg.lstsq(Ms / Ev.ravel()[:, None], (D / Ev).ravel(), rcond=None)[0]; ch[k] = (((D.ravel() - Ms @ bs) / Ev.ravel()) ** 2).sum()
    c0 = min(ch.values()); print('      single-source delta chi2:', {k: round(v - c0, 1) for k, v in ch.items()})
    OUT[s] = dict(dE=dE, dN=dN, eE=sx * Jc, eN=sy * Jc); np.save(os.path.join(wd, f'D_s{s}.npy'), np.array([D, E, O]))
Ea = np.array([[v['dE'], v['dN'], v['eE'], v['eN']] for v in OUT.values()]); wE = 1 / Ea[:, 2] ** 2; wN = 1 / Ea[:, 3] ** 2
print(f'weighted mean centroid offset: E {np.sum(Ea[:, 0] * wE) / wE.sum():+.2f} +- {wE.sum() ** -0.5:.2f}"  N {np.sum(Ea[:, 1] * wN) / wN.sum():+.2f} +- {wN.sum() ** -0.5:.2f}"; chi2 vs 0: {np.sum(Ea[:, 0] ** 2 * wE) + np.sum(Ea[:, 1] ** 2 * wN):.1f}/{2 * len(Ea)}')
json.dump(OUT, open(os.path.join(wd, 'diffimg.json'), 'w'), indent=1)
