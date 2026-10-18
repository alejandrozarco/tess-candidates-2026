"""Aperture photometry on the TESScut cutouts: circular apertures around the target and a local annulus background.

usage: python scripts/aperture_photometry.py <tic> [--work ./work]
input:  <work>/<tic>/cut/s<sector>.npz (from tesscut_download.py)
output: <work>/<tic>/cut/lc_s<sector>.npz with tm (BTJD), q (quality), bkg (annulus background per pixel, e/s),
        f1.0 ... f3.0 (summed background-subtracted flux in circles of radius 1.0-3.0 px around the target),
        cx, cy (flux-weighted centroid within r < 2.5 px) and fs2.0 (r < 2.0 px sum of the SPOC-background-subtracted flux).
The background is the median of the lowest 70 % of the pixels in the annulus 6 < r < 9.5 px (rejects stars).
The tables in data/ use the r = 2.0 px aperture.
"""
import argparse, glob, os
import numpy as np
from targets import workdir, add_common_args

args = add_common_args(argparse.ArgumentParser(description=__doc__)).parse_args()
wd = workdir(args, args.tic); N = 21; yy, xx = np.indices((N, N))
for fn in sorted(glob.glob(os.path.join(wd, 'cut', 's*.npz'))):
    s = int(os.path.basename(fn)[1:-4]); out = os.path.join(wd, 'cut', f'lc_s{s}.npz')
    if os.path.exists(out):
        continue
    z = np.load(fn); tm = z['tm']; f = z['f'].astype(float); q = z['q']; tx, ty = float(z['tx']), float(z['ty'])
    ok = np.isfinite(tm) & np.all(np.isfinite(f), axis=(1, 2)); tm, f, q = tm[ok], f[ok], q[ok]
    r = np.hypot(xx - tx, yy - ty); ann = (r > 6) & (r < 9.5)
    bkg = np.array([np.median(np.sort(fr[ann])[:int(ann.sum() * 0.7)]) for fr in f])
    g = f - bkg[:, None, None]
    res = {'tm': tm, 'q': q, 'bkg': bkg, 'tx': tx, 'ty': ty}
    for rad in [1.0, 1.5, 2.0, 2.5, 3.0]:
        res[f'f{rad}'] = g[:, r < rad].sum(1)
    m = r < 2.5; w = np.clip(g[:, m], 0, None)
    res['cx'] = (w * xx[m]).sum(1) / w.sum(1); res['cy'] = (w * yy[m]).sum(1) / w.sum(1)
    res['fs2.0'] = f[:, r < 2.0].sum(1)
    np.savez(out, **res); print(s, len(tm), 'cadences; median annulus background', round(float(np.median(bkg)), 2), 'e/s/px')
