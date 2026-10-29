#!/usr/bin/env python3
"""TOI-421 (TIC 94986319): the fits of fit_toi421.py and planet_b_toi421.py repeated on own photometry of the TESS
full-frame images, which does not use SPOC PDC (the SPOC 2-min products of sectors 5, 6, 32 are spoc-5.0.20, those of
98 and 107 spoc-5.0.126).
Photometry: TESScut 25 x 25-pixel cutouts; raw = FLUX + FLUX_BKG (the SPOC FFI background is added back); own background
= median of the annulus 6-11 px (pixels above median + 3 MAD of the median image excluded); circular aperture of radius
R = 2.0, 2.5, 3.0 px on the flux-weighted centroid of the median image; QUALITY == 0; no co-trending (the fits use
per-transit quadratic baselines); errors scaled to the point-to-point scatter. One code path for every sector.
Cadences: 1800 s (sectors 5, 6), 600 s (32), 200 s (98, 107); the models are integrated over each exposure.
usage: python scripts/toi421/ffi_toi421.py   -> work/toi421/ffi/ (cutouts), work/toi421/data/s<sector>_ffiR<R>.npz,
       work/toi421/out/ffi.json
"""
import glob, json, os, re, sys
import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_toi421 as F
import planet_b_toi421 as PBF

CUT = os.path.join(F.WORK, 'ffi'); DATA = os.path.join(F.WORK, 'data'); OUT = os.path.join(F.WORK, 'out')
SECTORS = (5, 6, 32, 98, 107)
EXP = {5: 1800., 6: 1800., 32: 600., 98: 200., 107: 200.}
RADII = (2.0, 2.5, 3.0)


def download():
    """TESScut REST API at the target position of the SPOC light-curve header (RA_OBJ, DEC_OBJ)"""
    import io, urllib.request, zipfile
    h = fits.getheader(sorted(glob.glob(os.path.join(DATA, '*s0005*_lc.fits')))[0], 0)
    os.makedirs(CUT, exist_ok=True)
    for s in SECTORS:
        d = os.path.join(CUT, 's%03d' % s)
        if glob.glob(os.path.join(d, '*astrocut.fits')):
            continue
        os.makedirs(d, exist_ok=True)
        url = ('https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra=%.6f&dec=%.6f&y=25&x=25&sector=%d'
               % (h['RA_OBJ'], h['DEC_OBJ'], s))
        print('TESScut', url, flush=True)
        zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(url, timeout=600).read())).extractall(d)


def photometry(R):
    tags = []
    for s in SECTORS:
        f = sorted(glob.glob(os.path.join(CUT, 's%03d' % s, '*astrocut.fits')))[0]
        d = fits.open(f)[1].data
        q = (d['QUALITY'] == 0) & np.isfinite(d['TIME'])
        cube = (d['FLUX'] + d['FLUX_BKG'])[q]; t = d['TIME'][q]
        ok = np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1); cube, t = cube[ok], t[ok]
        med = np.nanmedian(cube, axis=0); ny, nx = med.shape; yy, xx = np.mgrid[:ny, :nx]
        box = (np.abs(yy - ny // 2) <= 2) & (np.abs(xx - nx // 2) <= 2)
        w = np.clip(med - np.median(med), 0, None) * box
        cy, cx = (w * yy).sum() / w.sum(), (w * xx).sum() / w.sum()
        r = np.hypot(yy - cy, xx - cx)
        ann = (r >= 6) & (r <= 11); mm = med[ann]; good = ann.copy()
        good[ann] = mm < np.median(mm) + 3 * 1.4826 * np.median(np.abs(mm - np.median(mm)))
        bkg = np.median(cube[:, good], axis=1)
        ap = r <= R; fl = (cube[:, ap] - bkg[:, None]).sum(axis=1)
        npix, nb = ap.sum(), good.sum()
        e = np.sqrt(np.clip(fl, 0, None) + npix * (1 + npix / nb) * np.var(cube[:, good] - bkg[:, None], axis=1))
        p2p = 1.4826 * np.median(np.abs(np.diff(fl))) / np.sqrt(2); e = e * p2p / np.median(e)
        tag = 's%03d_ffiR%.1f' % (s, R)
        np.savez(os.path.join(DATA, tag + '.npz'), t=t, pdc=fl, pdce=e, sap=fl, sape=e, q=np.zeros(len(t), int), exp=EXP[s],
                 crowd=1.0, flfr=1.0, cam=0, ccd=0, npix=int(npix))
        tags.append(tag)
    return tags


def main():
    download()
    res = {}
    for R in RADII:
        tags = photometry(R)
        ws = F.all_windows(tags=tags)
        A = F.per_transit(ws); B = F.per_sector(ws)
        res['R%.1f' % R] = dict(c_per_transit=F.short(A), c_slope=F.slope(A), c_per_sector=F.short(B), c_sector_slope=F.slope(B))
        if R == 2.5:
            res['R2.5']['b_per_sector'] = {}
            for tag in tags:
                r = PBF.fit_b(PBF.wins([tag]), 'batman'); res['R2.5']['b_per_sector'][tag] = r
            E = PBF.wins(tags[:3]); L = PBF.wins(tags[3:])
            res['R2.5']['b_epochs'] = {k: dict(early=PBF.fit_b(E, k), late=PBF.fit_b(L, k)) for k in ('batman', 'box')}
        print('R %.1f c slope %+.4f +- %.4f (z %.1f)' % (R, res['R%.1f' % R]['c_slope']['slope'], res['R%.1f' % R]['c_slope']['eslope_scaled'],
                                                       res['R%.1f' % R]['c_slope']['z']), flush=True)
    res['_meta'] = __doc__
    json.dump(res, open(os.path.join(OUT, 'ffi.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
