"""Build the TRICERATOPS target object: TIC field stars, TESScut images, per-sector apertures, calc_depths; pickle it.

usage: python scripts/triceratops_target.py <tic> --sectors 39,93,100,102 --depth 0.00168 [--work ./work]
output: <work>/<tic>/triceratops_target.pkl and <work>/<tic>/gaia_background_<tic>.csv (the Gaia DR3 background
        population drawn by TRICERATOPS; TRILEGAL was not reachable when this was run, so the Gaia option is used).
Apertures: for every sector the pixels within 2.0 px of the target position in the TRICERATOPS cutout (the same
aperture as aperture_photometry.py); for sector 13 of TIC 324429952 the SPOC 2-min pipeline aperture read from the
light-curve file. depth = fitted k^2 of the raw aperture flux (dilution_corrected=False).
Requires triceratops 1.1.0 (pip install triceratops) and network access to MAST/Gaia.
"""
import argparse, glob, os, pickle, sys, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
import triceratops.triceratops as tr
from triceratops.funcs import query_gaia_background
from astropy.io import fits
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--sectors', required=True); p.add_argument('--depth', type=float, required=True)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); tic = args.tic; sectors = [int(x) for x in args.sectors.split(',')]
t = tr.target(ID=tic, sectors=np.array(sectors))
for a in range(5):
    if isinstance(t.trilegal_fname, str):
        break
    print('background query failed; retry', a + 1, flush=True); time.sleep(20)
    t.trilegal_fname = query_gaia_background(t.stars.ra.values[0], t.stars.dec.values[0], tic, verbose=1)
if not isinstance(t.trilegal_fname, str):
    sys.exit('no background population')
dst = os.path.join(wd, f'gaia_background_{tic}.csv'); os.replace(t.trilegal_fname, dst); t.trilegal_fname = dst
print('background population:', dst, sum(1 for _ in open(dst)) - 1, 'rows; sectors with images:', len(t.pix_coords), 'of', len(sectors))
aps = []
for k, s in enumerate(sectors):
    X, Y = t.pix_coords[k][0]
    if s == c.get('spoc_2min'):
        h = fits.open(glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{s:04d}-{tic:016d}*', '*lc.fits'))[0])['APERTURE']
        a = h.data.astype(int); rr, cc = np.nonzero(a & 2); ap = np.column_stack([cc + h.header['CRVAL1P'], rr + h.header['CRVAL2P']]); src = 'SPOC'
    else:
        c0, r0 = int(np.floor(X)) - 4, int(np.floor(Y)) - 4; cc, rr = np.meshgrid(np.arange(c0, c0 + 10), np.arange(r0, r0 + 10)); cc, rr = cc.ravel(), rr.ravel()
        m = np.hypot(cc - X, rr - Y) < 2.0; ap = np.column_stack([cc[m], rr[m]]); src = 'r2.0'
    aps.append(ap); print(f'S{s}: target px ({X:.3f},{Y:.3f}) aperture {src} npix={len(ap)}')
t.calc_depths(tdepth=args.depth, all_ap_pixels=aps, dilution_corrected=False)
import pandas as pd; pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(t.stars[['ID', 'Tmag', 'ra', 'dec', 'mass', 'rad', 'Teff', 'plx', 'sep (arcsec)', 'fluxratio', 'tdepth', 'disposition']].head(25).to_string())
print('stars with tdepth > 0 (possible hosts):', int((t.stars.tdepth > 0).sum()))
t.aps = aps; pickle.dump(t, open(os.path.join(wd, 'triceratops_target.pkl'), 'wb'))
