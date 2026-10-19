"""Download 21x21 pixel TESScut cutouts (SPOC-calibrated full-frame images) for one target, one file per sector.

usage: python scripts/tesscut_download.py <tic> [--work ./work] [--sectors 12,100,103]
output: <work>/<tic>/cut/s<sector>.npz with time (BTJD), flux and error cubes (e/s), SPOC background cube,
        quality flags, position-correction columns, the target pixel position and the cutout WCS header.
The cutout is requested at the TIC v8.2 coordinates and the target pixel is that position, except for targets with
cut_center='gaia' in targets.py (TIC 89247735, 0.11 arcsec/yr): there the cutout is requested at the Gaia DR3 epoch
2016.0 position and the target pixel is the proper-motion-corrected Gaia position at the median time of the sector.
Requires astroquery and astropy; network access to MAST. Retries each sector four times.
"""
import argparse, os, time
import numpy as np
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from targets import TARGETS, workdir, add_common_args

ap = add_common_args(argparse.ArgumentParser(description=__doc__)); ap.add_argument('--sectors', default=None, help='comma-separated subset of the sectors in targets.py')
args = ap.parse_args()
c = TARGETS[args.tic]; wd = workdir(args, args.tic); N = 21
GAIA = c.get('cut_center') == 'gaia'
coord = SkyCoord(c['ra'], c['dec'], unit='deg') if GAIA else SkyCoord(c['tic_ra'], c['tic_dec'], unit='deg')
for sec in ([int(x) for x in args.sectors.split(',')] if args.sectors else c['sectors']):
    fn = os.path.join(wd, 'cut', f's{sec}.npz')
    if os.path.exists(fn):
        continue
    hd = None
    for a in range(4):
        try:
            hd = Tesscut.get_cutouts(coordinates=coord, size=N, sector=sec)[0]; break
        except Exception as e:
            print('retry', sec, str(e)[:100], flush=True); time.sleep(15 * (a + 1))
    if hd is None:
        print(sec, 'FAILED'); continue
    d = hd[1].data; w = WCS(hd[2].header); pos = coord; extra = {}
    if GAIA:
        yr = 2016.0 + (np.nanmedian(d['TIME']) + 2457000 - 2457388.5) / 365.25
        pos = SkyCoord(c['ra'] + c['pmra'] / 3.6e6 * (yr - 2016) / np.cos(np.radians(c['dec'])), c['dec'] + c['pmdec'] / 3.6e6 * (yr - 2016), unit='deg')
        extra = dict(ra=pos.ra.deg, dec=pos.dec.deg, epoch=yr)
    tx, ty = [float(v) for v in w.world_to_pixel(pos)]
    np.savez_compressed(fn, tm=d['TIME'], f=d['FLUX'].astype(np.float32), fe=d['FLUX_ERR'].astype(np.float32),
                        bkg=d['FLUX_BKG'].astype(np.float32), q=d['QUALITY'], pc1=d['POS_CORR1'], pc2=d['POS_CORR2'],
                        tx=tx, ty=ty, cam=hd[0].header.get('CAMERA'), ccd=hd[0].header.get('CCD'), wcs=hd[2].header.tostring(), **extra)
    print(sec, 'ok', len(d), 'cam', hd[0].header.get('CAMERA'), hd[0].header.get('CCD'), 'target px', round(tx, 2), round(ty, 2), flush=True)
