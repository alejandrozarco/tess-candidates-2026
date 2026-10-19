"""Download the public MAST products used for the other TESS results (README section 6).

usage: python scripts/other_download.py 11613676 [--work ./work]
       python scripts/other_download.py 124751008 [--work ./work]
       python scripts/other_download.py 200322593 [--controls] [--work ./work]
11613676 (CD-30 448): 21x21 TESScut cutouts of sectors 3, 30, 97, 104 and 106, requested at the Gaia DR3 position
    (<work>/11613676/cut/s<sector>.npz, same format as tesscut_download.py; target pixel = proper-motion-corrected Gaia
    position at the median time of the sector), and the SPOC 2-min light curves of sectors 3 and 30.
124751008 (2MASS J06162084-3509125): 21x21 TESScut cutouts of sectors 6, 7, 33, 87, 98, 108 and 109 in the same
    format (about 170 MB compressed); where TESScut returns more than one cutout for a sector, the one with the most
    complete cadences is kept. The star has no SPOC 2-min data.
200322593 (TOI-540): SPOC 2-min light curves (*-s_lc.fits) and target pixel files (*-s_tp.fits) of sectors 4, 5, 6,
    31, 32 and 98 (about 330 MB) and the 20-s light curves (*-a_fast-lc.fits) of sectors 31, 32 and 98 (about 50 MB). --controls: SPOC 2-min light curves of the control planets of targets.py for their
    first sector, the last two sectors before 98, and sectors 98 and later (about 150 MB).
MAST products go to <work>/<tic>/dl/mastDownload/TESS/ (control stars: <work>/controls/dl/). Existing files are not
downloaded again. Requires astroquery and astropy; network access to MAST.
"""
import argparse, glob, os, re, time
import numpy as np
from astroquery.mast import Observations, Tesscut
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from targets import OTHER

p = argparse.ArgumentParser(description=__doc__); p.add_argument('tic', type=int, choices=sorted(OTHER))
p.add_argument('--work', default='./work'); p.add_argument('--controls', action='store_true')
args = p.parse_args(); c = OTHER[args.tic]; wd = os.path.join(args.work, str(args.tic)); os.makedirs(os.path.join(wd, 'cut'), exist_ok=True)


def retry(f, *a, n=4, **k):
    for i in range(n):
        try:
            return f(*a, **k)
        except Exception as e:
            print('retry', i + 1, str(e)[:120], flush=True); time.sleep(20 * (i + 1))
    raise RuntimeError('MAST request failed')


def spoc_2min(tic, sectors, kinds, outdir, fast=False):
    """SPOC 2-min products of one TIC (fast=True: 20-s light curves): kinds from {'lc', 'tp'}; sectors = list of ints or a function of the available list"""
    sfx = 'a_fast' if fast else 's'
    obs = retry(Observations.query_criteria, obs_collection='TESS', target_name=str(tic), dataproduct_type='timeseries', provenance_name='SPOC')
    have = sorted(int(m.group(1)) for o in obs['obs_id'] if (m := re.match(rf'tess\d+-s(\d{{4}})-\d{{16}}-\d{{4}}-{sfx}$', str(o))))
    use = sectors(have) if callable(sectors) else [s for s in sectors if s in have]
    print(f'TIC {tic}: SPOC {"20-s" if fast else "2-min"} sectors {have}; downloading {use} ({",".join(kinds)})', flush=True)
    for s in use:
        o = obs[[bool(re.match(rf'tess\d+-s{s:04d}-{tic:016d}-\d{{4}}-{sfx}$', str(x))) for x in obs['obs_id']]]
        prods = retry(Observations.get_product_list, o)
        want = prods[[any(str(f).endswith(f'-{sfx}{"-" if fast else "_"}{k}.fits') for k in kinds) for f in prods['productFilename']]]
        todo = want[[not glob.glob(os.path.join(outdir, 'mastDownload', 'TESS', '*', str(f))) for f in want['productFilename']]]
        if len(todo):
            retry(Observations.download_products, todo, download_dir=outdir)
        print(f'   S{s}: {len(want)} files ({len(want) - len(todo)} present)', flush=True)


if args.tic in (11613676, 124751008):
    coord = SkyCoord(c['ra'], c['dec'], unit='deg'); N = 21
    for sec in c['sectors']:
        fn = os.path.join(wd, 'cut', f's{sec}.npz')
        if os.path.exists(fn):
            continue
        cuts = retry(Tesscut.get_cutouts, coordinates=coord, size=N, sector=sec)
        hd = cuts[0] if args.tic == 11613676 else max(cuts, key=lambda h: int(np.isfinite(h[1].data['FLUX']).all(axis=(1, 2)).sum()))
        d = hd[1].data; w = WCS(hd[2].header)
        yr = 2016.0 + (np.nanmedian(d['TIME']) + 2457000 - 2457388.5) / 365.25
        pos = SkyCoord(c['ra'] + c['pmra'] / 3.6e6 * (yr - 2016) / np.cos(np.radians(c['dec'])), c['dec'] + c['pmdec'] / 3.6e6 * (yr - 2016), unit='deg')
        tx, ty = [float(v) for v in w.world_to_pixel(pos)]
        np.savez_compressed(fn, tm=d['TIME'], f=d['FLUX'].astype(np.float32), fe=d['FLUX_ERR'].astype(np.float32),
                            bkg=d['FLUX_BKG'].astype(np.float32), q=d['QUALITY'], tx=tx, ty=ty, cam=hd[0].header.get('CAMERA'),
                            ccd=hd[0].header.get('CCD'), wcs=hd[2].header.tostring(), ra=pos.ra.deg, dec=pos.dec.deg, epoch=yr)
        print(sec, 'ok', len(d), 'cam', hd[0].header.get('CAMERA'), hd[0].header.get('CCD'), 'target px', round(tx, 2), round(ty, 2), flush=True)
    if c['spoc_2min']:
        spoc_2min(args.tic, c['spoc_2min'], ['lc'], os.path.join(wd, 'dl'))
else:
    spoc_2min(args.tic, c['sectors'], ['lc', 'tp'], os.path.join(wd, 'dl'))
    spoc_2min(args.tic, [31, 32, 98], ['lc'], os.path.join(wd, 'dl'), fast=True)
    if args.controls:
        cd = os.path.join(args.work, 'controls', 'dl')
        for tic in sorted(set(x[1] for x in c['controls'])):
            spoc_2min(tic, lambda h: sorted(set([h[0]] + [s for s in h if s < 98][-2:] + [s for s in h if s >= 98])), ['lc'], cd)
print('done')
