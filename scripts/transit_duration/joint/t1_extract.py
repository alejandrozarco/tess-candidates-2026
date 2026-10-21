"""T1: extract SPOC LC FITS (2-min and 20-s) for L 98-59 with telemetry. Own code.
Usage: t1_extract.py <fits_dir> <out_dir>
Stores all finite-time cadences with QUALITY so dumps/scattered-light flags can be audited later."""
import sys, os, glob, numpy as np
from astropy.io import fits
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for f in sorted(glob.glob(os.path.join(src, '*lc.fits'))):
    fast = 'fast-lc' in f
    with fits.open(f, memmap=False) as h:
        h0, h1 = h[0].header, h[1].header; d = h[1].data
        col = lambda c: np.array(d[c], float) if c in d.columns.names else np.full(len(d), np.nan)
        t = col('TIME'); ok = np.isfinite(t)
        arr = {c.lower(): col(c)[ok] for c in ['TIME', 'SAP_FLUX', 'SAP_FLUX_ERR', 'SAP_BKG', 'PDCSAP_FLUX', 'PDCSAP_FLUX_ERR',
                                                'MOM_CENTR1', 'MOM_CENTR2', 'POS_CORR1', 'POS_CORR2']}
        arr['quality'] = np.array(d['QUALITY'], int)[ok]
        meta = dict(sector=h0['SECTOR'], tic=h0['TICID'], fast=fast, cam=h0.get('CAMERA'), ccd=h0.get('CCD'),
                    crowd=h1.get('CROWDSAP', np.nan), flfr=h1.get('FLFRCSAP', np.nan), procver=str(h0.get('PROCVER', '')),
                    data_rel=h0.get('DATA_REL', -1), pdcmethd=str(h1.get('PDCMETHD', '')), pdcvar=h1.get('PDCVAR', np.nan),
                    cdpp1=h1.get('CDPP1_0', np.nan), exptime=h1.get('TIMEDEL', np.nan) * 86400, numcad=int(ok.sum()),
                    ap_npix=None)
        try:
            ap = h[2].data; meta['ap_npix'] = int(((ap & 2) > 0).sum())
        except Exception: pass
    tag = '%d_s%03d_%s' % (meta['tic'], meta['sector'], '20s' if fast else '2m')
    np.savez(os.path.join(out, tag + '.npz'), meta=np.array([meta], dtype=object), **arr)
    g = arr['quality'] == 0
    print(tag, ok.sum(), g.sum(), 'cam%s ccd%s' % (meta['cam'], meta['ccd']), '%.4f' % meta['crowd'], meta['procver'], meta['data_rel'], meta['pdcmethd'], meta['ap_npix'], flush=True)
