#!/usr/bin/env python3
"""TOI-421 (TIC 94986319): download the SPOC 2-min and 20-s light curves of sectors 5, 6, 32, 98 and 107 from MAST,
check their SHA-256 against the files used for section 6.12, and convert them to npz.

usage: python scripts/toi421/download_toi421.py
writes work/toi421/data/<fits> and work/toi421/data/s<sector>_<2m|20s>.npz
       (keys: t BTJD, pdc, pdce, sap, sape, bkg, q, exp in s, crowd, flfr, cam, ccd, npix),
       work/toi421/out/lc_header_info.json (CROWDSAP, FLFRCSAP, camera, CCD, aperture size, SPOC version per file)
"""
import hashlib, json, os, urllib.request
import numpy as np
from astropy.io import fits

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.path.join(ROOT, 'work', 'toi421'); DATA = os.path.join(WORK, 'data'); OUT = os.path.join(WORK, 'out')
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/'
FILES = {   # file name: sha256 of the version used (MAST, retrieved 2026-10-29)
    'tess2018319095959-s0005-0000000094986319-0125-s_lc.fits': '884b71658fe341e29a9c9eb5193cd1b6a98ddf2d296385eed726c41bd0ccfd93',
    'tess2018349182500-s0006-0000000094986319-0126-s_lc.fits': '11282571ce22e836fca68468c9cacd36e4a156662d216bde7ac1501f9a47da1d',
    'tess2020324010417-s0032-0000000094986319-0200-s_lc.fits': '56657e250f8a6b6543ae8bd3e4a97e2b986cf1f06c68ab58cc2e2d64b8d45fce',
    'tess2025312202959-s0098-0000000094986319-0298-a_fast-lc.fits': 'f99acf4c246488e8bec1420f94be77bd87ff33a559650843c147537e31779763',
    'tess2025312202959-s0098-0000000094986319-0298-s_lc.fits': '871909062004b65de8fe66932bc0af7c648aa72c5aa1bf9bc8c2ea8cfccc04be',
    'tess2026221103959-s0107-0000000094986319-0309-a_fast-lc.fits': 'c3f33366ec81529f26d1984a03b2373b19ae0f3d71cee0139637e543c6fde244',
    'tess2026221103959-s0107-0000000094986319-0309-s_lc.fits': '1c2fd8f0f1352f02a5d38989ca17e81bb595d40210647093c26d7b4b9e066b35',
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(DATA, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    info = {}
    for fn, ref in FILES.items():
        path = os.path.join(DATA, fn)
        if not os.path.exists(path):
            print('downloading', fn, flush=True)
            urllib.request.urlretrieve(URL + fn, path)
        h = sha256(path)
        print(fn, 'sha256', 'matches' if h == ref else 'DIFFERS (MAST product changed?) ' + h)
        h0 = fits.getheader(path, 0); h1 = fits.getheader(path, 1); d = fits.getdata(path, 1); ap = fits.getdata(path, 2)
        npix = int(((ap & 2) > 0).sum())
        fast = 'fast' in fn; tag = 's%03d_%s' % (h0['SECTOR'], '20s' if fast else '2m')
        ok = np.isfinite(d['TIME'])
        np.savez(os.path.join(DATA, tag + '.npz'), t=d['TIME'][ok], pdc=d['PDCSAP_FLUX'][ok], pdce=d['PDCSAP_FLUX_ERR'][ok],
                 sap=d['SAP_FLUX'][ok], sape=d['SAP_FLUX_ERR'][ok], bkg=d['SAP_BKG'][ok], q=d['QUALITY'][ok],
                 exp=20. if fast else 120., crowd=h1.get('CROWDSAP'), flfr=h1.get('FLFRCSAP'), cam=h0['CAMERA'],
                 ccd=h0['CCD'], npix=npix)
        info[tag] = dict(file=fn, sha256=h, crowdsap=h1.get('CROWDSAP'), flfrcsap=h1.get('FLFRCSAP'), cam=h0['CAMERA'],
                         ccd=h0['CCD'], npix=npix, procver=h0.get('PROCVER'), data_rel=h0.get('DATA_REL'))
    json.dump(info, open(os.path.join(OUT, 'lc_header_info.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
