"""Download the public data used for TOI-1117 c (TIC 295541511) in README section 6.10.

usage: python scripts/toi1117/download_toi1117.py [--work ./work]
output: <work>/295541511/dl/lc/   SPOC 2-min (sectors 39, 100-103) and 20-s (39) light curves; TESS-SPOC (13, 39), QLP (13, 39,
                                  100-103), eleanor-lite (13) and TARS (13, 39) light curves
        <work>/295541511/dl/tpf/  SPOC 2-min target pixel files (39, 100-103) and TESS-SPOC FFI target pixel files (13, 39)
        <work>/295541511/dl/dv/   SPOC and TESS-SPOC data-validation XML reports
        <work>/295541511/dl/products.txt   the MAST product list used (provenance, sector, URI)
        <work>/295541511/dl/gaia_150as.csv Gaia DR3 sources within 150" (VizieR I/355/gaiadr3)
        <work>/295541511/dl/tic_150as.csv  TIC v8.2 sources within 150" (MAST catalogs; TIC and Gaia DR2 identifiers)
        <work>/295541511/dl/arxiv/2506.05521.tar.gz  the arXiv source of Lockley et al. (2025), which contains the HARPS table
The product list comes from astroquery.mast (Observations.query_criteria on target_name, then get_product_list); files are
fetched from the MAST download service. Existing files are not downloaded again. The downloaded files are read as data only.
"""
import argparse, hashlib, os, sys, urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import TIC, RA, DEC, dldir

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); args = p.parse_args()


def want(f):
    if f.endswith('.jpeg') or 'diamante' in f:
        return None
    if f.endswith('-s_lc.fits') or f.endswith('a_fast-lc.fits') or f.endswith('llc.fits') or ('eleanor' in f and f.endswith('_lc.fits')) \
            or ('tars' in f and f.endswith('_lc.fits')) or ('tess-spoc' in f and f.endswith('_lc.fits')):
        return 'lc'
    if f.endswith('-s_tp.fits') or ('tess-spoc' in f and f.endswith('_tp.fits')):
        return 'tpf'
    if f.endswith('dvr.xml'):
        return 'dv'
    return None


def fetch(url, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return False
    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'python-urllib'}), timeout=900) as r, open(out + '.part', 'wb') as f:
        f.write(r.read())
    os.rename(out + '.part', out); print('downloaded', os.path.basename(out), os.path.getsize(out)); return True


from astroquery.mast import Observations
obs = Observations.query_criteria(target_name=str(TIC), obs_collection=['TESS', 'HLSP'])
prods = Observations.get_product_list(obs).to_pandas().drop_duplicates('productFilename')
lines = []
for r in prods.itertuples():
    sub = want(r.productFilename)
    if sub is None:
        continue
    lines.append(f'{r.obs_collection} {r.project} {r.productFilename} {r.dataURI}')
    fetch('https://mast.stsci.edu/api/v0.1/Download/file?uri=' + urllib.parse.quote(r.dataURI, safe=''), os.path.join(dldir(args.work, sub), r.productFilename))
with open(os.path.join(dldir(args.work), 'products.txt'), 'w') as f:
    f.write('\n'.join(sorted(set(lines))) + '\n')
print(len(set(lines)), 'products listed')

g = os.path.join(dldir(args.work), 'gaia_150as.csv')
if not os.path.exists(g):
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    V = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'BPmag', 'RPmag', 'Plx', 'pmRA', 'pmDE', 'RUWE', '_r'], row_limit=-1)
    t = V.query_region(SkyCoord(RA, DEC, unit='deg'), radius=150 * u.arcsec, catalog='I/355/gaiadr3')[0].to_pandas().sort_values('_r')
    if len(t) < 10:
        sys.exit('Gaia DR3 query returned too few rows; not written')
    t.to_csv(g, index=False); print('Gaia DR3 sources within 150":', len(t))

g = os.path.join(dldir(args.work), 'tic_150as.csv')
if not os.path.exists(g):
    from astroquery.mast import Catalogs
    t = Catalogs.query_region(f'{RA} {DEC}', radius=150 / 3600, catalog='TIC').to_pandas()[['ID', 'GAIA', 'ra', 'dec', 'Tmag', 'GAIAmag', 'dstArcSec']]
    if len(t) < 10:
        sys.exit('TIC query returned too few rows; not written')
    t.to_csv(g, index=False); print('TIC sources within 150":', len(t))

a = os.path.join(dldir(args.work, 'arxiv'), '2506.05521.tar.gz')
fetch('https://arxiv.org/e-print/2506.05521', a)
print('arXiv source sha256', hashlib.sha256(open(a, 'rb').read()).hexdigest())
