"""Download the public TESS products of TOI-2638.01 (TIC 458342769) used in README section 6.9.

usage: python scripts/toi2638/download_toi2638.py [--work ./work] [--tpf]
output: <work>/458342769/dl/lc/   SPOC 2-min light curves (sectors 36, 37, 63, 64, 90, 99), TESS-SPOC (36, 37, 63, 64),
                                  QLP (10, 36, 37, 63, 64, 90, 99) and TGLC (10) light curves of the target
        <work>/458342769/dl/nb17/ QLP light curves of the 17" neighbour TIC 458342741 (TYC 8605-330-1)
        <work>/458342769/dl/tpf/  SPOC 2-min target pixel files (with --tpf; about 280 MB)
        <work>/458342769/dl/products.txt  the MAST product list (provenance, sector, URI)
The product list comes from the MAST CAOM service (Mast.Caom.Filtered on target_name); files are fetched from the
MAST download service. Existing files are not downloaded again.
"""
import argparse, json, os, sys, urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import TIC, TIC_NB17, dldir

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--tpf', action='store_true')
args = p.parse_args()
WANT = {'SPOC': [36, 37, 63, 64, 90, 99], 'TESS-SPOC': [36, 37, 63, 64], 'QLP': [10, 36, 37, 63, 64, 90, 99], 'TGLC': [10]}


def mast(req):
    data = urllib.parse.urlencode({'request': json.dumps(req)}).encode()
    with urllib.request.urlopen(urllib.request.Request('https://mast.stsci.edu/api/v0/invoke', data=data), timeout=180) as r:
        return json.loads(r.read())


def products(tic):
    req = {'service': 'Mast.Caom.Filtered', 'format': 'json', 'params': {
        'columns': 'obs_collection,provenance_name,sequence_number,dataURL',
        'filters': [{'paramName': 'target_name', 'values': [str(tic), f'TIC {tic}']}]}}
    rows = mast(req)['data']
    if tic != TIC:
        return rows
    req['params']['filters'] = [{'paramName': 'target_name', 'values': ['5353676910584543872', 'gaiaid-5353676910584543872']}]
    try:
        rows += mast(req)['data']          # TGLC products are listed under the Gaia DR3 source_id
    except Exception as e:
        print('TGLC query failed:', e)
    return rows


def fetch(uri, out):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return
    url = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=' + urllib.parse.quote(uri, safe='')
    with urllib.request.urlopen(url, timeout=900) as r, open(out + '.part', 'wb') as f:
        f.write(r.read())
    os.rename(out + '.part', out); print('downloaded', os.path.basename(out), os.path.getsize(out))


dl = dldir(args.work)
for sub in ('lc', 'nb17') + (('tpf',) if args.tpf else ()):
    os.makedirs(os.path.join(dl, sub), exist_ok=True)
lines = []
for tic, sub in ((TIC, 'lc'), (TIC_NB17, 'nb17')):
    for r in products(tic):
        prov, sec, uri = r.get('provenance_name'), r.get('sequence_number'), r.get('dataURL') or ''
        if not uri.endswith('.fits') or prov is None:
            continue
        prov = prov.upper()
        if tic == TIC_NB17 and prov != 'QLP':
            continue
        if prov not in WANT or sec not in WANT[prov]:
            continue
        if prov == 'SPOC' and not uri.endswith('-s_lc.fits'):
            continue
        lines.append(f'{tic} {prov} {sec} {uri}')
        fetch(uri, os.path.join(dl, sub, os.path.basename(uri)))
        if prov == 'SPOC' and args.tpf and tic == TIC:
            tp = uri.replace('_lc.fits', '_tp.fits'); fetch(tp, os.path.join(dl, 'tpf', os.path.basename(tp)))
with open(os.path.join(dl, 'products.txt'), 'w') as f:
    f.write('\n'.join(sorted(set(lines))) + '\n')
print(len(set(lines)), 'products listed in', os.path.join(dl, 'products.txt'))
