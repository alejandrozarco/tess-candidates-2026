#!/usr/bin/env python3
"""TOI-421 (TIC 94986319): CHEOPS inputs for section 6.12 -> work/toi421/cheops/data/
  - PIPE PSF photometry of the six CH_PR100024 visits published by Krenn et al. (2024), VizieR J/A+A/686/A301/cheops1
    (tab-separated, all columns); the CDS server is tried first, then the IDIA mirror;
  - the DACE CHEOPS archive index (public mode; the six visits and their start/end times) and the DRP 14.1.2
    DEFAULT-aperture light curves (archive revision V0300) of the six visits (DACE target name TESS_94986319).
usage: python scripts/toi421/cheops/download_cheops.py   (needs dace-query and pandas)"""
import hashlib, os, urllib.request
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA = os.path.join(ROOT, 'work', 'toi421', 'cheops', 'data')
KEYS = ['CH_PR100024_TG013601_V0300', 'CH_PR100024_TG014701_V0300', 'CH_PR100024_TG014702_V0300',
        'CH_PR100024_TG014703_V0300', 'CH_PR100024_TG015301_V0300', 'CH_PR100024_TG015302_V0300']


def vizier(fn):
    q = 'viz-bin/asu-tsv?-source=J/A%2BA/686/A301/cheops1&-out.max=unlimited&-out.all'
    for host in ('https://vizier.cds.unistra.fr/', 'https://vizier.idia.ac.za/'):
        try:
            txt = urllib.request.urlopen(host + q, timeout=300).read()
            if b'NFlux' in txt:
                open(fn, 'wb').write(txt); print('VizieR', host, len(txt), 'bytes'); return
        except Exception as e:
            print('VizieR', host, 'failed:', repr(e)[:120])
    raise SystemExit('VizieR table not retrieved')


def main():
    os.makedirs(DATA, exist_ok=True)
    fn = os.path.join(DATA, 'VizieR_J_A+A_686_A301_cheops1.tsv')
    if not os.path.exists(fn):
        vizier(fn)
    from dace_query.cheops import Cheops
    idx = os.path.join(DATA, 'dace_cheops_full_index.csv')
    if not os.path.exists(idx):
        Cheops.query_database(limit=200000, output_format='pandas').to_csv(idx)
    for k in KEYS:
        out = os.path.join(DATA, f'{k}_DRP_default.csv')
        if not os.path.exists(out):
            Cheops.get_lightcurve(target='TESS_94986319', aperture='default', filters={'file_key': {'equal': [k]}},
                                  output_format='pandas').to_csv(out, index=False)
    for f in sorted(os.listdir(DATA)):
        print(hashlib.sha256(open(os.path.join(DATA, f), 'rb').read()).hexdigest(), f)


if __name__ == '__main__':
    main()
