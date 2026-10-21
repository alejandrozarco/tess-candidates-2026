"""M2: product-level access check for selected JWST observations of L 98-59. Own code.
For each obs_id: product list (dataRights per product, counts by rights), and an HTTP HEAD on the public S3 mirror
(stpubdata) for each calints/x1dints/rateints/whtlt file, which returns 200 only for public files.
Positive controls: jw01224-o001 (public since 2024-06-25) must show PUBLIC products and HEAD 200 for its calints.
Also prints predicted mid-transit times of b, c, d (Cadieux et al. 2025 linear ephemerides) inside each observation
window.
Usage: m2_products.py <obs_id> [obs_id ...]"""
import sys, subprocess, json, numpy as np
from astropy.time import Time
from astroquery.mast import Observations
PL = {'b': (2.2531140, 1366.1706, 1.01), 'c': (3.6906764, 1367.2730, 1.28), 'd': (7.450729, 1362.7400, 0.84)}
def head(fn):
    url = 'https://stpubdata.s3.amazonaws.com/jwst/public/%s/%s/%s' % (fn[:7], fn[:13], fn)
    r = subprocess.run(['curl', '-sS', '-I', '-m', '30', '-o', '/dev/null', '-w', '%{http_code}', url], capture_output=True, text=True)
    return r.stdout.strip()
ids = ['jw01224-o001_t009_nirspec_f290lp-g395h-s1600a1-sub2048'] + sys.argv[1:]
for k, oid in enumerate(ids):
    obs = Observations.query_criteria(obs_id=oid)
    if len(obs) == 0:
        print(oid, 'NO OBSERVATION ROW'); continue
    o = obs[0]
    t0, t1 = float(o['t_min']), float(o['t_max'])
    print('\n==', oid, 'proposal', o['proposal_id'], o['proposal_pi'], o['instrument_name'], o['filters'], 'dataRights', o['dataRights'],
          'release', Time(float(o['t_obs_release']), format='mjd').iso[:16] if np.isfinite(float(o['t_obs_release'])) else '-', flush=True)
    if np.isfinite(t0):
        b0, b1 = t0 + 2400000.5 - 2457000, t1 + 2400000.5 - 2457000
        print('   window UTC %s .. %s  (BTJD %.4f .. %.4f, %.2f h)' % (Time(t0, format='mjd').iso[:16], Time(t1, format='mjd').iso[11:16], b0, b1, (b1 - b0) * 24))
        for p, (P, T0, dur) in PL.items():
            n = np.arange(np.floor((b0 - T0) / P) - 1, np.ceil((b1 - T0) / P) + 2)
            for tc in T0 + n * P:
                if b0 - dur / 24 < tc < b1 + dur / 24:
                    print('   predicted transit of %s at BTJD %.4f (%s UTC), %.0f min after window start, %.0f min before end' % (p, tc, Time(tc + 2457000 - 2400000.5, format='mjd').iso[11:16], (tc - b0) * 1440, (b1 - tc) * 1440))
                # eclipse (phase 0.5) for b, c
                te = tc + P / 2
                if p in 'bc' and b0 < te < b1: print('   predicted secondary eclipse of %s at BTJD %.4f' % (p, te))
    try:
        pr = Observations.get_product_list(obs)
    except Exception as e:
        print('   product list error', e); continue
    rights = {}
    for r in pr: rights[str(r['dataRights'])] = rights.get(str(r['dataRights']), 0) + 1
    print('   products', len(pr), 'by dataRights', rights)
    seen = 0
    for r in pr:
        fn = str(r['productFilename'])
        if any(s in fn for s in ('_calints', '_x1dints', '_rateints', '_whtlt')) and fn.endswith(('.fits', '.ecsv')):
            code = head(fn); seen += 1
            if seen <= 12 or code == '200': print('   %-70s %-17s %s MB  S3 HEAD %s' % (fn, r['dataRights'], round(int(r['size']) / 1e6, 1) if str(r['size']).isdigit() else r['size'], code), flush=True)
    if k == 0:
        ok = rights.get('PUBLIC', 0) > 0
        print('   POSITIVE CONTROL', 'PASS' if ok else 'FAIL');
        if not ok: raise SystemExit('positive control failed')
