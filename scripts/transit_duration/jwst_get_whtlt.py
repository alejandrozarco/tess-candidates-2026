import sys, os, subprocess
from astroquery.mast import Observations
outdir = sys.argv[1]
for oid in sys.argv[2:]:
    obs = Observations.query_criteria(obs_id=oid)
    pr = Observations.get_product_list(obs)
    for r in pr:
        fn = str(r['productFilename'])
        if fn.endswith('_whtlt.ecsv'):
            p = os.path.join(outdir, fn)
            if not os.path.exists(p):
                subprocess.run(['curl', '-sS', '-f', '-L', '--max-time', '300', '-o', p, 'https://mast.stsci.edu/api/v0.1/Download/file?uri=' + str(r['dataURI'])])
            print(fn, os.path.getsize(p))
