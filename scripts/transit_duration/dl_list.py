import sys, json
from astroquery.mast import Observations
tic = sys.argv[1]
obs = Observations.query_criteria(target_name=tic, obs_collection="TESS", dataproduct_type="timeseries")
prods = Observations.get_product_list(obs)
sel = [r for r in prods if str(r['productFilename']).endswith(('_lc.fits','_fast-lc.fits','_tp.fits','fast-tp.fits'))]
out = [{'fn':str(r['productFilename']),'uri':str(r['dataURI']),'size':int(r['size'])} for r in sel]
json.dump(out, open(sys.argv[2],'w'), indent=0)
print(len(out), sum(o['size'] for o in out)/1e6, 'MB')
