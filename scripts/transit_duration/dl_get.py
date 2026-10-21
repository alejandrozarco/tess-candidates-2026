import sys, json, os, subprocess
prods = json.load(open(sys.argv[1])); outdir = sys.argv[2]; kinds = sys.argv[3].split(',')
os.makedirs(outdir, exist_ok=True)
for o in prods:
    fn=o['fn']
    k = 'fasttp' if fn.endswith('fast-tp.fits') else 'tp' if fn.endswith('_tp.fits') else 'fastlc' if fn.endswith('fast-lc.fits') else 'lc'
    if k not in kinds: continue
    p=os.path.join(outdir,fn)
    if os.path.exists(p) and os.path.getsize(p)==o['size']: continue
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+o['uri']
    r=subprocess.run(['curl','-sS','-f','-L','--max-time','900','-o',p+'.part',url])
    if r.returncode==0 and os.path.getsize(p+'.part')==o['size']: os.rename(p+'.part',p); print('ok',fn,flush=True)
    else: print('FAIL',fn,r.returncode,flush=True)
