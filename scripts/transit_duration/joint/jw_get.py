"""Download named JWST products from the public STScI S3 mirror (stpubdata) into an untrusted dir.
Usage: jw_get.py <outdir> <file1> [file2 ...]"""
import sys, os, subprocess
outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
for fn in sys.argv[2:]:
    p = os.path.join(outdir, fn)
    if os.path.exists(p): print('have', fn, flush=True); continue
    url = 'https://stpubdata.s3.amazonaws.com/jwst/public/%s/%s/%s' % (fn[:7], fn[:13], fn)
    r = subprocess.run(['curl', '-sS', '-f', '-L', '--retry', '3', '--max-time', '5400', '-o', p + '.part', url])
    if r.returncode == 0: os.rename(p + '.part', p); print('ok', fn, os.path.getsize(p), flush=True)
    else: print('FAIL', fn, r.returncode, flush=True)
