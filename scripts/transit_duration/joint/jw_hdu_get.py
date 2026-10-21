"""Fetch selected HDUs (by EXTNAME) of a public JWST FITS product from the STScI S3 mirror with HTTP range requests,
and write them as a smaller valid FITS file (primary + chosen HDUs). Saves bandwidth (e.g. calints SCI only).
Usage: jw_hdu_get.py <outdir> <EXT1,EXT2,...> <file1> [file2 ...]"""
import sys, os, subprocess
outdir, want = sys.argv[1], sys.argv[2].split(','); os.makedirs(outdir, exist_ok=True)
def rng(url, a, b, dest=None):
    """Range fetch; large ranges go in 8 MB chunks, each retried, aborting connections slower than 30 kB/s for 20 s."""
    CH = 8 * 1024 * 1024
    if dest and b - a + 1 > CH:
        with open(dest, 'wb') as fo:
            for x in range(a, b + 1, CH):
                y = min(b, x + CH - 1); part = dest + '.c'
                rng(url, x, y, part); fo.write(open(part, 'rb').read()); os.remove(part)
        return True
    for k in range(12):
        args = ['curl', '-sS', '-f', '-L', '--max-time', '900', '--speed-limit', '30000', '--speed-time', '20', '-r', '%d-%d' % (a, b)]
        if dest: r = subprocess.run(args + ['-o', dest, url])
        else: r = subprocess.run(args + [url], capture_output=True)
        if r.returncode == 0:
            if dest:
                if os.path.getsize(dest) == b - a + 1: return True
            elif len(r.stdout) == b - a + 1: return r.stdout
    raise RuntimeError('range fetch failed %s %d-%d' % (url, a, b))
def header_at(url, pos):
    raw = b''
    while True:
        blk = rng(url, pos + len(raw), pos + len(raw) + 2879); raw += blk
        cards = [blk[i * 80:(i + 1) * 80].decode('ascii', 'replace') for i in range(36)]
        if any(c.rstrip() == 'END' for c in cards): break
    kv = {}
    for i in range(len(raw) // 80):
        c = raw[i * 80:(i + 1) * 80].decode('ascii', 'replace')
        if c[8:10] == '= ': kv[c[:8].strip()] = c[10:].split('/')[0].strip().strip("'").strip()
    return raw, kv
for fn in sys.argv[3:]:
    out = os.path.join(outdir, fn.replace('.fits', '_sub.fits'))
    if os.path.exists(out): print('have', out, flush=True); continue
    url = 'https://stpubdata.s3.amazonaws.com/jwst/public/%s/%s/%s' % (fn[:7], fn[:13], fn)
    pos = 0; hdus = []; first = True
    while True:
        try: raw, kv = header_at(url, pos)
        except RuntimeError: break
        n = abs(int(kv.get('BITPIX', 8))) // 8; nax = int(kv.get('NAXIS', 0))
        size = 0
        if nax > 0:
            size = n
            for a in range(1, nax + 1): size *= int(kv['NAXIS%d' % a])
        size += int(kv.get('PCOUNT', 0)) * (1 if nax > 0 else 0)
        dpos = pos + len(raw); padded = (size + 2879) // 2880 * 2880
        name = kv.get('EXTNAME', 'PRIMARY' if first else '?')
        hdus.append((name, pos, raw, dpos, size, padded)); print(fn, name, size, flush=True)
        pos = dpos + padded; first = False
    tmp = out + '.part'
    with open(tmp, 'wb') as fo:
        for name, hpos, raw, dpos, size, padded in hdus:
            if name != 'PRIMARY' and name not in want: continue
            fo.write(raw)
            if padded:
                part = tmp + '.blk'; rng(url, dpos, dpos + padded - 1, part)
                fo.write(open(part, 'rb').read()); os.remove(part)
    os.rename(tmp, out); print('ok', out, os.path.getsize(out), flush=True)
