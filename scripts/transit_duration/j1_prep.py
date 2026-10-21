"""Prepare JWST pipeline white-light curves: pick column, normalise, clip, trim. Usage: j1_prep.py <ecsv> <col> <out.npz> [trim_min]"""
import sys, numpy as np
from astropy.table import Table
f, col, out = sys.argv[1:4]; trim = float(sys.argv[4]) if len(sys.argv) > 4 else 30.
t = Table.read(f, format='ascii.ecsv')
tcol = [c for c in t.colnames if c.startswith('BJD_TDB')]
tc = 'BJD_TDB_' + col.split('_')[-1] if ('BJD_TDB_' + col.split('_')[-1]) in t.colnames else tcol[0]
x = np.array(t[col], float); tt = np.array(t[tc], float) + 2400000.5 - 2457000.0  # MJD-like BJD -> BTJD
good = np.isfinite(x) & np.isfinite(tt)
med = np.nanmedian(x[good]); y = x / med
# running-median clip
from scipy.ndimage import median_filter
rm = median_filter(np.where(good, y, 1.0), size=31, mode='nearest')
res = y - rm; s = 1.4826 * np.nanmedian(np.abs(res[good]))
good &= np.abs(res) < 5 * s
good &= tt > tt[good].min() + trim / 1440.
print(f.split('/')[-1], col, 'n', len(x), 'kept', good.sum(), 'scatter_ppm %.0f' % (s * 1e6), 'span_h %.2f' % ((tt[good].max() - tt[good].min()) * 24), 'cad_s %.1f' % (np.median(np.diff(tt)) * 86400))
np.savez(out, t=tt[good], f=y[good] / np.median(y[good]), e=np.full(good.sum(), s))
