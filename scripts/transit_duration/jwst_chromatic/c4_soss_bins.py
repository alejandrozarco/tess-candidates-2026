"""C4: NIRISS SOSS spectroscopic light curves in wavelength bins from the stage-2 x1dints (orders 1 and 2), same
cleaning as ../joint/j1b_soss.py (per-column normalisation, temporal 5-sigma clip against a 21-integration
running median, rejection of columns with scatter > 3x the median, inverse-variance weighted mean per bin).
Regressors as in j1b: background level and the order-1 red/blue 'tilt'. Own code.
Usage: c4_soss_bins.py <out.npz> <x1dints1> [x1dints2 ...]"""
import sys, numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from binsdef import BINS_UM
outp = sys.argv[1]
D = {1: dict(F=[], W=None), 2: dict(F=[], W=None)}; T, B = [], []
for fn in sys.argv[2:]:
    h = fits.open(fn, memmap=True)
    for x in h:
        if x.name != 'EXTRACT1D': continue
        o = x.header.get('SPORDER')
        if o not in D: continue
        d = x.data; D[o]['F'].append(np.array(d['FLUX'], float)); D[o]['W'] = np.array(d['WAVELENGTH'], float)[0]
        if o == 1: T.append(np.array(d['TDB-MID'], float)); B.append(np.nanmedian(np.array(d['BACKGROUND'], float), 1))
T = np.concatenate(T) + 2400000.5 - 2457000.0; B = np.concatenate(B); srt = np.argsort(T)
out = dict(t=None, bkg=None)
def clean(F, W, lo, hi):
    med = np.nanmedian(F, 0); ok = np.isfinite(med) & (med > 0) & (W >= lo) & (W < hi)
    R = F[:, ok] / med[ok]
    rm = median_filter(np.nan_to_num(R, nan=1.0), size=(21, 1), mode='nearest'); d = R - rm; sc = 1.4826 * np.nanmedian(np.abs(d), 0)
    R = np.where(np.abs(d) < 5 * sc[None, :], R, np.nan)
    good = np.isfinite(sc) & (sc < 3 * np.nanmedian(sc)); R = R[:, good]; w = 1 / sc[good] ** 2
    return np.nansum(R * w, 1) / np.nansum(np.isfinite(R) * w, 1), int(good.sum())
F1 = np.vstack(D[1]['F'])[srt]; W1 = D[1]['W']
blue, _ = clean(F1, W1, 0.90, 1.5); red, _ = clean(F1, W1, 1.5, 2.80)
wl, _ = clean(F1, W1, 0.90, 2.80)
fin = np.isfinite(wl) & np.isfinite(blue) & np.isfinite(red)
out = dict(t=T[srt][fin], bkg=B[srt][fin], tilt=(red / blue)[fin], f_white=wl[fin] / np.nanmedian(wl[fin]))
for key, o in (('SOSS1', 1), ('SOSS2', 2)):
    F = np.vstack(D[o]['F'])[srt]; W = D[o]['W']
    for lo, hi in BINS_UM[key]:
        f, n = clean(F, W, lo, hi); f = f[fin]
        f = np.where(np.isfinite(f), f, np.nanmedian(f))
        out['f_%.3f_%.3f' % (lo, hi)] = f / np.nanmedian(f); out['n_%.3f_%.3f' % (lo, hi)] = n
        s = 1.4826 * np.nanmedian(np.abs(np.diff(f))) / np.sqrt(2)
        print(key, '%.3f-%.3f' % (lo, hi), 'columns', n, 'pt-to-pt %.0f ppm' % (s * 1e6))
np.savez(outp, **out); print('saved', outp, len(out['t']))
