"""J1b: own white light from NIRISS SOSS x1dints (order 1, 0.90-2.80 um, same band for both visits). Per-column
normalisation, rejection of noisy columns, temporal (per-column running-median) and spatial outlier rejection,
inverse-variance weighted mean. Regressor proxies: background level and the spectral 'tilt' (red/blue flux ratio).
Own code. Usage: j1b_soss.py <out.npz> <x1dints1> [x1dints2 ...]"""
import sys, numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
outp = sys.argv[1]; F, T, B, W = [], [], [], None
for fn in sys.argv[2:]:
    h = fits.open(fn, memmap=True)
    for x in h:
        if x.name != 'EXTRACT1D' or x.header.get('SPORDER') != 1: continue
        d = x.data
        F.append(np.array(d['FLUX'], float)); T.append(np.array(d['TDB-MID'], float)); B.append(np.nanmedian(np.array(d['BACKGROUND'], float), 1))
        W = np.array(d['WAVELENGTH'], float)[0]
F = np.vstack(F); T = np.concatenate(T) + 2400000.5 - 2457000.0; B = np.concatenate(B)
o = np.argsort(T); F, T, B = F[o], T[o], B[o]
print('cube', F.shape, 'span %.2f h' % ((T[-1] - T[0]) * 24))
med = np.nanmedian(F, 0); ok = np.isfinite(med) & (med > 0) & (W > 0.90) & (W < 2.80)
R = F[:, ok] / med[ok]
# temporal outliers per column against a 21-integration running median
rm = median_filter(np.nan_to_num(R, nan=1.0), size=(21, 1), mode='nearest')
d = R - rm; sc = 1.4826 * np.nanmedian(np.abs(d), 0)
R = np.where(np.abs(d) < 5 * sc[None, :], R, np.nan)
good = np.isfinite(sc) & (sc < 3 * np.nanmedian(sc))
R = R[:, good]; w = 1 / sc[good] ** 2
wl = np.nansum(R * w, 1) / np.nansum(np.isfinite(R) * w, 1)
blue = W[ok][good] < 1.5
tilt = np.nanmean(R[:, ~blue], 1) / np.nanmean(R[:, blue], 1)
s = 1.4826 * np.nanmedian(np.abs(np.diff(wl))) / np.sqrt(2)
print('columns used', good.sum(), 'of', ok.sum(), 'white pt-to-pt %.0f ppm' % (s * 1e6))
fin = np.isfinite(wl)
np.savez(outp, t=T[fin], f=wl[fin] / np.nanmedian(wl[fin]), bkg=B[fin], tilt=tilt[fin], yshift=np.zeros(fin.sum()), xshift=np.zeros(fin.sum()), width=np.zeros(fin.sum()))
