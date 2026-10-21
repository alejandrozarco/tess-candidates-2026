"""C1: copy of ../joint/j1_extract.py (section 6.6.10); the only change is that the per-column aperture fluxes (Fcol,
one column per selected detector column, wavelengths in lam) are also saved, so that spectroscopic light curves can be
binned later (c3_binfit.py). The white light f is identical to that of j1_extract.py.
J1: own white-light extraction from NIRSpec G395H BOTS calints (SCI, WAVELENGTH, INT_TIMES subsets). Own code.
Per segment: static bad pixels masked (excluded from the aperture sum), temporal cosmic-ray/outlier replacement per pixel, per-column background (1/f) from rows >= 5 px off the
trace, trace centre y_c(x) from the median frame (poly deg 3), box aperture |y - y_c| <= HW, other NaN pixels filled from
the scaled median profile. Regressors per integration: trace y-shift and width (flux-weighted moments), dispersion
x-shift (cross-correlation of the 1D spectrum with the median spectrum), background level.
Usage: c1_extract_cols.py <out.npz> <lam_min> <lam_max> <HW> <seg1_sub.fits> [seg2_sub.fits ...]"""
import sys, numpy as np
from astropy.io import fits
outp, lmin, lmax, HW = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
files = sys.argv[5:]
# reference trace and median frame from the first segment
def load(f):
    h = fits.open(f, memmap=True)
    return h, h['SCI'].data, np.array(h['WAVELENGTH'].data, float), np.array(h['INT_TIMES'].data['int_mid_BJD_TDB'], float)
h, S, W, T = load(files[0])
cube0 = np.array(S[:min(300, S.shape[0])], float)
med = np.nanmedian(cube0, 0); tsd = 1.4826 * np.nanmedian(np.abs(cube0 - med), 0); del cube0
ny, nx = med.shape; yy = np.arange(ny)
# static bad pixels: outliers against an along-dispersion running median, or with anomalous temporal scatter
from scipy.ndimage import median_filter
sm = median_filter(np.nan_to_num(med, nan=0.), size=(1, 9), mode='nearest')
rs = med - sm; rsig = 1.4826 * np.nanmedian(np.abs(rs), 1, keepdims=True)
tsm = median_filter(np.nan_to_num(tsd, nan=0.), size=(1, 9), mode='nearest')
BAD = np.isfinite(med) & ((np.abs(rs) > 10 * rsig + 0.1 * np.abs(sm)) | (tsd > 5 * tsm + 1e-30))
print('static bad pixels', int(BAD.sum()), flush=True)
med[BAD] = np.nan
cols = np.where(np.sum(np.isfinite(med), 0) > 8)[0]
yc_raw = np.full(nx, np.nan)
for c in cols:
    p = np.nan_to_num(med[:, c], nan=0.); j = np.argmax(p); a, b = max(0, j - 3), min(ny, j + 4)
    if p[a:b].sum() > 0: yc_raw[c] = np.sum(p[a:b] * yy[a:b]) / p[a:b].sum()
ok = np.isfinite(yc_raw)
cc = np.polyfit(np.where(ok)[0], yc_raw[ok], 3)
for it in range(3):
    rr = yc_raw[ok] - np.polyval(cc, np.where(ok)[0]); good = np.abs(rr) < 4 * np.std(rr) + 0.05
    cc = np.polyfit(np.where(ok)[0][good], yc_raw[ok][good], 3)
yc = np.polyval(cc, np.arange(nx))
lam = np.array([W[int(round(yc[c])), c] if 0 <= round(yc[c]) < ny else np.nan for c in range(nx)])
use = (lam >= lmin) & (lam <= lmax) & np.isfinite(lam)
Y = yy[:, None] - yc[None, :]
apm = (np.abs(Y) <= HW) & use[None, :]
bkm = (np.abs(Y) >= HW + 2) & use[None, :]
print('cols used', use.sum(), 'lam %.3f-%.3f' % (np.nanmin(lam[use]), np.nanmax(lam[use])), 'ap pix/col %.1f' % (apm.sum() / use.sum()), flush=True)
medap = np.where(apm, med, np.nan)
T_all, F, Fcol, yshift, width, bk = [], [], [], [], [], []
spec_med = None
for f in files:
    h, S, W, T = load(f)
    n = S.shape[0]
    for a in range(0, n, 100):
        blk = np.array(S[a:a + 100], float); blk[:, BAD] = np.nan
        # temporal outlier replacement against the block median (5-sigma, robust per pixel)
        bm = np.nanmedian(blk, 0); bs = 1.4826 * np.nanmedian(np.abs(blk - bm), 0) + 1e-30
        bad = np.abs(blk - bm) > 6 * bs
        blk[bad] = np.broadcast_to(bm, blk.shape)[bad]
        for k in range(blk.shape[0]):
            im = blk[k]
            bg = np.nanmedian(np.where(bkm, im, np.nan), 0)            # per-column background (1/f)
            im = im - np.nan_to_num(bg)[None, :]
            ap = np.where(apm, im, np.nan)
            # fill NaNs in the aperture from the scaled median profile, per column
            nanap = apm & ~np.isfinite(ap)
            if nanap.any():
                sc = np.nansum(np.where(apm & np.isfinite(ap), ap, 0), 0) / np.nansum(np.where(apm & np.isfinite(ap), medap, 0), 0)
                ap[nanap] = (medap * sc[None, :])[nanap]
            col = np.nansum(ap, 0)
            wcol = np.nansum(np.where(apm, np.clip(np.nan_to_num(ap), 0, None), 0), 0)
            ym = np.nansum(np.where(apm, np.clip(np.nan_to_num(ap), 0, None) * Y, 0), 0) / np.where(wcol > 0, wcol, np.nan)
            y2 = np.nansum(np.where(apm, np.clip(np.nan_to_num(ap), 0, None) * Y ** 2, 0), 0) / np.where(wcol > 0, wcol, np.nan)
            wt = wcol[use]
            yshift.append(np.nansum(ym[use] * wt) / np.nansum(wt)); width.append(np.sqrt(np.nansum((y2[use] - ym[use] ** 2) * wt) / np.nansum(wt)))
            bk.append(np.nanmedian(bg[use])); F.append(np.sum(col[use])); Fcol.append(col[use].astype(np.float32))
        T_all.append(T[a:a + 100])
    print(f.split('/')[-1], n, flush=True)
T_all = np.concatenate(T_all); F = np.array(F); Fcol = np.array(Fcol)
# dispersion-direction shift: cross-correlate high-passed spectra with the median spectrum (parabolic peak)
ms = np.median(Fcol, 0)
def hp(s): return s / np.convolve(s, np.ones(31) / 31, mode='same') - 1
ref = hp(ms)[20:-20]; xs = []
for s in Fcol:
    v = hp(s)[20:-20]; lags = np.arange(-3, 4); cc_ = [np.sum(ref[3:-3] * np.roll(v, l)[3:-3]) for l in lags]
    j = int(np.argmax(cc_)); j = min(max(j, 1), len(lags) - 2); y0, y1, y2_ = cc_[j - 1], cc_[j], cc_[j + 1]
    xs.append(lags[j] + 0.5 * (y0 - y2_) / (y0 - 2 * y1 + y2_))
F = F / np.median(F)
np.savez(outp, t=T_all + 2400000.5 - 2457000.0,   # INT_TIMES BJD_TDB columns are MJD-style (BJD - 2400000.5)
         f=F, yshift=np.array(yshift), width=np.array(width), xshift=np.array(xs), bkg=np.array(bk),
         lam_min=np.nanmin(lam[use]), lam_max=np.nanmax(lam[use]), ncol=use.sum(), hw=HW, spec_med=ms, lam=lam[use], Fcol=Fcol.astype(np.float32), Fwl_raw=np.array(F))
d = np.diff(F); print('saved', outp, 'n', len(F), 'pt-to-pt rms %.0f ppm' % (1e6 * np.std(d) / np.sqrt(2)), 'yshift rms %.4f width rms %.4f xshift rms %.4f' % (np.std(yshift), np.std(width), np.std(xs)))
