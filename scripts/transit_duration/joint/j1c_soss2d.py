"""J1c: own order-1 box extraction from NIRISS SOSS SUBSTRIP96 calints (SCI + INT_TIMES subsets). Own code.
Static bad pixels (along-dispersion outliers / anomalous temporal scatter), per-pixel temporal outlier replacement,
per-column 1/f background from rows > HW+4 px off the order-1 trace, trace y_c(x) from the median frame (poly deg 5),
box |y - y_c| <= HW. Columns selected by wavelength (0.90-2.80 um) using a column->wavelength map from an x1dints file
of the same mode. Regressors: trace y-shift and width, dispersion x-shift, background.
Usage: j1c_soss2d.py <out.npz> <HW> <x1dints_for_wavelength> <seg1_sub.fits> [...]"""
import sys, numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
outp, HW, x1d = sys.argv[1], float(sys.argv[2]), sys.argv[3]; files = sys.argv[4:]
hx = fits.open(x1d, memmap=True); lamcol = None
for x in hx:
    if x.name == 'EXTRACT1D' and x.header.get('SPORDER') == 1:
        lamcol = np.array(x.data['WAVELENGTH'], float); lamcol = lamcol[0] if lamcol.ndim > 1 else lamcol; break
h = fits.open(files[0], memmap=True); S = h['SCI'].data
cube0 = np.array(S[:min(200, S.shape[0])], float); med = np.nanmedian(cube0, 0); tsd = 1.4826 * np.nanmedian(np.abs(cube0 - med), 0); del cube0
ny, nx = med.shape; yy = np.arange(ny)
sm = median_filter(np.nan_to_num(med), size=(1, 9), mode='nearest'); rs = med - sm; rsig = 1.4826 * np.nanmedian(np.abs(rs), 1, keepdims=True)
tsm = median_filter(np.nan_to_num(tsd), size=(1, 9), mode='nearest')
BAD = ~np.isfinite(med) | (np.abs(rs) > 10 * rsig + 0.1 * np.abs(sm)) | (tsd > 5 * tsm + 1e-30)
print('static bad', int(BAD.sum()), 'frame', med.shape, flush=True)
m2 = np.where(BAD, 0, med)
# order-1 trace: brightest smoothed row per column, refined by a flux-weighted centroid, polynomial fit
yc_raw = np.full(nx, np.nan)
for c in range(4, nx - 4):
    p = median_filter(m2[:, c], 5); j = int(np.argmax(p)); a, b = max(0, j - 8), min(ny, j + 9)
    w = np.clip(m2[a:b, c], 0, None)
    if w.sum() > 0: yc_raw[c] = np.sum(w * yy[a:b]) / w.sum()
lam = lamcol[:nx] if lamcol is not None and len(lamcol) >= nx else np.full(nx, np.nan)
use = (lam >= 0.90) & (lam <= 2.80) & np.isfinite(yc_raw)
xs_ = np.where(use)[0]; cc = np.polyfit(xs_, yc_raw[use], 5)
for it in range(3):
    rr = yc_raw[use] - np.polyval(cc, xs_); g = np.abs(rr) < 3 * np.std(rr) + 0.2; cc = np.polyfit(xs_[g], yc_raw[use][g], 5)
yc = np.polyval(cc, np.arange(nx)); Y = yy[:, None] - yc[None, :]
use &= (yc - HW >= 0) & (yc + HW <= ny - 1)
apm = (np.abs(Y) <= HW) & use[None, :]; bkm = (np.abs(Y) > HW + 4) & use[None, :] & ~BAD
print('columns used', use.sum(), 'lam %.3f-%.3f' % (np.nanmin(lam[use]), np.nanmax(lam[use])), 'ap pix/col %.1f' % (apm.sum() / use.sum()), flush=True)
medap = np.where(apm & ~BAD, med, np.nan)
T_all, F, Fcol, ysh, wid, bk = [], [], [], [], [], []
for f in files:
    h = fits.open(f, memmap=True); S = h['SCI'].data; T = np.array(h['INT_TIMES'].data['int_mid_BJD_TDB'], float)
    n = S.shape[0]
    for a in range(0, n, 50):
        blk = np.array(S[a:a + 50], float); blk[:, BAD] = np.nan
        bm = np.nanmedian(blk, 0); bs = 1.4826 * np.nanmedian(np.abs(blk - bm), 0) + 1e-30
        bad = np.abs(blk - bm) > 6 * bs; blk[bad] = np.broadcast_to(bm, blk.shape)[bad]
        for k in range(blk.shape[0]):
            im = blk[k]
            bg = np.nanmedian(np.where(bkm, im, np.nan), 0); im = im - np.nan_to_num(bg)[None, :]
            ap = np.where(apm, im, np.nan); nanap = apm & ~np.isfinite(ap)
            if nanap.any():
                okp = apm & np.isfinite(ap)
                sc = np.nansum(np.where(okp, ap, 0), 0) / np.nansum(np.where(okp, medap, 0), 0); ap[nanap] = (medap * sc[None, :])[nanap]
            col = np.nansum(ap, 0); wp = np.clip(np.nan_to_num(ap), 0, None)
            wcol = wp.sum(0); ym = (wp * Y).sum(0) / np.where(wcol > 0, wcol, np.nan); y2 = (wp * Y ** 2).sum(0) / np.where(wcol > 0, wcol, np.nan)
            wt = wcol[use]
            ysh.append(np.nansum(ym[use] * wt) / np.nansum(wt)); wid.append(np.sqrt(np.nansum((y2[use] - ym[use] ** 2) * wt) / np.nansum(wt)))
            bk.append(np.nanmedian(bg[use])); F.append(np.sum(col[use])); Fcol.append(col[use].astype(np.float32))
        T_all.append(T[a:a + 50])
    print(f.split('/')[-1], n, flush=True)
T_all = np.concatenate(T_all); F = np.array(F); Fcol = np.array(Fcol)
ms = np.median(Fcol, 0)
def hp(s): return s / np.convolve(s, np.ones(31) / 31, mode='same') - 1
ref = hp(ms)[20:-20]; xsh = []
for s in Fcol:
    v = hp(s)[20:-20]; lags = np.arange(-3, 4); c_ = [np.sum(ref[3:-3] * np.roll(v, l)[3:-3]) for l in lags]
    j = int(np.argmax(c_)); j = min(max(j, 1), len(lags) - 2); y0, y1, y2_ = c_[j - 1], c_[j], c_[j + 1]
    xsh.append(lags[j] + 0.5 * (y0 - y2_) / (y0 - 2 * y1 + y2_))
F = F / np.median(F)
np.savez(outp, t=T_all + 2400000.5 - 2457000.0, f=F, yshift=np.array(ysh), width=np.array(wid), xshift=np.array(xsh), bkg=np.array(bk), tilt=np.zeros(len(F)),
         lam_min=np.nanmin(lam[use]), lam_max=np.nanmax(lam[use]), ncol=use.sum(), hw=HW)
print('saved', outp, 'n', len(F), 'pt-to-pt rms %.0f ppm' % (1e6 * np.std(np.diff(F)) / np.sqrt(2)))
