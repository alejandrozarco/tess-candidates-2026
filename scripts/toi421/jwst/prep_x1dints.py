"""TOI-421 b (TIC 94986319): build a white-light light curve from MAST JWST x1dints (program 1935)
by summing the 1D spectra per integration over a fixed wavelength range, using only wavelength
columns that are finite and DQ==0 in >=99% of integrations (same column set for all integrations).
Data are read as FITS tables only. Time = TDB-MID (BMJD_TDB, mid-integration).
Usage: python -I prep_x1dints.py <out.npz> <label> <ld_band> <order> <wmin_um> <wmax_um> <meta> <x1dints...>
"""
import sys
import numpy as np
from astropy.io import fits

out, label, ld_band, order, wmin, wmax, meta = sys.argv[1:8]
order, wmin, wmax = int(order), float(wmin), float(wmax)
files = sys.argv[8:]
T, F, E, W = [], [], [], None
effint = None
for fn in sorted(files):
    with fits.open(fn, memmap=True) as h:
        effint = h[0].header.get("EFFINTTM")
        for ext in h:
            if ext.name != "EXTRACT1D": continue
            sp = ext.header.get("SPORDER", 1)
            if sp != order: continue
            d = ext.data
            T.append(np.array(d["TDB-MID"], float))
            F.append(np.array(d["FLUX"], float))
            E.append(np.array(d["FLUX_ERROR"], float))
            dq = np.array(d["DQ"])
            F[-1][dq != 0] = np.nan
            if W is None: W = np.nanmedian(np.array(d["WAVELENGTH"], float), axis=0)
t = np.concatenate(T); fl = np.concatenate(F); er = np.concatenate(E)
o = np.argsort(t); t, fl, er = t[o], fl[o], er[o]
inr = (W >= wmin) & (W <= wmax)
good_col = inr & (np.mean(np.isfinite(fl) & np.isfinite(er), axis=0) >= 0.99)
# drop columns with non-positive median flux or robust relative scatter > 5x the typical column
cmed = np.nanmedian(fl, axis=0)
cmad = 1.4826 * np.nanmedian(np.abs(fl - cmed), axis=0)
with np.errstate(invalid="ignore", divide="ignore"):
    rel = cmad / cmed
good_col &= (cmed > 0)
good_col &= rel < 5 * np.nanmedian(rel[good_col])
fs = fl[:, good_col]; es = er[:, good_col]
ok_row = np.all(np.isfinite(fs), axis=1)
# replace residual NaNs in a row by column median scaled to row (only rows with <1% missing)
frac_bad = np.mean(~np.isfinite(fs), axis=1)
keep = frac_bad < 0.01
colmed = np.nanmedian(fs, axis=0)
for i in np.where(keep & ~ok_row)[0]:
    b = ~np.isfinite(fs[i]); s = np.nanmedian(fs[i][~b] / colmed[~b])
    fs[i][b] = colmed[b] * s; es[i][b] = np.nanmedian(es[:, b], axis=0)
wl = np.sum(fs, axis=1); we = np.sqrt(np.sum(es ** 2, axis=1))
t, wl, we = t[keep], wl[keep], we[keep]
# 5-sigma running-median outlier clip (31-point window)
from scipy.ndimage import median_filter
rm = median_filter(wl, 31, mode="nearest")
dev = wl - rm; mad = 1.4826 * np.median(np.abs(dev - np.median(dev)))
ok = np.abs(dev) < 5 * mad
nclip = int((~ok).sum())
t, wl, we = t[ok], wl[ok], we[ok]
med = np.median(wl)
dt = float(np.median(np.diff(t)))
ptp = 1.4826 * np.median(np.abs(np.diff(wl) - np.median(np.diff(wl)))) / np.sqrt(2)
err_note = "pipeline FLUX_ERROR propagated"
if np.median(we) > 5 * ptp or not np.all(np.isfinite(we)):
    we = np.full_like(wl, ptp)
    err_note = "pipeline FLUX_ERROR unusable (median >5x point-to-point scatter); constant point-to-point MAD scatter used"
meta = meta + "; errors: " + err_note
meta_full = (meta + f"; built from {len(files)} x1dints files (order {order}); wavelength {W[good_col].min():.3f}-"
             f"{W[good_col].max():.3f} um, {good_col.sum()} columns; rows dropped (>1% bad px) {int((~keep).sum())}; "
             f"5-sigma clipped {nclip}; time=BMJD_TDB TDB-MID; EFFINTTM={effint} s")
np.savez(out, t=t, f=wl / med, e=we / med, exptime=(effint / 86400.0 if effint else dt), meta=meta_full,
         label=label, ld_band=ld_band, wave_range=np.array([W[good_col].min(), W[good_col].max()]))
print(meta_full); print(len(t), "points; cadence", dt * 86400, "s; span", (t[-1] - t[0]) * 24, "h; err ppm", np.median(we / med) * 1e6)
