"""Own TPF photometry + difference images for TOI-620. Usage: p1_tpf.py <tpf.fits> <outprefix> <P> <T0>"""
import sys, numpy as np, json
from astropy.io import fits
from astropy.wcs import WCS
f, outp, P, T0 = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
h = fits.open(f, memmap=False)
d = h[1].data; hd = h[1].header; ap = h[2].data
t = np.array(d['TIME'], float); fl = np.array(d['FLUX'], float); bk = np.array(d['FLUX_BKG'], float); q = np.array(d['QUALITY'])
sec = h[0].header['SECTOR']
good = np.isfinite(t) & (q == 0) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
t, fl = t[good], fl[good]
w = WCS(h[2].header)
ra0, de0 = h[0].header['RA_OBJ'], h[0].header['DEC_OBJ']
x0, y0 = w.all_world2pix(ra0, de0, 0)
# neighbour TIC 296739875 position from TIC is not in header; locate brightest other peak in median image
med = np.nanmedian(fl, axis=0)
ny, nx = med.shape
yy, xx = np.mgrid[:ny, :nx]
# local background: median of faintest 30% pixels per cadence
bgp = np.nanpercentile(fl.reshape(len(t), -1), 20, axis=1)
fl2 = fl - bgp[:, None, None]
pipe = (ap & 2) > 0
phase = (t - T0 + 0.5 * P) % P - 0.5 * P
intr = np.abs(phase) < 0.018; oot = (np.abs(phase) > 0.045) & (np.abs(phase) < 0.15)
# difference image per transit (local), averaged
ep = np.round((t - T0) / P)
diffs = []
for e in np.unique(ep[intr]):
    a = intr & (ep == e); b_ = oot & (ep == e)
    if a.sum() > 10 and b_.sum() > 30: diffs.append(np.nanmean(fl2[b_], 0) - np.nanmean(fl2[a], 0))
dimg = np.nanmean(diffs, 0) if diffs else np.full(med.shape, np.nan)
def cen(img):
    img = np.where(np.isfinite(img), img, 0); im = np.clip(img, 0, None)
    return float((im * xx).sum() / im.sum()), float((im * yy).sum() / im.sum())
res = dict(sector=int(sec), target_xy=[float(x0), float(y0)], ntransits=len(diffs), crowd=h[1].header.get('CROWDSAP'), pipe_npix=int(pipe.sum()),
           diff_centroid=cen(dimg), med_centroid_in_pipe=cen(np.where(pipe, med, 0)))
# brightest pixels
order = np.argsort(med.ravel())[::-1][:6]
res['brightest_px'] = [[int(i % nx), int(i // nx), float(med.ravel()[i])] for i in order]
# circular apertures centred on target
lcs = {}
for r in (1.0, 1.5, 2.0, 2.5, 3.0):
    m = (xx - x0) ** 2 + (yy - y0) ** 2 <= r * r
    lcs['r%.1f' % r] = (fl2[:, m].sum(1), int(m.sum()))
lcs['pipe'] = (fl2[:, pipe].sum(1), int(pipe.sum()))
np.savez(outp + '_lc.npz', t=t, **{k: v[0] / np.nanmedian(v[0]) for k, v in lcs.items()})
np.save(outp + '_diff.npy', dimg); np.save(outp + '_med.npy', med)
res['ap_npix'] = {k: v[1] for k, v in lcs.items()}
# fraction of median-image flux in each aperture coming from beyond 2 px of target (crude contamination proxy)
json.dump(res, open(outp + '.json', 'w'), indent=1)
print(json.dumps(res))
