"""Build the light curve given to TRICERATOPS: TTV-corrected fold of the jointly detrended epochs, 10-min bins.

usage: python scripts/triceratops_prep.py <tic> --fit <fit json> [--window 1.0] [--order 2] [--drop 93,100] [--tag _J1.0] [--work ./work]
input:  <work>/<tic>/cut/lc_s<sector>.npz (r = 2.0 px), the SPOC 2-min SAP light curve for S13 of TIC 324429952,
        per-epoch mid-transit times from <work>/<tic>/<fit json>
output: <work>/<tic>/lc_<tic><tag>.npz with t (days from mid-transit, binned), f (binned flux), sig (out-of-transit
        scatter of the bins, the single flux error given to TRICERATOPS), the unbinned points and the fitted shape.
Per epoch, a polynomial baseline (order --order over +-window d) is solved linearly together with a common batman
transit shape (k, b, a/R* fitted here, mid-transit times fixed); detrended flux = flux / baseline. Points within
+-0.5 d of mid-transit are folded on the per-epoch times and averaged in 10-min bins with inverse-variance weights
(per-epoch residual rms); 5-sigma outliers are clipped. The nominal runs use window 1.0 d for TIC 402710280,
0.6 d for TIC 324429952 and 0.8 d for TIC 89247735 (README). Raw aperture flux, not dilution corrected (TRICERATOPS corrects for the aperture).
"""
import argparse, glob, json, os
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from targets import TARGETS, workdir, add_common_args, import_batman, lab, sector_of

batman = import_batman()
p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--fit', required=True); p.add_argument('--window', type=float, default=1.0)
p.add_argument('--order', type=int, default=2); p.add_argument('--drop', default=''); p.add_argument('--tag', default='')
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); tic = args.tic; DROP = [lab(x) for x in args.drop.split(',') if x]
BIN = 10 / 1440.; HALF = 0.5
tcs = {lab(k): v['tc'] for k, v in json.load(open(os.path.join(wd, args.fit))).items() if lab(k) not in DROP}
E = {}
for s, tc in sorted(tcs.items(), key=lambda kv: c['ncyc'][kv[0]]):
    sec = sector_of(s)
    if sec == c.get('spoc_2min'):
        d = fits.open(glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{sec:04d}-{tic:016d}*', '*lc.fits'))[0])[1].data
        t = np.array(d['TIME'], float); f = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY']); src = 'SPOC 2-min SAP'
    else:
        z = np.load(os.path.join(wd, 'cut', f'lc_s{sec}.npz')); t = z['tm']; f = z['f2.0']; q = z['q']; src = 'TESScut r=2.0 px'
    m = (np.abs(t - tc) < args.window) & np.isfinite(f) & (q == 0); t, f = t[m], f[m]; f = f / np.median(f); cad = np.median(np.diff(t))
    pp = batman.TransitParams(); pp.t0 = tc; pp.per = c['P']; pp.rp = 0.04; pp.a = 50; pp.inc = 89.5; pp.ecc = 0; pp.w = 90; pp.u = c['u']; pp.limb_dark = 'quadratic'
    E[s] = dict(t=t, f=f, dt=t - tc, p=pp, mod=batman.TransitModel(pp, t, supersample_factor=max(1, int(round(cad * 1440 / 2))), exp_time=cad), src=src, cad=cad * 1440)


def ep(s, k, b, aR):
    e = E[s]; pp = e['p']; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(min(b / aR, 1))); mm = e['mod'].light_curve(pp)
    A = np.c_[[mm * e['dt'] ** j for j in range(args.order + 1)]].T; cc = np.linalg.lstsq(A, e['f'], rcond=None)[0]
    return e['f'] - A @ cc, np.c_[[e['dt'] ** j for j in range(args.order + 1)]].T @ cc


sig = {s: 1.0 for s in E}; res = lambda x: np.concatenate([ep(s, *x)[0] / sig[s] for s in E])
x = least_squares(res, [0.04, 0.4, 50], bounds=([0.01, 0, 5], [0.2, 1.0, 300])).x
for s in E:
    sig[s] = np.std(ep(s, *x)[0])
k, b, aR = least_squares(res, x, bounds=([0.01, 0, 5], [0.2, 1.0, 300])).x
T14 = c['P'] / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / aR / np.sin(np.arccos(b / aR))) * 24
T = []; F = []; W = []; info = {}
for s, e in E.items():
    r, base = ep(s, k, b, aR); fd = e['f'] / base; keep = (np.abs(e['dt']) < HALF) & (np.abs(r) < 5 * sig[s])
    T.append(e['dt'][keep]); F.append(fd[keep]); W.append(np.full(keep.sum(), 1 / sig[s] ** 2)); info[s] = dict(src=e['src'], cadence_min=round(e['cad'], 1), rms_ppm=round(sig[s] * 1e6), clipped=int((~keep & (np.abs(e['dt']) < HALF)).sum()))
T = np.concatenate(T); F = np.concatenate(F); W = np.concatenate(W); edges = np.arange(-HALF, HALF + 1e-9, BIN); idx = np.digitize(T, edges) - 1; tb = []; fb = []
for i in range(len(edges) - 1):
    w = idx == i
    if w.sum() >= 3:
        tb.append(np.average(T[w], weights=W[w])); fb.append(np.average(F[w], weights=W[w]))
tb = np.array(tb); fb = np.array(fb); o = np.abs(tb) > T14 / 24 / 2 + 0.02; sg = np.std(fb[o] - 1, ddof=1); r = fb[o] - 1
print(f'TIC {tic}{args.tag}: epochs {list(E)}; k={k:.4f} (k^2={k * k * 1e6:.0f} ppm) b={b:.2f} a/R*={aR:.1f} T14={T14:.2f} h; {len(tb)} bins, out-of-transit bin scatter {sg * 1e6:.0f} ppm, lag-1 autocorrelation {np.corrcoef(r[:-1], r[1:])[0, 1]:.2f}\n   {json.dumps(info)}')
np.savez(os.path.join(wd, f'lc_{tic}{args.tag}.npz'), t=tb, f=fb, sig=sg, dep=k * k, P=c['P'], t_unb=T, f_unb=F, k=k, b=b, aR=aR, T14h=T14,
         note=f'TIC {tic}; {"sectors" if all(isinstance(x, int) for x in E) else "epochs"} {list(E)}; TTV-corrected fold on per-epoch tc from {args.fit}; per-epoch order-{args.order} baseline over +-{args.window} d solved jointly with common transit shape; 10-min inverse-variance-weighted phase bins within +-{HALF} d; raw aperture flux (not dilution corrected)')
