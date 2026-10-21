"""Convert the own TPF aperture light curves of p1_tpf.py into the light-curve format of v1_extract.py, so that
v2_windows.py and the later steps can be run on them unchanged (pdc = sap = aperture flux / median; errors = the
point-to-point scatter). Usage: p3_aplc.py <tic> <aperture, e.g. r1.5> <out_dir> <p1 prefix> [<p1 prefix> ...]"""
import sys, os, numpy as np
tic, ap, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
for pre in sys.argv[4:]:
    z = np.load(pre + '_lc.npz'); t = z['t']; f = z[ap]; ok = np.isfinite(f) & np.isfinite(t)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok]); e = np.full(len(f), np.std(np.diff(f)) / np.sqrt(2))
    sec = int(os.path.basename(pre).lstrip('s'))
    meta = dict(sector=sec, tic=tic, fast=False, crowd=np.nan, flfr=np.nan, procver='own_tpf_' + ap, dataver=-1, exptime=120.0, cam=0, ccd=0)
    np.savez(os.path.join(out, '%d_s%03d_2m.npz' % (tic, sec)), t=t, pdc=f, pdce=e, sap=f, sape=e, bkg=np.zeros(len(f)), meta=np.array([meta], dtype=object))
    print(sec, ap, len(t), 'rms_ppm %.0f' % (e[0] * 1e6))
