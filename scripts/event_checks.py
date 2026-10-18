"""Checks of individual non-transit box events: depth per aperture, background, centroid, flagged cadences, QLP SAP.

usage: python scripts/event_checks.py <tic> --events 13:1682.134,93:3849.745 [--work ./work]
input:  <work>/<tic>/cut/lc_s<sector>.npz (aperture_photometry.py); optionally the QLP light curve of the sector in
        <work>/<tic>/dl/mastDownload/HLSP/hlsp_qlp_tess_ffi_s<sector>-*/*.fits (MAST, SAP_FLUX, quality 0)
output: printed table and <work>/<tic>/event_checks.json
For every event time a box of width 0.8 x T14 is compared with a quadratic baseline fitted to the cadences
|dt| > T14/2 + 0.04 d within +-1 d (quality 0), for the r = 1.5, 2.0 and 3.0 px apertures, the annulus background
and the flux-weighted centroid (x, y; relative units). A dip on the target has the same depth in every aperture,
no background or centroid excursion and appears in QLP SAP; a dip that grows with the aperture radius, moves the
centroid or coincides with flagged cadences or a background excursion is not on the target.
"""
import argparse, glob, json, os
import numpy as np
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--events', required=True)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); HW = c['T14h'] / 48; out = []


def box(t, f, tc, hw=HW * 0.8):
    m = np.isfinite(f) & (np.abs(t - tc) < 1.0); t, f = t[m], f[m]; dt = t - tc; win = np.abs(dt) < hw; oot = np.abs(dt) > HW + 0.04
    if win.sum() < 3 or oot.sum() < 20:
        return None, None
    cc = np.polyfit(dt[oot], f[oot], 2); r = f / np.polyval(cc, dt)
    return round(float((1 - np.mean(r[win])) * 1e6)), round(float(np.std(r[oot]) / np.sqrt(win.sum()) * 1e6))


for ev in args.events.split(','):
    s, tc = int(ev.split(':')[0]), float(ev.split(':')[1]); z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); q = z['q']; t = z['tm']; m = q == 0; row = dict(sector=s, tc=tc)
    for ap in ('f1.5', 'f2.0', 'f3.0'):
        row[ap] = box(t[m], z[ap][m], tc)
    row['bkg_rel'] = box(t[m], z['bkg'][m] / np.median(z['bkg'][m]), tc)
    for k in ('cx', 'cy'):
        v = z[k][m]; ok = np.isfinite(v); row[k + '_rel'] = box(t[m][ok], v[ok] / np.nanmedian(v), tc)
    row['n_q0_in_box'] = int((m & (np.abs(t - tc) < HW * 0.8)).sum()); row['n_flagged_in_box'] = int(((q != 0) & (np.abs(t - tc) < HW * 0.8)).sum())
    qf = glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'HLSP', f'hlsp_qlp_tess_ffi_s{s:04d}-*', '*.fits'))
    if qf:
        from astropy.io import fits
        d = fits.open(qf[0])[1].data; tq = np.array(d['TIME'], float); fq = np.array(d['SAP_FLUX'], float); mq = (np.array(d['QUALITY']) == 0) & np.isfinite(fq)
        row['qlp_sap'] = box(tq[mq], fq[mq], tc)
    print(row); out.append(row)
json.dump(out, open(os.path.join(wd, 'event_checks.json'), 'w'), indent=1)
