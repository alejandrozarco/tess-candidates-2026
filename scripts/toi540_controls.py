"""Control planets for the TOI-540 b sector-98 test: the same depth measurement on known TOIs observed at 2-min cadence.

usage: python scripts/toi540_controls.py [--work ./work] [--out data]
input:  SPOC 2-min light curves of the control stars in <work>/controls/dl/mastDownload/TESS/ (other_download.py 200322593 --controls)
output: data/control_planets_s98.csv (format in the file header)
For each control (ExoFOP TOI-table values of 2026-10-18 in targets.py): PDCSAP flux, quality 0, biweight detrending
with a window of max(0.25 d, 4 T14), a transit template of the TOI-table duration and depth (b = 0.5, quadratic limb
darkening (0.25, 0.35)), and the peak depth fitted at the TOI-table ephemeris, with the error scaled by the red-noise
factor at T14. Also the best depth within +-max(60 min, 5 sigma of the propagated ephemeris), to allow for ephemeris
drift. Sectors: the first, the last two before sector 98, and sectors 98 and later. TOI-540.01 is measured in the
same way at its TOI-table ephemeris (sectors 4-98) for comparison, and SPOC light-curve diagnostics are listed.
"""
import argparse, csv, glob, os
import numpy as np
from astropy.io import fits
from wotan import flatten
from targets import OTHER, import_batman

batman = import_batman()
p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
args = p.parse_args(); DL = os.path.join(args.work, 'controls', 'dl', 'mastDownload', 'TESS')
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')


def shape(P, T14_min, depth):
    k = np.sqrt(depth * 1e-6); aRs = np.linspace(1.5, 400, 20000); b = 0.5
    t14 = P / np.pi * np.arcsin(np.clip(np.sqrt((1 + k)**2 - b**2) / (aRs * np.sin(np.arccos(b / aRs))), 0, 1)) * 1440
    return k, aRs[np.argmin(np.abs(t14 - T14_min))], b


GRID = {}
def template(t, tc, P, k, aR, b, texp):
    key = (P, k, aR, b, round(texp * 86400))
    if key not in GRID:
        g = np.arange(-0.2, 0.2, 10 / 86400.); pp = batman.TransitParams(); pp.t0 = 0; pp.per = P; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(b / aR)); pp.ecc = 0; pp.w = 90
        pp.u = [0.25, 0.35]; pp.limb_dark = 'quadratic'; sh = 1 - batman.TransitModel(pp, g, supersample_factor=5, exp_time=texp).light_curve(pp); GRID[key] = (g, sh / sh.max())
    g, sh = GRID[key]
    return np.interp((((t - tc) / P + 0.5) % 1 - 0.5) * P, g, sh, left=0, right=0)


def beta_factor(t, r, cad, tbin):
    n = max(int(round(tbin / cad)), 1); rr = r[np.argsort(t)]; nb = len(rr) // n
    return 1.0 if nb < 10 else max(1.0, np.std(rr[:nb * n].reshape(nb, n).mean(1)) / (np.std(r) / np.sqrt(n) * np.sqrt(nb / (nb - 1))))


rows = []; c = OTHER[200322593]; T540 = ('540.01', 200322593, c['toi'][1], 4.8953157e-06, c['toi'][0], 0.0009524355, 0.488774921389765, 2318.66)
for toi, tic, P, eP, E0, eE, T14h, dep in c['controls'] + [T540]:
    E0 -= 2457000; T14 = T14h * 60; k, aR, b = shape(P, T14, dep)
    dl = DL if tic != 200322593 else os.path.join(args.work, str(tic), 'dl', 'mastDownload', 'TESS')
    for fn in sorted(glob.glob(os.path.join(dl, f'*-{tic:016d}-*-s', '*-s_lc.fits')), key=lambda f: int(fits.getheader(f)['SECTOR'])):
        h0 = fits.getheader(fn); h1 = fits.getheader(fn, 1); d = fits.open(fn)[1].data; t = np.array(d['TIME'], float); f = np.array(d['PDCSAP_FLUX'], float); q = np.array(d['QUALITY'])
        m = np.isfinite(t) & np.isfinite(f) & (q == 0)
        diag = [f"{np.nanmedian(np.array(d['SAP_FLUX'], float)[m]):.0f}", f"{np.mean(np.isfinite(np.array(d['MOM_CENTR1'], float)[m])):.2f}", f"{h1.get('CDPP1_0', np.nan):.0f}", h0.get('PROCVER', '')]
        if m.sum() < 1000:
            rows.append([toi, tic, h0['SECTOR'], f"{h0['CAMERA']}/{h0['CCD']}", int(m.sum()), '', '', '', '', '', '', '', '', f'{dep:.0f}'] + diag); continue
        t, f = t[m], f[m] / np.nanmedian(f[m]); cad = float(np.median(np.diff(t)))
        tr = flatten(t, f, method='biweight', window_length=max(0.25, 4 * T14 / 1440), return_trend=True, break_tolerance=0.3, edge_cutoff=0.)[1]
        f = f / tr; ok = np.isfinite(f); sg = 1.4826 * np.median(np.abs(f[ok] - 1)); ok &= (f - 1 < 4 * sg) & (f - 1 > -8 * sg); t, f = t[ok], f[ok]
        En = np.round((np.median(t) - E0) / P); tc = E0 + En * P; stc = np.hypot(eE, En * eP) * 1440; beta = beta_factor(t, f - 1, cad, T14 / 1440)
        def fit(o):
            s = template(t, tc + o / 1440, P, k, aR, b, cad); ss = np.sum(s * s); return -np.sum((f - 1) * s) / ss, sg / np.sqrt(ss) * beta
        A, eA = fit(0.0); W = max(60, 5 * stc); offs = np.arange(-W, W + 0.01, 1.0); fits_ = np.array([fit(o) for o in offs]); i = np.argmax(fits_[:, 0] / fits_[:, 1])
        ntr = len(np.unique(np.round((t[np.abs(((t - tc) / P + 0.5) % 1 - 0.5) * P < T14 / 2880] - tc) / P)))
        rows.append([toi, tic, h0['SECTOR'], f"{h0['CAMERA']}/{h0['CCD']}", int(m.sum()), ntr, f'{stc:.1f}', f'{A*1e6:.0f}', f'{eA*1e6:.0f}', f'{A/eA:.1f}', f'{offs[i]:+.0f}', f'{fits_[i,0]*1e6:.0f}', f'{fits_[i,0]/fits_[i,1]:.1f}', f'{dep:.0f}'] + diag)
        print(rows[-1], flush=True)
with open(os.path.join(args.out, 'control_planets_s98.csv'), 'w') as fo:
    fo.write(''.join('# ' + l + '\n' for l in [
        f'Control planets for the TOI-540 b sector-98 test: SPOC 2-min PDCSAP light curves, same template-depth measurement as for TOI-540 b (scripts/toi540_controls.py), {DATE}',
        'ephemerides, durations and depths: ExoFOP TOI table of 2026-10-18 (targets.py); cam_ccd = TESS camera/CCD; n_q0 = quality-0 cadences; n_transits = predicted transits with data',
        'sigma_tc_min = propagated TOI-table ephemeris error; depth_ppm, err, snr = peak depth of the template at the TOI-table ephemeris (error scaled by beta at T14)',
        'best_offset_min, best_depth_ppm, best_snr = highest-S/N template position within +-max(60 min, 5 sigma_tc) of the ephemeris (ephemeris drift); toi_depth_ppm = TOI-table depth',
        'SPOC light-curve diagnostics: sap_median = median SAP_FLUX (e/s, quality 0); centroid_frac = fraction of quality-0 cadences with a MOM_CENTR1 value; cdpp1_ppm = CDPP1_0 header keyword; procver = SPOC version',
        'last rows: TOI-540.01 itself, measured in the same way at its TOI-table ephemeris (ExoFOP, 2026-10-18)']))
    w = csv.writer(fo); w.writerow(['toi', 'tic', 'sector', 'cam_ccd', 'n_q0', 'n_transits', 'sigma_tc_min', 'depth_ppm', 'depth_err_ppm', 'snr', 'best_offset_min', 'best_depth_ppm', 'best_snr', 'toi_depth_ppm', 'sap_median', 'centroid_frac', 'cdpp1_ppm', 'procver'])
    w.writerows(rows)
print('written', os.path.join(args.out, 'control_planets_s98.csv'))
