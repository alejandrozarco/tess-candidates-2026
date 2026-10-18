"""Write the data/ tables for one target from the fit output and the photometry.

usage: python scripts/make_tables.py <tic> --fit <fit json> [--extra-fit <json> --extra-sector 66] [--window 1.0] [--work ./work] [--out data]
input:  <work>/<tic>/<fit json> (fit_transits.py), <work>/<tic>/cut/lc_s<sector>.npz and s<sector>.npz, the SPOC 2-min
        light curve for sector 13 of TIC 324429952, and <work>/<tic>/D_s<epoch>.npy, diffimg.json (difference_images.py)
output: data/transit_times_<tic>.csv, folded_<tic>.csv, model_<tic>.csv, sectors_<tic>.csv, diffimage_<tic>.csv,
        diffimage_sources_<tic>.csv, diffimage_centroids_<tic>.csv (formats described in the comment headers of the files).
Epochs are labelled by sector, or '<sector>a/b' when a sector holds two transits; the tables carry the sector and the
cycle number (epoch). --window = half-width (d) of the light-curve windows (the fit window of fit_transits.py).
The folded light curve divides each epoch by the quadratic baseline solved jointly with the common transit model
(same model as fit_transits.py; for TIC 402710280 the shape is refitted with the mid-transit times fixed).
"""
import argparse, csv, glob, json, os
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord
from targets import TARGETS, workdir, add_common_args, import_batman, lab, sector_of

batman = import_batman()
p = add_common_args(argparse.ArgumentParser(description=__doc__))
p.add_argument('--fit', required=True); p.add_argument('--extra-fit', default=None); p.add_argument('--extra-sector', type=int, default=None)
p.add_argument('--out', default='data'); p.add_argument('--date', default=''); p.add_argument('--window', type=float, default=1.0)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); tic = args.tic; os.makedirs(args.out, exist_ok=True)
PS = f'{c["P"]:.10g}'; EPS = f'{c["eP"]:.{len(PS.split(".")[1])}f}'   # period and error with the decimals given in targets.py
DATE = args.date or __import__('datetime').datetime.utcnow().strftime('%Y-%m-%d')


def hdr(lines):
    return ''.join('# ' + l + '\n' for l in lines)


def load_epoch(s, tc, win=args.window):
    s = sector_of(s)
    if s == c.get('spoc_2min'):
        fn = glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{s:04d}-{tic:016d}*', '*lc.fits'))[0]
        d = fits.open(fn)[1].data; t = np.array(d['TIME'], float); f = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY']); prod = 'SPOC 2-min SAP'
    else:
        z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); t = z['tm']; f = z['f2.0']; q = z['q']; prod = 'TESScut r=2.0 px'
    m = (np.abs(t - tc) < win) & np.isfinite(f) & (q == 0); t, f = t[m], f[m]
    return t, f / np.median(f), prod


def model(t, tc, k, b, aR, cad):
    pp = batman.TransitParams(); pp.t0 = tc; pp.per = c['P']; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(min(b / aR, 1)))
    pp.ecc = 0; pp.w = 90; pp.u = c['u']; pp.limb_dark = 'quadratic'
    return batman.TransitModel(pp, t, supersample_factor=max(1, int(round(cad * 1440 / 2))), exp_time=cad).light_curve(pp)


def resid_epoch(E, s, k, b, aR):
    e = E[s]; mm = model(e['t'], e['tc'], k, b, aR, e['cad']); x = e['t'] - e['tc']
    A = np.c_[[mm * x ** j for j in range(3)]].T; cc = np.linalg.lstsq(A, e['f'], rcond=None)[0]
    return e['f'] - A @ cc, np.c_[[x ** j for j in range(3)]].T @ cc, mm


J = json.load(open(os.path.join(wd, args.fit))); tcs = {lab(k): v for k, v in J.items()}
ORDER = lambda d: sorted(d.items(), key=lambda kv: c['ncyc'][kv[0]])
# ---- transit times
rows = []
for s, v in ORDER(tcs):
    n = c['ncyc'][s]; prod = 'SPOC 2-min SAP' if sector_of(s) == c.get('spoc_2min') else 'TESScut r=2.0 px'
    rows.append((n, sector_of(s), prod, v, (v['tc'] - c['T0'] - n * c['P']) * 1440, 'used in ephemeris'))
if args.extra_fit:
    v = json.load(open(os.path.join(wd, args.extra_fit)))[str(args.extra_sector)]; n = c['ncyc'][args.extra_sector]
    rows.append((n, args.extra_sector, 'TESScut r=2.0 px', v, (v['tc'] - c['T0'] - n * c['P']) * 1440, 'scattered light; not used in ephemeris'))
rows.sort(key=lambda r: r[0])
with open(os.path.join(args.out, f'transit_times_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic} = Gaia DR3 {c["gaia_dr3"]}: mid-transit times from the common-shape fit (scripts/fit_transits.py), {DATE}',
                  f'Linear ephemeris used for O-C: T0 = BTJD {c["T0"]:.5f} +- {c["eT0"]:.5f}, P = {PS} +- {EPS} d; epoch = cycle number relative to T0',
                  'tc_btjd = BJD_TDB - 2457000. tc_err = max(formal error x beta, prayer-bead error), in days; tc_err_min in minutes',
                  'depth_ppm = per-epoch depth with the common shape (k free, b and a/R* fixed); rms_ppm = per-epoch residual rms',
                  'product: TESScut 21x21 cutout, circular aperture r = 2.0 px, annulus background; or SPOC 2-min SAP_FLUX (S13 of TIC 324429952)',
                  'Not used: QLP light curves (QLP timing error in S99-S104; see README).']))
    w = csv.writer(fo); w.writerow(['epoch', 'sector', 'product', 'tc_btjd', 'tc_bjd_tdb', 'tc_err_d', 'tc_err_min', 'oc_min', 'depth_ppm', 'depth_err_ppm', 'rms_ppm', 'beta', 'prayer_bead_err_min', 'note'])
    for n, s, prod, v, oc, note in rows:
        w.writerow([n, s, prod, f'{v["tc"]:.6f}', f'{v["tc"] + 2457000:.6f}', f'{v["etot"]:.6f}', f'{v["etot"] * 1440:.1f}', f'{oc:+.2f}', f'{v["dep"]:.0f}', f'{v["edep"]:.0f}', f'{v["rms"]:.0f}', f'{v["beta"]:.2f}', f'{v["epb"] * 1440:.1f}', note])
# ---- folded light curves and model
E = {}
for s, v in ORDER(tcs):
    t, f, prod = load_epoch(s, v['tc']); E[s] = dict(t=t, f=f, cad=np.median(np.diff(t)), tc=v['tc'], prod=prod, used=1)
FITSEC = list(E)
if args.extra_fit:
    v = json.load(open(os.path.join(wd, args.extra_fit)))[str(args.extra_sector)]; t, f, prod = load_epoch(args.extra_sector, v['tc'])
    E[args.extra_sector] = dict(t=t, f=f, cad=np.median(np.diff(t)), tc=v['tc'], prod=prod, used=0)
if c['fix']:
    k, b, aR = c['fix']
else:
    sig = {s: 1.0 for s in FITSEC}
    res = lambda x: np.concatenate([resid_epoch(E, s, *x)[0] / sig[s] for s in FITSEC])
    x = least_squares(res, [c['k0'], c['b0'], c['aR0']], bounds=([0.01, 0, 5], [0.2, 1.0, 300]), x_scale=[0.003, 0.1, 5]).x
    for s in FITSEC:
        sig[s] = np.std(resid_epoch(E, s, *x)[0])
    k, b, aR = least_squares(res, x, bounds=([0.01, 0, 5], [0.2, 1.0, 300]), x_scale=[0.003, 0.1, 5]).x
print(f'shape: k={k:.4f} b={b:.2f} a/R*={aR:.1f}')
with open(os.path.join(args.out, f'folded_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic} = Gaia DR3 {c["gaia_dr3"]}: per-epoch light curves within +-{args.window:.1f} d of each fitted mid-transit time, {DATE}',
                  'flux = aperture flux divided by the per-epoch quadratic baseline that was solved jointly with the common transit model (scripts/fit_transits.py); not dilution corrected',
                  f'dt_days = time - tc(epoch) (TTV-corrected); phase = ((time - T0)/P + 0.5) mod 1 - 0.5 with T0 = BTJD {c["T0"]:.5f}, P = {PS} d',
                  f'common shape used for the model: k = Rp/R* = {k:.4f}, b = {b:.2f}, a/R* = {aR:.1f}, quadratic limb darkening u = {c["u"]} (TESS band)',
                  'flux_err = per-epoch residual rms; quality flag 0 cadences only; in_ephemeris = 0 marks the sector-66 event of TIC 402710280 (scattered light), which is not used in the ephemeris or the fold']))
    w = csv.writer(fo); w.writerow(['sector', 'product', 'epoch', 'in_ephemeris', 'time_btjd', 'dt_days', 'phase', 'flux', 'flux_err', 'model'])
    for s in sorted(E, key=lambda x: c['ncyc'][x]):
        r, base, mm = resid_epoch(E, s, k, b, aR); fd = E[s]['f'] / base; rms = np.std(r); ph = ((E[s]['t'] - c['T0']) / c['P'] + 0.5) % 1 - 0.5
        for ti, phi, fi, mi in zip(E[s]['t'], ph, fd, mm):
            w.writerow([sector_of(s), E[s]['prod'], c['ncyc'][s], E[s]['used'], f'{ti:.6f}', f'{ti - E[s]["tc"]:.6f}', f'{phi:.6f}', f'{fi:.6f}', f'{rms:.6f}', f'{mi:.6f}'])
tt = np.linspace(-0.5, 0.5, 2001)
with open(os.path.join(args.out, f'model_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic}: transit model used in the figures (batman, quadratic limb darkening u = {c["u"]}), k = {k:.4f}, b = {b:.2f}, a/R* = {aR:.1f}, P = {PS} d; no exposure-time smearing', 'dt_days = time from mid-transit']))
    w = csv.writer(fo); w.writerow(['dt_days', 'flux'])
    for a, bb in zip(tt, model(tt, 0.0, k, b, aR, 1e-4)):
        w.writerow([f'{a:.5f}', f'{bb:.7f}'])
# ---- sector coverage
with open(os.path.join(args.out, f'sectors_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic}: TESS sectors with full-frame-image cutouts used here, {DATE}',
                  'span = first/last cadence (BTJD) in the 21x21 TESScut; cadence_min = median cadence; predicted = mid-transit times of the linear ephemeris inside the span',
                  'covered = number of quality-0 cadences within +-T14/2 of each predicted time (0 = in a data gap)']))
    w = csv.writer(fo); w.writerow(['sector', 'camera', 'ccd', 'cadence_min', 'span_start_btjd', 'span_end_btjd', 'n_cadences_q0', 'predicted_transits_btjd', 'covered_cadences', 'transit_fitted'])
    for s in c['sectors']:
        z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); cc = np.load(os.path.join(wd, 'cut', f's{s}.npz'))
        t = z['tm']; q = z['q']; ok = np.isfinite(t); t = t[ok]; q = q[ok]; cad = np.median(np.diff(t)) * 1440
        n0 = int(np.floor((t.min() - c['T0']) / c['P'])) - 1; preds = []; cov = []
        for n in range(n0, n0 + 4):
            tp = c['T0'] + n * c['P']
            if t.min() - 0.2 < tp < t.max() + 0.2:
                preds.append(f'{tp:.3f}'); cov.append(int(((np.abs(t - tp) < c['T14h'] / 48) & (q == 0)).sum()))
        fitted = s in [sector_of(x) for x in tcs] or s == args.extra_sector
        w.writerow([s, int(cc['cam']), int(cc['ccd']), f'{cad:.1f}', f'{t.min():.3f}', f'{t.max():.3f}', int((q == 0).sum()), ';'.join(preds), ';'.join(map(str, cov)), 'yes' if fitted else 'no'])
# ---- difference images
srcrows = []
with open(os.path.join(args.out, f'diffimage_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic}: per-epoch difference images (out-of-transit model minus in-transit mean, e/s per cadence) from the 21x21 TESScut, quality-0 cadences (scripts/difference_images.py), {DATE}',
                  'sector = TESS sector; epoch = cycle number relative to T0 (transit_times_<tic>.csv)',
                  'x, y = cutout pixel indices (0-based); diff = difference image; oot = mean out-of-transit image; err = per-pixel error of the difference (empty where it was not computed)',
                  'pixel scale 21 arcsec/px; source positions in diffimage_sources_<tic>.csv']))
    w = csv.writer(fo); w.writerow(['sector', 'epoch', 'x', 'y', 'diff', 'oot', 'err'])
    for s, _ in ORDER(tcs):
        fn = os.path.join(wd, f'D_s{s}.npy')
        if not os.path.exists(fn):
            continue
        D = np.load(fn); Dd, Ee, Oo = D[0], D[1], D[2]; sec = sector_of(s); n = c['ncyc'][s]
        for y in range(21):
            for x in range(21):
                w.writerow([sec, n, x, y, f'{Dd[y, x]:.4f}', f'{Oo[y, x]:.3f}', f'{Ee[y, x]:.4f}'])
        cc = np.load(os.path.join(wd, 'cut', f's{sec}.npz')); wcs = WCS(Header.fromstring(str(cc['wcs']))); yr = 2016 + (tcs[s]['tc'] + 2457000 - 2457389) / 365.25
        for name, (ra, de, pa, pd) in c['sources'].items():
            dra = pa * (yr - 2016) / 3.6e6 / np.cos(np.radians(de)); dde = pd * (yr - 2016) / 3.6e6
            px, py = [float(v) for v in wcs.world_to_pixel(SkyCoord(ra + dra, de + dde, unit='deg'))]
            srcrows.append([sec, n, name, f'{ra:.7f}', f'{de:.7f}', f'{px:.3f}', f'{py:.3f}'])
with open(os.path.join(args.out, f'diffimage_sources_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic}: cutout pixel positions (from the TESScut WCS, proper motion applied where listed) of the target and the brighter Gaia DR3 sources near it, per epoch, {DATE}', 'ra, dec = Gaia DR3 epoch 2016.0 (deg); sector, epoch as in diffimage_<tic>.csv']))
    w = csv.writer(fo); w.writerow(['sector', 'epoch', 'source', 'ra_deg', 'dec_deg', 'x_px', 'y_px'])
    for r in srcrows:
        w.writerow(r)
fj = os.path.join(wd, 'diffimg.json')
if os.path.exists(fj):
    DJ = json.load(open(fj))
    with open(os.path.join(args.out, f'diffimage_centroids_{tic}.csv'), 'w') as fo:
        fo.write(hdr([f'TIC {tic}: per-epoch offset of the difference-image centroid from the out-of-transit centroid, with errors from fake transit windows (scripts/difference_images.py), {DATE}',
                      'dE, dN = offset in arcsec (East, North); dx, dy = offset in cutout pixels (21 arcsec/px); psf_sigma_px = Gaussian sigma fitted to the out-of-transit image']))
        w = csv.writer(fo); w.writerow(['sector', 'epoch', 'dE_arcsec', 'dN_arcsec', 'err_E_arcsec', 'err_N_arcsec', 'dx_px', 'dy_px', 'err_x_px', 'err_y_px', 'psf_sigma_px'])
        for s, _ in ORDER(tcs):
            v = DJ.get(str(s))
            if v:
                w.writerow([sector_of(s), c['ncyc'][s]] + [f'{v[x]:.2f}' for x in ('dE', 'dN', 'eE', 'eN')] + [f'{v[x]:.3f}' for x in ('dx_px', 'dy_px', 'ex_px', 'ey_px', 'psf_sigma_px')])
print('tables written to', args.out)
