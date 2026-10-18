"""Per-epoch difference images from the TESScut cutouts: centroid offsets and a joint PSF fit on the listed sources.

usage: python scripts/difference_images.py <tic> --fit <fit json> [--work ./work]
input:  <work>/<tic>/cut/s<sector>.npz (tesscut_download.py), mid-transit times from <work>/<tic>/<fit json>
output: <work>/<tic>/D_s<epoch>.npy = [difference image, per-pixel error, out-of-transit mean image] (21x21 each),
        <work>/<tic>/diffimg.json with the per-epoch centroid offsets (pixels, and arcsec E/N) and their errors, the
        Gaussian PSF sigma of the out-of-transit image and, when requested, the Gaia PSF-position test.
Difference image = (out-of-transit model of every pixel, polynomial in time of order oot_order (default 2) fitted to
|dt| > oot_min (default 0.15 d) within 0.7 d) minus the in-transit mean (|dt| < hw). A free 2-D Gaussian is fitted
to the difference image and to the out-of-transit image; their centroid difference is the offset. Its error is the
scatter of the same measurement over 80 fake windows at random out-of-transit times, at least fake_excl (default
0.8 d) from any predicted transit and with at least fake_min_frac of the in-transit cadence count (the real
difference image is added so that the signal-to-noise is the same). A joint linear fit of Gaussian PSFs placed at the
catalogue positions of the sources in targets.py (sigma from the out-of-transit fit) gives the share of the dimming
per source; single-source fits give delta chi2 between sources.
Per-target settings (`diff` in targets.py; TIC 89247735): hw, oot_min, oot_order, the Gaussian-fit radius and sigma
start/minimum, the random seed and fake-window settings, and gaia_radius/gaia_gmax for the Gaia PSF-position test:
for every Gaia DR3 source with G < gaia_gmax within gaia_radius arcsec of the TIC position (positions propagated with
proper motion to the epoch), the difference image within 6 px of the target is fitted with amplitude x PSF(source
position) + constant, and chi2 (noise = scatter of the difference image beyond 6 px) is compared with the fit at the
target position; the implied eclipse depth of each source is the fitted amplitude over the out-of-transit peak of the
target scaled by the Gaia G-band flux ratio. The Gaia query (all sources, no magnitude cut) is cached in
<work>/<tic>/gaia_sources.json and is also used by dilution_correction.py.
"""
import argparse, json, os
import numpy as np
from scipy.optimize import least_squares
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord
from targets import TARGETS, workdir, add_common_args, lab, sector_of

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--fit', required=True); p.add_argument('--hw', type=float, default=None, help='in-transit half-window (d); default T14/2 - 0.03')
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); DS = c.get('diff', {})
J = json.load(open(os.path.join(wd, args.fit))); hw = args.hw or DS.get('hw') or max(c['T14h'] / 48 - 0.03, 0.05)
OOTMIN = DS.get('oot_min', 0.15); OORD = DS.get('oot_order', 2); CRAD = DS.get('cen_radius', 3.5); S0 = DS.get('sigma0', 0.9); SMIN = DS.get('sigma_min', 0.4)
N = 21; yy, xx = np.indices((N, N)); src = c['sources']


def g2(p):
    a, x0, y0, s, cst = p; return a * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s * s)) + cst


def cen(img, x0, y0, rad=CRAD):
    w = np.hypot(xx - x0, yy - y0) < rad
    return least_squares(lambda p: (g2(p) - img)[w], [max(img[int(round(y0)), int(round(x0))], 1), x0, y0, S0, 0],
                         bounds=([-np.inf, x0 - 3, y0 - 3, SMIN, -np.inf], [np.inf, x0 + 3, y0 + 3, 3, np.inf])).x


def dimg(tm, f, tc):
    m = np.abs(tm - tc) < 0.7; t = tm[m] - tc; F = f[m]; win = np.abs(t) < hw; oot = np.abs(t) > OOTMIN
    A = np.c_[[t[oot] ** j for j in range(OORD + 1)]].T; cc = np.linalg.lstsq(A, F[oot].reshape(oot.sum(), -1), rcond=None)[0]
    pred = (np.c_[[t[win] ** j for j in range(OORD + 1)]].T @ cc).reshape(-1, N, N); r = F[oot] - (A @ cc).reshape(-1, N, N)
    return (pred - F[win]).mean(0), F[oot].mean(0), r.std(0) / np.sqrt(win.sum())


def near_transit(t, dmin=0.8):
    n = np.round((t - c['T0']) / c['P']); return abs(t - c['T0'] - n * c['P']) < dmin


GAIA = []
if DS.get('gaia_radius'):
    gf = os.path.join(wd, 'gaia_sources.json')
    if not os.path.exists(gf):
        from astroquery.gaia import Gaia
        r0 = DS['gaia_radius'] / 3600; q0 = (f"select source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, DISTANCE(POINT({c['tic_ra']},{c['tic_dec']}),POINT(ra,dec))*3600 as sep "
                                             f"from gaiadr3.gaia_source where 1=CONTAINS(POINT({c['tic_ra']},{c['tic_dec']}),CIRCLE(ra,dec,{r0})) order by sep")
        tb = Gaia.launch_job(q0).get_results()
        rows = [{k: (None if np.ma.is_masked(r[k]) else (str(r[k]) if k == 'source_id' else float(r[k]))) for k in tb.colnames} for r in tb]
        json.dump(rows, open(gf, 'w'), indent=1)
    GAIA = [r for r in json.load(open(gf)) if r['phot_g_mean_mag'] is not None and r['phot_g_mean_mag'] < DS.get('gaia_gmax', 19)]
    tg = [r for r in GAIA if str(r['source_id']) == str(c['gaia_dr3'])][0]
    print(f'Gaia PSF-position test: {len(GAIA)} Gaia DR3 sources with G < {DS.get("gaia_gmax", 19)} within {DS["gaia_radius"]:.0f}"')

rng = np.random.default_rng(DS.get('seed', 2)); OUT = {}; STACK = {}; FEX = DS.get('fake_excl', 0.8); FMIN = DS.get('fake_min_frac', 0.0)
for s_, v in sorted(J.items(), key=lambda kv: c['ncyc'][lab(kv[0])]):
    s = lab(s_); sec = sector_of(s); tc = v['tc']; z = np.load(os.path.join(wd, 'cut', f's{sec}.npz')); tm = z['tm']; f = z['f'].astype(float); q = z['q']
    ok = np.isfinite(tm) & np.all(np.isfinite(f), axis=(1, 2)) & (q == 0); tm, f = tm[ok], f[ok]
    w = WCS(Header.fromstring(str(z['wcs']))); yr = 2016 + (tc + 2457000 - 2457389) / 365.25; pos = {}
    for k, (ra, de, pa, pd) in src.items():
        dra = pa * (yr - 2016) / 3.6e6 / np.cos(np.radians(de)); dde = pd * (yr - 2016) / 3.6e6
        pos[k] = [float(x) for x in w.world_to_pixel(SkyCoord(ra + dra, de + dde, unit='deg'))]
    tx, ty = pos['target']; D, O, E = dimg(tm, f, tc); pD = cen(D, tx, ty); pO = cen(O, tx, ty)
    offs = []
    for _ in range(80):
        tf = rng.uniform(tm.min() + 0.8, tm.max() - 0.8)
        if abs(tf - tc) < FEX or near_transit(tf, FEX) or (np.abs(tm - tf) < hw).sum() < max(3, FMIN * (np.abs(tm - tc) < hw).sum()):
            continue
        try:
            Df, _, _ = dimg(tm, f, tf); offs.append(cen(Df + D, tx, ty)[1:3])
        except Exception:
            pass
    offs = np.array(offs); sx, sy = offs[:, 0].std(), offs[:, 1].std()
    sky = lambda x, y: (lambda cc: (cc.ra.deg, cc.dec.deg))(w.pixel_to_world(x, y))
    rD, dD = sky(pD[1], pD[2]); rO, dO = sky(pO[1], pO[2]); rT, dT = sky(tx, ty)
    dE = (rD - rO) * np.cos(np.radians(dT)) * 3600; dN = (dD - dO) * 3600; Jc = np.abs(np.linalg.det(w.pixel_scale_matrix)) ** 0.5 * 3600
    print(f'S{s} tc={tc:.3f}: difference-image centroid minus out-of-transit centroid: dE={dE:+.1f}" dN={dN:+.1f}" (+-{sx * Jc:.1f}, {sy * Jc:.1f}" from {len(offs)} fake windows); OOT centroid - target = ({(rO - rT) * np.cos(np.radians(dT)) * 3600:+.1f},{(dO - dT) * 3600:+.1f})"')
    print(f'      in pixels: dx={pD[1] - pO[1]:+.3f} dy={pD[2] - pO[2]:+.3f} px (+-{sx:.3f}, {sy:.3f}); OOT PSF sigma {pO[3]:.2f} px; {int((np.abs(tm - tc) < hw).sum())} in-transit cadences')
    sg = pO[3]; Ev = np.where(np.isfinite(E) & (E > 0), E, np.nanmedian(E))
    psf = lambda x0, y0: (lambda G: (G / G.sum()).ravel())(np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sg * sg)))
    M = np.array([psf(*pos[k]) for k in src] + [np.ones(N * N)]).T; Mw = M / Ev.ravel()[:, None]
    b = np.linalg.lstsq(Mw, (D / Ev).ravel(), rcond=None)[0]; e = np.sqrt(np.diag(np.linalg.inv(Mw.T @ Mw))); bo = np.linalg.lstsq(M, O.ravel(), rcond=None)[0]
    for i, k in enumerate(src):
        print(f'      joint PSF fit: {k:52s} dflux {b[i]:+8.1f} +- {e[i]:5.1f} ({b[i] / e[i]:+.1f} sigma); OOT flux {bo[i]:9.0f}; implied depth {b[i] / bo[i] * 100 if bo[i] > 0 else np.nan:7.3f} %')
    ch = {}
    for i, k in enumerate(src):
        Ms = np.c_[M[:, i], M[:, -1]]; bs = np.linalg.lstsq(Ms / Ev.ravel()[:, None], (D / Ev).ravel(), rcond=None)[0]; ch[k] = (((D.ravel() - Ms @ bs) / Ev.ravel()) ** 2).sum()
    c0 = min(ch.values()); print('      single-source delta chi2:', {k: round(v - c0, 1) for k, v in ch.items()})
    OUT[str(s)] = dict(sector=sec, epoch=c['ncyc'][s], tc=tc, dE=dE, dN=dN, eE=sx * Jc, eN=sy * Jc, dx_px=pD[1] - pO[1], dy_px=pD[2] - pO[2], ex_px=sx, ey_px=sy,
                       psf_sigma_px=pO[3], tgt_px=[tx, ty], n_fake=len(offs)); np.save(os.path.join(wd, f'D_s{s}.npy'), np.array([D, E, O]))
    if GAIA:
        def gpx(r):
            ra = r['ra'] + (r['pmra'] or 0) / 3.6e6 * (yr - 2016) / np.cos(np.radians(r['dec'])); de = r['dec'] + (r['pmdec'] or 0) / 3.6e6 * (yr - 2016)
            return [float(x) for x in w.world_to_pixel(SkyCoord(ra, de, unit='deg'))]
        gx, gy = gpx(tg); mk = np.hypot(xx - gx, yy - gy) < 6; noise = np.std((D - np.median(D))[~mk]); peak = O[int(round(gy)), int(round(gx))] - np.median(O); rows = []
        for r in GAIA:
            x0, y0 = gpx(r); G = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sg * sg)); A_ = np.c_[G[mk], np.ones(mk.sum())]; cf = np.linalg.lstsq(A_, D[mk], rcond=None)[0]
            fr = 10 ** (-0.4 * (r['phot_g_mean_mag'] - tg['phot_g_mean_mag']))
            rows.append(dict(source_id=str(r['source_id']), G=r['phot_g_mean_mag'], sep=r['sep'], chi2=float(np.sum((D[mk] - A_ @ cf) ** 2) / noise ** 2), implied_depth=float(cf[0] / (peak * fr))))
        c0 = [r for r in rows if r['source_id'] == str(c['gaia_dr3'])][0]['chi2']
        for r in rows:
            r['dchi2'] = r['chi2'] - c0; STACK.setdefault(r['source_id'], [r['G'], r['sep'], 0.0]); STACK[r['source_id']][2] += r['dchi2']
        rows.sort(key=lambda r: r['chi2']); OUT[str(s)]['gaia_psf_test'] = rows
        print('      Gaia PSF-position test (lowest chi2 first):', '; '.join(f'{r["source_id"]} G={r["G"]:.2f} {r["sep"]:.1f}" dchi2={r["dchi2"]:+.0f} depth={r["implied_depth"] * 100:.1f}%' for r in rows[:4]))
if STACK:
    print('Gaia PSF-position test, delta chi2 against the target summed over epochs:', '; '.join(f'{k} G={g:.2f} {sp:.1f}" {d:+.0f}' for k, (g, sp, d) in sorted(STACK.items(), key=lambda kv: kv[1][2])[:8]))
Ea = np.array([[v['dE'], v['dN'], v['eE'], v['eN']] for v in OUT.values()]); wE = 1 / Ea[:, 2] ** 2; wN = 1 / Ea[:, 3] ** 2
print(f'weighted mean centroid offset: E {np.sum(Ea[:, 0] * wE) / wE.sum():+.2f} +- {wE.sum() ** -0.5:.2f}"  N {np.sum(Ea[:, 1] * wN) / wN.sum():+.2f} +- {wN.sum() ** -0.5:.2f}"; chi2 vs 0: {np.sum(Ea[:, 0] ** 2 * wE) + np.sum(Ea[:, 1] ** 2 * wN):.1f}/{2 * len(Ea)}')
if STACK:
    OUT['gaia_psf_test_stacked'] = {k: dict(G=g, sep=sp, dchi2_sum=d) for k, (g, sp, d) in STACK.items()}
json.dump(OUT, open(os.path.join(wd, 'diffimg.json'), 'w'), indent=1)
