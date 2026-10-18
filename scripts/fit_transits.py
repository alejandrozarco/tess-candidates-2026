"""Common-shape transit fit with per-epoch mid-transit times and quadratic baselines; timing errors; linear ephemeris.

usage: python scripts/fit_transits.py <tic> [--work ./work] [--aperture 2.0] [--epochs 39,93,100,102] [--window 1.0]
                                        [--order 2] [--fix k,b,aR] [--bfix b] [--arfix a/R*] [--sig-from <fit json>]
input:  <work>/<tic>/cut/lc_s<sector>.npz (aperture_photometry.py); for TIC 324429952 sector 13 the SPOC 2-min
        light curve <work>/<tic>/dl/mastDownload/TESS/*s0013-0000000324429952*/*lc.fits (SAP_FLUX, quality 0)
output: <work>/<tic>/fit_r<aperture>_<epochs>[_fix][_w<window>][_o<order>][_b<b>][_aR<a/R*>].json (window, order only if not 1.0, 2) with, per epoch (labelled by sector, or
        '<sector>a/b' for two transits in one sector): tc (BTJD), formal error ef, beta (time-averaging factor from
        20-90 min bins), prayer-bead error epb, etot = max(ef*beta, epb), per-epoch depth (k free, shape fixed) and
        residual rms; and <...>_shape.json with the common shape (k, b, a/R*, errors, rho*, T14, T23), chi2 and the
        linear ephemeris.
Model: batman (quadratic limb darkening), circular orbit, common k = Rp/R*, b, a/R* for all epochs, one quadratic
baseline per epoch multiplied into the model and solved linearly at every function evaluation, exposure-time
supersampling. Quality-0 cadences within +-window d of each mid-transit time. For TIC 324429952 the shape is fixed
(--fix, default from targets.py) at the solution constrained to the TIC stellar density; for TIC 402710280 and
TIC 89247735 it is free. --bfix / --arfix fix b or a/R* alone (chi2 profiles); --sig-from takes the per-epoch noise
from an earlier fit instead of re-estimating it, so that chi2 values of a profile are comparable.
"""
import argparse, glob, json, os
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as C2
from targets import TARGETS, workdir, add_common_args, import_batman, lab, sector_of

batman = import_batman()
p = add_common_args(argparse.ArgumentParser(description=__doc__))
p.add_argument('--aperture', default='2.0'); p.add_argument('--epochs', default=None); p.add_argument('--window', type=float, default=1.0)
p.add_argument('--order', type=int, default=2); p.add_argument('--fix', default='auto', help="'k,b,aR', 'none' or 'auto' (targets.py)")
p.add_argument('--allq', action='store_true', help='keep flagged cadences (used only for the S66 check of TIC 402710280)')
p.add_argument('--bfix', type=float, default=None); p.add_argument('--arfix', type=float, default=None); p.add_argument('--sig-from', default=None)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic)
USE = [lab(x) for x in args.epochs.split(',')] if args.epochs else c['epochs_fit']
FIX = c['fix'] if args.fix == 'auto' else (None if args.fix == 'none' else tuple(float(v) for v in args.fix.split(',')))
P = c['P']; tc0 = c['tc_guess']; ncyc = c['ncyc']; u = c['u']
D = {}
for s in USE:
    sec = sector_of(s)
    if sec == c.get('spoc_2min'):
        from astropy.io import fits
        fn = glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{sec:04d}-{args.tic:016d}*', '*lc.fits'))[0]
        d = fits.open(fn)[1].data; t = np.array(d['TIME'], float); f = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY'])
    else:
        z = np.load(os.path.join(wd, 'cut', f'lc_s{sec}.npz')); t = z['tm']; f = z[f'f{args.aperture}']; q = z['q']
    m = (np.abs(t - tc0[s]) < args.window) & np.isfinite(f)
    if not args.allq:
        m &= (q == 0)
    t, f = t[m], f[m]; f = f / np.median(f); cad = np.median(np.diff(t))
    pp = batman.TransitParams(); pp.t0 = tc0[s]; pp.per = P; pp.rp = c['k0']; pp.a = c['aR0']; pp.inc = 89.6; pp.ecc = 0; pp.w = 90; pp.u = u; pp.limb_dark = 'quadratic'
    D[s] = dict(t=t, f=f, cad=cad, p=pp, mod=batman.TransitModel(pp, t, supersample_factor=max(1, int(round(cad * 1440 / 2))), exp_time=cad))


def lc(s, tc, k, b, aR):
    d = D[s]; pp = d['p']; pp.t0 = tc; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(min(b / aR, 1))); return d['mod'].light_curve(pp)


def resid_epoch(s, tc, k, b, aR, ret=False):
    d = D[s]; m = lc(s, tc, k, b, aR); x = d['t'] - tc0[s]
    A = np.c_[[m * x ** j for j in range(args.order + 1)]].T; cc = np.linalg.lstsq(A, d['f'], rcond=None)[0]; r = d['f'] - A @ cc
    return (r, A @ cc) if ret else r


def res(pp):
    k, b, aR = pp[:3]; tcs = pp[3:]
    return np.concatenate([resid_epoch(s, tcs[i], k, b, aR) for i, s in enumerate(USE)])


p0 = [c['k0'], c['b0'], c['aR0']] + [tc0[s] for s in USE]
lo = [0.01, 0, 5] + [tc0[s] - 0.05 for s in USE]; hi = [0.2, 1.0, 200] + [tc0[s] + 0.05 for s in USE]
if FIX:
    p0[:3] = list(FIX); lo[:3] = [v - 1e-6 for v in FIX]; hi[:3] = [v + 1e-6 for v in FIX]
for j, v in ((1, args.bfix), (2, args.arfix)):
    if v is not None:
        p0[j] = v; lo[j] = v - 1e-6; hi[j] = v + 1e-6
fit = least_squares(res, p0, bounds=(lo, hi), x_scale=[0.003, 0.1, 5] + [0.002] * len(USE))
if args.sig_from:
    SF = json.load(open(os.path.join(wd, args.sig_from))); sig = {s: SF[str(s)]['rms'] * 1e-6 for s in USE}
else:
    sig = {s: np.std(resid_epoch(s, fit.x[3 + i], *fit.x[:3])) for i, s in enumerate(USE)}


def resw(pp):
    k, b, aR = pp[:3]; tcs = pp[3:]
    return np.concatenate([resid_epoch(s, tcs[i], k, b, aR) / sig[s] for i, s in enumerate(USE)])


fit = least_squares(resw, fit.x, bounds=(lo, hi), x_scale=[0.003, 0.1, 5] + [0.002] * len(USE))
J = fit.jac; cov = np.linalg.pinv(J.T @ J); e = np.sqrt(np.abs(np.diag(cov))); k, b, aR = fit.x[:3]
rho = 3 * np.pi * aR ** 3 / (6.674e-8 * (P * 86400) ** 2)
T14 = P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / aR / np.sin(np.arccos(b / aR))) * 24
T23 = P / np.pi * np.arcsin(np.sqrt(max((1 - k) ** 2 - b ** 2, 0)) / aR / np.sin(np.arccos(b / aR))) * 24
chi2 = float(np.sum(fit.fun ** 2))
print(f'aperture r={args.aperture} epochs {USE}: k={k:.4f}+-{e[0]:.4f} (k^2={k * k * 1e6:.0f} ppm) b={b:.2f}+-{e[1]:.2f} a/R*={aR:.1f}+-{e[2]:.1f} '
      f'rho*={rho:.2f} g/cc (circular) T14={T14:.2f} h T23={T23:.2f} h; chi2={chi2:.1f} for {len(fit.fun)} points; chi2_red={chi2 / (len(fit.fun) - len(fit.x)):.3f}')
out = {}
for i, s in enumerate(USE):
    tc = fit.x[3 + i]; ef = e[3 + i]; r, base = resid_epoch(s, tc, k, b, aR, ret=True)
    betas = []
    for nb in [int(round(x / (D[s]['cad'] * 1440))) for x in (20, 30, 45, 60, 90)]:
        if nb < 2:
            continue
        M = len(r) // nb; rb = r[:M * nb].reshape(M, nb).mean(1); betas.append(np.std(rb) / (np.std(r) / np.sqrt(nb) * np.sqrt(M / (M - 1))))
    beta = max(1, np.median(betas))
    tcs = []; n = len(r)
    for sh in np.linspace(0, n - 1, min(n, 150)).astype(int):      # prayer bead: cyclic residual shifts, refit tc only
        d0 = D[s]['f']; D[s]['f'] = base + np.roll(r, sh)
        tcs.append(least_squares(lambda x: resid_epoch(s, x[0], k, b, aR), [tc], bounds=([tc - 0.05], [tc + 0.05]), x_scale=[0.002]).x[0]); D[s]['f'] = d0
    epb = np.std(tcs)
    f2 = least_squares(lambda x: resid_epoch(s, tc, x[0], b, aR) / sig[s], [k], bounds=([0.0], [0.2])); ek = np.sqrt(np.linalg.inv(f2.jac.T @ f2.jac))[0, 0]
    out[s] = dict(tc=tc, ef=ef, beta=beta, eb=ef * beta, epb=epb, etot=max(ef * beta, epb), dep=f2.x[0] ** 2 * 1e6, edep=2 * f2.x[0] * ek * 1e6, rms=sig[s] * 1e6)
    print(f'  S{s}: tc={tc:.5f} formal {ef * 1440:.1f} min, beta={beta:.2f} -> {ef * beta * 1440:.1f} min, prayer-bead {epb * 1440:.1f} min; depth {out[s]["dep"]:.0f}+-{out[s]["edep"]:.0f} ppm; rms {sig[s] * 1e6:.0f} ppm')
n = np.array([ncyc[s] for s in USE]); T = np.array([out[s]['tc'] for s in USE]); E = np.array([out[s]['etot'] for s in USE])
A = np.c_[np.ones(len(n)), n]; W = 1 / E ** 2; C = np.linalg.inv(A.T @ (A * W[:, None])); bb = C @ (A.T @ (W * T)); oc = T - A @ bb; chi = np.sum(oc ** 2 * W); dof = len(n) - 2
print(f'  linear ephemeris: T0={bb[0]:.5f}+-{np.sqrt(C[0, 0]):.5f} P={bb[1]:.5f}+-{np.sqrt(C[1, 1]):.5f}; O-C (min) {np.round(oc * 1440, 1)} +- {np.round(E * 1440, 1)}; chi2={chi:.1f}/{dof} p={C2.sf(chi, dof):.3g}')
tag = (f'fit_r{args.aperture}_{"_".join(map(str, USE))}{"_fix" if FIX else ""}{"_allq" if args.allq else ""}'
       f'{f"_w{args.window:g}" if args.window != 1.0 else ""}{f"_o{args.order}" if args.order != 2 else ""}'
       f'{f"_b{args.bfix:g}" if args.bfix is not None else ""}{f"_aR{args.arfix:g}" if args.arfix is not None else ""}.json')
json.dump({str(s): v for s, v in out.items()}, open(os.path.join(wd, tag), 'w'), indent=1)
json.dump(dict(k=k, ek=e[0], b=b, eb=e[1], aR=aR, eaR=e[2], rho_circ=rho, T14_h=T14, T23_h=T23, chi2=chi2, npts=len(fit.fun), aperture=args.aperture,
               window=args.window, order=args.order, u=u, T0=bb[0], eT0=float(np.sqrt(C[0, 0])), P=bb[1], eP=float(np.sqrt(C[1, 1])), chi2_ephem=chi, dof_ephem=dof,
               oc_min={str(s): float(o * 1440) for s, o in zip(USE, oc)}), open(os.path.join(wd, tag[:-5] + '_shape.json'), 'w'), indent=1)
print('written', tag)
