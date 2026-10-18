"""Companion-mass limits at the transit period from the Gaia DR3 radial-velocity scatter, and reference amplitudes.

usage: python scripts/gaia_companion_limits.py <tic> [--shape <fit>_shape.json] [--work ./work]
input:  Gaia DR3 gaia_source row of the target (queried once, cached in <work>/<tic>/gaia_target.json), the period of
        targets.py, the TIC v8.2 stellar mass and distance (`star` in targets.py), optionally a/R* from a fit_transits.py
        shape file for the ellipsoidal amplitudes
output: printed limits and <work>/<tic>/gaia_limits.json
Radial velocity: rv_amplitude_robust (peak-to-peak of the epoch radial velocities after outlier rejection) is taken
as an upper bound on 2K, so K < rv_amplitude_robust / 2 and, for a circular orbit with M2 << M*,
M2 sin i < K / K_1MJ with K_1MJ = 28.4329 m/s (M*/Msun)^(-2/3) (P/yr)^(-1/3). Astrometry: semi-major axis of the
star's orbit about the centre of mass, a1 = a M2 / (M* + M2), divided by the distance (companion light neglected).
Reference amplitudes: ellipsoidal ~ 1.5 (M2/M*) (a/R*)^-3, Doppler beaming ~ 4 K / c, with K from the full two-body
expression.
"""
import argparse, json, os
import numpy as np
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--shape', default=None)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); st = c['star']
gf = os.path.join(wd, 'gaia_target.json')
if not os.path.exists(gf):
    from astroquery.gaia import Gaia
    cols = ('source_id, parallax, ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, astrometric_excess_noise, astrometric_excess_noise_sig, non_single_star, '
            'radial_velocity, radial_velocity_error, rv_nb_transits, rv_nb_deblended_transits, rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof')
    tb = Gaia.launch_job(f"select {cols} from gaiadr3.gaia_source where source_id = {c['gaia_dr3']}").get_results()
    json.dump({k: (None if np.ma.is_masked(tb[k][0]) else (str(tb[k][0]) if k == 'source_id' else float(tb[k][0]))) for k in tb.colnames}, open(gf, 'w'), indent=1)
g = json.load(open(gf)); P = c['P']; Pyr = P / 365.25; Ms = st['M']; MJ = 9.546e-4
print('Gaia DR3', g['source_id'], {k: v for k, v in g.items() if k != 'source_id'})
Kper = 28.4329 * Ms ** (-2 / 3) * Pyr ** (-1 / 3); Kmax = g['rv_amplitude_robust'] / 2 * 1e3
out = dict(K_per_MJ_ms=Kper, K_max_kms=Kmax / 1e3, M2sini_max_MJ=Kmax / Kper, rv_amplitude_robust_kms=g['rv_amplitude_robust'],
           rv_nb_transits=g['rv_nb_transits'], rv_chisq_pvalue=g['rv_chisq_pvalue'])
print(f'K for 1 MJ (M2 << M*, circular, sin i = 1): {Kper:.1f} m/s; K < rv_amplitude_robust/2 = {Kmax / 1e3:.2f} km/s -> M2 sin i < {Kmax / Kper:.1f} MJ '
      f'(M* = {Ms} Msun, P = {P} d; {g["rv_nb_transits"]:.0f} Gaia RV transits, rv_chisq_pvalue {g["rv_chisq_pvalue"]:.2f})')
aR = json.load(open(os.path.join(wd, args.shape)))['aR'] if args.shape else None
out['reference'] = {}
for M2 in (30 * MJ, 0.03, 0.08, 0.3):
    a_au = ((Ms + M2) * Pyr ** 2) ** (1 / 3); a1_mas = a_au * M2 / (Ms + M2) / st['d_pc'] * 1e3
    K = 28.4329 * (M2 / MJ) * (Ms + M2) ** (-2 / 3) * Pyr ** (-1 / 3); r = dict(photocentre_semi_amplitude_mas=a1_mas, K_kms=K / 1e3, beaming_ppm=4 * K / 2.998e8 * 1e6)
    if aR:
        r['ellipsoidal_ppm'] = 1.5 * M2 / Ms / aR ** 3 * 1e6
    out['reference'][f'{M2:.4f}'] = r
    print(f'M2 = {M2:.4f} Msun ({M2 / MJ:.0f} MJ): astrometric semi-amplitude {a1_mas:.3f} mas, K = {K / 1e3:.2f} km/s, beaming {r["beaming_ppm"]:.0f} ppm'
          + (f', ellipsoidal {r["ellipsoidal_ppm"]:.1f} ppm (a/R* = {aR:.1f})' if aR else ''))
json.dump(out, open(os.path.join(wd, 'gaia_limits.json'), 'w'), indent=1)
