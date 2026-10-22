"""TOI-2638.01 (TIC 458342769): calibration of the two shape tests by simulation (README section 6.9).

usage: python scripts/toi2638/sims_toi2638.py [--work ./work] [--out data] [--nsim 500] [--ntrap 500] [--seed 1]
input:  the light curves of download_toi2638.py; data/ephemeris_458342769.csv, data/jointfit_458342769.csv and
        data/shape_458342769.csv written by fit_toi2638.py (read from --out)
output: <out>/simulations_458342769.csv (format in the file header)
Test 1 (joint models): the observed delta chi2 = chi2(M5) - chi2(M2) of the variant 'SPOC PDCSAP + QLP sector 10'
  (M5: one shape whose depth scales linearly in time, as a change of dilution would give; M2: impact parameter linear in
  time; same number of parameters). Simulated light curves = best-fit M5 or best-fit M2 model + the observed residuals
  of the M2 fit of each sector shifted cyclically by a random amount (keeps correlated noise); both models are refitted
  (two starting points each) and delta chi2 recorded.
Test 2 (trapezoid duration): the slope of the per-sector trapezoid T14 (SPOC PDCSAP, weights from the bootstrap errors
  of shape_458342769.csv). Simulations as in test 1 for the six SPOC sectors; trapezoids refitted per sector.
Requires numpy, scipy, astropy, batman-package.
"""
import argparse, csv, os, sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load, windows, Model, trapezoid, wlinfit, hdr, REF_YEAR, DATE

ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--work', default='./work'); ap.add_argument('--out', default='data')
ap.add_argument('--nsim', type=int, default=500); ap.add_argument('--ntrap', type=int, default=500); ap.add_argument('--seed', type=int, default=1)
args = ap.parse_args(); rng = np.random.default_rng(args.seed)


def rows(fn):
    return [r for r in csv.reader(l for l in open(os.path.join(args.out, fn)) if not l.startswith('#'))]


eph = [r for r in rows('ephemeris_458342769.csv')[1:] if r[0] == 'ephemeris'][0]; T0, P = float(eph[10]), float(eph[12])
J = {}
for r in rows('jointfit_458342769.csv')[1:]:
    if r[0] == 'SPOC PDCSAP + QLP sector 10':
        J.setdefault(r[1], dict(chi2=float(r[4]), x=[]))['x'].append(float(r[6]))
sh = {int(r[0]): r for r in rows('shape_458342769.csv')[1:]}
pdc = load(args.work, 'SPOC'); d = dict(pdc); d[10] = load(args.work, 'QLP')[10]
W = windows(d, P, T0); secs = sorted(W); yr = {s: W[s]['year'] - REF_YEAR for s in secs}
M = {s: Model(W[s]['t'], W[s]['cad'], P, T0) for s in secs}
m2 = lambda p: {s: 1 - M[s].depth(p[0], p[2] + p[3] * yr[s], p[1]) for s in secs}
m5 = lambda p: {s: 1 - (1 + p[3] * yr[s]) * M[s].depth(*p[:3]) for s in secs}
p2, p5 = J['M2']['x'], J['M5']['x']; mod2, mod5 = m2(p2), m5(p5)
resid = {s: W[s]['f'] - mod2[s] for s in secs}
obs1 = J['M5']['chi2'] - J['M2']['chi2']


def fitboth(F):
    f2 = lambda p: np.concatenate([(F[s] - v) / W[s]['e'] for s, v in m2(p).items()])
    f5 = lambda p: np.concatenate([(F[s] - v) / W[s]['e'] for s, v in m5(p).items()])
    c2 = min(2 * least_squares(f2, x0, bounds=([0.01, 2, 0, -0.5], [0.5, 40, 1.5, 0.5]), x_scale='jac', max_nfev=3000).cost for x0 in (p2, [0.055, 11, 0.92, -0.015]))
    c5 = min(2 * least_squares(f5, x0, bounds=([0.01, 0, 2, -0.3], [0.5, 1.5, 40, 0.3]), x_scale='jac', max_nfev=3000).cost for x0 in (p5, [0.065, 0.97, 9.5, 0.15]))
    return c5 - c2


T1 = {'M5': [], 'M2': []}
for i in range(args.nsim):
    for name, mod in (('M5', mod5), ('M2', mod2)):
        F = {s: mod[s] + np.roll(resid[s], rng.integers(len(resid[s]))) for s in secs}
        T1[name].append(float(fitboth(F)))
    if i % 10 == 0:
        print('test 1', i, T1['M5'][-1], T1['M2'][-1], flush=True)

# test 2: trapezoid T14 slope on the SPOC sectors
ts = [s for s in secs if s in sh and sh[s][5] != '']
eT = np.array([float(sh[s][6]) for s in ts]); T14obs = np.array([float(sh[s][5]) for s in ts]); x_ = np.array([yr[s] for s in ts])
obs2 = wlinfit(x_, T14obs, eT)[0][1]


def t14fit(s, y):
    ph = (W[s]['t'] - T0 + 0.5 * P) % P - 0.5 * P; e = W[s]['e']; best = None
    for x0 in ([0.003, 0.058, 0.9, 0], [0.003, 0.075, 0.6, 0]):
        r = least_squares(lambda p: (y - trapezoid(p, ph, W[s]['cad'])) / e, x0, bounds=([0, 0.02, 0.05, -0.01], [0.02, 0.12, 1.0, 0.01]))
        if best is None or r.cost < best.cost:
            best = r
    return best.x[1] * 1440


T2 = {'M5': [], 'M2': []}; noiseless = {n: [t14fit(s, 1 - m[s]) for s in ts] for n, m in (('M5', mod5), ('M2', mod2))}
for i in range(args.ntrap):
    for name, mod in (('M5', mod5), ('M2', mod2)):
        T = np.array([t14fit(s, 1 - (mod[s] + np.roll(resid[s], rng.integers(len(resid[s]))))) for s in ts])
        T2[name].append(float(wlinfit(x_, T, eT)[0][1]))
    if i % 20 == 0:
        print('test 2', i, T2['M5'][-1], T2['M2'][-1], flush=True)

with open(os.path.join(args.out, 'simulations_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: calibration of the shape tests by simulation; scripts/toi2638/sims_toi2638.py, {DATE} (seed {args.seed})',
            'truth M5 = one transit shape whose depth scales linearly in time (what a change of dilution gives); truth M2 = impact parameter linear in time; both at the best fits of '
            'jointfit_458342769.csv (variant SPOC PDCSAP + QLP sector 10); noise = observed M2 residuals of each sector shifted cyclically by a random amount',
            'test 1: delta chi2 = chi2(M5) - chi2(M2) after refitting both models; test 2: weighted slope of the per-sector trapezoid T14 (SPOC sectors, bootstrap errors of shape_458342769.csv)',
            'n_ge_observed = number of simulations with a value >= the observed one; z_equiv = (observed - mean) / std; noiseless_T14_min = trapezoid T14 per SPOC sector fitted to the noise-free truth',
            'test,truth,n_sim,observed,mean,std,median,p16,p84,maximum,n_ge_observed,z_equiv,noiseless_T14_min'])
    for test, D, obs in (('delta chi2 M5 - M2', T1, obs1), ('trapezoid T14 slope (min/yr)', T2, obs2)):
        for truth in ('M5', 'M2'):
            a = np.array(D[truth]); nl = ' '.join(f'{v:.1f}' for v in noiseless[truth]) if test.startswith('trap') else ''
            f.write(f'{test},{truth},{len(a)},{obs:.2f},{a.mean():.2f},{a.std():.2f},{np.median(a):.2f},{np.percentile(a, 16):.2f},{np.percentile(a, 84):.2f},'
                    f'{a.max():.2f},{int(np.sum(a >= obs))},{(obs - a.mean()) / a.std():.2f},{nl}\n')
            print(test, truth, f'mean {a.mean():.2f} std {a.std():.2f} n>=obs {np.sum(a >= obs)}/{len(a)}')
