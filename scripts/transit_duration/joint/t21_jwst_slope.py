"""T21: JWST-only duration slope from the joint JWST posteriors (shared a/R*): per posterior sample, weighted
least-squares slope of T14 against epoch over the visits; also the TESS-predicted T14 change between the epochs.
Usage: t21_jwst_slope.py <fitname> [<fitname> ...]"""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
EP = {}
for v, f in (('1201', 'wl_1201_soss'), ('1224', 'wl_1224_nrs1'), ('4098', 'wl_4098_nrs1'), ('4098o1', 'wl_4098o1_soss'), ('c2512a', None), ('c2512b', None)):
    if f: EP[v] = float(year(np.median(np.load(ROOT + '/out/jwst/%s.npz' % f)['t'])))
tess = np.load(ROOT + '/out/mcmc/main_lin_b_samples.npz')['dT14dt_min_per_yr']
out = {}
for fn in sys.argv[1:]:
    d = json.load(open(ROOT + '/out/jwst/fits/%s.json' % fn)); z = np.load(ROOT + '/out/jwst/fits/%s_samples.npz' % fn)
    names = [str(x) for x in z['names']]; S = z['samples']; ix = {n: i for i, n in enumerate(names)}
    vis = [v for v in d['T14'] if v in EP and (len(sys.argv) < 0 or True)]
    use = [v for v in vis if d['T14'][v][1] < 1.0]          # drop visits with T14 error > 1 min (4098o1 x1dints)
    P = PL['d'][0]
    aR = aR_from_lnrho(S[:, ix['lnrho']], P) if 'lnrho' in ix else np.full(len(S), d['cfg'].get('aR_value', np.nan))
    dets = {v: sorted(set(l['det'] for l in d['cfg']['lcs'] if l['visit'] == v)) for v in use}
    T = []
    for v in use:
        k = np.mean([S[:, ix['k_' + dd]] for dd in dets[v]], 0); b = S[:, ix['b_' + v]]
        x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(aR ** 2 - b ** 2); T.append(P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440)
    T = np.array(T).T; y = np.array([EP[v] for v in use]); w = 1 / T.std(0) ** 2
    A = np.vstack([np.ones_like(y), y - y.mean()]).T
    sl = np.array([np.linalg.lstsq(A * np.sqrt(w)[:, None], t * np.sqrt(w), rcond=None)[0][1] for t in T])
    pred = np.median(tess)
    out[fn] = dict(visits=use, epochs=y.tolist(), T14=[[float(np.median(T[:, j])), float(T[:, j].std())] for j in range(len(use))],
                   slope=[float(np.median(sl)), float(sl.std())], z=float(np.median(sl) / sl.std()), tess_slope=float(pred),
                   tension_sigma=float((np.median(sl) - pred) / np.hypot(sl.std(), tess.std())))
    print('%-20s visits %s  JWST slope %.2f +- %.2f min/yr (z %.1f); TESS %.2f; JWST-TESS %.1f sigma' % (fn, ','.join(use), np.median(sl), sl.std(), np.median(sl) / sl.std(), pred, out[fn]['tension_sigma']))
json.dump(out, open(ROOT + '/out/jwst_slopes.json', 'w'), indent=1)
