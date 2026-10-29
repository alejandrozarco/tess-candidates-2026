"""TOI-421 b (TIC 94986319; P 5.19759807 d, T0 BTJD 1441.283952): early (sectors 5, 6, 32) and late (98, 107) epoch
folds of the SPOC 2-min PDCSAP light curves, and per-sector fits. c transits masked (+-0.12 d). One tc offset per epoch,
quadratic baseline per window.
Models: (1) batman k 0.0279, a/R* 13.76, u 0.401/0.222, b free; (2) box depth in +-30 min; (3) trapezoid (T14, T12, depth).
usage: python scripts/toi421/planet_b_toi421.py   -> work/toi421/out/planet_b.json"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import fit_toi421 as adv
from scipy.optimize import least_squares
PB, T0B = 5.19759807, 1441.283952; HW = 0.15
def wins(tags):
    out = []
    for tag in tags:
        t, f, e, exp = adv.load(tag); sec = int(tag[1:4])
        for n in range(int((t.min() - T0B) / PB) - 1, int((t.max() - T0B) / PB) + 2):
            tc = T0B + n * PB; m = np.abs(t - tc) < HW
            nc = np.round((t - adv.T0) / adv.P); m &= np.abs(t - (adv.T0 + nc * adv.P)) > 0.12
            if m.sum() < 0.7 * 2 * HW / exp or (np.abs(t[m] - tc) < 0.022).sum() < 0.7 * 0.044 / exp: continue
            out.append(dict(tc0=tc, t=t[m], f=f[m], e=e[m], exp=exp, sec=sec, year=adv.year(tc)))
    return out
def fit_b(ws, kind='batman'):
    def mdl(t, tc, x, exp):
        if kind == 'batman':
            p = adv.batman.TransitParams(); p.t0 = tc; p.per = PB; p.rp = 0.0279; p.a = 13.76
            p.inc = np.degrees(np.arccos(min(x[0] / 13.76, 1))); p.ecc = 0; p.w = 90; p.u = [0.401, 0.222]; p.limb_dark = 'quadratic'
            return adv.batman.TransitModel(p, t, supersample_factor=5, exp_time=exp).light_curve(p)
        if kind == 'box':
            return 1 - x[0] * (np.abs(t - tc) < 30 / 1440.)
        return adv.trap(t, tc, x[0], x[1], x[2], exp)
    g0 = {'batman': [0.95], 'box': [5e-4], 'trap': [5e-4, 0.044, 0.008]}[kind]; ng = len(g0) + 1
    x0 = g0 + [0.0] + [1, 0, 0] * len(ws)
    def res(x):
        r = []
        for i, w in enumerate(ws):
            c = x[ng + 3 * i: ng + 3 * i + 3]; dt = w['t'] - w['tc0']
            r.append((w['f'] - mdl(w['t'], w['tc0'] + x[ng - 1], x[:ng - 1], w['exp']) * (c[0] + c[1] * dt + c[2] * dt * dt)) / w['e'])
        return np.concatenate(r)
    lo = {'batman': [0.5], 'box': [-0.003], 'trap': [0, 0.01, 0.001]}[kind] + [-0.02] + [0.9, -1, -10] * len(ws)
    hi = {'batman': [1.1], 'box': [0.003], 'trap': [0.003, 0.08, 0.03]}[kind] + [0.02] + [1.1, 1, 10] * len(ws)
    s = least_squares(res, x0, bounds=(lo, hi), x_scale='jac'); chi2r = np.sum(s.fun ** 2) / (s.fun.size - len(x0))
    C = np.linalg.pinv(s.jac.T @ s.jac) * max(1, chi2r); e = np.sqrt(np.abs(np.diag(C)))
    return dict(x=[float(v) for v in s.x[:ng]], ex=[float(v) for v in e[:ng]], chi2r=float(chi2r), ntr=len(ws),
                dt_min=float(s.x[ng - 1] * 1440), year=float(np.mean([w['year'] for w in ws])))
if __name__ == '__main__':
    E = wins(['s005_2m', 's006_2m', 's032_2m']); L = wins(['s098_2m', 's107_2m']); out = {}
    for kind in ('batman', 'box', 'trap'):
        out[kind] = dict(early=fit_b(E, kind), late=fit_b(L, kind)); print(kind, json.dumps(out[kind]))
    # per sector batman
    out['per_sector'] = {}
    for tag in ('s005_2m', 's006_2m', 's032_2m', 's098_2m', 's107_2m'):
        r = fit_b(wins([tag]), 'batman'); out['per_sector'][tag] = r; print(tag, r['ntr'], round(r['x'][0], 4), round(r['ex'][0], 4), round(r['dt_min'], 1))
    out['_meta'] = __doc__
    json.dump(out, open(os.path.join(adv.WORK, 'out', 'planet_b.json'), 'w'), indent=1)
