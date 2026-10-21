"""Fit JWST white-light transits (own code). Usage: j2_fit.py <out.json> <P> <mode> <name1> [name2 ...]
mode: single = each visit free (k,aR,b,LD,tc,ramp); joint = common k,aR,LD; per-visit b, tc, ramp.
Input: $TD_WORK/jw_<name>.npz (j1_prep.py); the first letter of <name> selects the planet (b, c, d).
Optional env AR_FIX=value fixes aR; NWALK, NSTEP set the emcee walkers and steps (defaults 40, 4000; 40, 3000 were used)."""
import sys, os, json, numpy as np, emcee
sys.path.insert(0, __file__.rsplit('/',1)[0])
from lib import ld_q2u, t14, t23, year, WORK
import batman
out, P, mode = sys.argv[1], float(sys.argv[2]), sys.argv[3]; names = sys.argv[4:]
arfix = float(os.environ['AR_FIX']) if 'AR_FIX' in os.environ else None
nwalk = int(os.environ.get('NWALK', 40)); nstep = int(os.environ.get('NSTEP', 4000))
V = []
for n in names:
    z = np.load(os.path.join(WORK, 'jw_%s.npz' % n)); t, f = z['t'], z['f']
    p = batman.TransitParams(); p.t0 = 0; p.per = P; p.rp = .045; p.a = 34; p.inc = 88.4; p.ecc = 0; p.w = 90; p.limb_dark = 'quadratic'; p.u = [.1, .2]
    V.append(dict(name=n, t=t, f=f, t0=t[0], p=p, m=None, s=np.std(np.diff(f)) / np.sqrt(2)))
eph = {'b':(2.2531140,1366.17056),'c':(3.6906764,1367.27303),'d':(7.450729,1362.74002)}
for v in V:
    Pp, T0 = eph[v['name'][0]]; E = np.round((v['t'].mean() - T0) / Pp); v['tc0'] = T0 + E * Pp
    v['m'] = batman.TransitModel(v['p'], v['t'] - v['tc0'])
def lc(v, dtc, k, aR, b, u1, u2):
    p = v['p']; p.t0 = dtc; p.rp = k; p.a = aR; p.inc = np.degrees(np.arccos(b / aR)); p.u = [u1, u2]
    return v['m'].light_curve(p)
def vis_chi2(v, dtc, k, aR, b, u1, u2, A, ltau, ret=False):
    m = lc(v, dtc, k, aR, b, u1, u2)
    x = (v['t'] - v['t'][0]); ramp = 1 + A * np.exp(-x / np.exp(ltau))
    X = np.vstack([np.ones_like(x), x, x * x]).T * (m * ramp)[:, None]
    c, *_ = np.linalg.lstsq(X, v['f'], rcond=None)
    r = v['f'] - X @ c
    if ret: return r, X @ c, m
    return np.sum(r * r) / v['s'] ** 2
nv = len(V)
# parameter layout
if mode == 'single':
    assert nv == 1
    lab = ['dtc', 'k', 'aR', 'b', 'q1', 'q2', 'A', 'ltau']
    def lnp(th):
        dtc, k, aR, b, q1, q2, A, ltau = th
        if arfix: aR = arfix
        if not (abs(dtc) < 0.02 and 0.01 < k < 0.2 and 5 < aR < 100 and 0 <= b < 1 + k and 0 < q1 < 1 and 0 < q2 < 1 and abs(A) < 0.01 and np.log(0.001) < ltau < np.log(0.1)): return -np.inf
        u1, u2 = ld_q2u(q1, q2)
        return -0.5 * vis_chi2(V[0], dtc, k, aR, b, u1, u2, A, ltau)
    x0 = [0, .046, 34., .92, .2, .3, 0., np.log(0.01)]
else:
    lab = ['k', 'aR', 'q1', 'q2'] + sum([['dtc%d' % i, 'b%d' % i, 'A%d' % i, 'ltau%d' % i] for i in range(nv)], [])
    def lnp(th):
        k, aR, q1, q2 = th[:4]
        if arfix: aR = arfix
        if not (0.01 < k < 0.2 and 5 < aR < 100 and 0 < q1 < 1 and 0 < q2 < 1): return -np.inf
        u1, u2 = ld_q2u(q1, q2); c = 0
        for i, v in enumerate(V):
            dtc, b, A, ltau = th[4 + 4 * i: 8 + 4 * i]
            if not (abs(dtc) < 0.02 and 0 <= b < 1 + k and abs(A) < 0.01 and np.log(0.001) < ltau < np.log(0.1)): return -np.inf
            c += vis_chi2(v, dtc, k, aR, b, u1, u2, A, ltau)
        return -0.5 * c
    x0 = [.046, 34., .2, .3] + sum([[0, .92, 0., np.log(0.01)] for _ in range(nv)], [])
from scipy.optimize import minimize
best = None
for b0 in ([0.90, 0.93, 0.95] if V[0]['name'][0] == 'd' else [0.3, 0.5]):
    for aR0 in ([30., 34., 38.] if V[0]['name'][0] == 'd' else ([16., 19., 22.] if V[0]['name'][0] == 'b' else [24., 26., 28.])):
        x = list(x0)
        if mode == 'single': x[2] = aR0; x[3] = b0
        else:
            x[1] = aR0
            for i in range(nv): x[5 + 4 * i] = b0
        r = minimize(lambda th: -lnp(th), x, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-6, adaptive=True))
        r = minimize(lambda th: -lnp(th), r.x, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-6, adaptive=True))
        if best is None or r.fun < best.fun: best = r
# rescale per-visit noise to chi2_red=1 and red-noise beta
for i, v in enumerate(V):
    th = best.x
    if mode == 'single': dtc, k, aR, b, q1, q2, A, ltau = th
    else:
        k, aR, q1, q2 = th[:4]; dtc, b, A, ltau = th[4 + 4 * i: 8 + 4 * i]
    if arfix: aR = arfix
    r, mod, m = vis_chi2(v, dtc, k, aR, b, *ld_q2u(q1, q2), A, ltau, ret=True)
    s_white = np.std(r)
    betas = []
    for nb in (5, 10, 20, 40, 80):
        nn = len(r) // nb
        rb = r[:nn * nb].reshape(nn, nb).mean(1)
        betas.append(np.std(rb) / (s_white / np.sqrt(nb)) * np.sqrt(nn / (nn - 1)))
    v['beta'] = max(1.0, float(np.median(betas))); v['s'] = s_white * v['beta']
    v['s_white'] = s_white
ndim = len(best.x)
rng = np.random.default_rng(1)
scale = np.maximum(np.abs(best.x) * 1e-3, 1e-5)
p0 = []
while len(p0) < nwalk:
    p = best.x + scale * rng.standard_normal(ndim)
    if np.isfinite(lnp(p)): p0.append(p)
p0 = np.array(p0)
sam = emcee.EnsembleSampler(nwalk, ndim, lnp)
sam.run_mcmc(p0, nstep, progress=False, skip_initial_state_check=True)
ch = sam.get_chain(discard=nstep // 2, thin=10, flat=True)
res = dict(names=names, mode=mode, arfix=arfix, labels=lab, best=best.x.tolist(), chi2=float(best.fun * 2),
           beta=[v['beta'] for v in V], s_white_ppm=[v['s_white'] * 1e6 for v in V], acc=float(np.mean(sam.acceptance_fraction)))
def q(x): return [float(np.percentile(x, 50)), float(np.percentile(x, 84) - np.percentile(x, 50)), float(np.percentile(x, 50) - np.percentile(x, 16))]
for i, v in enumerate(V):
    if mode == 'single':
        dtc, k, aR, b, q1, q2 = ch[:, 0], ch[:, 1], ch[:, 2], ch[:, 3], ch[:, 4], ch[:, 5]
    else:
        k, aR, q1, q2 = ch[:, 0], ch[:, 1], ch[:, 2], ch[:, 3]; dtc, b = ch[:, 4 + 4 * i], ch[:, 5 + 4 * i]
    if arfix: aR = np.full(len(k), arfix)
    T14 = t14(P, aR, b, k) * 1440; T23 = t23(P, aR, b, k) * 1440
    u1, u2 = ld_q2u(q1, q2)
    res[v['name']] = dict(year=float(year(v['tc0'])), tc=q(v['tc0'] + dtc), T14_min=q(T14), T23_min=q(T23), b=q(b), k=q(k), aR=q(aR), u1=q(u1), u2=q(u2), depth_ppm=q(k * k * 1e6))
    print(v['name'], 'yr %.3f' % res[v['name']]['year'], 'T14 %.2f +%.2f -%.2f' % tuple(res[v['name']]['T14_min']), 'b %.4f +%.4f -%.4f' % tuple(res[v['name']]['b']),
          'k %.4f' % res[v['name']]['k'][0], 'aR %.2f +%.2f -%.2f' % tuple(res[v['name']]['aR']), 'T23 %.1f' % res[v['name']]['T23_min'][0], 'beta %.2f' % v['beta'], 'tc %.5f' % res[v['name']]['tc'][0])
if mode != 'single' and nv > 1:
    for i in range(1, nv):
        dT = (t14(P, ch[:, 1] if not arfix else arfix, ch[:, 5 + 4 * i], ch[:, 0]) - t14(P, ch[:, 1] if not arfix else arfix, ch[:, 5], ch[:, 0])) * 1440
        db = ch[:, 5 + 4 * i] - ch[:, 5]
        res['dT14_%s_minus_%s' % (names[i], names[0])] = q(dT); res['db_%s_minus_%s' % (names[i], names[0])] = q(db)
        print('dT14 %s - %s = %.2f +%.2f -%.2f min ; db = %.4f +%.4f -%.4f' % (names[i], names[0], *q(dT), *q(db)))
json.dump(res, open(out, 'w'), indent=1)
