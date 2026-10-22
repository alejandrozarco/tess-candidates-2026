"""J2: joint fit of JWST white-light curves of L 98-59 d (own extractions, j1). Shared a/R* (stellar-density prior
variants), shared k and LD per detector across visits, per-visit mid-time and impact parameter (or one common b).
Systematics per light curve: polynomial in time + detector regressors (trace x/y shift, width) [+ exponential ramp],
linear coefficients profiled; white noise with per-curve sigma x beta (time-averaging), or a celerite Matern-3/2 GP
(linear coefficients by generalised least squares). emcee posterior; reports T14 per visit and the difference.
Own code. Usage: j2_fit.py <config.json> <outbase>"""
import sys, os, json, numpy as np, batman, emcee, time
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
cfg = json.load(open(sys.argv[1])); outb = sys.argv[2]
PLN = cfg.get('planet', 'd'); P = PL[PLN][0]
LDMODE = cfg.get('ld', 'free'); LAW = cfg.get('law', 'quadratic'); AR = cfg.get('aR', 'free'); GEOM = cfg.get('geom', 'diff')
REGS = cfg.get('regs', ['xshift', 'yshift', 'width']); RAMP = cfg.get('ramp', 0); POLY = cfg.get('poly', 2); NOISE = cfg.get('noise', 'white')
BIN = cfg.get('bin', 5); TRIM = cfg.get('trim_min', 10.)
lcs = []
for item in cfg['lcs']:        # {"file":..., "visit":..., "det":...}
    fn = item['file'] if os.path.isabs(item['file']) else os.path.join(ROOT, item['file'])   # config paths are relative to ROOT
    z = np.load(fn); t, f = z['t'], z['f']
    keep = (t - t.min()) * 1440 > TRIM
    if 'tmax_rel' in item: keep &= (t - t.min()) < item['tmax_rel']
    if 'tmin_btjd' in item: keep &= t >= item['tmin_btjd']
    R = {r: z[r] for r in REGS}
    t, f = t[keep], f[keep]; R = {r: v[keep] for r, v in R.items()}
    # bin
    n = len(t) // BIN * BIN
    tb = t[:n].reshape(-1, BIN).mean(1); fb = f[:n].reshape(-1, BIN).mean(1); Rb = {r: v[:n].reshape(-1, BIN).mean(1) for r, v in R.items()}
    # robust outlier clip against a running median (5 sigma)
    from scipy.ndimage import median_filter
    rm = fb - median_filter(fb, 15, mode='nearest'); s = 1.4826 * np.median(np.abs(rm)); ok = np.abs(rm) < 5 * s
    tb, fb = tb[ok], fb[ok]; Rb = {r: v[ok] for r, v in Rb.items()}
    tau = (tb - tb.mean()) / (tb.max() - tb.min())
    cols = [tau ** j for j in range(POLY + 1)] + [(v - v.mean()) / (v.std() + 1e-12) for v in Rb.values()]
    if RAMP: cols.append(np.exp(-(tb - tb.min()) / (RAMP / 1440.)))
    X = np.vstack(cols).T
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = 0.045; pm.a = 34; pm.inc = 88.4; pm.ecc = 0; pm.w = 90
    pm.limb_dark = LAW; pm.u = [0.1, 0.2]
    tc0 = item.get('tc_guess', None)
    if tc0 is None:   # predicted mid-time from the linear ephemeris
        nn = np.round((tb.mean() - PL[PLN][1]) / P); tc0 = PL[PLN][1] + nn * P
    m = batman.TransitModel(pm, tb, supersample_factor=3, exp_time=np.median(np.diff(tb)) * 2 / 3)   # midpoint sub-exposures
    lcs.append(dict(t=tb, f=fb, X=X, visit=item['visit'], det=item['det'], tc0=tc0, pm=pm, m=m, sig=s / np.sqrt(1.), beta=1.))
visits = sorted(set(l['visit'] for l in lcs)); dets = sorted(set(l['det'] for l in lcs))
names = ['lnrho' if AR in ('rho', 'free') else None] + ['k_' + d for d in dets]
names = [n for n in names if n]
if LDMODE in ('free', 'prior'):
    for d in dets: names += ['l1_' + d, 'l2_' + d]
if GEOM == 'diff': names += ['b_' + v for v in visits]
else: names += ['b']
names += ['dtc_' + v for v in visits]
ix = {n: i for i, n in enumerate(names)}
LDFIX = cfg.get('ld_fixed', None)     # e.g. {"q1":0.117,"q2":0.0}
def aRof(th):
    if AR == 'fix': return cfg['aR_value']
    return aR_from_lnrho(th[ix['lnrho']], P)
def ldof(th, d):
    if LDMODE in ('free', 'prior'): a, b = th[ix['l1_' + d]], th[ix['l2_' + d]]
    else: a, b = LDFIX['q1'], LDFIX['q2']
    if LAW == 'quadratic' or LAW == 'squareroot':
        s = np.sqrt(a); return [2 * s * b, s * (1 - 2 * b)] if LAW == 'quadratic' else [s * (1 - 2 * b), 2 * s * b]   # Kipping q -> (u1,u2) / sqrt-law
    if LAW == 'power2': return [a, b * 3]   # c in (0,1), alpha in (0,3)
def bof(th, v): return th[ix['b_' + v]] if GEOM == 'diff' else th[ix['b']]
def lnprior(th):
    lp = 0.
    if AR == 'rho': lp += -0.5 * ((th[ix['lnrho']] - LNRHO[0]) / LNRHO[1]) ** 2
    if AR == 'free' and not (8.0 < th[ix['lnrho']] < 11.0): return -np.inf
    for d in dets:
        if not (0.01 < th[ix['k_' + d]] < 0.1): return -np.inf
        if LDMODE in ('free', 'prior'):
            if not (0 < th[ix['l1_' + d]] < 1 and 0 < th[ix['l2_' + d]] < 1): return -np.inf
    for v in visits:
        if not (0 <= bof(th, v) < 1.1): return -np.inf
        if abs(th[ix['dtc_' + v]]) > 30: return -np.inf
    return lp
def lc_model(th, l):
    pm = l['pm']; aR = aRof(th); b = bof(th, l['visit'])
    if b >= aR: return None
    pm.t0 = l['tc0'] + th[ix['dtc_' + l['visit']]] / 1440.; pm.rp = th[ix['k_' + l['det']]]; pm.a = aR
    pm.inc = np.degrees(np.arccos(b / aR)); pm.u = ldof(th, l['det'])
    return l['m'].light_curve(pm)
def gp_for(l, h):
    import celerite
    from celerite import terms
    k = terms.Matern32Term(log_sigma=h[0], log_rho=h[1]); gp = celerite.GP(k); gp.compute(l['t'], np.exp(h[2]) if len(h) > 2 else l['sig'] * l['beta']); return gp
def chi2_lc(m, l, ret=False):
    A = m[:, None] * l['X']
    if NOISE == 'gp':
        gp = l['gp']; KiA = gp.apply_inverse(A); Kif = gp.apply_inverse(l['f'])[:, 0]
        c = np.linalg.solve(A.T @ KiA, A.T @ Kif); r = l['f'] - A @ c
        ll = gp.log_likelihood(r)
        return (-2 * ll, r) if ret else -2 * ll
    w = 1 / (l['sig'] * l['beta']) ** 2
    c = np.linalg.lstsq(A * np.sqrt(w), l['f'] * np.sqrt(w), rcond=None)[0]; r = l['f'] - A @ c
    return (np.sum(r * r * w), r) if ret else np.sum(r * r * w)
def lnpost(th, retchi=False):
    lp = lnprior(th)
    if not np.isfinite(lp): return (-np.inf, np.inf) if retchi else -np.inf
    c2 = 0.
    for l in lcs:
        m = lc_model(th, l)
        if m is None: return (-np.inf, np.inf) if retchi else -np.inf
        c2 += chi2_lc(m, l)
    return (lp - 0.5 * c2, c2) if retchi else lp - 0.5 * c2
x0 = []
for n in names:
    x0.append(dict(lnrho=9.5).get(n, cfg.get('k_guess', 0.045) if n.startswith('k_') else 0.3 if n.startswith('l1') else 0.2 if n.startswith('l2') else cfg.get('b_guess', 0.93) if n.startswith('b') else 0.))
x0 = np.array(x0)
if LDMODE == 'fixed' and LDFIX is None: raise SystemExit('ld_fixed needed')
def mapfit(x):
    best = None
    for i in range(4):
        r = minimize(lambda th: -lnpost(th), x if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=30000, maxfev=30000, xatol=1e-8, fatol=1e-6, adaptive=True))
        r = minimize(lambda th: -lnpost(th), r.x, method='Powell', options=dict(maxiter=30000, xtol=1e-8, ftol=1e-10))
        if best is None or r.fun < best.fun - 1e-3: best = r
        else: break
    return best
# timing guess: scan mid-time offsets per visit with the starting shape
if NOISE == 'gp':
    for l in lcs: l['gp'] = gp_for(l, [np.log(0.5 * l['sig']), np.log(15 / 1440.)])
for v in visits:
    best = None
    for dt in np.arange(-20, 20.1, 1.0):
        x = x0.copy(); x[ix['dtc_' + v]] = dt; val = lnpost(x)
        if best is None or val > best[0]: best = (val, dt)
    x0[ix['dtc_' + v]] = best[1]
r = mapfit(x0)
# noise: per-curve sigma from residual rms, beta from time-averaging (5-30 min), then (GP) hyper-parameters
for l in lcs:
    m = lc_model(r.x, l); _, res = chi2_lc(m, l, ret=True); l['sig'] = np.std(res)
    dtm = np.median(np.diff(l['t'])) * 1440; bs = []
    for nb in (5, 10, 15, 20, 30):
        N = max(1, int(round(nb / dtm))); M = len(res) // N
        if M < 4: continue
        bm = res[:M * N].reshape(M, N).mean(1); bs.append(np.std(bm) / (l['sig'] / np.sqrt(N) * np.sqrt(M / (M - 1))))
    l['beta'] = max(1., float(np.mean(bs))) if NOISE == 'white' else 1.
    if NOISE == 'gp':
        def nll(h):
            if not (-14 < h[0] < -4 and np.log(2 / 1440.) < h[1] < np.log(0.2) and -14 < h[2] < -4): return 1e30
            l['gp'] = gp_for(l, h); return chi2_lc(m, l)
        rh = minimize(nll, [np.log(l['sig'] * 0.3), np.log(15 / 1440.), np.log(l['sig'] * 0.9)], method='Nelder-Mead', options=dict(maxiter=4000)); l['gp'] = gp_for(l, rh.x); l['gp_h'] = rh.x.tolist()
    print(l['visit'], l['det'], 'n %d sig %.0f ppm beta %.2f' % (len(l['t']), l['sig'] * 1e6, l['beta']), l.get('gp_h', ''), flush=True)
r = mapfit(r.x)
lpm, chim = lnpost(r.x, retchi=True)
print('MAP', dict(zip(names, np.round(r.x, 5))), 'chi2 %.1f' % chim, flush=True)
nd = len(names); nw = max(32, 2 * nd + 2); ns = cfg.get('nsteps', 4000)
rng = np.random.default_rng(5)
p0 = []
while len(p0) < nw:
    x = r.x + rng.standard_normal(nd) * np.array([0.01 if n == 'lnrho' else 2e-4 if n.startswith('k') else 0.02 if n.startswith('l') else 2e-3 if n.startswith('b') else 0.05 for n in names])
    if np.isfinite(lnpost(x)): p0.append(x)
sam = emcee.EnsembleSampler(nw, nd, lnpost, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
sam.run_mcmc(np.array(p0), ns, progress=False)
burn = ns // 3; ch = sam.get_chain()[burn:]; flat = ch.reshape(-1, nd)
tau = emcee.autocorr.integrated_time(ch, tol=0, quiet=True)
# split R-hat across two halves of the walkers
import arviz as az
half = nw // 2; rh = {n: float(az.rhat(np.array([ch[:, :half, i].ravel(), ch[:, half:2 * half, i].ravel()]))) for i, n in enumerate(names)}
def t14(P, aR, b, k):
    x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(aR ** 2 - b ** 2); return P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
th = flat[::max(1, len(flat) // 20000)]
aRs = np.array([aRof(x) for x in th])
kv = {v: np.mean([th[:, ix['k_' + d]] for d in sorted(set(l['det'] for l in lcs if l['visit'] == v))], 0) for v in visits}
T = {v: t14(P, aRs, np.array([bof(x, v) for x in th]), kv[v]) for v in visits}
B = {v: np.array([bof(x, v) for x in th]) for v in visits}
def q(x): return [float(np.percentile(x, 50)), float(np.std(x)), float(np.percentile(x, 16)), float(np.percentile(x, 84))]
out = dict(cfg=cfg, names=names, map=dict(zip(names, map(float, r.x))), chi2_map=float(chim), npts=int(sum(len(l['t']) for l in lcs)),
           noise={l['visit'] + '_' + l['det']: dict(sig=float(l['sig']), beta=float(l['beta']), gp=l.get('gp_h')) for l in lcs},
           post={n: q(flat[:, i]) for i, n in enumerate(names)}, tau_max=float(tau.max()), n_over_tau=float((ns - burn) / tau.max()), rhat=rh,
           T14={v: q(T[v]) for v in visits}, b={v: q(B[v]) for v in visits}, aR=q(aRs))
if len(visits) == 2 and GEOM == 'diff':
    v1, v2 = visits; out['dT14'] = q(T[v2] - T[v1]); out['db'] = q(B[v2] - B[v1])
if GEOM == 'diff' and len(visits) >= 2:
    # same noise model, all b_v tied to one value: MAP chi2 for the constant-geometry hypothesis
    bi = [ix['b_' + v] for v in visits]
    def full(xr):
        th = np.insert(np.delete(r.x.copy(), bi), 0, 0)[1:]   # placeholder
        x = r.x.copy(); j = 0; keep = [i for i in range(nd) if i not in bi]
        for i in keep: x[i] = xr[j]; j += 1
        for i in bi: x[i] = xr[-1]
        return x
    keep = [i for i in range(nd) if i not in bi]
    xr0 = np.r_[r.x[keep], np.mean(r.x[bi])]
    best = None
    for it in range(4):
        rr = minimize(lambda xr: -lnpost(full(xr)), xr0 if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=30000, maxfev=30000, xatol=1e-8, fatol=1e-6, adaptive=True))
        rr = minimize(lambda xr: -lnpost(full(xr)), rr.x, method='Powell', options=dict(maxiter=30000, xtol=1e-8, ftol=1e-10))
        if best is None or rr.fun < best.fun - 1e-3: best = rr
        else: break
    lpc, chic = lnpost(full(best.x), retchi=True)
    out['const_same_noise'] = dict(chi2=float(chic), dchi2=float(chic - chim), ddof=len(visits) - 1)
    print('constant-geometry MAP with the same noise model: chi2 %.2f, dchi2 %.2f for %d dof' % (chic, chic - chim, len(visits) - 1))
json.dump(out, open(outb + '.json', 'w'), indent=1)
np.savez_compressed(outb + '_samples.npz', samples=flat[::max(1, len(flat) // 30000)], names=np.array(names))
print(json.dumps({k: out[k] for k in ('T14', 'b', 'aR', 'tau_max', 'n_over_tau')}), flush=True)
if 'dT14' in out: print('dT14 %.3f +- %.3f min, db %.4f +- %.4f' % (out['dT14'][0], out['dT14'][1], out['db'][0], out['db'][1]))
