"""C3: per-wavelength-bin transit fits of one JWST visit (chromatic duration test). Own code.
Light curve per bin from the per-column extraction (c1_extract_cols.py) or a SOSS bin file (c4_soss_bins.py); same
preparation as ../joint/j2_fit.py (first 10 min trimmed, 5-integration bins, 5-sigma running-median clip,
systematics = quadratic in time + detector regressors, linear coefficients profiled, white noise sigma x beta).
Fixed per visit (shared by all bins): period, mid-time and a/R* (from the white-light joint fit). Free per bin:
k, b, quadratic u1, u2. LD modes: 'prior' (Gaussian, sigma = LDSIG, centred on the Kurucz values of c2_ld.py) or
'free' (uniform, -0.5 < u1, u2 < 1, 0 < u1 + u2 < 1). Optional 'ramp_min': exponential ramp exp(-(t - t_start)/tau)
from the exposure start as an extra profiled regressor.
Outputs per bin (JSON lines, resumable): MAP, posterior (emcee) of k, b, u1, u2 and T14 (first to fourth contact),
noise sigma and beta, and a profile chi2(b) on a grid (k, u1, u2 re-optimised at each b) for the shared-b
(achromatic) model. Usage: c3_binfit.py <config.json>"""
import sys, json, os, numpy as np, batman, emcee
from scipy.optimize import minimize
from scipy.ndimage import median_filter
cfg = json.load(open(sys.argv[1]))
P, TC, AR = cfg['P'], cfg['tc_btjd'], cfg['aR']
LDMODE = cfg.get('ld', 'prior'); LDSIG = cfg.get('ld_sigma', 0.1); REGS = cfg.get('regs', ['xshift', 'yshift', 'width'])
BIN = cfg.get('bin', 5); TRIM = cfg.get('trim_min', 10.); POLY = cfg.get('poly', 2); NS = cfg.get('nsteps', 3000)
BGRID = np.arange(*cfg.get('bgrid', [0.890, 0.9705, 0.0015]))
ld = json.load(open(cfg['ld_file']))['bins'][cfg['ld_key']]
outp = cfg['out']; done = set()
if os.path.exists(outp):
    for s in open(outp): done.add(json.loads(s)['bin'])
z = np.load(cfg['file'])
t_all = z['t']
keep = (t_all - t_all.min()) * 1440 > TRIM
if 'tmax_rel' in cfg: keep &= (t_all - t_all.min()) < cfg['tmax_rel']
if 'tmin_btjd' in cfg: keep &= t_all >= cfg['tmin_btjd']
def binned(v): v = v[keep]; n = len(v) // BIN * BIN; return v[:n].reshape(-1, BIN).mean(1)
tb = binned(t_all); Rb = {r: binned(z[r]) for r in REGS}
tau = (tb - tb.mean()) / (tb.max() - tb.min())
pm = batman.TransitParams(); pm.t0 = TC; pm.per = P; pm.rp = 0.045; pm.a = AR; pm.inc = 88.4; pm.ecc = 0; pm.w = 90; pm.limb_dark = 'quadratic'; pm.u = [0.1, 0.2]
def t14(b, k): x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(AR ** 2 - b ** 2); return P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
def bin_flux(lo, hi):
    if 'Fcol' in z.files:
        sel = (z['lam'] >= lo) & (z['lam'] < hi); f = z['Fcol'][:, sel].sum(1).astype(float)
        return f / np.median(f), int(sel.sum())
    key = 'f_%.3f_%.3f' % (lo, hi); return z[key] / np.median(z[key]), int(z['n_' + key[2:]])
for ib, row in enumerate(ld):
    lo, hi = row['lo'], row['hi']; name = '%.3f-%.3f' % (lo, hi)
    if name in done: continue
    f, ncol = bin_flux(lo, hi); fb = binned(f)
    ok = np.isfinite(fb)
    rm = fb - median_filter(np.where(ok, fb, np.nanmedian(fb)), 15, mode='nearest'); s = 1.4826 * np.nanmedian(np.abs(rm)); ok &= np.abs(rm) < 5 * s
    t, y = tb[ok], fb[ok]
    cols_ = [tau[ok] ** j for j in range(POLY + 1)] + [(v[ok] - v[ok].mean()) / (v[ok].std() + 1e-12) for v in Rb.values()]
    if cfg.get('ramp_min'): cols_.append(np.exp(-(tb[ok] - t_all.min()) * 1440. / cfg['ramp_min']))   # exponential ramp from the exposure start
    X = np.vstack(cols_).T
    m = batman.TransitModel(pm, t, supersample_factor=3, exp_time=np.median(np.diff(t)))
    U0 = np.array([row['u1'], row['u2']]); state = dict(sig=s, beta=1.)
    def model(k, b, u1, u2):
        pm.rp = k; pm.inc = np.degrees(np.arccos(b / AR)); pm.u = [u1, u2]; return m.light_curve(pm)
    def chi2(th):
        k, b, u1, u2 = th
        if not (0.01 < k < 0.1 and 0 <= b < 1.1): return 1e30
        if not (-0.5 < u1 < 1 and -0.5 < u2 < 1 and 0 < u1 + u2 < 1): return 1e30
        mm = model(k, b, u1, u2); A = mm[:, None] * X
        c = np.linalg.lstsq(A, y, rcond=None)[0]; r = y - A @ c
        c2 = np.sum(r * r) / (state['sig'] * state['beta']) ** 2
        if LDMODE == 'prior': c2 += np.sum(((np.array([u1, u2]) - U0) / LDSIG) ** 2)
        return c2
    def resid(th):
        mm = model(*th); A = mm[:, None] * X; c = np.linalg.lstsq(A, y, rcond=None)[0]; return y - A @ c, mm, A @ c
    def fit(x0):
        best = None
        for it in range(3):
            r = minimize(chi2, x0 if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=8000, maxfev=8000, xatol=1e-7, fatol=1e-5, adaptive=True))
            if best is None or r.fun < best.fun - 1e-3: best = r
            else: break
        return best
    x0 = np.array([cfg.get('k_guess', 0.045), cfg.get('b_guess', 0.93), U0[0], U0[1]])
    best = None
    for bb in cfg.get('b_starts', [0.91, 0.93, 0.95]):
        x0[1] = bb; r = fit(x0)
        if best is None or r.fun < best.fun: best = r
    # noise: sigma from residuals, beta from time-averaging over 5-30 min
    res, _, _ = resid(best.x); state['sig'] = float(np.std(res))
    dtm = np.median(np.diff(t)) * 1440; bs = []
    for nb in (5, 10, 15, 20, 30):
        N = max(1, int(round(nb / dtm))); M = len(res) // N
        if M < 4: continue
        bm = res[:M * N].reshape(M, N).mean(1); bs.append(np.std(bm) / (state['sig'] / np.sqrt(N) * np.sqrt(M / (M - 1))))
    state['beta'] = max(1., float(np.mean(bs)))
    best = fit(best.x); c2min = float(best.fun)
    # profile chi2(b): k, u1, u2 re-optimised at each b (warm start along the grid, both directions)
    prof = np.full(len(BGRID), np.nan); j0 = int(np.argmin(np.abs(BGRID - best.x[1])))
    for order in (range(j0, len(BGRID)), range(j0, -1, -1)):
        xs = best.x[[0, 2, 3]].copy()
        for j in order:
            bj = BGRID[j]
            g = lambda q: chi2([q[0], bj, q[1], q[2]])
            r = minimize(g, xs, method='Nelder-Mead', options=dict(maxiter=3000, maxfev=3000, xatol=1e-7, fatol=1e-5, adaptive=True))
            r = minimize(g, r.x, method='Nelder-Mead', options=dict(maxiter=3000, maxfev=3000, xatol=1e-7, fatol=1e-5, adaptive=True))
            prof[j] = min(prof[j], r.fun) if np.isfinite(prof[j]) else r.fun; xs = r.x
    # posterior
    nw = 20; rng = np.random.default_rng(11 + ib)
    p0 = []
    while len(p0) < nw:
        x = best.x + rng.standard_normal(4) * np.array([3e-4, 3e-3, 0.02, 0.02])
        if chi2(x) < 1e29: p0.append(x)
    sam = emcee.EnsembleSampler(nw, 4, lambda th: -0.5 * chi2(th) if chi2(th) < 1e29 else -np.inf, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    sam.run_mcmc(np.array(p0), NS, progress=False)
    burn = NS // 3; ch = sam.get_chain()[burn:]; fl = ch.reshape(-1, 4)
    try: taus = emcee.autocorr.integrated_time(ch, tol=0, quiet=True); tmax = float(np.max(taus))
    except Exception: tmax = float('nan')
    T = t14(fl[:, 1], fl[:, 0])
    def q(x): return [float(np.median(x)), float(np.std(x)), float(np.percentile(x, 16)), float(np.percentile(x, 84))]
    resb, mmb, sysb = resid(best.x)
    out = dict(bin=name, lo=lo, hi=hi, ncol=ncol, npts=int(len(t)), sig_ppm=state['sig'] * 1e6, beta=state['beta'], ld_mode=LDMODE, ld_model=U0.tolist(),
               map=dict(k=best.x[0], b=best.x[1], u1=best.x[2], u2=best.x[3], T14=float(t14(best.x[1], best.x[0]))), chi2_map=c2min,
               post=dict(k=q(fl[:, 0]), b=q(fl[:, 1]), u1=q(fl[:, 2]), u2=q(fl[:, 3]), T14=q(T)), corr_b_u=[float(np.corrcoef(fl[:, 1], fl[:, 2])[0, 1]), float(np.corrcoef(fl[:, 1], fl[:, 3])[0, 1])],
               tau_max=tmax, n_over_tau=float((NS - burn) / tmax) if np.isfinite(tmax) else None, bgrid=BGRID.tolist(), prof_chi2=prof.tolist(),
               t=t.tolist(), resid=resb.tolist(), model=mmb.tolist())
    with open(outp, 'a') as fo: fo.write(json.dumps(out) + '\n')
    print('%s %s ncol %d sig %.0f ppm beta %.2f  T14 %.2f +- %.2f  b %.4f +- %.4f  k %.4f  u %.2f %.2f  N/tau %.0f' % (cfg.get('name', ''), name, ncol, state['sig'] * 1e6, state['beta'], out['post']['T14'][0], out['post']['T14'][1], out['post']['b'][0], out['post']['b'][1], out['post']['k'][0], out['post']['u1'][0], out['post']['u2'][0], out['n_over_tau'] or -1), flush=True)
