"""T6: emcee posterior for one model variant of the joint TESS inference. Several independent ensembles (for R-hat),
HDF5 checkpoint/resume, autocorrelation times, MAP polish, BIC, importance-sampling evidence. Own code.
Usage: t6_mcmc.py <prep.npz> <variant> <nens> <nwalk> <nsteps> <outbase>"""
import sys, os, json, numpy as np, time, emcee
from scipy.optimize import minimize
from scipy import stats
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *; from t4_core import *
prep, var, nens, nwalk, nsteps, outb = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
pr = np.load(prep, allow_pickle=True)
segs = list(np.load(str(pr['win']), allow_pickle=True))
D = TData(segs, nb=int(pr['nb'])); D.build(list(pr['keep'])); D.set_sig(pr['sig'])
LD = json.load(open(ROOT + '/out/claret_ld.json')); ld0 = (LD['u1'], LD['u2'], LD['sigma_prior'])
VARS = {
 'const': dict(geom='const'), 'lin_b': dict(geom='lin_b'), 'lin_i': dict(geom='lin_i'),
 'camp': dict(geom='camp'), 'camp_slope': dict(geom='camp_slope'),
 'lin_b_kcamp': dict(geom='lin_b', kmode='camp'), 'const_kcamp': dict(geom='const', kmode='camp'),
 'lin_b_ldfree': dict(geom='lin_b', ld='free'), 'const_ldfree': dict(geom='const', ld='free'),
 'lin_b_ecc': dict(geom='lin_b', ecc='lit'), 'const_ecc': dict(geom='const', ecc='lit'),
 'lin_b_fixtc': dict(geom='lin_b', timing='fixed'), 'const_fixtc': dict(geom='const', timing='fixed'),
 'lin_b_rhofree': dict(geom='lin_b', rho='wide'), 'const_rhofree': dict(geom='const', rho='wide'),
}
kw = VARS[var]
S = Spec(D=D, ld0=ld0, tc_fixed=pr['tc_pertr'], **kw)
# starting point from the lin_b MAP
m0 = dict(zip([str(x) for x in pr['names']], pr['map_lin_b']))
aR0 = aR_from_lnrho(m0['lnrho'], P_D)
def start():
    x = []
    for n in S.names:
        if n in m0: x.append(m0[n])
        elif n.startswith('k'): x.append(m0['k'])
        elif n.startswith('bc'):
            c = int(n[2]); x.append(m0['b0'] + m0['bdot'] * (S.camp_mid[c] - TREF))
        elif n == 'i0': x.append(np.degrees(np.arccos(m0['b0'] / aR0)))
        elif n == 'idot': x.append(-np.degrees(m0['bdot'] / aR0))
        elif n == 'secw': x.append(0.05)
        elif n == 'sesw': x.append(-0.03)
        else: raise ValueError(n)
    x = np.array(x)
    if 'bdot' in S.ix and S.geom == 'camp_slope': x[S.ix['bdot']] = m0['bdot']
    return x
x0 = start()
negl = lambda th: -S.lnpost(th)
def polish(x):
    best = None
    for i in range(4):
        r = minimize(negl, x if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-5, adaptive=True))
        r = minimize(negl, r.x, method='Powell', options=dict(maxiter=40000, xtol=1e-7, ftol=1e-10))
        if best is None or r.fun < best.fun - 1e-4: best = r
        else: break
    return best
t0 = time.time()
r = polish(x0); xmap = r.x
print(var, 'initial MAP lnpost %.2f (%.0f s)' % (-r.fun, time.time() - t0), flush=True)
# ensembles
ndim = len(S.names); nwalk = max(nwalk, 2 * ndim + 2); scale = np.array([1e-4 if n.startswith('k') else 0.01 if n == 'lnrho' else 0.01 if n.startswith('q') else 1e-3 if n.startswith('b') else 1e-4 if n == 'bdot' else 0.05 if n in ('i0',) else 1e-3 if n == 'idot' else 0.01 if n.startswith('se') else 0.1 for n in S.names])
samples_all, lnp_all, taus, accs = [], [], [], []
for e in range(nens):
    h5 = '%s_ens%d.h5' % (outb, e)
    be = emcee.backends.HDFBackend(h5)
    rng = np.random.default_rng(1000 + e)
    done = be.iteration if os.path.exists(h5) else 0
    if done == 0:
        be.reset(nwalk, ndim)
        p0 = []
        while len(p0) < nwalk:
            x = xmap + scale * rng.standard_normal(ndim) * (3 if e % 2 else 1)
            if np.isfinite(S.lnpost(x)): p0.append(x)
        p0 = np.array(p0)
    else: p0 = None
    sam = emcee.EnsembleSampler(nwalk, ndim, S.lnpost, backend=be, moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    left = nsteps - done
    if left > 0:
        sam.run_mcmc(p0, left, progress=False, store=True)
    ch = be.get_chain(); lp = be.get_log_prob()
    burn = nsteps // 3
    tau = emcee.autocorr.integrated_time(ch[burn:], tol=0, quiet=True)
    taus.append(tau); accs.append(float(np.mean(be.accepted / be.iteration)))
    samples_all.append(ch[burn:]); lnp_all.append(lp[burn:])
    print(var, 'ens', e, 'steps', be.iteration, 'acc %.2f' % accs[-1], 'tau max %.0f (N/tau = %.0f)' % (tau.max(), (nsteps - burn) / tau.max()), '%.0f s' % (time.time() - t0), flush=True)
# R-hat across ensembles (rank-normalised split R-hat, Vehtari et al. 2021) on per-ensemble step-major sequences
import arviz as az
arr = np.array([s.reshape(-1, ndim) for s in samples_all])           # (nens, draws, ndim)
rhat = {n: float(az.rhat(arr[:, :, i])) for i, n in enumerate(S.names)} if nens > 1 else {}
flat = arr.reshape(-1, ndim); lpf = np.array(lnp_all).reshape(-1)
# MAP polish from best sample
jb = np.argmax(lpf); r2 = polish(flat[jb])
if -r2.fun > -r.fun: xmap = r2.x
lpmap, chimap = S.lnpost(xmap, retchi=True)
nfree = ndim; npts = D.npts
bic = chimap + nfree * np.log(npts)
# derived quantities (thinned)
th = flat[::max(1, len(flat) // 20000)]
yrs = np.array([np.mean(D.yr[D.camp == c]) for c in range(4)])
der = {}
if S.geom in ('lin_b', 'const', 'lin_i'):
    T = np.array([S.t14_at(x, np.r_[yrs, TREF - 0.5, TREF + 0.5]) for x in th])
    der['T14_camp_min'] = T[:, :4]; der['dT14dt_min_per_yr'] = T[:, 5] - T[:, 4]
    if S.geom == 'lin_i':
        bb = []
        for x in th:
            p = S.unpack(x); g = lambda n: x[S.ix[n]]
            bb.append([p['aR'] * np.cos(np.radians(g('i0') + g('idot') * (y - TREF))) for y in (TREF - 0.5, TREF + 0.5)])
        bb = np.array(bb); der['dbdt_per_yr'] = bb[:, 1] - bb[:, 0]
if 'bdot' in S.ix: der['dbdt_per_yr'] = th[:, S.ix['bdot']]
der['aR'] = aR_from_lnrho(th[:, S.ix['lnrho']], P_D)
# importance-sampling evidence with a multivariate-t proposal fitted to the posterior: integrates the profiled likelihood
# (baselines not marginalised) times the unnormalised prior weights of lnpost; t20_summary.py adds the uniform-prior
# width of bdot only (the i0 prior of lin_i is not normalised) -> profile-evidence diagnostics, not normalised Bayes factors
mu = flat.mean(0); cov = np.cov(flat.T) * 1.3
qd = stats.multivariate_t(loc=mu, shape=cov, df=5)
rng = np.random.default_rng(7); NIS = 6000
xs = qd.rvs(size=NIS, random_state=rng); lq = qd.logpdf(xs)
lps = np.array([S.lnpost(x) for x in xs]); wts = lps - lq
good = np.isfinite(wts); mx = wts[good].max()
lnZ = mx + np.log(np.sum(np.exp(wts[good] - mx)) / NIS)
ess = np.sum(np.exp(wts[good] - mx)) ** 2 / np.sum(np.exp(2 * (wts[good] - mx)))
# Laplace evidence at MAP (numerical Hessian of lnpost via posterior covariance)
lnZ_lap = lpmap + 0.5 * ndim * np.log(2 * np.pi) + 0.5 * np.linalg.slogdet(np.cov(flat.T))[1]
def q(x): return [float(np.percentile(x, 50)), float(np.percentile(x, 16)), float(np.percentile(x, 84))]
out = dict(variant=var, spec=kw, names=S.names, npts=int(npts), ntr=int(D.n), nens=nens, nwalk=nwalk, nsteps=nsteps,
           map=dict(zip(S.names, map(float, xmap))), chi2_map=float(chimap), lnpost_map=float(lpmap), bic=float(bic),
           lnZ_is=float(lnZ), is_ess=float(ess), lnZ_laplace=float(lnZ_lap),
           tau=[t.tolist() for t in taus], acc=accs, rhat=rhat, post={n: q(flat[:, i]) for i, n in enumerate(S.names)},
           derived={k: (q(v) if v.ndim == 1 else [q(v[:, j]) for j in range(v.shape[1])]) for k, v in der.items()},
           camp_years=yrs.tolist(), runtime_s=time.time() - t0)
json.dump(out, open(outb + '.json', 'w'), indent=1)
np.savez_compressed(outb + '_samples.npz', samples=flat[::max(1, len(flat) // 40000)], names=np.array(S.names), **{k: v for k, v in der.items()})
print(json.dumps({k: out[k] for k in ('map', 'chi2_map', 'bic', 'lnZ_is', 'is_ess', 'lnZ_laplace', 'rhat')}, indent=None), flush=True)
print('post', json.dumps(out['post']), flush=True); print('derived', json.dumps(out['derived']), flush=True)
