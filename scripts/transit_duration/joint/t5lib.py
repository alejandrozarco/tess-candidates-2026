"""Fast preparation + least-squares fitting used by bootstraps, injections and per-transit work. Own code.
Same recipe as t5_prep.py (clip 4 up / 6 down after a first fit, per-transit MAD sigma, per-sector beta) but with
trust-region least squares instead of simplex."""
import numpy as np, json, sys
from scipy.optimize import least_squares
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *; from t4_core import *
LD = json.load(open(ROOT + '/out/claret_ld.json')); LD0 = (LD['u1'], LD['u2'], LD['sigma_prior'])
def u2q(u1, u2): return (u1 + u2) ** 2, u1 / (2 * (u1 + u2))
XS = dict(k=1e-3, lnrho=0.05, q1=0.05, q2=0.05, b0=0.005, bdot=0.001, dT0=0.3, dP=0.1, At=0.3, Bt=0.3, i0=0.02, idot=0.002,
          bc0=0.005, bc1=0.005, bc2=0.005, bc3=0.005, k0=1e-3, k1=1e-3, k2=1e-3, k3=1e-3, secw=0.02, sesw=0.02)
def lsq(S, x0):
    xs = np.array([XS[n] for n in S.names])
    r = least_squares(S.resid, x0, x_scale=xs, method='trf', diff_step=1e-4, max_nfev=200)
    J = r.jac; cov = None
    try: cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError: pass
    lp, c2 = S.lnpost(r.x, retchi=True)
    return r.x, c2, lp, cov
def default_x0(S, b0=0.924, bdot=0.0064, k=0.0441, lnrho=9.56):
    q1, q2 = u2q(LD0[0], LD0[1]); x = []
    for n in S.names:
        x.append(dict(k=k, lnrho=lnrho, q1=q1, q2=q2, b0=b0, bdot=bdot, dT0=0., dP=0., At=-0.8, Bt=1.5, secw=0.05, sesw=-0.03).get(n, k if n.startswith('k') else b0))
    return np.array(x)
def prepare(segs, mult=None, geom='lin_b', nb=3, verbose=False, regs=False):
    D = TData(segs, nb=nb, mult=mult, regs=regs)
    S = Spec(geom, D=D, ld0=LD0)
    x, c2, lp, cov = lsq(S, default_x0(S))
    p = S.unpack(x); m = D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2']); res = D.chi2(m, ret=True)[0]
    keep = []; sig = np.zeros(D.n)
    for i in range(D.n):
        sel = D.ti == i; rr = res[sel]; s = 1.4826 * np.median(np.abs(rr - np.median(rr)))
        bad = (rr > 4 * s) | (rr < -6 * s); bad = bad | np.roll(bad, 1) | np.roll(bad, -1)
        keep.append(~bad); sig[i] = np.std(rr[~bad])
    D.build(keep); D.set_sig(sig)
    m = D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2']); res = D.chi2(m, ret=True)[0]
    beta = np.ones(D.n)
    for s in np.unique(D.sector):
        bs = []
        for nbin in (10, 15, 20, 30):
            num, den = 0., 0.
            for i in np.where(D.sector == s)[0]:
                rr = res[D.ti == i] / sig[i]; N = int(round(nbin * 60 / D.exptime[i])); M = len(rr) // N
                if M < 3: continue
                bm = rr[:M * N].reshape(M, N).mean(1); num += np.sum(bm ** 2); den += M * (1. / N) * M / (M - 1)
            if den > 0: bs.append(np.sqrt(num / den))
        beta[D.sector == s] = max(1.0, np.mean(bs)) if bs else 1.0
    D.set_sig(sig * beta)
    return D
def fit_pair(D, extra=None):
    """LSQ fits of const and lin_b on prepared data; returns dict with bdot, its error, dchi2, dT14/dt."""
    out = {}
    for g in ('const', 'lin_b'):
        S = Spec(g, D=D, ld0=LD0)
        x, c2, lp, cov = lsq(S, default_x0(S))
        out[g] = dict(x=x, chi2=c2, lp=lp, cov=cov, S=S)
    S = out['lin_b']['S']; x = out['lin_b']['x']; cov = out['lin_b']['cov']
    T = S.t14_at(x, np.array([TREF - 0.5, TREF + 0.5])); dT = T[1] - T[0]
    # error of dT14/dt by linear propagation
    eps = 1e-6; gsl = np.zeros(len(x))
    for j in range(len(x)):
        xx = x.copy(); xx[j] += eps * max(1, abs(x[j])); Tj = S.t14_at(xx, np.array([TREF - 0.5, TREF + 0.5])); gsl[j] = ((Tj[1] - Tj[0]) - dT) / (eps * max(1, abs(x[j])))
    edT = float(np.sqrt(gsl @ cov @ gsl)) if cov is not None else np.nan
    ib = S.ix['bdot']
    return dict(bdot=float(x[ib]), bdot_err=float(np.sqrt(cov[ib, ib])) if cov is not None else np.nan, dT14=float(dT), dT14_err=edT,
                dchi2=float(out['const']['chi2'] - out['lin_b']['chi2']), dlp=float(out['lin_b']['lp'] - out['const']['lp']),
                x=x.tolist(), names=S.names, b0=float(x[S.ix['b0']]), k=float(x[S.ix['k']]), lnrho=float(x[S.ix['lnrho']]))
