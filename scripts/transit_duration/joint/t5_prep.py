"""T5: prepare a window set for the joint inference: MAP (evolving b, TTV sinusoid), residual clipping (4 sigma up,
6 down), per-transit white sigma from MAD, red-noise beta per sector (time-averaging, 10-30 min bins), per-transit
mid-times (for the fixed-times variant). Own code. Usage: t5_prep.py <win.npy> <out_prep.npz> [nb]"""
import sys, json, numpy as np, time
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *; from t4_core import *
win, outp = sys.argv[1], sys.argv[2]; nb = int(sys.argv[3]) if len(sys.argv) > 3 else 3
LD = json.load(open(ROOT + '/out/claret_ld.json')); ld0 = (LD['u1'], LD['u2'], LD['sigma_prior'])
segs = list(np.load(win, allow_pickle=True))
D = TData(segs, nb=nb)
def u2q(u1, u2): return (u1 + u2) ** 2, u1 / (2 * (u1 + u2))
q1, q2 = u2q(ld0[0], ld0[1])
def mapfit(S, x0, n=3):
    f = lambda th: -S.lnpost(th)
    best = None
    for i in range(n):
        r = minimize(f, x0 if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-7, fatol=1e-4, adaptive=True))
        r = minimize(f, r.x, method='Powell', options=dict(maxiter=20000, xtol=1e-6, ftol=1e-8))
        if best is None or r.fun < best.fun - 1e-6: best = r
        else: break
    return best
S = Spec('lin_b', D=D, ld0=ld0)
x0 = np.array([0.046, 9.5, q1, q2, 0.93, 0.005, 0., 0., 0., 0.])
t = time.time(); [S.lnpost(x0) for _ in range(20)]; print('lnpost eval %.2f ms, npts %d, ntr %d' % ((time.time() - t) / 20 * 1e3, D.npts, D.n), flush=True)
r = mapfit(S, x0); print('MAP1', dict(zip(S.names, np.round(r.x, 5))), 'chi2 %.1f' % (-2 * (r.fun + S.logprior(r.x))), flush=True)
# clipping and sigma
p = S.unpack(r.x); m = D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2']); res, mod, c = D.chi2(m, ret=True)
keep = []; sig = np.zeros(D.n); nclip = 0
for i in range(D.n):
    sel = D.ti == i; rr = res[sel]; s = 1.4826 * np.median(np.abs(rr - np.median(rr)))
    bad = (rr > 4 * s) | (rr < -6 * s); bad = bad | np.roll(bad, 1) | np.roll(bad, -1)
    keep.append(~bad); nclip += bad.sum(); sig[i] = np.std(rr[~bad])
D.build(keep); D.set_sig(sig)
# red-noise beta per sector
p = S.unpack(r.x); m = D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2']); res, mod, c = D.chi2(m, ret=True)
beta_sec = {}
for s in np.unique(D.sector):
    bs = []
    for nbin_min in (10, 15, 20, 30):
        num, den = [], []
        for i in np.where(D.sector == s)[0]:
            rr = res[D.ti == i] / sig[i]; N = int(round(nbin_min * 60 / D.exptime[i])); M = len(rr) // N
            if M < 3: continue
            bm = rr[:M * N].reshape(M, N).mean(1); num.append(np.sum(bm ** 2)); den.append(M * (1. / N) * M / (M - 1))
        if den: bs.append(np.sqrt(np.sum(num) / np.sum(den)))
    beta_sec[int(s)] = max(1.0, float(np.mean(bs)))
beta = np.array([beta_sec[int(s)] for s in D.sector])
D.set_sig(sig * beta)
print('clipped', nclip, 'beta per sector', {k: round(v, 2) for k, v in beta_sec.items()}, flush=True)
r = mapfit(S, r.x); thS = r.x
print('MAP2', dict(zip(S.names, np.round(r.x, 5))), 'chi2 %.1f' % (-2 * (r.fun + S.logprior(r.x))), flush=True)
# per-transit mid-times at the MAP shape (grid + parabola), for the fixed-times variant
p = S.unpack(thS); tcs = p['tc'].copy(); tce = np.zeros(D.n)
for i in range(D.n):
    grid = np.arange(-12, 12.01, 0.25) / 1440.; ch = []
    for g in grid:
        tc = p['tc'].copy(); tc[i] += g
        m = D.model(tc, p['b'], p['k'], p['aR'], p['u1'], p['u2']); ch.append(D.chi2(m))
    ch = np.array(ch); j = np.argmin(ch); j = min(max(j, 2), len(ch) - 3)
    a = np.polyfit(grid[j - 2:j + 3], ch[j - 2:j + 3], 2); tcs[i] = p['tc'][i] - a[1] / (2 * a[0]); tce[i] = 1 / np.sqrt(max(a[0], 1e-9))
oc = (tcs - p['tc']) * 1440
print('per-transit O-C about TTV model: rms %.2f min, median err %.2f min, chi2 %.1f / %d' % (np.std(oc), np.median(tce) * 1440, np.sum((oc / (tce * 1440)) ** 2), D.n), flush=True)
np.savez(outp, keep=np.array(keep, dtype=object), sig=sig * beta, sig_white=sig, beta=beta, map_lin_b=thS, names=np.array(S.names), tc_pertr=tcs, tc_pertr_err=tce, win=win, nb=nb)
print('saved', outp)
