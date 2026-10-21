"""Per-sector (or per-transit) measurements with a fixed global template.
(i) b-only fit (k, aR, LD fixed) -> T14 ; (ii) stretch w + amplitude A on the template (shape-agnostic width).
Usage: v6_persector.py <win.npy> <tim.npz> <P> <exptime> <tcmode pt|ttv|lin> <unit sector|transit> <out.json> k aR b u1 u2"""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/',1)[0])
from v3_core import Data
from lib import t14, year, Model, profile_baseline
from scipy.optimize import minimize_scalar, minimize
a = sys.argv
win, tim, P, exp, tcmode, unit, out = a[1], a[2], float(a[3]), float(a[4]), a[5], a[6], a[7]
k, aR, b0, u1, u2 = [float(x) for x in a[8:13]]
segs = list(np.load(win, allow_pickle=True)); T = np.load(tim, allow_pickle=True)
rows = T['rows']; masks = T['masks']; sig = T['sig']
ep = rows[:, 1]; tc_pt = rows[:, 2]; etc = rows[:, 3]; ok = np.isfinite(etc) & (rows[:, 6] == 0)
wt = 1 / etc[ok] ** 2
Alin = np.vstack([np.ones(ok.sum()), ep[ok]]).T
lin = np.linalg.solve(Alin.T @ (Alin * wt[:, None]), Alin.T @ (wt * tc_pt[ok])); tc_lin = lin[0] + lin[1] * ep
def ttvfit(Pt):
    B = np.vstack([np.ones(ok.sum()), ep[ok], np.sin(2*np.pi*tc_pt[ok]/Pt), np.cos(2*np.pi*tc_pt[ok]/Pt)]).T
    c = np.linalg.solve(B.T @ (B * wt[:, None]), B.T @ (wt * tc_pt[ok])); r = tc_pt[ok] - B @ c; return np.sum(r*r*wt), c
grid = np.linspace(300, 500, 401); Pt = grid[int(np.argmin([ttvfit(g)[0] for g in grid]))]; cs = ttvfit(Pt)[1]
tc_ttv = cs[0] + cs[1]*ep + cs[2]*np.sin(2*np.pi*tc_pt/Pt) + cs[3]*np.cos(2*np.pi*tc_pt/Pt)
sec = np.array([s['sector'] for s in segs])
tc_sec = tc_lin.copy()
for ss in np.unique(sec):
    m = (sec == ss) & ok
    if m.sum(): tc_sec[sec == ss] = tc_lin[sec == ss] + np.sum((tc_pt[m] - tc_lin[m]) / etc[m] ** 2) / np.sum(1 / etc[m] ** 2)
tcs = dict(pt=np.where(ok, tc_pt, tc_ttv), ttv=tc_ttv, lin=tc_lin, sec=tc_sec)[tcmode]
units = np.unique(sec) if unit == 'sector' else np.arange(len(segs))
res = []
for uu in units:
    idx = np.where(sec == uu)[0] if unit == 'sector' else np.array([uu])
    D = Data([segs[i] for i in idx], P, exp, tc=tcs[idx]); D.mask = [masks[i] for i in idx]; D.sig = sig[idx].copy(); D.build()
    # (i) b-only
    import os; bmin = float(os.environ.get('BMIN', 0.80)); bg = np.linspace(bmin, min(1 + k, 1.04), 241 if bmin > 0.5 else 521)
    ch = np.array([D.chi2(k, aR, bb, u1, u2) for bb in bg]); j = np.argmin(ch)
    chi_flat = np.sum(((D.f - profile_baseline(D.f, D.w, np.ones(len(D.f)), D.st, D.X)[0]) ** 2) * D.w)
    dof = len(D.f) - 3 * len(idx) - 1
    s = np.sqrt(max(1, ch[j] / dof))
    inside = bg[(ch - ch[j]) / s ** 2 <= 1]
    T = t14(P, aR, bg, k) * 1440
    Tin = T[(ch - ch[j]) / s ** 2 <= 1]
    # (ii) stretch + amplitude on template (template b = b0)
    base = Model(D.dt, P, exp).lc(k, aR, b0, u1, u2)
    def chi_w(w):
        mt = Model(D.dt / w, P / 1.0, exp / w).lc(k, aR, b0, u1, u2)  # time-stretched template, exposure scaled consistently
        dep = 1 - mt
        # amplitude profiled jointly with baseline: f = (1 - A dep) * poly -> nonlinear in A; do 1D inner min
        def chiA(A):
            return np.sum((D.f - profile_baseline(D.f, D.w, 1 - A * dep, D.st, D.X)[0]) ** 2 * D.w)
        r = minimize_scalar(chiA, bounds=(0.2, 2.5), method='bounded', options=dict(xatol=1e-4))
        return r.fun, r.x
    wg = np.linspace(0.6, 1.5, 91)
    cw = np.array([chi_w(w)[0] for w in wg]); jw = np.argmin(cw)
    win_ = wg[(cw - cw[jw]) / s ** 2 <= 1]
    Abest = chi_w(wg[jw])[1]
    yr = float(np.mean(year(tcs[idx])))
    res.append(dict(unit=int(uu), sectors=sorted(set(sec[idx].tolist())), year=yr, ntr=len(idx), b=float(bg[j]), b_lo=float(inside.min()), b_hi=float(inside.max()),
                    T14=float(T[j]), T14_lo=float(Tin.min()), T14_hi=float(Tin.max()), w=float(wg[jw]), w_lo=float(win_.min()), w_hi=float(win_.max()), A=float(Abest),
                    dchi2_transit=float(chi_flat - ch[j]), chi2red=float(ch[j] / dof), edge_b=bool(j in (0, len(bg) - 1))))
    r_ = res[-1]
    print('%s %4s yr %.2f ntr %d  b %.3f [%.3f,%.3f]  T14 %.1f [%.1f,%.1f]  w %.3f [%.3f,%.3f] A %.2f  dchi2 %.0f chi2r %.2f' % (unit, uu, yr, len(idx), r_['b'], r_['b_lo'], r_['b_hi'], r_['T14'], r_['T14_lo'], r_['T14_hi'], r_['w'], r_['w_lo'], r_['w_hi'], r_['A'], r_['dchi2_transit'], r_['chi2red']), flush=True)
# regressions with empirical scaling
def reg(x, y, e):
    A = np.vstack([np.ones_like(x), x - 2022]).T; W = 1 / e ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None])); c = cov @ (A.T @ (W * y)); r = y - A @ c
    chi = np.sum(r * r * W); dof = len(y) - 2; sc = np.sqrt(max(1, chi / dof))
    chic = np.sum((y - np.sum(W * y) / np.sum(W)) ** 2 * W)
    return dict(slope=float(c[1]), err_formal=float(np.sqrt(cov[1, 1])), err_scaled=float(np.sqrt(cov[1, 1]) * sc), chi2_lin=float(chi), chi2_const=float(chic), dof=dof,
                z_formal=float(c[1] / np.sqrt(cov[1, 1])), z_scaled=float(c[1] / np.sqrt(cov[1, 1]) / sc))
x = np.array([r['year'] for r in res])
good = np.array([not r['edge_b'] for r in res])
summ = {}
for lab, y, e in [('T14', [r['T14'] for r in res], [max(0.05, (r['T14_hi'] - r['T14_lo']) / 2) for r in res]),
                  ('b', [r['b'] for r in res], [max(1e-4, (r['b_hi'] - r['b_lo']) / 2) for r in res]),
                  ('w', [r['w'] for r in res], [max(1e-3, (r['w_hi'] - r['w_lo']) / 2) for r in res])]:
    y = np.array(y); e = np.array(e)
    summ[lab] = reg(x[good], y[good], e[good])
    # robust: unweighted OLS with bootstrap over units
    rng = np.random.default_rng(0); sl = []
    xx, yy = x[good], y[good]
    for _ in range(2000):
        ii = rng.integers(0, len(xx), len(xx))
        if np.ptp(xx[ii]) < 1: continue
        sl.append(np.polyfit(xx[ii], yy[ii], 1)[0])
    summ[lab]['ols_slope'] = float(np.polyfit(xx, yy, 1)[0]); summ[lab]['boot_err'] = float(np.std(sl)); summ[lab]['boot_z'] = float(summ[lab]['ols_slope'] / np.std(sl))
    print(lab, json.dumps(summ[lab]))
json.dump(dict(config=a[1:], template=[k, aR, b0, u1, u2], Pttv=float(Pt), units=res, regress=summ), open(out, 'w'), indent=1)
