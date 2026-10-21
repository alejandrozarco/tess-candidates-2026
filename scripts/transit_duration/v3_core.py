"""Core fitting machinery: data container, per-transit timing, group fits. Own code."""
import numpy as np, sys
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/',1)[0]); from lib import Model, profile_baseline, ld_q2u, t14, t23, year

class Data:
    def __init__(self, segs, P, exptime_s, tc=None, nb=3, inwin=None):
        self.segs = segs; self.P = P; self.exptime = exptime_s; self.nb = nb
        self.tc = np.array([s['tc_lin'] for s in segs]) if tc is None else np.array(tc)
        self.mask = [np.ones(len(s['t']), bool) for s in segs]
        self.sig = np.array([np.std(np.diff(s['f'])) / np.sqrt(2) for s in segs])
        self.build()
    def build(self, idx=None):
        idx = range(len(self.segs)) if idx is None else idx
        self.idx = list(idx)
        dt, f, w, st, X = [], [], [], [], []
        n = 0
        for i in self.idx:
            s = self.segs[i]; m = self.mask[i]
            d = s['t'][m] - self.tc[i]
            st.append(n); n += m.sum()
            dt.append(d); f.append(s['f'][m]); w.append(np.full(m.sum(), 1 / self.sig[i] ** 2))
            X.append(np.vstack([d ** j for j in range(self.nb)]).T)
        self.dt = np.concatenate(dt); self.f = np.concatenate(f); self.w = np.concatenate(w)
        self.st = np.array(st); self.X = np.concatenate(X)
        self.model = Model(self.dt, self.P, self.exptime)
    def chi2(self, k, aR, b, u1, u2, ret=False):
        m = self.model.lc(k, aR, b, u1, u2)
        if m is None: return 1e30
        mod, c = profile_baseline(self.f, self.w, m, self.st, self.X)
        r = self.f - mod
        if ret: return r, mod, m
        return np.sum(r * r * self.w)

def ld_prior(u1, u2, u0, su):
    if u0 is None: return 0.
    return ((u1 - u0[0]) / su) ** 2 + ((u2 - u0[1]) / su) ** 2

def fit_common(D, x0, u0=None, su=0.1, fix_ld=None, fix_aR=None, aR_prior=None):
    """theta = k, aR, b, q1, q2 (q's ignored if fix_ld)."""
    def fun(th):
        k, aR, b = th[0], th[1], th[2]
        if fix_aR is not None: aR = fix_aR
        if fix_ld is not None: u1, u2 = fix_ld
        else:
            q1, q2 = th[3], th[4]
            if not (0 < q1 < 1 and 0 < q2 < 1): return 1e30
            u1, u2 = ld_q2u(q1, q2)
        if not (0.005 < k < 0.3 and 2 < aR < 300 and 0 <= b < 1 + k): return 1e30
        pen = ld_prior(u1, u2, u0, su) if fix_ld is None else 0.
        if aR_prior is not None: pen += ((aR - aR_prior[0]) / aR_prior[1]) ** 2
        return D.chi2(k, aR, b, u1, u2) + pen
    best = None
    for start in x0:
        r = minimize(fun, start, method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-6, fatol=1e-4))
        r = minimize(fun, r.x, method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-7, fatol=1e-5))
        if best is None or r.fun < best.fun: best = r
    return best

def clip(D, pars, nsig_up=4.0, nsig_dn=6.0, D14=None):
    k, aR, b, u1, u2 = pars
    m = D.model.lc(k, aR, b, u1, u2)
    mod, c = profile_baseline(D.f, D.w, m, D.st, D.X)
    r = D.f - mod
    nclip = 0
    for j, i in enumerate(D.idx):
        a = D.st[j]; e = D.st[j + 1] if j + 1 < len(D.st) else len(D.f)
        rr = r[a:e]; s = 1.4826 * np.median(np.abs(rr - np.median(rr)))
        D.sig[i] = s
        bad = (rr > nsig_up * s) | (rr < -nsig_dn * s)
        # grow by one cadence each side
        bad = bad | np.roll(bad, 1) | np.roll(bad, -1)
        ii = np.where(D.mask[i])[0]
        nclip += bad.sum()
        D.mask[i][ii[bad]] = False
    D.build(D.idx)
    return nclip

def time_grid_fit(D, i, pars, half_min=20., step_min=0.25):
    """Per-transit mid-time: chi2 profile over shifts with fixed shape and per-transit quadratic baseline."""
    k, aR, b, u1, u2 = pars
    s = D.segs[i]; m_ = D.mask[i]
    t = s['t'][m_]; f = s['f'][m_]; w = np.full(len(t), 1 / D.sig[i] ** 2)
    shifts = np.arange(-half_min, half_min + 1e-9, step_min) / 1440.
    chi = []
    for sh in shifts:
        dt = t - (D.tc[i] + sh)
        mdl = Model(dt, D.P, D.exptime).lc(k, aR, b, u1, u2)
        X = np.vstack([dt ** j for j in range(D.nb)]).T
        mod, _ = profile_baseline(f, w, mdl, np.array([0]), X)
        chi.append(np.sum((f - mod) ** 2 * w))
    chi = np.array(chi); j = np.argmin(chi)
    # flat model chi2 (no transit)
    X = np.vstack([(t - D.tc[i]) ** j for j in range(D.nb)]).T
    mod0, _ = profile_baseline(f, w, np.ones(len(t)), np.array([0]), X)
    chi0 = np.sum((f - mod0) ** 2 * w)
    if 2 <= j <= len(chi) - 3:
        a = np.polyfit(shifts[j-2:j+3], chi[j-2:j+3], 2)
        best = -a[1] / (2 * a[0]); err = 1 / np.sqrt(a[0]) if a[0] > 0 else np.nan
    else:
        best, err = shifts[j], np.nan
    dof = len(t) - D.nb - 1
    return D.tc[i] + best, err * np.sqrt(max(1, chi[j] / dof)), chi[j], chi0, dof, (j in (0, len(chi) - 1))
