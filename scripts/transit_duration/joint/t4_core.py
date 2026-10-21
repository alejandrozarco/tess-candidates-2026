"""T4: joint transit model of all TESS transits of L 98-59 d. Own code.
Geometry per transit (b_i), optional per-campaign k, Kipping q1,q2, stellar density -> a/R*, optional e/omega (local
approximation near conjunction), TTV sinusoid or fixed per-transit times. Per-transit quadratic baselines are profiled
(linear least squares). batman's C quadratic-LD routine is called on our own sky-separation array."""
import numpy as np, sys
from batman import _quadratic_ld
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
P_D = PL['d'][0]; T0_D = PL['d'][1]
TREF = 2022.0    # reference epoch (yr) for b(t) and T14(t)

class TData:
    def __init__(self, segs, ss=5, near=0.03, nb=3, sig=None, mult=None, regs=False):
        self.segs = segs; self.n = len(segs); self.nb = nb
        self.tc_lin = np.array([s['tc_lin'] for s in segs]); self.ep = np.array([s['epoch'] for s in segs])
        self.yr = year(self.tc_lin); self.camp = np.array([campaign(s['sector']) for s in segs]); self.sector = np.array([s['sector'] for s in segs])
        self.exptime = np.array([s['exptime'] for s in segs])
        self.sig = np.array([np.std(np.diff(s['f'])) / np.sqrt(2) for s in segs]) if sig is None else np.array(sig)
        self.ss = ss; self.nearw = near; self.regs = regs
        self.mult = np.ones(self.n) if mult is None else np.asarray(mult, float)   # bootstrap multiplicities
        self.build()
    def build(self, keep=None):
        """keep: optional list of per-segment boolean masks (e.g. after clipping)."""
        t, f, ti, X = [], [], [], []
        for i, s in enumerate(self.segs):
            m = np.ones(len(s['t']), bool) if keep is None else keep[i]
            d = s['t'][m] - self.tc_lin[i]
            t.append(s['t'][m]); f.append(s['f'][m]); ti.append(np.full(m.sum(), i))
            Xi = np.vstack([(d / 0.16) ** j for j in range(self.nb)]).T
            if self.regs and 'reg' in s:
                R = np.nan_to_num(s['reg'][m] - np.nanmean(s['reg'][m], 0)); sd = R.std(0); sd[sd == 0] = 1
                Xi = np.hstack([Xi, R / sd])
            X.append(Xi)
        self.keep = keep
        self.t = np.concatenate(t); self.f = np.concatenate(f); self.ti = np.concatenate(ti); self.X = np.concatenate(X)
        self.st = np.r_[0, np.where(np.diff(self.ti) != 0)[0] + 1]
        self.w = self.mult[self.ti] / self.sig[self.ti] ** 2
        nr = np.abs(self.t - self.tc_lin[self.ti]) < self.nearw
        self.inear = np.where(nr)[0]
        off = (np.arange(self.ss) - (self.ss - 1) / 2.) / self.ss
        self.tss = self.t[self.inear][:, None] + off[None, :] * (self.exptime[self.ti[self.inear]][:, None] / 86400.)
        self.ti_near = self.ti[self.inear]
        self.npts = len(self.t)
        self._precompute()
    def set_sig(self, sig):
        self.sig = np.array(sig); self.w = self.mult[self.ti] / self.sig[self.ti] ** 2; self._precompute()
    def model(self, tc, b, k, aR, u1, u2, e=0., w=np.pi / 2):
        """tc, b: per-transit arrays; k: scalar or per-transit array."""
        fac = (1 + e * np.sin(w)) ; rc = aR * (1 - e * e) / fac            # distance at conjunction (R*)
        neff = 2 * np.pi / P_D * fac ** 2 / (1 - e * e) ** 1.5              # local angular rate
        dt = self.tss - tc[self.ti_near][:, None]
        th = neff * dt
        bi = b[self.ti_near][:, None]
        z = np.sqrt((rc * np.sin(th)) ** 2 + (bi * np.cos(th)) ** 2)
        m = np.ones(self.npts)
        if np.ndim(k) == 0:
            fl = _quadratic_ld._quadratic_ld(np.ascontiguousarray(z.ravel()), float(k), float(u1), float(u2), 1).reshape(z.shape)
        else:
            kk = k[self.ti_near]; fl = np.empty_like(z)
            for kv in np.unique(kk):
                s = kk == kv
                fl[s] = _quadratic_ld._quadratic_ld(np.ascontiguousarray(z[s].ravel()), float(kv), float(u1), float(u2), 1).reshape(z[s].shape)
        m[self.inear] = fl.mean(1)
        return m
    def _precompute(self):
        """Constant (far-from-transit, m = 1) contributions to the per-transit normal equations."""
        nb = self.X.shape[1]; n = self.n; self.nbx = nb
        far = np.ones(self.npts, bool); far[self.inear] = False
        self.pairs = [(j, l) for j in range(nb) for l in range(j, nb)]
        wx = self.w * far
        self.Mfar = np.zeros((n, nb, nb)); self.Vfar = np.zeros((n, nb))
        for (j, l) in self.pairs:
            v = np.bincount(self.ti, wx * self.X[:, j] * self.X[:, l], minlength=n); self.Mfar[:, j, l] = v; self.Mfar[:, l, j] = v
        for j in range(nb): self.Vfar[:, j] = np.bincount(self.ti, wx * self.f * self.X[:, j], minlength=n)
        self.Sff = np.bincount(self.ti, self.w * self.f ** 2, minlength=n)
        self.Xn = self.X[self.inear]; self.wn = self.w[self.inear]; self.fn = self.f[self.inear]
    def chi2(self, m, ret=False):
        if ret is True:
            A = m[:, None] * self.X; Aw = A * self.w[:, None]
            M = np.einsum('ni,nj->nij', Aw, A); V = Aw * self.f[:, None]
            Ms = np.add.reduceat(M, self.st, axis=0); Vs = np.add.reduceat(V, self.st, axis=0)
            c = np.linalg.solve(Ms, Vs[..., None])[..., 0]
            mod = np.sum(A * c[self.ti], axis=1)
            return self.f - mod, mod, c
        mn = m[self.inear]; n = self.n; tn = self.ti_near
        M = self.Mfar.copy(); V = self.Vfar.copy()
        wm2 = self.wn * mn * mn; wmf = self.wn * mn * self.fn
        for (j, l) in self.pairs:
            v = np.bincount(tn, wm2 * self.Xn[:, j] * self.Xn[:, l], minlength=n); M[:, j, l] += v
            if j != l: M[:, l, j] += v
        for j in range(self.nbx): V[:, j] += np.bincount(tn, wmf * self.Xn[:, j], minlength=n)
        ok = self.mult > 0
        c = np.zeros_like(V); c[ok] = np.linalg.solve(M[ok], V[ok][..., None])[..., 0]
        if ret == 'resid':
            mod = m * np.sum(self.X * c[self.ti], axis=1)
            return (self.f - mod) * np.sqrt(self.w)
        return float(np.sum(self.Sff - np.sum(V * c, axis=1)))

def t14_from(P, aR, b, k, e=0., w=np.pi / 2):
    fac = 1 + e * np.sin(w); rc = aR * (1 - e * e) / fac; neff = 2 * np.pi / P * fac ** 2 / (1 - e * e) ** 1.5
    x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(rc ** 2 - b ** 2 + 1e-12)
    return 2 / neff * np.arcsin(np.clip(x, 0, 1))

def ld_q2u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)

class Spec:
    """Parameter layout and priors for one model variant."""
    def __init__(self, geom='lin_b', kmode='common', ld='prior', ecc='zero', timing='ttv', rho='prior', tc_fixed=None, D=None, ld0=None):
        self.geom, self.kmode, self.ld, self.ecc, self.timing, self.rho = geom, kmode, ld, ecc, timing, rho
        self.tc_fixed = tc_fixed; self.D = D; self.ld0 = ld0
        names = []
        names += ['k%d' % c for c in range(4)] if kmode == 'camp' else ['k']
        names += ['lnrho', 'q1', 'q2']
        if geom == 'const': names += ['b0']
        elif geom == 'lin_b': names += ['b0', 'bdot']
        elif geom == 'lin_i': names += ['i0', 'idot']
        elif geom == 'camp': names += ['bc%d' % c for c in range(4)]
        elif geom == 'camp_slope': names += ['bc%d' % c for c in range(4)] + ['bdot']
        if ecc == 'lit': names += ['secw', 'sesw']
        if timing == 'ttv': names += ['dT0', 'dP', 'At', 'Bt']
        self.names = names; self.ix = {n: i for i, n in enumerate(names)}
        self.camp_mid = np.array([np.mean(D.yr[D.camp == c]) for c in range(4)]) if D is not None else None
    def unpack(self, th):
        g = lambda n: th[self.ix[n]]
        D = self.D
        k = np.array([g('k%d' % c) for c in range(4)])[D.camp] if self.kmode == 'camp' else g('k')
        lnrho = g('lnrho'); aR = aR_from_lnrho(lnrho, P_D)
        u1, u2 = ld_q2u(g('q1'), g('q2'))
        if self.ecc == 'lit':
            e = g('secw') ** 2 + g('sesw') ** 2; w = np.arctan2(g('sesw'), g('secw'))
        else: e, w = 0., np.pi / 2
        dy = D.yr - TREF
        if self.geom == 'const': b = np.full(D.n, g('b0'))
        elif self.geom == 'lin_b': b = g('b0') + g('bdot') * dy
        elif self.geom == 'lin_i':
            inc = np.radians(g('i0') + g('idot') * dy); fac = 1 + e * np.sin(w)
            b = aR * np.cos(inc) * (1 - e * e) / fac
        elif self.geom in ('camp', 'camp_slope'):
            b = np.array([g('bc%d' % c) for c in range(4)])[D.camp]
            if self.geom == 'camp_slope': b = b + g('bdot') * (D.yr - self.camp_mid[D.camp])
        if self.timing == 'ttv':
            tc = D.tc_lin + g('dT0') / 1440. + D.ep * (g('dP') / 86400.) + (g('At') * np.sin(2 * np.pi * D.tc_lin / PTTV) + g('Bt') * np.cos(2 * np.pi * D.tc_lin / PTTV)) / 1440.
        else: tc = self.tc_fixed
        return dict(k=k, aR=aR, u1=u1, u2=u2, e=e, w=w, b=b, tc=tc, lnrho=lnrho)
    def logprior(self, th):
        g = lambda n: th[self.ix[n]]
        lp = 0.
        for n in self.names:
            if n.startswith('k') and not (0.01 < g(n) < 0.1): return -np.inf
        if not (0 < g('q1') < 1 and 0 < g('q2') < 1): return -np.inf
        if self.rho == 'prior': lp += -0.5 * ((g('lnrho') - LNRHO[0]) / LNRHO[1]) ** 2
        elif not (8.0 < g('lnrho') < 11.0): return -np.inf
        if self.ld == 'prior':
            u1, u2 = ld_q2u(g('q1'), g('q2')); lp += -0.5 * (((u1 - self.ld0[0]) / self.ld0[2]) ** 2 + ((u2 - self.ld0[1]) / self.ld0[2]) ** 2)
        if 'bdot' in self.ix and not (-0.05 < g('bdot') < 0.05): return -np.inf
        if 'idot' in self.ix and not (-0.5 < g('idot') < 0.5): return -np.inf
        if 'i0' in self.ix and not (80 < g('i0') < 90): return -np.inf
        if self.ecc == 'lit':
            if g('secw') ** 2 + g('sesw') ** 2 > 0.3: return -np.inf
            lp += -0.5 * (((g('secw') - 0.05) / 0.04) ** 2 + ((g('sesw') + 0.03) / 0.055) ** 2)   # Cadieux+2025 Table (joint RV-TTV)
        if self.timing == 'ttv':
            for n, lim in (('dT0', 30), ('dP', 30), ('At', 15), ('Bt', 15)):
                if abs(g(n)) > lim: return -np.inf
        return lp
    def lnpost(self, th, retchi=False):
        lp = self.logprior(th)
        if not np.isfinite(lp): return (-np.inf, np.inf) if retchi else -np.inf
        p = self.unpack(th)
        kmax = np.max(p['k'])
        if np.any(p['b'] < 0) or np.any(p['b'] > 1 + kmax - 0.002) or np.any(p['b'] > 0.99 * p['aR']): return (-np.inf, np.inf) if retchi else -np.inf
        m = self.D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2'], p['e'], p['w'])
        c2 = self.D.chi2(m)
        return (lp - 0.5 * c2, c2) if retchi else lp - 0.5 * c2
    def prior_resid(self, th):
        g = lambda n: th[self.ix[n]]; r = []
        if self.rho == 'prior': r.append((g('lnrho') - LNRHO[0]) / LNRHO[1])
        if self.ld == 'prior':
            u1, u2 = ld_q2u(g('q1'), g('q2')); r += [(u1 - self.ld0[0]) / self.ld0[2], (u2 - self.ld0[1]) / self.ld0[2]]
        if self.ecc == 'lit': r += [(g('secw') - 0.05) / 0.04, (g('sesw') + 0.03) / 0.055]
        return np.array(r)
    def resid(self, th):
        """Weighted residual vector with Gaussian-prior residuals appended (for least-squares fitting)."""
        lp = self.logprior(th)
        npr = len(self.prior_resid(self.x_ok)) if hasattr(self, 'x_ok') else None
        p = self.unpack(th) if np.isfinite(lp) else None
        bad = p is None or np.any(p['b'] < 0) or np.any(p['b'] > 1 + np.max(p['k']) - 0.002)
        if bad:
            return np.full(self.D.npts + (npr if npr is not None else 0), 1e3)
        self.x_ok = th.copy()
        m = self.D.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2'], p['e'], p['w'])
        r = self.D.chi2(m, ret='resid')
        return np.r_[r, self.prior_resid(th)]
    def t14_at(self, th, yr):
        """T14 (min) at epoch(s) yr for parameter vector th."""
        p = self.unpack(th); g = lambda n: th[self.ix[n]]
        k = p['k'] if np.ndim(p['k']) == 0 else np.mean(p['k'])
        dy = np.atleast_1d(yr) - TREF
        if self.geom == 'lin_b': b = g('b0') + g('bdot') * dy
        elif self.geom == 'const': b = np.full(len(dy), g('b0'))
        elif self.geom == 'lin_i':
            fac = 1 + p['e'] * np.sin(p['w']); b = p['aR'] * np.cos(np.radians(g('i0') + g('idot') * dy)) * (1 - p['e'] ** 2) / fac
        else: return None
        return t14_from(P_D, p['aR'], b, k, p['e'], p['w']) * 1440.
