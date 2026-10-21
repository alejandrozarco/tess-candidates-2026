"""Own transit-analysis helpers for the independent verification (not the survey's code)."""
import numpy as np, glob, os, batman
WORK = os.environ.get('TD_WORK', os.path.join('work', 'transit_duration'))   # work directory, relative to the repository root
LCDIR = os.environ.get('LCDIR', os.path.join(WORK, 'lc'))                     # extracted light curves (v1_extract.py)

def load(tic, cad='2m', flux='pdc'):
    out = []
    for f in sorted(glob.glob(os.path.join(LCDIR, '%d_s*_%s.npz' % (tic, cad)))):
        z = np.load(f, allow_pickle=True)
        m = z['meta'][0]
        out.append(dict(sector=m['sector'], t=z['t'], f=z[flux], e=z[flux + 'e'], meta=m))
    return out

def year(btjd):
    return 2000.0 + (btjd + 2457000.0 - 2451544.5) / 365.25

def ld_q2u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)

def t14(P, aR, b, k):
    si = np.sqrt(1 - (b / aR) ** 2)
    x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * si)
    return P / np.pi * np.arcsin(np.clip(x, 0, 1))

def t23(P, aR, b, k):
    si = np.sqrt(1 - (b / aR) ** 2)
    x = np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (aR * si)
    return P / np.pi * np.arcsin(np.clip(x, 0, 1))

class Model:
    """batman model on concatenated, per-transit-centred times (dt = t - tc_i)."""
    def __init__(self, dt, P, exptime_s, ss=None):
        self.p = batman.TransitParams()
        p = self.p
        p.t0 = 0.; p.per = P; p.rp = 0.04; p.a = 30.; p.inc = 89.; p.ecc = 0.; p.w = 90.; p.limb_dark = 'quadratic'; p.u = [0.2, 0.4]
        if ss is None: ss = max(1, int(round(exptime_s / 20.)))
        if exptime_s > 30:
            self.m = batman.TransitModel(p, dt, supersample_factor=ss, exp_time=exptime_s / 86400.)
        else:
            self.m = batman.TransitModel(p, dt)
    def lc(self, k, aR, b, u1, u2):
        p = self.p
        if b >= aR: return None
        p.rp = k; p.a = aR; p.inc = np.degrees(np.arccos(b / aR)); p.u = [u1, u2]
        return self.m.light_curve(p)

def profile_baseline(f, w, m, seg_starts, X):
    """For each segment fit f ~ m * (X @ c) by weighted LSQ; return model and chi2. X: (n, nb) basis."""
    A = m[:, None] * X
    nb = X.shape[1]
    Aw = A * w[:, None]
    # per segment normal equations
    M = np.einsum('ni,nj->nij', Aw, A)
    V = Aw * f[:, None]
    Ms = np.add.reduceat(M, seg_starts, axis=0)
    Vs = np.add.reduceat(V, seg_starts, axis=0)
    c = np.linalg.solve(Ms, Vs[..., None])[..., 0]
    seglen = np.diff(np.append(seg_starts, len(f)))
    cfull = np.repeat(c, seglen, axis=0)
    mod = np.sum(A * cfull, axis=1)
    return mod, c
