"""Stage: template fit + per-transit mid-times + flare clipping. Usage: v4_timing.py <win.npy> <P> <exptime_s> <out.npz> k aR b u1 u2"""
import sys, numpy as np, json, time
sys.path.insert(0, __file__.rsplit('/',1)[0])
from v3_core import Data, fit_common, clip, time_grid_fit
from lib import ld_q2u, t14, year
win, P, exp, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
pars = [float(x) for x in sys.argv[5:10]]
U0 = (pars[3], pars[4])
segs = list(np.load(win, allow_pickle=True))
D = Data(segs, P, exp)
t0 = time.time()
for it in range(3):
    # timing with current template
    for i in range(len(segs)):
        tc, e, c, c0, dof, edge = time_grid_fit(D, i, pars)
        if np.isfinite(e) and not edge: D.tc[i] = tc
    D.build()
    nc = clip(D, pars)
    # global common-shape fit (LD with Claret prior sigma 0.1)
    def q_from_u(u1, u2):
        q1 = (u1 + u2) ** 2; q2 = u1 / (2 * (u1 + u2)); return q1, q2
    q1, q2 = q_from_u(*pars[3:5])
    r = fit_common(D, [[pars[0], pars[1], pars[2], q1, q2]], u0=U0, su=0.1)
    k, aR, b, q1, q2 = r.x; u1, u2 = ld_q2u(q1, q2)
    pars = [k, aR, b, u1, u2]
    print('iter', it, 'clipped', nc, 'chi2', round(r.fun, 1), 'N', len(D.f), 'k %.4f aR %.2f b %.4f u %.3f %.3f T14 %.2f min' % (k, aR, b, u1, u2, t14(P, aR, b, k) * 1440), round(time.time() - t0), 's', flush=True)
rows = []
for i in range(len(segs)):
    tc, e, c, c0, dof, edge = time_grid_fit(D, i, pars)
    rows.append((segs[i]['sector'], segs[i]['epoch'], tc, e, c0 - c, dof, edge))
rows = np.array(rows, dtype=float)
np.savez(out, rows=rows, pars=np.array(pars), tc=D.tc, sig=D.sig, masks=np.array(D.mask, dtype=object), allow_pickle=True)
# linear ephemeris from per-transit times
ok = np.isfinite(rows[:, 3]) & (rows[:, 6] == 0)
A = np.vstack([np.ones(ok.sum()), rows[ok, 1]]).T
wt = 1 / rows[ok, 3] ** 2
cov = np.linalg.inv(A.T @ (A * wt[:, None])); sol = cov @ (A.T @ (wt * rows[ok, 2]))
oc = (rows[ok, 2] - A @ sol) * 1440
chi2r = np.sum((oc / (rows[ok, 3] * 1440)) ** 2) / (ok.sum() - 2)
print('linear ephem T0 %.6f P %.7f ; O-C rms %.2f min ; median err %.2f min ; chi2_red %.2f' % (sol[0], sol[1], np.std(oc), np.median(rows[ok, 3]) * 1440, chi2r))
for r_, o in zip(rows[ok], oc):
    print('  S%02d E%4d tc %.5f err %.2f min  O-C %+6.2f min  dchi2 %.0f' % (r_[0], r_[1], r_[2], r_[3] * 1440, o, r_[4]))
