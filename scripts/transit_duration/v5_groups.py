"""Group fits. Usage: v5_groups.py <win.npy> <tim.npz> <P> <exptime> <tcmode: pt|ttv|lin> <ldmode: prior|fixed|free> <grouping: year|sector|half> <out.json> [u1 u2] [aRprior mu sig]"""
import sys, numpy as np, json, time
sys.path.insert(0, __file__.rsplit('/',1)[0])
from v3_core import Data, ld_prior
from lib import ld_q2u, t14, t23, year
from scipy.optimize import minimize
a = sys.argv
win, tim, P, exp, tcmode, ldmode, grouping, out = a[1], a[2], float(a[3]), float(a[4]), a[5], a[6], a[7], a[8]
U0 = (float(a[9]), float(a[10])) if len(a) > 10 else (0.17, 0.44)
aRp = (float(a[11]), float(a[12])) if len(a) > 12 else None
segs = list(np.load(win, allow_pickle=True)); T = np.load(tim, allow_pickle=True)
rows = T['rows']; pars0 = T['pars']; masks = T['masks']; sig = T['sig']
ep = rows[:, 1]; tc_pt = rows[:, 2]; etc = rows[:, 3]
ok = np.isfinite(etc) & (rows[:, 6] == 0)
# linear ephemeris
A = np.vstack([np.ones(ok.sum()), ep[ok]]).T; wt = 1 / etc[ok] ** 2
lin = np.linalg.solve(A.T @ (A * wt[:, None]), A.T @ (wt * tc_pt[ok]))
tc_lin = lin[0] + lin[1] * ep
# sinusoidal TTV (super-period free around 396 d) on top of linear
def ttvfit(Pttv):
    B = np.vstack([np.ones(ok.sum()), ep[ok], np.sin(2*np.pi*tc_pt[ok]/Pttv), np.cos(2*np.pi*tc_pt[ok]/Pttv)]).T
    c = np.linalg.solve(B.T @ (B * wt[:, None]), B.T @ (wt * tc_pt[ok]))
    r = tc_pt[ok] - B @ c; return np.sum(r*r*wt), c
grid = np.linspace(300, 500, 401); chis = [ttvfit(g)[0] for g in grid]; Pttv = grid[int(np.argmin(chis))]
_, cs = ttvfit(Pttv)
tc_ttv = cs[0] + cs[1]*ep + cs[2]*np.sin(2*np.pi*tc_pt/Pttv) + cs[3]*np.cos(2*np.pi*tc_pt/Pttv)
ttv_amp = np.hypot(cs[2], cs[3]) * 1440
chi_lin = np.sum(((tc_pt[ok]-tc_lin[ok])/etc[ok])**2); chi_ttv = min(chis)
tcs = dict(pt=np.where(ok, tc_pt, tc_ttv), ttv=tc_ttv, lin=tc_lin)[tcmode]
yrs = year(tcs); sec = np.array([s['sector'] for s in segs])
if grouping == 'year':
    edges = [0, 2020, 2022, 2024, 2100]; g = np.digitize(yrs, edges) - 1
elif grouping == 'half':
    g = (yrs > np.median(yrs)).astype(int)
elif grouping == 'sector':
    us = np.unique(sec); g = np.searchsorted(us, sec)
elif grouping == 'single':
    g = np.arange(len(segs))
G = np.unique(g)
Ds = []
for gi in G:
    idx = np.where(g == gi)[0]
    D = Data([segs[i] for i in idx], P, exp, tc=tcs[idx])
    D.mask = [masks[i] for i in idx]; D.sig = sig[idx].copy(); D.build()
    Ds.append(D)
gy = np.array([np.mean(yrs[g == gi]) for gi in G]); ntr = np.array([(g == gi).sum() for gi in G])
ng = len(G)
def unpack_ld(th):
    if ldmode == 'fixed': return U0, 0.
    q1, q2 = th
    if not (0 < q1 < 1 and 0 < q2 < 1): return None, 0
    u = ld_q2u(q1, q2)
    return u, (ld_prior(u[0], u[1], U0, 0.1) if ldmode == 'prior' else 0.)
nld = 0 if ldmode == 'fixed' else 2
def q_from_u(u1, u2): return [(u1 + u2) ** 2, u1 / (2 * (u1 + u2))]
def tot(th, kind):
    # kind C: [k,aR,b,(ld)] ; A: [k,aR,b_1..b_ng,(ld)] ; L: [k,aR,b0,bdot,(ld)] ; F: [k_g,aR_g,b_g]*ng + (ld)
    if kind == 'F':
        u, pen = unpack_ld(th[3*ng:3*ng+nld])
    else:
        u, pen = unpack_ld(th[len(th)-nld:] if nld else [])
    if u is None: return 1e30
    c = pen
    if kind != 'F':
        k, aR = th[0], th[1]
        if not (0.01 < k < 0.3 and 5 < aR < 200): return 1e30
        if aRp: c += ((aR - aRp[0]) / aRp[1]) ** 2
    for j, D in enumerate(Ds):
        if kind == 'C': b = th[2]
        elif kind == 'A': b = th[2 + j]
        elif kind == 'L': b = th[2] + th[3] * (gy[j] - 2022.0)
        elif kind == 'F':
            k, aR, b = th[3*j:3*j+3]
            if not (0.01 < k < 0.3 and 5 < aR < 200): return 1e30
            if aRp: c += ((aR - aRp[0]) / aRp[1]) ** 2
        if not (0 <= b < 1 + k): return 1e30
        c += D.chi2(k, aR, b, u[0], u[1])
    return c
def opt(f, x0, n=3):
    r = None
    for it in range(n):
        r2 = minimize(f, x0 if r is None else r.x, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-7, fatol=1e-5, adaptive=True))
        if r is None or r2.fun <= r.fun: r = r2
    return r
k0, aR0, b0 = pars0[0], pars0[1], pars0[2]
ld0 = q_from_u(*U0) if nld else []
t0 = time.time()
rC = opt(lambda th: tot(th, 'C'), [k0, aR0, b0] + ld0)
rA = opt(lambda th: tot(th, 'A'), list(rC.x[:2]) + [rC.x[2]] * ng + list(rC.x[3:]))
rL = opt(lambda th: tot(th, 'L'), list(rC.x[:3]) + [0.0] + list(rC.x[3:]))
xF0 = []
for j in range(ng): xF0 += [rA.x[0], rA.x[1], rA.x[2 + j]]
rF = opt(lambda th: tot(th, 'F'), xF0 + list(rA.x[2 + ng:]))
def ldof(r, kind):
    if not nld: return U0
    return ld_q2u(*(r.x[-nld:]))
res = dict(config=dict(win=win, tcmode=tcmode, ldmode=ldmode, grouping=grouping, U0=U0, aRprior=aRp, P=P, exptime=exp),
           ttv=dict(Pttv=Pttv, amp_min=ttv_amp, chi2_lin=chi_lin, chi2_ttv=chi_ttv, n=int(ok.sum())),
           groups=dict(year=gy.tolist(), ntr=ntr.tolist(), sectors=[sorted(set(sec[g == gi].tolist())) for gi in G], npts=[len(D.f) for D in Ds]),
           C=dict(chi2=rC.fun, x=rC.x.tolist()), A=dict(chi2=rA.fun, x=rA.x.tolist()), L=dict(chi2=rL.fun, x=rL.x.tolist()), F=dict(chi2=rF.fun, x=rF.x.tolist()))
uA = ldof(rA, 'A'); uF = ldof(rF, 'F')
res['A']['T14_min'] = [t14(P, rA.x[1], rA.x[2 + j], rA.x[0]) * 1440 for j in range(ng)]
res['A']['u'] = list(uA)
res['F']['T14_min'] = [t14(P, rF.x[3*j+1], rF.x[3*j+2], rF.x[3*j]) * 1440 for j in range(ng)]
res['F']['T23_min'] = [t23(P, rF.x[3*j+1], rF.x[3*j+2], rF.x[3*j]) * 1440 for j in range(ng)]
res['F']['u'] = list(uF)
res['L']['T14_ends'] = [t14(P, rL.x[1], rL.x[2] + rL.x[3] * (y - 2022), rL.x[0]) * 1440 for y in (gy[0], gy[-1])]
res['C']['T14_min'] = t14(P, rC.x[1], rC.x[2], rC.x[0]) * 1440
res['dchi2_A_vs_C'] = rC.fun - rA.fun; res['dchi2_L_vs_C'] = rC.fun - rL.fun; res['dchi2_F_vs_A'] = rA.fun - rF.fun
res['ndata'] = int(sum(len(D.f) for D in Ds)); res['runtime_s'] = time.time() - t0
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps({k: res[k] for k in ['ttv', 'groups', 'dchi2_A_vs_C', 'dchi2_L_vs_C', 'dchi2_F_vs_A']}, indent=None))
print('C', np.round(rC.x, 4), 'T14 %.2f' % res['C']['T14_min'])
print('A', np.round(rA.x, 4), 'T14', np.round(res['A']['T14_min'], 2))
print('L', np.round(rL.x, 4), 'T14 ends', np.round(res['L']['T14_ends'], 2))
print('F', np.round(rF.x, 4), 'T14', np.round(res['F']['T14_min'], 2), 'T23', np.round(res['F']['T23_min'], 2))
