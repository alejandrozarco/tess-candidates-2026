"""C9: dependence of the per-detector white-light T14 of d on the systematics model. Own code.
Per NIRSpec visit and detector (own white light, $TDJ_WORK/out/jwst/wl_*.npz, and the pipeline NRS2 white light
wlpipe_*.npz, of ../joint/run_joint.sh), least-squares fit with fixed P, mid-time and a/R* (white-light joint MAP), free k and b, quadratic LD
with the Kurucz prior (sigma 0.1); systematics variants: quadratic only; quadratic + each regressor alone; + all three;
linear + all three; cubic + all three; + exponential ramp (15 min); trimmed to matched coverage (30 min before first
contact and after last contact). T14 error from the curvature of the profile chi2 in b. Usage: c9_sysvariants.py <out.json>"""
import sys, json, os, numpy as np, batman
from scipy.optimize import minimize
from scipy.ndimage import median_filter
W = os.path.join(os.environ.get('TDJ_WORK', 'work/transit_duration_joint'), 'out', 'jwst')   # white light of ../joint/run_joint.sh step 9
H = os.environ.get('TDC_WORK', 'work/transit_duration_chromatic')
LD = json.load(open(H + '/out/ld_kurucz.json'))['bins']
P, AR = 7.450729, 36.25; TC = {'1224': 3121.11245, '4098': 3381.88704}
def t14(b, k): x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(AR ** 2 - b ** 2); return P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
VAR = {'poly2': (2, [], 0, None), 'poly2+xshift': (2, ['xshift'], 0, None), 'poly2+yshift': (2, ['yshift'], 0, None), 'poly2+width': (2, ['width'], 0, None),
       'poly2+regs (section 6.6.10)': (2, ['xshift', 'yshift', 'width'], 0, None), 'poly1+regs': (1, ['xshift', 'yshift', 'width'], 0, None), 'poly3+regs': (3, ['xshift', 'yshift', 'width'], 0, None),
       'poly2+regs+ramp15': (2, ['xshift', 'yshift', 'width'], 15, None), 'poly2+regs, +-30 min of contacts': (2, ['xshift', 'yshift', 'width'], 0, 30.)}
out = {}
for v in ('1224', '4098'):
    for det, fn in (('NRS1', 'wl_%s_nrs1.npz' % v), ('NRS2', 'wl_%s_nrs2.npz' % v), ('NRS2pipe', 'wlpipe_%s_nrs2.npz' % v)):
        z = np.load(W + '/' + fn); t0 = z['t'].min(); keep = (z['t'] - t0) * 1440 > 10
        def bn(x): x = x[keep]; n = len(x) // 5 * 5; return x[:n].reshape(-1, 5).mean(1)
        regs_av = [r for r in ('xshift', 'yshift', 'width') if r in z.files]
        tb, fb = bn(z['t']), bn(z['f']); R = {r: bn(z[r]) for r in regs_av}
        rm = fb - median_filter(fb, 15, mode='nearest'); s0 = 1.4826 * np.median(np.abs(rm)); ok0 = np.abs(rm) < 5 * s0
        u0 = np.array([LD[det[:4] + '_white'][0]['u1'], LD[det[:4] + '_white'][0]['u2']])
        seen = set()
        for name, (poly, regs, ramp, win) in VAR.items():
            regs = [r for r in regs if r in R]          # the pipeline white light has no regressors
            if (poly, tuple(regs), ramp, win) in seen: continue
            seen.add((poly, tuple(regs), ramp, win))
            ok = ok0.copy()
            if win: ok &= np.abs(tb - TC[v]) * 1440 < 22.5 + win
            t, y = tb[ok], fb[ok]; tau = (t - t.mean()) / (t.max() - t.min())
            cols = [tau ** j for j in range(poly + 1)] + [(R[r][ok] - R[r][ok].mean()) / (R[r][ok].std() + 1e-12) for r in regs]
            if ramp: cols.append(np.exp(-(t - t0) * 1440 / ramp))
            X = np.vstack(cols).T
            pm = batman.TransitParams(); pm.t0 = TC[v]; pm.per = P; pm.rp = 0.045; pm.a = AR; pm.inc = 88.4; pm.ecc = 0; pm.w = 90; pm.limb_dark = 'quadratic'; pm.u = list(u0)
            m = batman.TransitModel(pm, t, supersample_factor=3, exp_time=np.median(np.diff(t)) * 2 / 3)   # midpoint sub-exposures
            st = dict(sig=s0)
            def chi(th):
                k, b, u1, u2 = th
                if not (0.01 < k < 0.1 and 0 < b < 1.05 and -0.5 < u1 < 1 and -0.5 < u2 < 1 and 0 < u1 + u2 < 1): return 1e30
                pm.rp = k; pm.inc = np.degrees(np.arccos(b / AR)); pm.u = [u1, u2]; mm = m.light_curve(pm); A = mm[:, None] * X
                c = np.linalg.lstsq(A, y, rcond=None)[0]; r = y - A @ c; st['r'] = r
                return np.sum(r * r) / st['sig'] ** 2 + np.sum(((np.array([u1, u2]) - u0) / 0.1) ** 2)
            best = None
            for b0 in (0.92, 0.93, 0.94):
                r = minimize(chi, [0.045, b0, u0[0], u0[1]], method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-6, adaptive=True))
                if best is None or r.fun < best.fun: best = r
            chi(best.x); st['sig'] = float(np.std(st['r'])); best = minimize(chi, best.x, method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-8, fatol=1e-6, adaptive=True))
            # profile in b for the error
            bg = best.x[1] + np.linspace(-0.006, 0.006, 13); pr = []
            for bb in bg:
                rr = minimize(lambda q: chi([q[0], bb, q[1], q[2]]), best.x[[0, 2, 3]], method='Nelder-Mead', options=dict(maxiter=6000, maxfev=6000, xatol=1e-8, fatol=1e-6, adaptive=True)); pr.append(rr.fun)
            pp = np.polyfit(bg, pr, 2); sb = 1 / np.sqrt(pp[0]) if pp[0] > 0 else np.nan
            T = t14(best.x[1], best.x[0]); sT = abs(t14(best.x[1] + sb, best.x[0]) - T)
            out.setdefault(v, {}).setdefault(det, {})[name] = dict(T14=[float(T), float(sT)], b=[float(best.x[1]), float(sb)], k=float(best.x[0]), u=[float(best.x[2]), float(best.x[3])], n=int(len(t)), sig_ppm=st['sig'] * 1e6)
            print('%s %-8s %-34s T14 %.2f +- %.2f  b %.4f  k %.4f  u %.2f %.2f  n %d' % (v, det, name, T, sT, best.x[1], best.x[0], best.x[2], best.x[3], len(t)), flush=True)
# differences
print('\nDelta T14 (4098 - 1224) per detector and variant')
for det in ('NRS1', 'NRS2', 'NRS2pipe'):
    for name in VAR:
        a, b = out['1224'][det].get(name), out['4098'][det].get(name)
        if a and b: print('  %-8s %-34s %.2f +- %.2f' % (det, name, b['T14'][0] - a['T14'][0], np.hypot(a['T14'][1], b['T14'][1])))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
