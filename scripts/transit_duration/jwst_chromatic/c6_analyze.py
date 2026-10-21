"""C6: chromatic duration test from the per-bin fits (c3_binfit.py outputs). Own code.
Per visit and LD mode:
  - T14, b per bin (posterior median, sd);
  - shared-b (achromatic) model from the summed profile chi2(b) of all bins, and per-detector shared b; T14 at those b
    with the per-bin MAP k averaged; Delta chi2 of free-b-per-bin vs shared b (achromaticity test);
  - step model (one b per detector) vs smooth model (b linear in wavelength across both detectors), both 2 parameters;
  - adjacent bins across the detector gap; within-detector wavelength slopes extrapolated across the gap;
Between NIRSpec visits: per-bin Delta T14 (4098 minus 1224), per-detector weighted means, step vs smooth.
Per visit, chi2_const and chi2_step are local quadratic minima of the summed profiles, while chi2_linear and
chi2_step_plus_linear are evaluated on linearly interpolated profiles. dchi2_linear_vs_const and
dchi2_linear_given_step mix the two approximations and can be negative for nested models; they and
dchi2_step_given_linear (the comparisons involving b linear in wavelength) are withheld from the published tables
(make_tables_chromatic.py). The between-visit regressions are not affected.
Usage: c6_analyze.py <fits_dir> <out.json>"""
import sys, json, glob, os, numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as CHI2
fd, outp = sys.argv[1], sys.argv[2]
def load(name):
    f = os.path.join(fd, name + '.jsonl')
    if not os.path.exists(f): return None
    return [json.loads(s) for s in open(f)]
AR = {'d': 36.25, 'c': 22.14}; P = {'d': 7.450729, 'c': 3.6906764}
def t14(b, k, pl): a = AR[pl]; x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(a ** 2 - b ** 2); return P[pl] / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
def wmean(v, e): w = 1 / np.asarray(e) ** 2; m = np.sum(w * v) / np.sum(w); return float(m), float(1 / np.sqrt(np.sum(w))), float(np.sum(w * (np.asarray(v) - m) ** 2))
def wline(x, v, e):
    w = 1 / np.asarray(e) ** 2; A = np.vstack([np.ones_like(x), x]).T; C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ (A.T @ (w * v))
    return p, C, float(np.sum(w * (v - A @ p) ** 2))
def prof_interp(rows):
    fs = []
    for r in rows:
        g = np.array(r['bgrid']); c = np.array(r['prof_chi2']); ok = np.isfinite(c)
        fs.append((g[ok], c[ok]))
    return fs
def shared_b(fs, sel=None):
    """minimum of the summed profile chi2 over a common b; 1-sigma from Delta chi2 = 1 (parabola near the minimum)."""
    sel = range(len(fs)) if sel is None else sel
    g = fs[sel[0]][0]; tot = np.zeros(len(g))
    for i in sel: tot += np.interp(g, fs[i][0], fs[i][1])
    j = int(np.argmin(tot)); lo, hi = max(0, j - 3), min(len(g), j + 4)
    pp = np.polyfit(g[lo:hi], tot[lo:hi], 2); b0 = -pp[1] / (2 * pp[0]) if pp[0] > 0 else g[j]
    cmin = np.polyval(pp, b0) if pp[0] > 0 else tot[j]; sb = 1 / np.sqrt(pp[0]) if pp[0] > 0 else np.nan
    return float(b0), float(sb), float(cmin), float(np.sum([fs[i][1].min() for i in sel]))
def chi_of_bfunc(fs, bvals):
    return float(np.sum([np.interp(b, f[0], f[1], left=1e6, right=1e6) for f, b in zip(fs, bvals)]))
res = {}
VIS = [('d1224', 'd'), ('d4098', 'd'), ('d1201', 'd'), ('c2512a', 'c'), ('c2512b', 'c')]
for ldm in ('prior', 'free', 'ramp'):
    for vis, pl in VIS:
        dets = [('NRS1', '%s_nrs1_%s' % (vis, ldm)), ('NRS2', '%s_nrs2_%s' % (vis, ldm))] if vis != 'd1201' else [('SOSS1', 'd1201_soss1_%s' % ldm), ('SOSS2', 'd1201_soss2_%s' % ldm)]
        rows, det_of = [], []
        for det, nm in dets:
            r = load(nm)
            if r: rows += r; det_of += [det] * len(r)
        if not rows: continue
        det_of = np.array(det_of); lam = np.array([0.5 * (r['lo'] + r['hi']) for r in rows])
        T = np.array([r['post']['T14'][0] for r in rows]); Te = np.array([r['post']['T14'][1] for r in rows])
        B = np.array([r['post']['b'][0] for r in rows]); Be = np.array([r['post']['b'][1] for r in rows])
        K = np.array([r['map']['k'] for r in rows])
        good = Te < 5.0      # drop unconstrained bins (e.g. SOSS order 2) from the combinations
        fs = prof_interp(rows)
        out = dict(bins=[dict(det=d, lam=float(l), lo=r['lo'], hi=r['hi'], T14=[float(a), float(e)], b=[float(bb), float(be)], k=r['post']['k'][:2], u1=r['post']['u1'][:2], u2=r['post']['u2'][:2],
                              sig_ppm=r['sig_ppm'], beta=r['beta'], n_over_tau=r['n_over_tau'], used=bool(g)) for d, l, r, a, e, bb, be, g in zip(det_of, lam, rows, T, Te, B, Be, good)])
        ig = np.where(good)[0]
        b_all, sb_all, c_sh, c_free = shared_b(fs, list(ig))
        out['shared_all'] = dict(b=[b_all, sb_all], T14=float(t14(b_all, np.mean(K[ig]), pl)), dchi2_free_vs_shared=c_sh - c_free, dof=len(ig) - 1, p=float(CHI2.sf(c_sh - c_free, len(ig) - 1)))
        out['per_det'] = {}
        for det in sorted(set(det_of[ig])):
            s = [i for i in ig if det_of[i] == det]
            bd, sbd, cd, cfd = shared_b(fs, s)
            Td = t14(bd, np.mean(K[s]), pl); dT = abs(t14(bd + sbd, np.mean(K[s]), pl) - Td)
            wm = wmean(T[s], Te[s])
            p, C, c2l = wline(lam[s] - lam[s].mean(), T[s], Te[s])
            out['per_det'][det] = dict(b_shared=[bd, sbd], T14_shared=[float(Td), float(dT)], dchi2_within=cd - cfd, dof_within=len(s) - 1, p_within=float(CHI2.sf(cd - cfd, len(s) - 1)),
                                       T14_wmean=[wm[0], wm[1]], chi2_wmean=wm[2], slope_T14_per_um=[float(p[1]), float(np.sqrt(C[1, 1]))], lam_mean=float(lam[s].mean()))
        dd = sorted(out['per_det'])
        if len(dd) == 2:
            a, b2 = out['per_det'][dd[0]], out['per_det'][dd[1]]
            out['det_diff_T14_shared'] = [a['T14_shared'][0] - b2['T14_shared'][0], float(np.hypot(a['T14_shared'][1], b2['T14_shared'][1]))]
            out['det_diff_T14_wmean'] = [a['T14_wmean'][0] - b2['T14_wmean'][0], float(np.hypot(a['T14_wmean'][1], b2['T14_wmean'][1]))]
            # step vs smooth: b = per-detector constants (2 par) vs b = b0 + b1 (lam - 3.77) (2 par), vs constant (1 par)
            s_ = list(ig); fss = [fs[i] for i in s_]; lm = lam[s_]; isd2 = det_of[s_] == dd[1]
            c_const = shared_b(fs, s_)[2]
            c_step = shared_b(fs, [i for i in s_ if det_of[i] == dd[0]])[2] + shared_b(fs, [i for i in s_ if det_of[i] == dd[1]])[2]
            lin = minimize(lambda q: chi_of_bfunc(fss, q[0] + q[1] * (lm - 3.77)), [b_all, 0.0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-6, maxiter=4000))
            both = minimize(lambda q: chi_of_bfunc(fss, q[0] + q[1] * (lm - 3.77) + q[2] * isd2), [b_all, 0.0, 0.0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-6, maxiter=6000))
            out['model_compare'] = dict(chi2_const=c_const, chi2_step=float(c_step), chi2_linear=float(lin.fun), chi2_step_plus_linear=float(both.fun),
                                        dchi2_step_vs_const=c_const - c_step, dchi2_linear_vs_const=c_const - float(lin.fun), linear_db_dum=float(lin.x[1]),
                                        dchi2_step_given_linear=float(lin.fun - both.fun), dchi2_linear_given_step=float(c_step - both.fun), step_db_given_linear=float(both.x[2]))
            # adjacent bins across the gap, and within-detector slopes extrapolated to the gap
            i1 = max([i for i in ig if det_of[i] == dd[0]], key=lambda i: lam[i]); i2 = min([i for i in ig if det_of[i] == dd[1]], key=lambda i: lam[i])
            out['adjacent'] = dict(bins=[rows[i1]['bin'], rows[i2]['bin']], dT14=[float(T[i1] - T[i2]), float(np.hypot(Te[i1], Te[i2]))])
            gap = 0.5 * (rows[i1]['hi'] + rows[i2]['lo'])
            ex = {}
            for det in dd:
                s = [i for i in ig if det_of[i] == det]; p, C, _ = wline(lam[s] - gap, T[s], Te[s]); ex[det] = [float(p[0]), float(np.sqrt(C[0, 0]))]
            out['extrapolated_to_gap'] = dict(lam_gap=gap, T14=ex, diff=[ex[dd[0]][0] - ex[dd[1]][0], float(np.hypot(ex[dd[0]][1], ex[dd[1]][1]))])
        res['%s_%s' % (vis, ldm)] = out
    # between NIRSpec visits of d, bin by bin
    a, b = res.get('d1224_%s' % ldm), res.get('d4098_%s' % ldm)
    if a and b:
        ba = {(x['det'], x['lo']): x for x in a['bins']}; bb = {(x['det'], x['lo']): x for x in b['bins']}
        rows = []
        for key in sorted(ba):
            if key in bb: rows.append(dict(det=key[0], lam=ba[key]['lam'], dT14=[bb[key]['T14'][0] - ba[key]['T14'][0], float(np.hypot(bb[key]['T14'][1], ba[key]['T14'][1]))], db=[bb[key]['b'][0] - ba[key]['b'][0], float(np.hypot(bb[key]['b'][1], ba[key]['b'][1]))]))
        D = dict(bins=rows, per_det={})
        for det in ('NRS1', 'NRS2'):
            s = [r for r in rows if r['det'] == det]; v = np.array([r['dT14'][0] for r in s]); e = np.array([r['dT14'][1] for r in s]); l = np.array([r['lam'] for r in s])
            m, sm, c2 = wmean(v, e); p, C, c2l = wline(l - l.mean(), v, e)
            sh = [a['per_det'][det]['T14_shared'], b['per_det'][det]['T14_shared']]
            D['per_det'][det] = dict(dT14_wmean=[m, sm], chi2=c2, dof=len(s) - 1, p=float(CHI2.sf(c2, len(s) - 1)), slope_per_um=[float(p[1]), float(np.sqrt(C[1, 1]))],
                                     dT14_shared=[sh[1][0] - sh[0][0], float(np.hypot(sh[0][1], sh[1][1]))])
        v = np.array([r['dT14'][0] for r in rows]); e = np.array([r['dT14'][1] for r in rows]); l = np.array([r['lam'] for r in rows]); isd2 = np.array([r['det'] == 'NRS2' for r in rows])
        m, sm, c2c = wmean(v, e); p, C, c2l = wline(l - 3.77, v, e)
        w = 1 / e ** 2; A = np.vstack([np.ones_like(l), isd2.astype(float)]).T; Cs = np.linalg.inv(A.T @ (A * w[:, None])); ps = Cs @ (A.T @ (w * v)); c2s = float(np.sum(w * (v - A @ ps) ** 2))
        A3 = np.vstack([np.ones_like(l), l - 3.77, isd2.astype(float)]).T; C3 = np.linalg.inv(A3.T @ (A3 * w[:, None])); p3 = C3 @ (A3.T @ (w * v)); c23 = float(np.sum(w * (v - A3 @ p3) ** 2))
        D['all'] = dict(dT14_wmean=[m, sm], chi2_const=c2c, dof=len(v) - 1, p_const=float(CHI2.sf(c2c, len(v) - 1)), chi2_linear=c2l, slope_per_um=[float(p[1]), float(np.sqrt(C[1, 1]))],
                        chi2_step=c2s, step_NRS2_minus_NRS1=[float(ps[1]), float(np.sqrt(Cs[1, 1]))], chi2_step_plus_linear=c23, step_given_linear=[float(p3[2]), float(np.sqrt(C3[2, 2]))], slope_given_step=[float(p3[1]), float(np.sqrt(C3[1, 1]))])
        res['between_%s' % ldm] = D
json.dump(res, open(outp, 'w'), indent=1)
# printout
for k, v in res.items():
    if k.startswith('between'):
        print('\n==', k)
        for det, d in v['per_det'].items(): print('  %s  dT14 wmean %.2f +- %.2f (chi2 %.1f/%d, p %.2f)  shared-b %.2f +- %.2f  slope %.2f +- %.2f min/um' % (det, *d['dT14_wmean'], d['chi2'], d['dof'], d['p'], *d['dT14_shared'], *d['slope_per_um']))
        a = v['all']; print('  all: wmean %.2f +- %.2f chi2_const %.1f/%d p %.3f | linear chi2 %.1f slope %.2f +- %.2f | step chi2 %.1f NRS2-NRS1 %.2f +- %.2f | step|linear %.2f +- %.2f' % (*a['dT14_wmean'], a['chi2_const'], a['dof'], a['p_const'], a['chi2_linear'], *a['slope_per_um'], a['chi2_step'], *a['step_NRS2_minus_NRS1'], *a['step_given_linear']))
        continue
    print('\n==', k)
    for x in v['bins']: print('  %s %.3f-%.3f  T14 %6.2f +- %.2f  b %.4f +- %.4f  k %.4f  u1 %.2f u2 %.2f  sig %4.0f beta %.2f %s' % (x['det'], x['lo'], x['hi'], *x['T14'], *x['b'], x['k'][0], x['u1'][0], x['u2'][0], x['sig_ppm'], x['beta'], '' if x['used'] else 'UNUSED'))
    s = v['shared_all']; print('  shared b (all bins) %.4f +- %.4f  T14 %.2f  dchi2(free b per bin) %.1f for %d dof, p %.3f' % (*s['b'], s['T14'], s['dchi2_free_vs_shared'], s['dof'], s['p']))
    for det, d in v['per_det'].items(): print('  %s shared b %.4f +- %.4f T14 %.2f +- %.2f | wmean %.2f +- %.2f (chi2 %.1f/%d) | within dchi2 %.1f/%d p %.2f | slope %.2f +- %.2f min/um' % (det, *d['b_shared'], *d['T14_shared'], *d['T14_wmean'], d['chi2_wmean'], d['dof_within'], d['dchi2_within'], d['dof_within'], d['p_within'], *d['slope_T14_per_um']))
    if 'det_diff_T14_shared' in v:
        m = v['model_compare']; print('  det diff (first - second) shared %.2f +- %.2f, wmean %.2f +- %.2f' % (*v['det_diff_T14_shared'], *v['det_diff_T14_wmean']))
        print('  chi2 const %.1f step %.1f linear %.1f step+lin %.1f | step vs const %.1f, linear vs const %.1f, step|linear %.1f, linear|step %.1f' % (m['chi2_const'], m['chi2_step'], m['chi2_linear'], m['chi2_step_plus_linear'], m['dchi2_step_vs_const'], m['dchi2_linear_vs_const'], m['dchi2_step_given_linear'], m['dchi2_linear_given_step']))
        print('  adjacent %s vs %s: dT14 %.2f +- %.2f | extrapolated to %.2f um: diff %.2f +- %.2f' % (*v['adjacent']['bins'], *v['adjacent']['dT14'], v['extrapolated_to_gap']['lam_gap'], *v['extrapolated_to_gap']['diff']))
