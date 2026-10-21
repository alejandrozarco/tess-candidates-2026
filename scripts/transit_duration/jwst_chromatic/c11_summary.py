"""C11: JWST duration summary with an empirical detector-systematics floor. Own code.
Inputs: out/chromatic.json (c6), out/sysvariants.json (c9). sigma_sys per detector-visit is estimated from the
within-visit NRS1 - NRS2 differences (shared-b T14, Kurucz-prior LD): sigma_sys = sqrt(mean(d_v^2 - s_v^2) / 2),
assuming both detectors carry independent systematics of equal size. Slopes are weighted least squares in time.
Usage: c11_summary.py <out.json>"""
import json, os, sys, numpy as np
H = os.environ.get('TDC_WORK', 'work/transit_duration_chromatic')
R = json.load(open(H + '/out/chromatic.json')); S = json.load(open(H + '/out/sysvariants.json'))
YR = {'1201': 2022.992, '1224': 2023.481, '4098': 2024.195}; TESS = (-1.185, 0.112)
out = {}
for ldm in ('prior', 'free'):
    T = {(v, d): R['d%s_%s' % (v, ldm)]['per_det'][d]['T14_shared'] for v in ('1224', '4098') for d in ('NRS1', 'NRS2')}
    T[('1201', 'SOSS1')] = R['d1201_%s' % ldm]['per_det']['SOSS1']['T14_shared']
    dd = [T[(v, 'NRS1')][0] - T[(v, 'NRS2')][0] for v in ('1224', '4098')]; ss = [np.hypot(T[(v, 'NRS1')][1], T[(v, 'NRS2')][1]) for v in ('1224', '4098')]
    ssys = float(np.sqrt(max(0, np.mean(np.array(dd) ** 2 - np.array(ss) ** 2)) / 2))
    o = dict(T14={'%s_%s' % k: v for k, v in T.items()}, det_diff=dd, det_diff_err=ss, sigma_sys_per_detector=ssys,
             interaction=[dd[1] - dd[0], float(np.hypot(*ss))])
    def visit_avg(v):
        a, b = T[(v, 'NRS1')], T[(v, 'NRS2')]; m = 0.5 * (a[0] + b[0]); st = 0.5 * np.hypot(a[1], b[1]); return m, st, float(np.hypot(st, ssys / np.sqrt(2)))
    va = {v: visit_avg(v) for v in ('1224', '4098')}
    o['visit_avg'] = va
    o['pair_dT14_avg'] = [va['4098'][0] - va['1224'][0], float(np.hypot(va['4098'][1], va['1224'][1])), float(np.hypot(va['4098'][2], va['1224'][2]))]
    def slope(pts):
        x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts]); e = np.array([p[2] for p in pts]); w = 1 / e ** 2
        A = np.vstack([np.ones_like(x), x - 2023.5]).T; C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ (A.T @ (w * y))
        chi_tess = None
        return [float(p[1]), float(np.sqrt(C[1, 1]))]
    sl = {}
    for lab, sysf in (('formal', 0.0), ('with sigma_sys', ssys)):
        pts_avg = [(YR['1224'], va['1224'][0], va['1224'][1] if not sysf else va['1224'][2]), (YR['4098'], va['4098'][0], va['4098'][1] if not sysf else va['4098'][2])]
        ni = T[('1201', 'SOSS1')]; pts3 = pts_avg + [(YR['1201'], ni[0], float(np.hypot(ni[1], sysf)))]
        sl[lab] = {'NIRSpec pair, detector average': slope(pts_avg), 'NIRISS 1201 + NIRSpec pair (avg)': slope(pts3)}
        for d in ('NRS1', 'NRS2'):
            sl[lab]['NIRSpec pair, %s only' % d] = slope([(YR[v], T[(v, d)][0], float(np.hypot(T[(v, d)][1], sysf))) for v in ('1224', '4098')])
            sl[lab]['NIRISS 1201 + %s pair' % d] = slope([(YR[v], T[(v, d)][0], float(np.hypot(T[(v, d)][1], sysf))) for v in ('1224', '4098')] + [(YR['1201'], ni[0], float(np.hypot(ni[1], sysf)))])
    o['slopes'] = sl
    out[ldm] = o
    print('== LD', ldm)
    print('  per-detector T14 (shared b): ' + ', '.join('%s %s %.2f +- %.2f' % (k[0], k[1], *v) for k, v in T.items()))
    print('  NRS1 - NRS2: 1224 %.2f +- %.2f, 4098 %.2f +- %.2f; change of the difference %.2f +- %.2f (%.1f sigma)' % (dd[0], ss[0], dd[1], ss[1], o['interaction'][0], o['interaction'][1], abs(o['interaction'][0] / o['interaction'][1])))
    print('  empirical sigma_sys per detector-visit %.2f min' % ssys)
    print('  detector-average T14: 1224 %.2f +- %.2f (stat) +- %.2f (with sys); 4098 %.2f +- %.2f / %.2f; pair dT14 %.2f +- %.2f (stat) / %.2f (with sys); TESS predicts -0.85 +- 0.08' % (*va['1224'], *va['4098'], *o['pair_dT14_avg']))
    for lab, d in sl.items():
        for k, v in d.items(): print('  slope [%s] %-34s %.2f +- %.2f min/yr  (TESS %.2f; diff %.1f sigma)' % (lab, k, v[0], v[1], TESS[0], (v[0] - TESS[0]) / np.hypot(v[1], TESS[1])))
# own vs pipeline NRS2 (c9; own: systematics model of section 6.6.10; pipeline: quadratic only, no regressors available)
for v in ('1224', '4098'):
    a = S[v]['NRS2']['poly2+regs (section 6.6.10)']['T14']; b = S[v]['NRS2pipe']['poly2']['T14']; c = S[v]['NRS1']['poly2+regs (section 6.6.10)']['T14']
    print('visit %s white light (Kurucz-prior LD, fixed a/R*): own NRS1 %.2f, own NRS2 %.2f, pipeline NRS2 %.2f (own - pipeline NRS2 %.2f)' % (v, c[0], a[0], b[0], a[0] - b[0]))
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=lambda x: list(x) if isinstance(x, tuple) else float(x))
