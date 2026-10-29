"""TOI-421 c and b: recomputes the predictions of section 6.12 (the frozen file data/predictions_94986319.json was written
by this computation on 2026-10-29 22:54 UTC) and compares them with the frozen values.
c: linear trend from the per-transit fit (baseline.json A rows; a/R* 31.91, k 0.05459, u 0.401/0.222, circular);
b: early/late epoch fold at fixed k 0.0279, a/R* 13.76 (planet_b.json) -> linear in time; 60-min box depth from the model.
Constant alternatives: c b = 2018-20 mean; b-planet b = 2018-20 mean.
usage: python scripts/toi421/predict_toi421.py   -> work/toi421/out/predictions_recomputed.json"""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import fit_toi421 as adv, batman
O = os.path.join(adv.WORK, 'out')
A = json.load(open(os.path.join(O, 'baseline.json')))['A']
x = np.array([r['year'] for r in A]); y = np.array([r['b'] for r in A]); e = np.array([r['eb'] for r in A]); w = 1 / e ** 2
M = np.vstack([np.ones_like(x), x - 2022]).T; C = np.linalg.inv(M.T @ (M * w[:, None])); p = C @ M.T @ (w * y)
pb = json.load(open(os.path.join(O, 'planet_b.json'))); be, bl = pb['batman']['early'], pb['batman']['late']
sl_b = (bl['x'][0] - be['x'][0]) / (bl['year'] - be['year']); esl_b = np.hypot(bl['ex'][0], be['ex'][0]) / (bl['year'] - be['year'])
def T14(b, aR, k, P): return P / np.pi * np.arcsin(np.sqrt(max((1 + k) ** 2 - b * b, 0)) / np.sqrt(aR * aR - b * b)) * 1440
def boxdepth(b, k=0.0279, aR=13.76, P=5.19759807):
    t = np.linspace(-30, 30, 601) / 1440.; q = batman.TransitParams(); q.t0 = 0; q.per = P; q.rp = k; q.a = aR
    q.inc = np.degrees(np.arccos(min(b / aR, 1))); q.ecc = 0; q.w = 90; q.u = [0.401, 0.222]; q.limb_dark = 'quadratic'
    return float((1 - batman.TransitModel(q, t).light_curve(q)).mean() * 1e6)
out = dict(meta=__doc__, c_fit=dict(b2022=float(p[0]), slope=float(p[1]), cov=C.tolist()),
           b_fit=dict(b_early=be['x'][0], year_early=be['year'], slope=sl_b, eslope=esl_b), rows=[])
for yr in (2026.79, 2026.87, 2027.0, 2027.5, 2028.0, 2029.0, 2030.0):
    v = np.array([1, yr - 2022]); bc = float(v @ p); ebc = float(np.sqrt(v @ C @ v))
    bbp = be['x'][0] + sl_b * (yr - be['year']); ebbp = float(np.hypot(be['ex'][0], esl_b * (yr - be['year'])))
    out['rows'].append(dict(year=yr, c_b=bc, c_eb=ebc, c_T14=T14(bc, 31.91, 0.05459, 16.06750764),
                            c_T14_const_early=T14(0.7654, 31.91, 0.05459, 16.06750764),
                            b_b=bbp, b_eb=ebbp, b_T14=T14(bbp, 13.76, 0.0279, 5.19759807), b_box60_depth=boxdepth(bbp),
                            b_box60_depth_const=boxdepth(be['x'][0]), b_T14_const=T14(be['x'][0], 13.76, 0.0279, 5.19759807)))
    r = out['rows'][-1]; print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
# planet-b box depth calibration: model vs measured early/late
print('box model early/late', boxdepth(be['x'][0]), boxdepth(bl['x'][0]), 'measured', pb['box']['early']['x'][0] * 1e6, pb['box']['late']['x'][0] * 1e6)
# next c transits (linear ephemeris from the TTV test; TTV rms ~3 min)
tt = json.load(open(os.path.join(O, 'ttv.json')))['linear_ephem']; T0c, Pc = tt['T0'], tt['P']
nxt = []
for n in range(int((4400 - T0c) / Pc), int((4560 - T0c) / Pc) + 1):
    tc = T0c + n * Pc; nxt.append(dict(n=n, BTJD=round(tc, 4), BJD_TDB=round(tc + 2457000, 4), year=round(adv.year(tc), 3)))
out['next_c_transits'] = nxt; print(nxt)
json.dump(out, open(os.path.join(O, 'predictions_recomputed.json'), 'w'), indent=1)
frozen = json.load(open(os.path.join(adv.ROOT, 'data', 'predictions_94986319.json')))
dmax = max(abs(a[k] - b[k]) for a, b in zip(out['rows'], frozen['rows']) for k in a if k != 'year')
dt = max(abs(a['BTJD'] - b['BTJD']) for a, b in zip(out['next_c_transits'], frozen['next_c_transits']))
print('largest difference from the frozen predictions: %.3g (rows), %.3g d (transit times)' % (dmax, dt))
