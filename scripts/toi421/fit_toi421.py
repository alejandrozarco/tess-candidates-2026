"""TOI-421 c (TIC 94986319; P 16.0675 d): impact parameter and transit duration per transit, per sector and per epoch
group, with the variants of README section 6.12.3-6.12.4.
Conventions unless a test says otherwise: circular, quadratic LD, model integrated over the exposure,
template k 0.05459, a/R* 31.910, u 0.401/0.222, P 16.06750764, T0 BTJD 2195.30635; planet-b transits (P 5.19759807,
T0 1441.283952) masked if within the window +- 0.06 d. Windows +-0.20 d around each predicted c transit, QUALITY == 0.
Per-transit tc and quadratic baseline always free.
usage: python scripts/toi421/fit_toi421.py <test>   (tests: baseline aR ld dilution duration cadence spots ttv inject
       inject_trend loso); reads work/toi421/data/*.npz (download_toi421.py), writes work/toi421/out/<test>.json
"""
import sys, json, os, glob
import numpy as np, batman
from scipy.optimize import least_squares

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.path.join(ROOT, 'work', 'toi421')
P, T0, K0, AR0, U0 = 16.06750764, 2195.30635, 0.05459, 31.910, (0.401, 0.222)
PB, T0B, HWB = 5.19759807, 1441.283952, 0.035 + 0.06
HALF = 0.20


def year(btjd): return 2014.0 + (btjd + 2457000.0 - 2456658.5) / 365.25


def load(tag, flux='pdc', crowd_override=None, bin_to=None):
    d = np.load(os.path.join(WORK, 'data', tag + '.npz'))
    # native-endian float64 copies: FITS arrays are big-endian and batman's C code reads the raw buffer
    t, q = np.ascontiguousarray(d['t'], dtype=np.float64), d['q']
    f, e = np.ascontiguousarray(d[flux], dtype=np.float64), np.ascontiguousarray(d[flux + 'e'], dtype=np.float64)
    ok = (q == 0) & np.isfinite(f) & np.isfinite(e) & np.isfinite(t)
    t, f, e = t[ok], f[ok], e[ok]
    med = np.nanmedian(f)
    if flux == 'sap':        # own crowding/flux-fraction correction (SPOC convention)
        c = float(d['crowd']) if crowd_override is None else crowd_override
        f = (f - (1 - c) * med) / float(d['flfr']); e = e / float(d['flfr'])
    elif crowd_override is not None:   # re-dilute PDCSAP to an alternative crowding value
        c0 = float(d['crowd']); m2 = np.nanmedian(f)
        # PDC flux = (raw - (1-c0) m)/flfr ; alternative: (raw - (1-c) m)/flfr  -> shift by (c0-c) m_raw/flfr ~ (c0-c)*m2/c0
        f = f + (c0 - crowd_override) * m2 / c0
    m = np.nanmedian(f); f, e = f / m, e / m
    exp = float(d['exp']) / 86400.
    if bin_to:
        nb = int(round(bin_to / (exp * 86400.)))
        if nb > 1:
            edges = np.arange(t.min(), t.max() + bin_to / 86400., bin_to / 86400.)
            idx = np.digitize(t, edges)
            tt, ff, ee = [], [], []
            for i in np.unique(idx):
                s = idx == i
                if s.sum() >= 0.8 * nb:
                    tt.append(t[s].mean()); ff.append(f[s].mean()); ee.append(np.sqrt(np.sum(e[s] ** 2)) / s.sum())
            t, f, e = np.array(tt), np.array(ff), np.array(ee); exp = bin_to / 86400.
    return t, f, e, exp


def windows_from(t, f, e, exp, sec, inj=None, offset=0.0):
    """inj: dict(b=..,aR=..,k=..) -> multiply a synthetic c transit at tc+offset into windows centred there
    (only windows that do not contain a real c or b transit)."""
    out = []
    for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2):
        tc = T0 + n * P + offset; m = np.abs(t - tc) < HALF
        if m.sum() < 0.8 * 2 * HALF / exp: continue
        if (np.abs(t[m] - tc) < 0.06).sum() < 0.8 * 0.12 / exp: continue
        nb = np.round((tc - T0B) / PB)
        dtb = tc - (T0B + nb * PB)
        if abs(dtb) < 0.058 + 0.035: continue          # planet-b transit overlaps c's transit itself -> skip window
        if abs(dtb) < HALF + 0.035:                     # planet-b transit in the baseline -> mask its points (+-0.035 d)
            m = m & (np.abs(t - (T0B + nb * PB)) > 0.035)
        if inj is not None:
            # skip if a real c transit is within the window
            nc = np.round((tc - T0) / P)
            if abs(tc - (T0 + nc * P)) < HALF + 0.08: continue
        tt, ff, ee = t[m], f[m].copy(), e[m]
        if inj is not None:
            ff = ff * model(tt, tc, inj['b'], inj.get('aR', AR0), inj.get('k', K0), exp)
        out.append(dict(sec=sec, n=n, tc0=tc, t=tt, f=ff, e=ee, exp=exp, year=year(tc)))
    return out


def model(t, tc, b, aR, k, exp, u=U0):
    ss = max(3, int(round(exp / (30 / 86400.))) | 1) if exp > 40 / 86400. else 1
    p = batman.TransitParams(); p.t0 = tc; p.per = P; p.rp = k; p.a = aR
    p.inc = np.degrees(np.arccos(min(b / aR, 1.0))); p.ecc = 0.; p.w = 90.; p.u = list(u); p.limb_dark = 'quadratic'
    if ss > 1: return batman.TransitModel(p, t, supersample_factor=ss, exp_time=exp).light_curve(p)
    return batman.TransitModel(p, t).light_curve(p)


def trap(t, tc, depth, T14, T12, exp):
    """trapezoid (duration-only, convention-free), integrated over exposure by 7-point supersampling"""
    ss = 7 if exp > 40 / 86400. else 1
    off = (np.arange(ss) - (ss - 1) / 2) * exp / ss if ss > 1 else np.zeros(1)
    tt = (t[:, None] + off[None, :]) - tc; x = np.abs(tt)
    T12 = min(T12, T14 / 2 * 0.999)
    y = np.where(x <= T14 / 2 - T12, depth, np.where(x < T14 / 2, depth * (T14 / 2 - x) / T12, 0.0))
    return 1 - y.mean(axis=1)


def fit(ws, b0=0.75, aR=None, k=None, u=U0, free_aR=False, free_k=False, free_u=False, prior_aR=None, clip=None):
    """shared b (+aR, +k, +u) over the windows ws; per-window tc + quadratic baseline."""
    aR = AR0 if aR is None else aR; k = K0 if k is None else k
    gnames = ['b'] + (['aR'] if free_aR else []) + (['k'] if free_k else []) + (['q1', 'q2'] if free_u else [])
    x0 = [b0] + ([aR] if free_aR else []) + ([k] if free_k else [])
    if free_u:
        q1 = (u[0] + u[1]) ** 2; q2 = u[0] / (2 * (u[0] + u[1])); x0 += [q1, q2]
    ng = len(x0)
    for w in ws: x0 += [w['tc0'], 1.0, 0.0, 0.0]
    masks = [np.ones(len(w['t']), bool) for w in ws]

    def unpack(x):
        g = dict(zip(gnames, x[:ng]))
        if free_u:
            sq = np.sqrt(g['q1']); uu = (2 * sq * g['q2'], sq * (1 - 2 * g['q2']))
        else: uu = u
        return g['b'], g.get('aR', aR), g.get('k', k), uu

    def resid(x):
        b, a, kk, uu = unpack(x); r = []
        for i, w in enumerate(ws):
            tc, c0, c1, c2 = x[ng + 4 * i: ng + 4 * i + 4]; dt = w['t'] - w['tc0']
            rr = (w['f'] - model(w['t'], tc, b, a, kk, w['exp'], uu) * (c0 + c1 * dt + c2 * dt * dt)) / w['e']
            r.append(rr[masks[i]])
        if prior_aR is not None and free_aR: r.append(np.array([(a - prior_aR[0]) / prior_aR[1]]))
        return np.concatenate(r)
    lo = [0.0] + ([8.] if free_aR else []) + ([0.02] if free_k else []) + ([0., 0.] if free_u else [])
    hi = [1.1] + ([80.] if free_aR else []) + ([0.1] if free_k else []) + ([1., 1.] if free_u else [])
    for w in ws: lo += [w['tc0'] - 0.04, 0.9, -1, -10]; hi += [w['tc0'] + 0.04, 1.1, 1, 10]
    x0 = np.clip(x0, np.array(lo) + 1e-9, np.array(hi) - 1e-9)
    s = least_squares(resid, x0, bounds=(lo, hi), x_scale='jac')
    if clip:
        for it in range(2):
            b, a, kk, uu = unpack(s.x)
            for i, w in enumerate(ws):
                tc, c0, c1, c2 = s.x[ng + 4 * i: ng + 4 * i + 4]; dt = w['t'] - w['tc0']
                rr = (w['f'] - model(w['t'], tc, b, a, kk, w['exp'], uu) * (c0 + c1 * dt + c2 * dt * dt)) / w['e']
                sd = 1.4826 * np.median(np.abs(rr - np.median(rr)))
                masks[i] = np.abs(rr - np.median(rr)) < clip * sd
            s = least_squares(resid, s.x, bounds=(lo, hi), x_scale='jac')
    npts = s.fun.size; dof = npts - len(x0); chi2r = float(np.sum(s.fun ** 2) / dof)
    try: C = np.linalg.inv(s.jac.T @ s.jac) * max(chi2r, 1.0)
    except np.linalg.LinAlgError: C = np.full((len(x0), len(x0)), np.nan)
    e = np.sqrt(np.abs(np.diag(C)))
    b, a, kk, uu = unpack(s.x)
    T14 = P / np.pi * np.arcsin(np.sqrt(max((1 + kk) ** 2 - b * b, 0)) / np.sqrt(a * a - b * b)) * 1440
    T23 = P / np.pi * np.arcsin(np.sqrt(max((1 - kk) ** 2 - b * b, 0)) / np.sqrt(a * a - b * b)) * 1440
    r = dict(b=float(b), eb=float(e[0]), chi2r=chi2r, ntr=len(ws), year=float(np.mean([w['year'] for w in ws])),
             T14=float(T14), T23=float(T23), tcs=[float(v) for v in s.x[ng::4]], etcs=[float(v) for v in e[ng::4]],
             secs=sorted(set(int(w['sec']) for w in ws)))
    gi = 1
    if free_aR: r.update(aR=float(a), eaR=float(e[gi]), corr_b_aR=float(C[0, gi] / (e[0] * e[gi]))); gi += 1
    if free_k: r.update(k=float(kk), ek=float(e[gi])); gi += 1
    if free_u: r.update(u1=float(uu[0]), u2=float(uu[1]))
    return r


def fit_trap(ws, T14=0.118, T12=0.012, depth=0.0034):
    """duration-only: shared depth, T14, T12 over ws; per-window tc + quadratic baseline. Returns minutes."""
    ng = 3; x0 = [depth, T14, T12]
    for w in ws: x0 += [w['tc0'], 1.0, 0.0, 0.0]

    def resid(x):
        dp, a, c = x[:3]; r = []
        for i, w in enumerate(ws):
            tc, c0, c1, c2 = x[ng + 4 * i: ng + 4 * i + 4]; dt = w['t'] - w['tc0']
            r.append((w['f'] - trap(w['t'], tc, dp, a, c, w['exp']) * (c0 + c1 * dt + c2 * dt * dt)) / w['e'])
        return np.concatenate(r)
    lo = [0.0005, 0.06, 0.002]; hi = [0.01, 0.2, 0.05]
    for w in ws: lo += [w['tc0'] - 0.04, 0.9, -1, -10]; hi += [w['tc0'] + 0.04, 1.1, 1, 10]
    s = least_squares(resid, x0, bounds=(lo, hi), x_scale='jac')
    dof = s.fun.size - len(x0); chi2r = float(np.sum(s.fun ** 2) / dof)
    C = np.linalg.inv(s.jac.T @ s.jac) * max(chi2r, 1.0); e = np.sqrt(np.abs(np.diag(C)))
    dp, a, c = s.x[:3]
    cT = C[1, 2]
    T23 = a - 2 * c; eT23 = np.sqrt(C[1, 1] + 4 * C[2, 2] - 4 * cT)
    return dict(depth_ppm=float(dp * 1e6), edepth=float(e[0] * 1e6), T14=float(a * 1440), eT14=float(e[1] * 1440),
                T12=float(c * 1440), eT12=float(e[2] * 1440), T23=float(T23 * 1440), eT23=float(eT23 * 1440),
                chi2r=chi2r, ntr=len(ws), year=float(np.mean([w['year'] for w in ws])),
                secs=sorted(set(int(w['sec']) for w in ws)))


def slope(rows, key='b', ekey='eb'):
    y = np.array([r[key] for r in rows]); e = np.array([r[ekey] for r in rows]); x = np.array([r['year'] for r in rows])
    w = 1 / e ** 2; A = np.vstack([np.ones_like(x), x - 2022]).T
    C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ A.T @ (w * y)
    chi2 = float(np.sum(w * (y - A @ p) ** 2)); dof = max(len(y) - 2, 1)
    es = float(np.sqrt(C[1, 1] * max(1, chi2 / dof)))
    early = x < 2023; late = ~early
    def wm(s): ww = w[s]; return float(np.sum(ww * y[s]) / ww.sum()), float(1 / np.sqrt(ww.sum()))
    out = dict(slope=float(p[1]), eslope=float(np.sqrt(C[1, 1])), eslope_scaled=es, chi2=chi2, dof=dof, z=float(p[1] / es))
    if early.any() and late.any():
        (m1, e1), (m2, e2) = wm(early), wm(late)
        out.update(early=m1, e_early=e1, late=m2, e_late=e2, step=m2 - m1, estep=float(np.hypot(e1, e2)),
                   zstep=float((m2 - m1) / np.hypot(e1, e2)))
    return out


SECS2M = ['s005_2m', 's006_2m', 's032_2m', 's098_2m', 's107_2m']


def all_windows(flux='pdc', tags=SECS2M, crowd=None, bin_to=None, inj=None, offset=0.0):
    ws = []
    for tag in tags:
        sec = int(tag[1:4]); c = crowd.get(sec) if isinstance(crowd, dict) else None
        t, f, e, exp = load(tag, flux, c, bin_to)
        ws += windows_from(t, f, e, exp, sec, inj, offset)
    return ws


def per_transit(ws, **kw):
    rows = []
    for w in ws:
        r = fit([w], **kw); r.update(sec=w['sec'], n=w['n']); rows.append(r)
    return rows


def per_sector(ws, **kw):
    rows = []
    for s in sorted(set(w['sec'] for w in ws)):
        g = [w for w in ws if w['sec'] == s]; r = fit(g, **kw); r['sec'] = s; rows.append(r)
    return rows


def short(rows, keys=('sec', 'n', 'year', 'b', 'eb', 'T14')):
    return [{k: (round(r[k], 4) if isinstance(r.get(k), float) else r.get(k)) for k in keys if k in r} for r in rows]


def T_baseline():
    ws = all_windows(); A = per_transit(ws); B = per_sector(ws)
    D = per_sector(ws, free_k=True)
    return dict(windows=[(w['sec'], w['n']) for w in ws], A=short(A), A_slope=slope(A), B=short(B), B_slope=slope(B),
                D=short(D, ('sec', 'b', 'eb', 'k', 'ek')), D_slope=slope(D), D_k_slope=slope(D, 'k', 'ek'))


def T_aR():
    ws = all_windows(); res = {}
    for aR in (27.0, 29.2, 31.91, 34.0, 36.0):
        A = per_transit(ws, aR=aR); res['fixed_aR_%.2f' % aR] = dict(slope=slope(A), T14_first=A[0]['T14'])
    # global a/R* free (shared by all transits), b per epoch: fit early and late jointly with shared aR
    early = [w for w in ws if w['year'] < 2023]; late = [w for w in ws if w['year'] > 2023]
    g = fit(ws, free_aR=True); res['global_free_aR_shared_b'] = {k: g[k] for k in ('b', 'eb', 'aR', 'eaR', 'chi2r')}
    # two-epoch joint with shared free aR: implement via scan over aR (profile): for each aR, fit early & late b
    prof = []
    for aR in np.arange(24.0, 40.01, 1.0):
        e1 = fit(early, aR=aR); l1 = fit(late, aR=aR)
        chi = (e1['chi2r'] * 1.0)  # placeholder (chi2r not additive) -> recompute total chi2 below
        prof.append(dict(aR=float(aR), b_early=e1['b'], eb_early=e1['eb'], b_late=l1['b'], eb_late=l1['eb'],
                         db=l1['b'] - e1['b'], z=(l1['b'] - e1['b']) / np.hypot(e1['eb'], l1['eb'])))
    res['profile_two_epoch'] = prof
    # per-sector free aR with density prior: rho* 0.85/0.87^3 -> circular a/R* 29.2; eccentric c (e~0.13-0.19) widens it
    for pr in ((29.2, 1.5), (31.9, 3.0)):
        B = per_sector(ws, free_aR=True, prior_aR=pr)
        res['per_sector_free_aR_prior_%.1f_%.1f' % pr] = dict(rows=short(B, ('sec', 'b', 'eb', 'aR', 'eaR')), slope=slope(B),
                                                          aR_slope=slope(B, 'aR', 'eaR'))
    B = per_sector(ws, free_aR=True)
    res['per_sector_free_aR_noprior'] = dict(rows=short(B, ('sec', 'b', 'eb', 'aR', 'eaR', 'corr_b_aR')), slope=slope(B),
                                             aR_slope=slope(B, 'aR', 'eaR'), T14_slope=None)
    # early/late with free aR each (2 epochs)
    e2 = fit(early, free_aR=True); l2 = fit(late, free_aR=True)
    res['epoch_free_aR_each'] = dict(early={k: e2[k] for k in ('b', 'eb', 'aR', 'eaR', 'T14', 'T23')},
                                     late={k: l2[k] for k in ('b', 'eb', 'aR', 'eaR', 'T14', 'T23')})
    return res


def T_ld():
    ws = all_windows(); res = {}
    for name, u in (('template_0.401_0.222', U0), ('claret_like_0.44_0.20', (0.44, 0.20)), ('low_0.30_0.25', (0.30, 0.25)),
                    ('high_0.50_0.15', (0.50, 0.15)), ('linear_0.60_0', (0.60, 0.0))):
        A = per_transit(ws, u=u); res[name] = slope(A)
    g = fit(ws, free_u=True); res['global_free_LD'] = {k: g[k] for k in ('b', 'u1', 'u2', 'chi2r')}
    A = per_transit(ws, u=(g['u1'], g['u2'])); res['per_transit_at_global_free_LD'] = slope(A)
    early = [w for w in ws if w['year'] < 2023]; late = [w for w in ws if w['year'] > 2023]
    e2 = fit(early, free_u=True); l2 = fit(late, free_u=True)
    res['epoch_free_LD_each'] = dict(early={k: e2[k] for k in ('b', 'eb', 'u1', 'u2')}, late={k: l2[k] for k in ('b', 'eb', 'u1', 'u2')})
    return res


def T_dilution():
    res = {}
    info = json.load(open(os.path.join(WORK, 'out', 'lc_header_info.json')))
    crowd = {int(k[1:4]): v['crowdsap'] for k, v in info.items() if k.endswith('2m')}
    res['crowdsap'] = crowd
    # SAP with own SPOC-convention correction (same crowd values) and SAP uncorrected (crowd=1)
    A = per_transit(all_windows('sap', crowd=crowd)); res['SAP_crowd_corrected'] = slope(A)
    A = per_transit(all_windows('sap', crowd={s: 1.0 for s in crowd})); res['SAP_no_crowd_corr'] = slope(A)
    # crowding shifted in opposite directions for the two epochs (early sectors more diluted than SPOC states, late less, and reverse)
    for d in (+0.01, -0.01, +0.02):
        cr = {s: (c - d if s < 50 else c + d) for s, c in crowd.items()}
        A = per_transit(all_windows('pdc', crowd=cr)); res['PDC_crowd_shift_early-%.2f_late+%.2f' % (d, d)] = slope(A)
    # k-free per transit (depth absorbs any dilution)
    A = per_transit(all_windows(), free_k=True); res['per_transit_free_k'] = dict(b=slope(A), k=slope(A, 'k', 'ek'))
    return res


def T_duration():
    res = {}
    ws = all_windows()
    A = []
    for w in ws:
        r = fit_trap([w]); r.update(sec=w['sec'], n=w['n']); A.append(r)
    res['per_transit'] = [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()} for r in A]
    for key in ('T14', 'T23', 'T12', 'depth_ppm'):
        res['slope_' + key] = slope(A, key, 'e' + ('depth' if key == 'depth_ppm' else key))
    B = []
    for s in sorted(set(w['sec'] for w in ws)):
        r = fit_trap([w for w in ws if w['sec'] == s]); r['sec'] = s; B.append(r)
    res['per_sector'] = [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()} for r in B]
    for key in ('T14', 'T23', 'T12', 'depth_ppm'):
        res['sector_slope_' + key] = slope(B, key, 'e' + ('depth' if key == 'depth_ppm' else key))
    early = [w for w in ws if w['year'] < 2023]; late = [w for w in ws if w['year'] > 2023]
    res['epoch'] = dict(early=fit_trap(early), late=fit_trap(late))
    # same on SAP
    wsap = all_windows('sap', crowd={5: 0.97656471, 6: 0.98091507, 32: 0.97440773, 98: 0.9725526, 107: 0.97939783})
    A2 = []
    for w in wsap:
        r = fit_trap([w]); r.update(sec=w['sec'], n=w['n']); A2.append(r)
    res['SAP_slope_T14'] = slope(A2, 'T14', 'eT14'); res['SAP_slope_T23'] = slope(A2, 'T23', 'eT23')
    return res


def T_cadence():
    res = {}
    # 20-s for S98/S107 replacing the 2-min
    tags = ['s005_2m', 's006_2m', 's032_2m', 's098_20s', 's107_20s']
    A = per_transit(all_windows(tags=tags)); res['20s_late'] = dict(slope=slope(A), rows=short([a for a in A if a['sec'] > 50]))
    B = per_sector(all_windows(tags=['s098_20s', 's107_20s'])); res['20s_late_per_sector'] = short(B)
    for bt in (600., 1800.):
        A = per_transit(all_windows(bin_to=bt)); res['binned_%ds_all' % bt] = slope(A)
    # late binned to 30 min only (mimic the early FFI cadence asymmetry)
    ws = []
    for tag in SECS2M:
        sec = int(tag[1:4]); t, f, e, exp = load(tag, 'pdc', None, 1800. if sec > 50 else None); ws += windows_from(t, f, e, exp, sec)
    res['late_binned_1800_early_2m'] = slope(per_transit(ws))
    ws = []
    for tag in SECS2M:
        sec = int(tag[1:4]); t, f, e, exp = load(tag, 'pdc', None, 1800. if sec < 50 else None); ws += windows_from(t, f, e, exp, sec)
    res['early_binned_1800_late_2m'] = slope(per_transit(ws))
    return res


def T_spots():
    res = {}
    ws = all_windows()
    A = per_transit(ws, clip=4.0); res['clipped_4sigma'] = slope(A)
    A0 = per_transit(ws)
    # leave-one-transit-out
    loo = []
    for i in range(len(A0)):
        s = slope([a for j, a in enumerate(A0) if j != i]); loo.append(dict(drop=(A0[i]['sec'], A0[i]['n']), slope=s['slope'], z=s['z']))
    res['leave_one_transit_out'] = loo
    # in-transit residual structure per transit (MAP of per-transit fit): binned 10-min residuals inside vs outside
    stats = []
    for w, a in zip(ws, A0):
        tc = a['tcs'][0]
        m = model(w['t'], tc, a['b'], AR0, K0, w['exp'])
        # baseline: refit quadratic on residual ratio
        ratio = w['f'] / m; dt = w['t'] - w['tc0']; cf = np.polyfit(dt, ratio, 2); r = ratio - np.polyval(cf, dt)
        edges = np.arange(w['t'].min(), w['t'].max() + 1e-6, 10 / 1440.)
        idx = np.digitize(w['t'], edges); rb, tb = [], []
        for k in np.unique(idx):
            s = idx == k
            if s.sum() >= 3: rb.append(r[s].mean() / (np.std(r) / np.sqrt(s.sum()))); tb.append(w['t'][s].mean())
        rb, tb = np.array(rb), np.array(tb); inn = np.abs(tb - tc) < 0.058
        stats.append(dict(sec=w['sec'], n=w['n'], in_max_pos_z=float(rb[inn].max()), in_mean_z2=float(np.mean(rb[inn] ** 2)),
                          out_mean_z2=float(np.mean(rb[~inn] ** 2)), rms_ppm=float(np.std(r) * 1e6)))
    res['residual_bins'] = stats
    # out-of-transit variability per sector (rms of 1-d binned PDCSAP and SAP, ppm)
    var = {}
    for tag in SECS2M:
        for fl in ('pdc', 'sap'):
            t, f, e, exp = load(tag, fl)
            edges = np.arange(t.min(), t.max(), 1.0); idx = np.digitize(t, edges)
            mb = [f[idx == k].mean() for k in np.unique(idx) if (idx == k).sum() > 100]
            var[tag + '_' + fl] = dict(rms_1d_ppm=float(np.std(mb) * 1e6), ptp_1d_ppm=float((np.max(mb) - np.min(mb)) * 1e6))
    res['variability'] = var
    return res


def T_ttv():
    res = {}
    ws = all_windows(); A0 = per_transit(ws)
    tcs = np.array([a['tcs'][0] for a in A0]); etc = np.array([a['etcs'][0] for a in A0]); n = np.array([a['n'] for a in A0])
    w = 1 / etc ** 2; M = np.vstack([np.ones_like(n), n]).T.astype(float)
    p = np.linalg.solve(M.T @ (M * w[:, None]), M.T @ (w * tcs)); oc = (tcs - M @ p) * 1440
    res['linear_ephem'] = dict(T0=float(p[0]), P=float(p[1]), oc_min=[round(float(v), 2) for v in oc], eoc_min=[round(float(v), 2) for v in etc * 1440],
                               chi2=float(np.sum((oc / (etc * 1440)) ** 2)), n=len(oc))
    # misfold test: fit per sector with ONE tc offset per sector (fold at linear ephemeris) vs per-transit tc
    def fixed_tc_fit(group):
        # shift each window so that predicted tc (from linear ephemeris) is the fixed centre; fit b + 1 common offset
        x0 = [0.75, 0.0] + [1.0, 0.0, 0.0] * len(group)
        def resid(x):
            b, dto = x[:2]; r = []
            for i, g in enumerate(group):
                tcp = p[0] + p[1] * g['n'] + dto; c0, c1, c2 = x[2 + 3 * i: 5 + 3 * i]; dt = g['t'] - g['tc0']
                r.append((g['f'] - model(g['t'], tcp, b, AR0, K0, g['exp']) * (c0 + c1 * dt + c2 * dt * dt)) / g['e'])
            return np.concatenate(r)
        lo = [0, -0.03] + [0.9, -1, -10] * len(group); hi = [1.1, 0.03] + [1.1, 1, 10] * len(group)
        s = least_squares(resid, x0, bounds=(lo, hi), x_scale='jac')
        chi2r = np.sum(s.fun ** 2) / (s.fun.size - len(x0)); C = np.linalg.inv(s.jac.T @ s.jac) * max(1, chi2r)
        return dict(b=float(s.x[0]), eb=float(np.sqrt(C[0, 0])), year=float(np.mean([g['year'] for g in group])))
    B_fixed = []; B_free = per_sector(ws)
    for s in sorted(set(x['sec'] for x in ws)):
        r = fixed_tc_fit([x for x in ws if x['sec'] == s]); r['sec'] = s; B_fixed.append(r)
    res['per_sector_linear_ephem_fold'] = dict(rows=short(B_fixed, ('sec', 'b', 'eb')), slope=slope(B_fixed))
    res['per_sector_free_tc'] = dict(rows=short(B_free, ('sec', 'b', 'eb')), slope=slope(B_free))
    # whole-data single fold at the linear ephemeris (worst case: one global b per epoch, global tc offset)
    early = [x for x in ws if x['year'] < 2023]; late = [x for x in ws if x['year'] > 2023]
    res['epoch_fold_fixed'] = dict(early=fixed_tc_fit(early), late=fixed_tc_fit(late))
    return res


def T_inject(step=0.25, mode='const'):
    """red-noise injection: synthetic c transits (constant b = 0.748, or trend) at offsets from the real ephemeris,
    into the real 2-min PDCSAP data; same per-transit fitter; null distribution of the slope z and the early/late step."""
    data = {}
    for tag in SECS2M:
        sec = int(tag[1:4]); data[sec] = load(tag)
    rows = []
    for off in np.arange(0.6, P - 0.6, step):
        ws = []
        for sec, (t, f, e, exp) in data.items():
            for w in windows_from(t, f, e, exp, sec, inj=None, offset=off):
                nc = np.round((w['tc0'] - T0) / P)
                if abs(w['tc0'] - (T0 + nc * P)) < HALF + 0.08: continue
                bt = 0.748 if mode == 'const' else 0.748 - 0.0051 * (w['year'] - 2022.4)
                w['f'] = w['f'] * model(w['t'], w['tc0'], bt, AR0, K0, w['exp']); ws.append(w)
        early = [w for w in ws if w['year'] < 2023]; late = [w for w in ws if w['year'] > 2023]
        if len(early) < 3 or len(late) < 3: continue
        A = per_transit(ws); s = slope(A)
        B = per_sector(ws); sb = slope(B)
        rows.append(dict(offset=float(off), n_early=len(early), n_late=len(late), slope=s['slope'], z=s['z'],
                         zstep=s.get('zstep'), chi2=s['chi2'], dof=s['dof'], sector_slope=sb['slope'], sector_z=sb['z'],
                         sector_b=[round(r['b'], 4) for r in B], sector_eb=[round(r['eb'], 4) for r in B], sector_secs=[r['sec'] for r in B]))
        print('off %.2f ne %d nl %d slope %+.4f z %+.2f zstep %+.2f | sector z %+.2f' % (off, len(early), len(late), s['slope'], s['z'], s.get('zstep', np.nan), sb['z']), flush=True)
    z = np.array([r['z'] for r in rows]); zs = np.array([r['sector_z'] for r in rows]); sl = np.array([r['slope'] for r in rows])
    summ = dict(n=len(rows), z_mean=float(z.mean()), z_std=float(z.std()), z_absmax=float(np.abs(z).max()),
                frac_absz_ge5=float(np.mean(np.abs(z) >= 5)), frac_absz_ge3=float(np.mean(np.abs(z) >= 3)),
                slope_mean=float(sl.mean()), slope_std=float(sl.std()),
                sector_z_mean=float(zs.mean()), sector_z_std=float(zs.std()), sector_z_absmax=float(np.abs(zs).max()))
    return dict(mode=mode, summary=summ, rows=rows)


def T_loso():
    """leave one sector out: per-transit slope (fit A) without each sector"""
    ws = all_windows(); A = per_transit(ws); res = {}
    for s in sorted(set(a['sec'] for a in A)):
        res['without_S%d' % s] = slope([a for a in A if a['sec'] != s])
    return res


TESTS = dict(baseline=T_baseline, aR=T_aR, ld=T_ld, dilution=T_dilution, duration=T_duration, cadence=T_cadence,
             spots=T_spots, ttv=T_ttv, inject=lambda: T_inject(0.25, 'const'), inject_trend=lambda: T_inject(0.5, 'trend'), loso=T_loso)

if __name__ == '__main__':
    name = sys.argv[1]
    res = TESTS[name]()
    res['_meta'] = dict(planet='TOI-421 c (TIC 94986319)', test=name, conventions=__doc__.split('\n')[2:6])
    json.dump(res, open(os.path.join(WORK, 'out', name + '.json'), 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k not in ('rows', 'per_transit', 'residual_bins', 'profile_two_epoch')}, default=float, indent=0)[:6000])
