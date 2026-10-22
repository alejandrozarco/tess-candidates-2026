"""C7: white-light NRS1/NRS2 diagnostics per NIRSpec visit of d. Own code.
(1) Joint least-squares fit of the simultaneous NRS1 and NRS2 white light curves ($TDJ_WORK/out/jwst/wl_*.npz of ../joint/run_joint.sh) with
    fixed P, mid-time, a/R* (MAP of the joint white-light fit ns_main of section 6.6.10: a/R* 36.25, mid-times as below), k and quadratic LD per detector (Kurucz prior, sigma 0.1), and
    either one b for both detectors or one b per detector; same systematics model as j2_fit.py. Delta chi2.
(2) Residuals per detector in contact-defined windows (pre-ingress, ingress, flat bottom, egress, post-egress).
(3) Ratio light curve NRS1/NRS2 (same integrations): transit nearly cancels; detrended with a quadratic in time and the
    regressors, compared with the ratio of the two best common-geometry models and the per-detector-b models.
(4) Correlation of the 5-min-binned residuals (each detector and their difference) with the regressors and time.
(5) Spot-crossing search: largest 10-min boxcar excess inside the flat bottom, both detectors.
Usage: c7_white_diag.py <out.json> <fig.png>"""
import sys, json, os, numpy as np, batman
from scipy.optimize import minimize
from scipy.ndimage import median_filter
W = os.path.join(os.environ.get('TDJ_WORK', 'work/transit_duration_joint'), 'out', 'jwst')   # white light of ../joint/run_joint.sh step 9
H = os.environ.get('TDC_WORK', 'work/transit_duration_chromatic')
LD = json.load(open(H + '/out/ld_kurucz.json'))['bins']
P, AR = 7.450729, 36.25; TC = {'1224': 3121.11245, '4098': 3381.88704}
REGS = ['xshift', 'yshift', 'width']
def prep(f):
    z = np.load(f); t = z['t']; keep = (t - t.min()) * 1440 > 10
    def bn(v): v = v[keep]; n = len(v) // 5 * 5; return v[:n].reshape(-1, 5).mean(1)
    tb, fb = bn(t), bn(z['f']); R = {r: bn(z[r]) for r in REGS + ['bkg']}
    rm = fb - median_filter(fb, 15, mode='nearest'); s = 1.4826 * np.median(np.abs(rm)); ok = np.abs(rm) < 5 * s
    return tb, fb, R, ok, t.min()
def t14(b, k): x = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / np.sqrt(AR ** 2 - b ** 2); return P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
def t23(b, k): x = np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / np.sqrt(AR ** 2 - b ** 2); return P / np.pi * np.arcsin(np.clip(x, 0, 1)) * 1440
out = {}; figdat = {}
for v in ('1224', '4098'):
    L = {}
    for det in ('nrs1', 'nrs2'):
        tb, fb, R, ok, t0 = prep(W + '/wl_%s_%s.npz' % (v, det))
        tau = (tb - tb.mean()) / (tb.max() - tb.min())
        X = np.vstack([tau ** j for j in range(3)] + [(R[r] - R[r].mean()) / (R[r].std() + 1e-12) for r in REGS]).T
        pm = batman.TransitParams(); pm.t0 = TC[v]; pm.per = P; pm.rp = 0.045; pm.a = AR; pm.inc = 88.4; pm.ecc = 0; pm.w = 90; pm.limb_dark = 'quadratic'; pm.u = [0, 0.2]
        m = batman.TransitModel(pm, tb, supersample_factor=3, exp_time=np.median(np.diff(tb)) * 2 / 3)   # midpoint sub-exposures
        u0 = np.array([LD[det.upper() + '_white'][0]['u1'], LD[det.upper() + '_white'][0]['u2']])
        L[det] = dict(t=tb, f=fb, X=X, ok=ok, R=R, pm=pm, m=m, u0=u0, sig=1.4826 * np.median(np.abs(np.diff(fb[ok]))) / np.sqrt(2), t0=t0)
    def lcmod(l, k, b, u1, u2):
        l['pm'].rp = k; l['pm'].inc = np.degrees(np.arccos(b / AR)); l['pm'].u = [u1, u2]; return l['m'].light_curve(l['pm'])
    def chi_det(l, k, b, u1, u2, ret=False):
        if not (0.01 < k < 0.1 and 0 < b < 1.05 and -0.5 < u1 < 1 and -0.5 < u2 < 1 and 0 < u1 + u2 < 1): return 1e30
        mm = lcmod(l, k, b, u1, u2); A = (mm[:, None] * l['X'])[l['ok']]; y = l['f'][l['ok']]
        c = np.linalg.lstsq(A, y, rcond=None)[0]; r = y - A @ c
        c2 = np.sum(r ** 2) / l['sig'] ** 2 + np.sum(((np.array([u1, u2]) - l['u0']) / 0.1) ** 2)
        if ret: return c2, mm, (mm[:, None] * l['X']) @ c, (l['X'] @ c)
        return c2
    def common(th): b = th[0]; return chi_det(L['nrs1'], th[1], b, th[2], th[3]) + chi_det(L['nrs2'], th[4], b, th[5], th[6])
    def separate(th): return chi_det(L['nrs1'], th[1], th[0], th[3], th[4]) + chi_det(L['nrs2'], th[5], th[2], th[6], th[7])
    def opt(f, x):
        best = None
        for i in range(5):
            r = minimize(f, x if best is None else best.x, method='Nelder-Mead', options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-6, adaptive=True))
            if best is None or r.fun < best.fun - 1e-4: best = r
            else: break
        return best
    rs = opt(separate, [0.93, 0.045, 0.93, 0.0, 0.2, 0.045, 0.0, 0.2])
    # re-estimate sigma from residuals of the separate fit, then refit both
    for det, (k, b, u1, u2) in (('nrs1', (rs.x[1], rs.x[0], rs.x[3], rs.x[4])), ('nrs2', (rs.x[5], rs.x[2], rs.x[6], rs.x[7]))):
        l = L[det]; _, mm, full, sysm = chi_det(l, k, b, u1, u2, ret=True); l['sig'] = float(np.std((l['f'] - full)[l['ok']]))
    rs = opt(separate, rs.x); rc = opt(common, [0.5 * (rs.x[0] + rs.x[2]), rs.x[1], rs.x[3], rs.x[4], rs.x[5], rs.x[6], rs.x[7]])
    o = dict(sig_ppm={d: L[d]['sig'] * 1e6 for d in L}, sep=dict(b1=rs.x[0], b2=rs.x[2], k1=rs.x[1], k2=rs.x[5], T14_1=float(t14(rs.x[0], rs.x[1])), T14_2=float(t14(rs.x[2], rs.x[5])), chi2=float(rs.fun)),
             com=dict(b=rc.x[0], k1=rc.x[1], k2=rc.x[4], T14_1=float(t14(rc.x[0], rc.x[1])), T14_2=float(t14(rc.x[0], rc.x[4])), chi2=float(rc.fun)))
    o['dchi2_common_minus_separate'] = o['com']['chi2'] - o['sep']['chi2']
    # windows (contacts of the common fit)
    T14c, T23c = t14(rc.x[0], 0.5 * (rc.x[1] + rc.x[4])), t23(rc.x[0], 0.5 * (rc.x[1] + rc.x[4]))
    tc = TC[v]; c1, c2, c3, c4 = tc - T14c / 2880, tc - T23c / 2880, tc + T23c / 2880, tc + T14c / 2880
    o['contacts_min_from_start'] = [float((c - L['nrs1']['t0']) * 1440) for c in (c1, c2, c3, c4)]
    res = {}; fd = {}
    for det, pars in (('nrs1', (rc.x[1], rc.x[0], rc.x[2], rc.x[3])), ('nrs2', (rc.x[4], rc.x[0], rc.x[5], rc.x[6]))):
        l = L[det]; _, mm, full, sysm = chi_det(l, *pars, ret=True); r = (l['f'] - full) / sysm
        r[~l['ok']] = np.nan; res[det] = r; fd[det] = dict(t=l['t'], f=l['f'] / sysm, model=mm, resid=r)
        # separate-b model for the same detector, for the figure
        ps = (rs.x[1], rs.x[0], rs.x[3], rs.x[4]) if det == 'nrs1' else (rs.x[5], rs.x[2], rs.x[6], rs.x[7])
        _, mms, fulls, syss = chi_det(l, *ps, ret=True); fd[det]['model_sep'] = mms; fd[det]['sys_sep'] = syss
    t = L['nrs1']['t']
    win = dict(pre=(t < c1), ingress=(t >= c1) & (t < c2), bottom=(t >= c2) & (t < c3), egress=(t >= c3) & (t < c4), post=(t >= c4))
    o['window_mean_resid_ppm'] = {det: {w: [float(np.nanmean(res[det][s]) * 1e6), float(np.nanstd(res[det][s]) / np.sqrt(np.sum(np.isfinite(res[det][s]))) * 1e6), int(np.sum(np.isfinite(res[det][s])))] for w, s in win.items()} for det in res}
    dr = res['nrs1'] - res['nrs2']
    o['window_mean_resid_diff_ppm'] = {w: [float(np.nanmean(dr[s]) * 1e6), float(np.nanstd(dr[s]) / np.sqrt(np.sum(np.isfinite(dr[s]))) * 1e6)] for w, s in win.items()}
    # correlations with regressors on 5-min bins (~10 points of 5 integrations)
    nb = 10; M = len(t) // nb
    def bb(x): return np.nanmean(x[:M * nb].reshape(M, nb), 1)
    cor = {}
    for name, series in (('nrs1', res['nrs1']), ('nrs2', res['nrs2']), ('nrs1_minus_nrs2', dr)):
        y = bb(series); cor[name] = {}
        for rn in REGS + ['bkg']:
            for det in ('nrs1', 'nrs2'):
                x = bb(L[det]['R'][rn]); ok = np.isfinite(x) & np.isfinite(y)
                rr = float(np.corrcoef(x[ok], y[ok])[0, 1]); rng = np.random.default_rng(1)
                perm = np.array([np.corrcoef(x[ok], rng.permutation(y[ok]))[0, 1] for _ in range(2000)])
                cor[name]['%s_%s' % (rn, det)] = [rr, float(np.mean(np.abs(perm) >= abs(rr)))]
        x = bb(t); ok = np.isfinite(y); cor[name]['time'] = [float(np.corrcoef(x[ok], y[ok])[0, 1]), None]
    o['resid_corr_5min'] = cor
    # regressor excursions during ingress/egress relative to their out-of-transit scatter
    o['regressor_in_window_sigma'] = {}
    for det in ('nrs1', 'nrs2'):
        o['regressor_in_window_sigma'][det] = {}
        for rn in REGS + ['bkg']:
            x = L[det]['R'][rn]; base = x[win['pre'] | win['post']]; mu, sd = np.mean(base), np.std(base)
            o['regressor_in_window_sigma'][det][rn] = {w: float((np.mean(x[s]) - mu) / (sd / np.sqrt(max(1, s.sum())))) for w, s in win.items() if s.sum() > 0}
    # spot crossing: max 10-min boxcar excess inside the flat bottom
    sc = {}
    for det in ('nrs1', 'nrs2'):
        r = res[det]; s = win['bottom']; rb = r[s]; tt = t[s]; n = max(2, int(round(10 / (np.median(np.diff(t)) * 1440))))
        if len(rb) > n:
            box = np.convolve(np.nan_to_num(rb), np.ones(n) / n, mode='valid'); j = int(np.argmax(box))
            sc[det] = dict(max_box_ppm=float(box[j] * 1e6), at_min_from_tc=float((tt[j + n // 2] - tc) * 1440), sigma_box_ppm=float(L[det]['sig'] / np.sqrt(n) * 1e6))
    o['spot_box'] = sc
    # ratio light curve NRS1 / NRS2: detrend with the same regressors on out-of-transit data only
    ratio = L['nrs1']['f'] / L['nrs2']['f']; okr = L['nrs1']['ok'] & L['nrs2']['ok']
    X = np.hstack([L['nrs1']['X'], L['nrs2']['X'][:, 3:]]); oot = (win['pre'] | win['post']) & okr
    c = np.linalg.lstsq(X[oot], ratio[oot], rcond=None)[0]; rr = ratio / (X @ c) - 1
    mc = fd['nrs1']['model'] / fd['nrs2']['model'] - 1; ms = fd['nrs1']['model_sep'] / fd['nrs2']['model_sep'] - 1
    o['ratio'] = dict(rms_ppm=float(np.std(rr[okr]) * 1e6), chi2_vs_common_model=float(np.sum(((rr - mc)[okr]) ** 2) / np.var((rr - mc)[oot])), chi2_vs_sep_model=float(np.sum(((rr - ms)[okr]) ** 2) / np.var((rr - mc)[oot])), npts=int(okr.sum()))
    o['ratio']['dchi2_common_minus_sep'] = o['ratio']['chi2_vs_common_model'] - o['ratio']['chi2_vs_sep_model']
    # linear decomposition of the ratio: baseline (X) + a * (common-model ratio) + beta * (separate - common); beta = 1 means the
    # ratio supports the full per-detector T14 difference of the separate fits, beta = 0 means common geometry
    Xr = np.hstack([X, mc[:, None], (ms - mc)[:, None]]); cr = np.linalg.lstsq(Xr[okr], ratio[okr], rcond=None)[0]
    rres = ratio[okr] - Xr[okr] @ cr; Cr = np.linalg.inv(Xr[okr].T @ Xr[okr]) * np.var(rres)
    nb2 = 10; M2 = len(rres) // nb2; bet = np.std(rres[:M2 * nb2].reshape(M2, nb2).mean(1)) / (np.std(rres) / np.sqrt(nb2))
    o['ratio']['beta_sep'] = [float(cr[-1] / np.mean(X @ c)), float(np.sqrt(Cr[-1, -1]) / np.mean(X @ c) * max(1, bet))]
    o['ratio']['a_common'] = [float(cr[-2] / np.mean(X @ c)), float(np.sqrt(Cr[-2, -2]) / np.mean(X @ c) * max(1, bet))]
    o['ratio']['red_beta'] = float(bet)
    o['ratio']['T14diff_sep'] = o['sep']['T14_1'] - o['sep']['T14_2']
    o['ratio']['T14diff_supported'] = [o['ratio']['beta_sep'][0] * o['ratio']['T14diff_sep'], abs(o['ratio']['beta_sep'][1] * o['ratio']['T14diff_sep'])]
    o['ratio']['window_mean_minus_common_ppm'] = {w: float(np.mean((rr - mc)[s & okr]) * 1e6) for w, s in win.items()}
    o['ratio']['window_mean_ppm'] = {w: [float(np.mean(rr[s & okr]) * 1e6), float(np.std(rr[oot]) / np.sqrt(max(1, (s & okr).sum())) * 1e6)] for w, s in win.items()}
    fd['ratio'] = dict(t=t, r=rr, mc=mc, ms=ms); fd['contacts'] = [c1, c2, c3, c4]; fd['tc'] = tc
    figdat[v] = fd; out[v] = o
    print('== visit', v, json.dumps({k: o[k] for k in ('sep', 'com', 'dchi2_common_minus_separate', 'contacts_min_from_start')}, default=float))
    print('   window resid (ppm):', json.dumps(o['window_mean_resid_ppm']))
    print('   NRS1-NRS2 resid by window:', json.dumps(o['window_mean_resid_diff_ppm']))
    print('   ratio:', json.dumps(o['ratio']))
    print('   spot box:', json.dumps(o['spot_box']))
    sig = {n: {k: x for k, x in d.items() if x[1] is not None and x[1] < 0.05} for n, d in cor.items()}
    print('   resid-regressor correlations with permutation p < 0.05:', json.dumps(sig))
    print('   regressor excursions (sigma of mean) by window:', json.dumps({d: {r: {w: round(x, 1) for w, x in q.items()} for r, q in rr_.items()} for d, rr_ in o['regressor_in_window_sigma'].items()}))
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(3, 2, figsize=(12, 10), sharex='col')
for j, v in enumerate(('1224', '4098')):
    fd = figdat[v]; tc = fd['tc']
    for det, col in (('nrs1', 'C0'), ('nrs2', 'C3')):
        x = (fd[det]['t'] - tc) * 1440
        ax[0, j].plot(x, fd[det]['f'] - (0.003 if det == 'nrs2' else 0), '.', ms=2, color=col, alpha=0.4)
        ax[0, j].plot(x, fd[det]['model'] - (0.003 if det == 'nrs2' else 0), '-', color='k', lw=1)
        r = fd[det]['resid']; n = 10; M = len(r) // n
        ax[1, j].plot(x[:M * n].reshape(M, n).mean(1), np.nanmean(r[:M * n].reshape(M, n), 1) * 1e6, 'o-', ms=3, color=col, label=det.upper())
    x = (fd['ratio']['t'] - tc) * 1440; r = fd['ratio']['r']; n = 10; M = len(r) // n
    ax[2, j].plot(x[:M * n].reshape(M, n).mean(1), r[:M * n].reshape(M, n).mean(1) * 1e6, 'ko', ms=3, label='NRS1/NRS2 - 1 (5-min bins)')
    ax[2, j].plot(x, fd['ratio']['mc'] * 1e6, 'C2-', label='common b')
    ax[2, j].plot(x, fd['ratio']['ms'] * 1e6, 'C1-', label='b per detector')
    for a in ax[:, j]:
        for c in fd['contacts']: a.axvline((c - tc) * 1440, color='0.6', lw=0.8, ls=':')
    ax[0, j].set_title('visit %s (white light, common-b model)' % v); ax[1, j].set_ylabel('resid (ppm)'); ax[2, j].set_ylabel('ratio (ppm)'); ax[2, j].set_xlabel('min from mid-transit')
    ax[1, j].legend(fontsize=8); ax[2, j].legend(fontsize=8)
plt.tight_layout(); plt.savefig(sys.argv[2], dpi=110)
