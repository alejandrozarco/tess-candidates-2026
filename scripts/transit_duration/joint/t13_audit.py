"""T13: telemetry and stellar-activity audit of the per-transit durations (t8 output). Own code.
- linear T14(t) fit; residuals tested against each telemetry covariate (Pearson correlation of the standardised covariate
  with the residuals divided by their individual uncertainties, permutation p-value),
  and a multivariate weighted regression with all covariates (+ camera/CCD dummies) to see if the time slope survives;
- activity: per-sector flare rate and SAP variability vs sector-mean residual; rotation-phase sinusoid (Prot 77.3 d,
  Cadieux+2025) and a 30-150 d periodogram of residuals; in-transit spot-crossing bump search vs out-of-transit controls.
Usage: t13_audit.py <pertransit.json> <prep.npz> <main_posterior.json> <flares.npy> <out.json>"""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from t5lib import *
pj, prep, postf, flf, outp = sys.argv[1:6]
R = json.load(open(pj))['rows']; rng = np.random.default_rng(3)
y = np.array([r['year'] for r in R]); T = np.array([r['T14'] for r in R]); e = np.array([r['T14_err'] for r in R]); W = 1 / e ** 2
def wls(Xm, yv, w):
    A = Xm * w[:, None]; C = np.linalg.inv(Xm.T @ A); c = C @ (A.T @ yv); res = yv - Xm @ c
    chi2 = np.sum(w * res ** 2); dof = len(yv) - Xm.shape[1]; s = max(1., chi2 / dof)
    return c, np.sqrt(np.diag(C) * s), res, chi2, dof
X0 = np.vstack([np.ones_like(y), y - TREF]).T
c0, e0, res, chi2, dof = wls(X0, T, W)
out = dict(base=dict(slope=c0[1], err=e0[1], chi2=chi2, dof=dof, n=len(T)))
print('per-transit T14 slope %.3f +- %.3f min/yr (chi2 %.1f/%d)' % (c0[1], e0[1], chi2, dof))
z = res / e
COV = ['bkg_mean', 'bkg_in_minus_out', 'bkg_slope', 'pc_rms', 'mc_rms', 'mc_in_minus_out', 'n_desat_win', 'dt_desat', 'dt_stray', 'dt_edge',
       'n_flare_masked', 'cov_in', 'rms_raw', 'crowd', 'flfr', 'data_rel', 'ap_npix', 'other_in_window', 'dtc_min', 'sig']
def val(r, c):
    v = r[c]; v = float(v) if v is not None else np.nan
    if c in ('dt_desat', 'dt_stray', 'dt_edge'): v = np.log10(min(max(v, 1e-3), 30.))
    return v
corr = {}
for c in COV:
    x = np.array([val(r, c) for r in R]); ok = np.isfinite(x)
    if ok.sum() < 10 or np.std(x[ok]) == 0: continue
    xs = (x[ok] - np.mean(x[ok])) / np.std(x[ok])
    rho = np.corrcoef(xs, z[ok])[0, 1]
    perm = np.array([np.corrcoef(xs, rng.permutation(z[ok]))[0, 1] for _ in range(4000)])
    p = np.mean(np.abs(perm) >= abs(rho))
    rt = np.corrcoef(xs, y[ok])[0, 1]
    corr[c] = dict(r_resid=float(rho), p_perm=float(p), r_with_time=float(rt), n=int(ok.sum()))
    print('%-18s r(resid)=%+.3f p=%.3f   r(time)=%+.3f' % (c, rho, p, rt))
out['univariate'] = corr
# multivariate: time + all continuous covariates (standardised), then + camera/CCD dummies
cont = [c for c in corr]
Xc = []
for c in cont:
    x = np.array([val(r, c) for r in R]); x = np.where(np.isfinite(x), x, np.nanmean(x)); Xc.append((x - x.mean()) / (x.std() + 1e-12))
Xm = np.vstack([np.ones_like(y), y - TREF] + Xc).T
c1, e1, _, chi2m, dofm = wls(Xm, T, W)
camccd = np.array(['%d.%d' % (r['cam'], r['ccd']) for r in R]); u = sorted(set(camccd))
D = [(camccd == k).astype(float) for k in u[1:]]
Xm2 = np.vstack([np.ones_like(y), y - TREF] + Xc + D).T
c2, e2, _, chi2m2, dofm2 = wls(Xm2, T, W)
Xm3 = np.vstack([np.ones_like(y), y - TREF] + D).T
c3, e3, _, _, _ = wls(Xm3, T, W)
cont_nd = [c for c in cont if c != 'data_rel']; Xnd = np.vstack([np.ones_like(y), y - TREF] + [Xc[cont.index(c)] for c in cont_nd]).T
c4, e4, *_ = wls(Xnd, T, W)
dr = np.array([float(r['data_rel']) for r in R]); cp = np.array([r['camp'] for r in R]); drw = dr - np.array([dr[cp == k].mean() for k in cp])
Xdr = np.vstack([np.ones_like(y), y - TREF, (drw - drw.mean()) / (drw.std() + 1e-12)]).T; c5, e5, *_ = wls(Xdr, T, W)
print('slope with all covariates except data_rel: %.3f +- %.3f; with within-campaign data_rel: slope %.3f +- %.3f, data_rel coeff %.3f +- %.3f min/sd' % (c4[1], e4[1], c5[1], e5[1], c5[2], e5[2]))
out['multivariate_no_datarel'] = [c4[1], e4[1]]; out['within_campaign_datarel'] = dict(slope=[c5[1], e5[1]], coeff=[c5[2], e5[2]])
out['multivariate'] = dict(slope_all_cov=[c1[1], e1[1]], chi2=chi2m, dof=dofm, slope_all_cov_camccd=[c2[1], e2[1]], slope_camccd_only=[c3[1], e3[1]], n_camccd=len(u), covariates=cont)
print('slope with all %d covariates: %.3f +- %.3f; + %d cam/CCD dummies: %.3f +- %.3f; cam/CCD only: %.3f +- %.3f' % (len(cont), c1[1], e1[1], len(u) - 1, c2[1], e2[1], c3[1], e3[1]))
# per camera/CCD slopes
pcs = {}
for k in u:
    s = camccd == k
    if s.sum() >= 6 and np.ptp(y[s]) > 2:
        c, ee, *_ = wls(X0[s], T[s], W[s]); pcs[k] = [float(c[1]), float(ee[1]), int(s.sum())]
out['per_camccd'] = pcs; print('per cam/CCD slopes', pcs)
# ---- activity ----
secs = np.array([r['sector'] for r in R]); us = np.unique(secs)
fl = np.load(flf); lcs = {int(l['meta']['sector']): l for l in load_lc('2m')}
act = {}
for s in us:
    l = lcs[int(s)]; g = (l['quality'] == 0) & np.isfinite(l['sap_flux']); t = l['time'][g]; sap = l['sap_flux'][g] / np.nanmedian(l['sap_flux'][g])
    days = np.ptp(t); nfl = np.sum((fl[:, 0] == s) & (fl[:, 2] > 5)) if len(fl) else 0
    # SAP variability: rms of 6-h bin means after a linear trend per orbit-half (gaps > 1 d split)
    seg = np.split(np.arange(len(t)), np.where(np.diff(t) > 1.0)[0] + 1); bm = []
    for sg in seg:
        if len(sg) < 360: continue
        tt, ff = t[sg], sap[sg]; ff = ff - np.polyval(np.polyfit(tt - tt[0], ff, 1), tt - tt[0])
        nb = int(np.ptp(tt) / 0.25)
        for j in range(nb):
            m = (tt >= tt[0] + j * 0.25) & (tt < tt[0] + (j + 1) * 0.25)
            if m.sum() > 50: bm.append(ff[m].mean())
    act[int(s)] = dict(flare_rate=float(nfl / days), sap_var_ppm=float(np.std(bm) * 1e6) if bm else np.nan, mean_resid=float(np.average(res[secs == s], weights=W[secs == s])))
fr = np.array([act[int(s)]['flare_rate'] for s in us]); sv = np.array([act[int(s)]['sap_var_ppm'] for s in us]); mr = np.array([act[int(s)]['mean_resid'] for s in us])
def pc(a, b):
    ok = np.isfinite(a) & np.isfinite(b); r = np.corrcoef(a[ok], b[ok])[0, 1]
    pp = np.mean([abs(np.corrcoef(a[ok], rng.permutation(b[ok]))[0, 1]) >= abs(r) for _ in range(4000)]); return float(r), float(pp)
out['activity'] = dict(per_sector=act, flare_rate_vs_resid=pc(fr, mr), sapvar_vs_resid=pc(sv, mr), flare_rate_vs_time=pc(fr, np.array([y[secs == s].mean() for s in us])))
print('flare rate vs sector residual r=%.3f p=%.3f; SAP variability vs residual r=%.3f p=%.3f' % (*out['activity']['flare_rate_vs_resid'], *out['activity']['sapvar_vs_resid']))
# rotation phase
tb = np.array([lcs[int(r['sector'])]['time'][0] for r in R])  # placeholder not used
tt = (y - 2000.0) * 365.25
def sinfit(Pr):
    Xs = np.vstack([np.ones_like(y), y - TREF, np.sin(2 * np.pi * tt / Pr), np.cos(2 * np.pi * tt / Pr)]).T
    c, ee, _, ch, _ = wls(Xs, T, W); return c, ee, ch
c, ee, chs = sinfit(77.3)
amp = np.hypot(c[2], c[3]); out['rotation'] = dict(Prot=77.3, amp_min=float(amp), dchi2=float(chi2 - chs), slope_with_rot=[float(c[1]), float(ee[1])])
grid = np.linspace(30, 150, 200); dch = np.array([chi2 - sinfit(p)[2] for p in grid])
# false-alarm: same scan on permuted residuals
mx = []
for _ in range(200):
    Tp = X0 @ c0 + rng.permutation(res)
    def ch_p(Pr):
        Xs = np.vstack([np.ones_like(y), y - TREF, np.sin(2 * np.pi * tt / Pr), np.cos(2 * np.pi * tt / Pr)]).T
        return wls(Xs, Tp, W)[3]
    base = wls(X0, Tp, W)[3]; mx.append(max(base - ch_p(p) for p in grid))
out['rotation'].update(scan_max_dchi2=float(dch.max()), scan_best_P=float(grid[np.argmax(dch)]), scan_fap=float(np.mean(np.array(mx) >= dch.max())))
print('rotation 77.3 d: amp %.2f min dchi2 %.2f, slope %.3f +- %.3f; scan max dchi2 %.2f at %.1f d (FAP %.2f)' % (amp, chi2 - chs, c[1], ee[1], dch.max(), grid[np.argmax(dch)], out['rotation']['scan_fap']))
# spot-crossing bumps: main-model residuals in transit vs out-of-transit controls
pr = np.load(prep, allow_pickle=True); segs = list(np.load(str(pr['win']), allow_pickle=True))
Dd = TData(segs); Dd.build(list(pr['keep'])); Dd.set_sig(pr['sig'])
post = json.load(open(postf)); S = Spec('lin_b', D=Dd, ld0=LD0); th = np.array([post['map'][n] for n in S.names]); p = S.unpack(th)
m = Dd.model(p['tc'], p['b'], p['k'], p['aR'], p['u1'], p['u2']); rr = Dd.chi2(m, ret=True)[0]
bump_in, bump_out = [], []
for i in range(Dd.n):
    s = Dd.ti == i; t_ = Dd.t[s] - p['tc'][i]; r_ = rr[s] / Dd.sig[i]
    k5 = 5   # 10-min box (5 x 2-min)
    sm = np.convolve(r_, np.ones(k5) / np.sqrt(k5), mode='same')
    half = t14_from(P_D, p['aR'], p['b'][i], p['k']) / 2
    inn = np.abs(t_) < half * 0.9
    bump_in.append(sm[inn].max() if inn.sum() > k5 else np.nan)
    for off in (-0.11, 0.11):   # control stretches of equal length well outside transit
        oo = np.abs(t_ - off) < half * 0.9
        if oo.sum() > k5: bump_out.append(sm[oo].max())
bump_in = np.array(bump_in); bump_out = np.array(bump_out)
from scipy.stats import ks_2samp
ks = ks_2samp(bump_in[np.isfinite(bump_in)], bump_out)
rz = np.corrcoef(bump_in[np.isfinite(bump_in)], z[np.isfinite(bump_in)])[0, 1]
out['spots'] = dict(n_in_gt3=int(np.sum(bump_in > 3)), frac_in_gt3=float(np.mean(bump_in > 3)), frac_out_gt3=float(np.mean(bump_out > 3)),
                    median_in=float(np.nanmedian(bump_in)), median_out=float(np.median(bump_out)), ks_p=float(ks.pvalue), r_bump_vs_T14resid=float(rz),
                    chi2_in_per_pt=None)
print('spot bumps: in-transit max-box S/N >3 in %d/%d (%.2f) vs controls %.2f; KS p=%.3f; r(bump, T14 resid)=%.3f' % (np.sum(bump_in > 3), np.isfinite(bump_in).sum(), np.mean(bump_in > 3), np.mean(bump_out > 3), ks.pvalue, rz))
json.dump(out, open(outp, 'w'), indent=1, default=float)
