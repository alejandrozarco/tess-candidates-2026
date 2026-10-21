"""T22: significance diagnostics of the TESS duration trend. Combines: the posterior / nominal LR statistic, the
bootstrap widths, the null-injection distribution (width of z, exceedances of the observed statistic), and an
illustrative Gaussian look-elsewhere adjustment (independent trials assumed) (801 planets with a screening statistic out of the 802 searched). Own code.
Usage: t22_significance.py [survey_summary.csv]  (the survey table is not part of this repository; without it the
survey rank block is skipped, the look-elsewhere factors use the fixed trial numbers 801 and 3204)"""
import json, numpy as np, sys, os
from scipy import stats
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import ROOT
S = json.load(open(ROOT + '/out/summary.json'))
main = S['mcmc']['main_lin_b']; q = main['derived']['dT14dt_min_per_yr']; slope = q[0]; sig_post = 0.5 * (q[2] - q[1])
dchi2 = S['comparisons']['main_lin_b_vs_main_const']['dchi2']
res = dict(slope=slope, sig_post=sig_post, z_post=slope / sig_post, dchi2=dchi2, z_LR=float(np.sqrt(dchi2)))
for m in ('transit', 'sector'):
    if m in S.get('bootstrap', {}): res['z_boot_' + m] = slope / S['bootstrap'][m]['dT14_sd']; res['sd_boot_' + m] = S['bootstrap'][m]['dT14_sd']
inj = {}
for m in ('lc', 'pix2'):
    f = ROOT + '/out/inj_%s.jsonl' % m
    if not os.path.exists(f): continue
    R = [json.loads(l) for l in open(f)]
    dt = np.array([r['dT14'] for r in R]); de = np.array([r['dT14_err'] for r in R]); z = dt / de; dc = np.array([r['dchi2'] for r in R])
    zsd = float(np.std(z)); zsd_lo, zsd_hi = [float(x) for x in np.percentile([np.std(np.random.default_rng(i).choice(z, len(z))) for i in range(2000)], [16, 84])]
    inj[m] = dict(n=len(R), bias_min_per_yr=float(dt.mean()), bias_err=float(dt.std() / np.sqrt(len(dt))), sd_slope=float(dt.std()), z_mean=float(z.mean()), z_sd=zsd, z_sd_68=[zsd_lo, zsd_hi],
                  n_abs_slope_ge_obs=int(np.sum(np.abs(dt) >= abs(slope))), n_dchi2_ge_obs=int(np.sum(dc >= dchi2)), max_abs_z=float(np.max(np.abs(z))), dchi2_max=float(dc.max()))
    inj[m]['z_calibrated'] = float(slope / (sig_post * max(1.0, zsd)))
    inj[m]['z_calibrated_conservative'] = float(slope / (sig_post * max(1.0, zsd_hi)))
    inj[m]['z_vs_injection_slope_sd'] = float((slope - dt.mean()) / dt.std())
res['injections'] = inj
# campaign structure
c = S['comparisons']
res['within_campaign'] = dict(dchi2=c['v_camp_slope_vs_v_camp']['dchi2'], z=c['v_camp_slope_vs_v_camp']['z_LR'],
                              bdot=S['mcmc']['v_camp_slope']['post']['bdot'])
res['campaign_means_vs_line'] = dict(dchi2=c['v_camp_vs_main_lin_b']['dchi2'], ddof=2)
# look-elsewhere: survey of 801 analysed planets; 4 flag statistics -> up to 3204 trials
Ntr = {'planets (801)': 801, 'planets x 4 flag types (3204)': 3204}
zs = {'posterior': res['z_post'], 'LR': res['z_LR']}
if 'lc' in inj: zs['injection-calibrated'] = inj['lc']['z_calibrated']; zs['injection-calibrated (84th pct width)'] = inj['lc']['z_calibrated_conservative']
if 'sector' in S.get('bootstrap', {}): zs['sector-bootstrap'] = abs(res['z_boot_sector'])
le = {}
for zn, z in zs.items():
    z = abs(z); lp = 2 * stats.norm.sf(z)
    for tn, N in Ntr.items():
        pg = -np.expm1(N * np.log1p(-lp)) if lp > 0 else N * lp
        pg = max(pg, N * lp if lp < 1e-300 else pg)
        zg = stats.norm.isf(pg / 2) if pg > 0 else np.inf
        le['%s | %s' % (zn, tn)] = dict(z_local=z, p_local=lp, p_global=pg, z_global=zg)
res['look_elsewhere'] = le
# survey empirical: rank of L 98-59 d's screening statistic among the stable KP references
import csv
if len(sys.argv) > 1:
    rows = list(csv.DictReader(open(sys.argv[1])))
else: rows = []
zw = np.array([abs(float(r['zmax_width'])) for r in rows if r['status'] == 'ok' and r['zmax_width'] not in ('', 'nan')])
zst = np.array([abs(float(r['zmax_width'])) for r in rows if r['status'] == 'ok' and r['stable_ref'] == 'True' and r['zmax_width'] not in ('', 'nan')])
mine = [abs(float(r['zmax_width'])) for r in rows if r['pid'] == 'TOI-175.02'][0] if rows else None
if rows: res['survey'] = dict(n_planets=len(zw), screening_stat=mine, n_ge_all=int(np.sum(zw >= mine)), n_stable=len(zst), n_ge_stable=int(np.sum(zst >= mine)),
                     stable_mad_sigma=float(1.4826 * np.median(zst)), stable_99pct=float(np.percentile(zst, 99)))
json.dump(res, open(ROOT + '/out/significance.json', 'w'), indent=1, default=float)
for k, v in res.items():
    if k != 'look_elsewhere': print(k, v)
for k, v in le.items(): print('LEE %-62s z_local %.2f p_local %.2e -> p_global %.2e (z %.2f)' % (k, v['z_local'], v['p_local'], v['p_global'], v['z_global']))
