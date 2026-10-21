"""Write the tables of section 6.6.10 (joint analysis of L 98-59 d) from the outputs of run_joint.sh:
data/jointfit_307210830.csv (TESS posterior summaries and least-squares variants), data/significance_307210830.csv
(model comparison, bootstraps, injections, look-elsewhere, telemetry and activity audit) and
data/jwst_joint_307210830.csv (joint fits of the own JWST white-light extractions).
Usage (repository root): python scripts/transit_duration/joint/make_tables_joint.py   [$TDJ_WORK as in run_joint.sh]"""
import json, os, csv, datetime
import numpy as np
from scipy import stats
W = os.environ.get('TDJ_WORK', os.path.join('work', 'transit_duration_joint'))
O = os.path.join(W, 'out')
S = json.load(open(os.path.join(O, 'summary.json')))
G = json.load(open(os.path.join(O, 'significance.json')))
JS = json.load(open(os.path.join(O, 'jwst_slopes.json')))
DATE = os.environ.get('TABLE_DATE', datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d'))
f6 = lambda x: '' if x is None else '%.6g' % x

def write(fn, header, cols, rows):
    with open(fn, 'w', newline='') as fo:
        for h in header: fo.write('# ' + h + '\n')
        w = csv.writer(fo, lineterminator='\n'); w.writerow(cols)
        for r in rows: w.writerow(r)
    print('wrote', fn, len(rows), 'rows')

# ---------------------------------------------------------------- joint TESS fits
DESC = {'main_lin_b': 'main: b linear in time', 'main_const': 'main: constant b', 'v_lin_i': 'inclination linear in time',
        'v_lin_b_kcamp': 'b linear, one k per campaign', 'v_const_kcamp': 'constant b, one k per campaign',
        'v_lin_b_ldfree': 'b linear, limb darkening free', 'v_const_ldfree': 'constant b, limb darkening free',
        'v_lin_b_ecc': 'b linear, e and omega with priors', 'v_const_ecc': 'constant b, e and omega with priors',
        'v_lin_b_fixtc': 'b linear, per-transit mid-times fixed', 'v_const_fixtc': 'constant b, per-transit mid-times fixed',
        'v_lin_b_rhofree': 'b linear, uniform prior on ln rho*: U(8, 11)', 'v_camp': 'one b per campaign',
        'v_camp_slope': 'one b per campaign plus a common within-campaign slope of b'}
rows = []
for v, e in S['mcmc'].items():
    base = [v, DESC.get(v, v), 'emcee', 72, e['npar']]
    for st in ('chi2', 'bic', 'lnZ', 'lnZ_lap', 'ess', 'rhat_max', 'n_over_tau', 'nens', 'nsteps'):
        rows.append(base + ['fit:' + st, f6(e[st]), '', ''])
    for n, q in e['post'].items():
        rows.append(base + [n, f6(q[0]), f6(q[0] - q[1]), f6(q[2] - q[0])])
    for n, q in e['derived'].items():
        if n == 'T14_camp_min':
            if e['derived'].get('dT14dt_min_per_yr', [0])[0] == 0 and v.find('const') >= 0:
                q = q[:1]; labels = ['T14_min (constant)']
            else: labels = ['T14_min at %.2f' % y for y in (2019.13, 2021.09, 2023.28, 2025.44)]
            for lab, qq in zip(labels, q): rows.append(base + [lab, f6(qq[0]), f6(qq[0] - qq[1]), f6(qq[2] - qq[0])])
        elif n == 'dT14dt_min_per_yr' and q[0] == 0: continue
        else: rows.append(base + [n, f6(q[0]), f6(q[0] - q[1]), f6(q[2] - q[0])])
ALT = {'pdc2m': 'PDCSAP 2 min', 'pdc2m_regs': 'PDCSAP 2 min, pointing and background regressors', 'sap2m': 'SAP 2 min',
       'pdc2m_ebmodel': 'PDCSAP 2 min, eclipses modelled instead of masked', 'pdc20s': 'PDCSAP 20 s', 'sap20s': 'SAP 20 s',
       'own_spoc': 'own TPF photometry, SPOC aperture', 'own_spoc_regs': 'own TPF photometry, SPOC aperture, regressors',
       'own_r15': 'own TPF photometry, r = 1.5 px', 'own_r15_regs': 'own TPF photometry, r = 1.5 px, regressors',
       'own_r25': 'own TPF photometry, r = 2.5 px', 'own_r25_regs': 'own TPF photometry, r = 2.5 px, regressors'}
for v, a in S['alt'].items():
    base = ['alt_' + v, ALT.get(v, v), 'least squares', a['ntr'], 10]
    rows.append(base + ['dT14dt_min_per_yr', f6(a['dT14']), f6(a['dT14_err']), f6(a['dT14_err'])])
    rows.append(base + ['dbdt_per_yr', f6(a['bdot']), f6(a['bdot_err']), f6(a['bdot_err'])])
    rows.append(base + ['dchi2 (constant minus linear b)', f6(a['dchi2']), '', ''])
    for n in ('k', 'b0', 'lnrho'): rows.append(base + [n, f6(a[n]), '', ''])
write('data/jointfit_307210830.csv', [
    'L 98-59 d (TIC 307210830): joint fits of all usable TESS transits (scripts/transit_duration/joint/, run_joint.sh), %s' % DATE,
    'data: SPOC 2-min PDCSAP, quality 0, 72 transits in 26 sectors (2 to 98); planets b and c, the 1.0492-d eclipses (+-0.95 h) and flares masked',
    'model: batman quadratic limb darkening on an own sky-separation array, 5x supersampling, per-transit quadratic baseline profiled; noise: per-transit white sigma times per-sector beta',
    'priors: ln rho* [kg m-3] ~ N(9.50, 0.09) (Cadieux et al. 2025); u1, u2 ~ N(0.162, 0.1), N(0.441, 0.1) via Kipping q1, q2 (Claret 2017); TTV sinusoid at 396.32 d; e = 0 unless stated',
    'parameters: k = Rp/R*; lnrho; q1, q2; b0 = b at 2022.0; bdot = db/dt (per yr); i0, idot (deg, deg/yr); bc0-bc3 = b per campaign; k0-k3 = k per campaign; dT0 (min), dP (s), At, Bt (min) = ephemeris and TTV terms; secw, sesw = sqrt(e) cos/sin omega',
    'derived: T14 (first to fourth contact, min) at the four campaign mean epochs; dT14dt_min_per_yr = T14(2022.5) - T14(2021.5); aR = a/R*',
    'method emcee: median, err_minus = median - 16th percentile, err_plus = 84th percentile - median; fit:* rows = MAP chi2, BIC, ln profile evidence (importance sampling of the profiled likelihood times the unnormalised prior weights, baselines not marginalised; lnZ_lap = Laplace), importance-sampling ESS, max rank-normalised split R-hat, N/tau, ensembles, steps per ensemble',
    'method least squares (alt_*): same windows and noise recipe, symmetric errors from the covariance; n_transits differs for the 20-s data (sectors 28-98 only)'],
    ['variant', 'description', 'method', 'n_transits', 'npar', 'quantity', 'median', 'err_minus', 'err_plus'], rows)

# ---------------------------------------------------------------- significance
NONNESTED = ('v_camp_vs_main_lin_b',)   # one b per campaign vs b linear: not nested, no chi2 likelihood-ratio calibration
rows = []
CMP = {'main_lin_b_vs_main_const': 'b linear vs constant (main)', 'v_lin_i_vs_main_const': 'i linear vs constant b',
       'v_lin_b_kcamp_vs_v_const_kcamp': 'b linear vs constant, k per campaign', 'v_lin_b_ldfree_vs_v_const_ldfree': 'b linear vs constant, LD free',
       'v_lin_b_ecc_vs_v_const_ecc': 'b linear vs constant, e and omega', 'v_lin_b_fixtc_vs_v_const_fixtc': 'b linear vs constant, mid-times fixed',
       'v_camp_slope_vs_v_camp': 'within-campaign slope of b vs one b per campaign', 'v_camp_vs_main_lin_b': 'one b per campaign vs b linear'}
for k, c in S['comparisons'].items():
    g = 'model comparison: ' + CMP.get(k, k)
    rows += [[g, 'delta chi2', f6(c['dchi2']), '', 'extra parameters %d' % c['ddof']], [g, 'delta BIC', f6(c['dBIC']), '', ''],
             [g, 'profile-evidence log contrast', f6(c['lnBF']), '', 'importance sampling of the profiled likelihood times the implemented prior weights; uniform prior width of the slope parameter included; baselines not marginalised, prior normalisations incomplete (i0 prior of the linear-i model); not a fully normalised Bayes factor' if c['lnBF'] is not None else 'not computed'],
             [g, 'likelihood-ratio z', '', '', 'not applicable: non-nested models'] if k in NONNESTED or c.get('z_LR') is None else [g, 'likelihood-ratio z', f6(c['z_LR']), '', 'p = %.3g' % c['p_LR']]]
rows.append(['posterior', 'dT14/dt (min/yr)', f6(G['slope']), f6(G['sig_post']), 'z = %.2f' % G['z_post']])
rows.append(['within campaigns', 'db/dt (per yr), common slope within campaigns', f6(S['mcmc']['v_camp_slope']['post']['bdot'][0]),
             f6(0.5 * (S['mcmc']['v_camp_slope']['post']['bdot'][2] - S['mcmc']['v_camp_slope']['post']['bdot'][1])), 'delta chi2 %.2f, z %.2f' % (G['within_campaign']['dchi2'], G['within_campaign']['z'])])
for m in ('transit', 'sector'):
    b = S['bootstrap'][m]
    rows.append(['block bootstrap: whole %ss' % m, 'dT14/dt (min/yr)', f6(b['dT14_med']), f6(b['dT14_sd']), 'n = %d; fraction >= 0: %g' % (b['n'], b['frac_dT14_ge0'])])
    rows.append(['block bootstrap: whole %ss' % m, 'db/dt (per yr)', f6(b['bdot_med']), f6(b['bdot_sd']), 'n = %d' % b['n']])
j = S['bootstrap']['jack']
for lab, d in zip(j['labels'], j['dT14_all']):
    rows.append(['jackknife', 'dT14/dt (min/yr), %s left out' % lab.replace('camp0', 'campaign 2019').replace('camp1', 'campaign 2021').replace('camp2', 'campaign 2023').replace('camp3', 'campaign 2025').replace('sec', 'sector '), f6(d), '', ''])
for m, lab in (('lc', 'null injections into the PDCSAP flux'), ('pix2', 'null injections into the target pixels')):
    i = G['injections'][m]; s2 = S['inj_' + m]
    rows += [[lab, 'n', i['n'], '', ''], [lab, 'mean dT14/dt (min/yr)', f6(i['bias_min_per_yr']), f6(i['bias_err']), ''],
             [lab, 'sd of dT14/dt (min/yr)', f6(i['sd_slope']), '', ''], [lab, 'z = slope/error: mean', f6(i['z_mean']), '', ''],
             [lab, 'z = slope/error: sd', f6(i['z_sd']), '', '68%% interval %.3f-%.3f' % tuple(i['z_sd_68'])],
             [lab, 'max delta chi2', f6(i['dchi2_max']), '', 'fraction > 9: %g; KS p vs chi2(1): %.3g' % (s2['frac_dchi2_gt9'], s2['ks_dchi2_vs_chi2_1'])],
             [lab, 'injections with |dT14/dt| >= observed', i['n_abs_slope_ge_obs'], '', ''], [lab, 'injections with delta chi2 >= observed', i['n_dchi2_ge_obs'], '', ''],
             [lab, 'z: observed slope / injection sd', f6(i['z_vs_injection_slope_sd']), '', '']]
# look-elsewhere: the entries of significance.json, and the slope / injection-scatter statistics
LE = dict(G['look_elsewhere'])
for m, lab in (('lc', 'slope / injection sd, light-curve level'), ('pix2', 'slope / injection sd, pixel level')):
    z = abs(G['injections'][m]['z_vs_injection_slope_sd']); lp = 2 * stats.norm.sf(z)
    for tn, N in (('planets (801)', 801), ('planets x 4 flag types (3204)', 3204)):
        pg = -np.expm1(N * np.log1p(-lp)); LE['%s | %s' % (lab, tn)] = dict(z_local=z, p_local=lp, p_global=pg, z_global=stats.norm.isf(pg / 2))
for k, v in LE.items():
    st, tn = k.split(' | ')
    st = {'injection-calibrated': 'posterior scaled by the injection scatter of slope/error', 'injection-calibrated (84th pct width)': 'posterior scaled by the injection scatter of slope/error (84th percentile)', 'LR': 'likelihood ratio (nominal)', 'posterior': 'posterior slope/width', 'sector-bootstrap': 'sector-bootstrap slope/width'}.get(st, st)
    rows.append(['look-elsewhere: ' + st, 'z_global, trials = ' + tn, f6(v['z_global']), '', 'z_local %.2f, p_local %.3g, p_global %.3g' % (v['z_local'], v['p_local'], v['p_global'])])
if 'survey' in G:
    s = G['survey']
    rows += [['survey', 'planets with a screening statistic', s['n_planets'], '', ''], ['survey', 'screening statistic of L 98-59 d (duration change)', f6(s['screening_stat']), '', ''],
             ['survey', 'planets with a statistic >= that of L 98-59 d', s['n_ge_all'], '', ''], ['survey', 'stable reference planets', s['n_stable'], '', ''],
             ['survey', 'stable reference planets with a statistic >= that of L 98-59 d', s['n_ge_stable'], '', '']]
A = S['audit']
rows.append(['telemetry audit', 'per-transit T14 slope (min/yr)', f6(A['base']['slope']), f6(A['base']['err']), 'chi2 %.1f for %d dof' % (A['base']['chi2'], A['base']['dof'])])
for k, u in A['univariate'].items():
    rows.append(['telemetry audit: correlation of T14 residual with covariate', k, f6(u['r_resid']), '', 'permutation p = %.3f; correlation with time %.2f' % (u['p_perm'], u['r_with_time'])])
rows.append(['telemetry audit', 'slope with all covariates except data release (min/yr)', f6(A['multivariate_no_datarel'][0]), f6(A['multivariate_no_datarel'][1]), ''])
rows.append(['telemetry audit', 'slope with camera/CCD offsets (min/yr)', f6(A['multivariate']['slope_camccd_only'][0]), f6(A['multivariate']['slope_camccd_only'][1]), ''])
rows.append(['telemetry audit', 'slope with all covariates including data release (min/yr)', f6(A['multivariate']['slope_all_cov'][0]), f6(A['multivariate']['slope_all_cov'][1]), 'data release correlates with time (r = 0.97)'])
rows.append(['telemetry audit', 'data release coefficient within campaigns (min per sd)', f6(A['within_campaign_datarel']['coeff'][0]), f6(A['within_campaign_datarel']['coeff'][1]), ''])
for cc, (sl, er, n) in A['per_camccd'].items():
    rows.append(['telemetry audit: per camera.CCD', 'slope (min/yr), camera.CCD %s' % cc, f6(sl), f6(er), 'n = %d transits' % n])
ac = A['activity']; ro = A['rotation']; sp = A['spots']
rows += [['activity', 'flare rate vs sector-mean T14 residual: r', f6(ac['flare_rate_vs_resid'][0]), '', 'p = %.3f' % ac['flare_rate_vs_resid'][1]],
         ['activity', '6-h SAP variability vs sector-mean T14 residual: r', f6(ac['sapvar_vs_resid'][0]), '', 'p = %.3f' % ac['sapvar_vs_resid'][1]],
         ['activity', 'flare rate vs time: r', f6(ac['flare_rate_vs_time'][0]), '', 'p = %.3f' % ac['flare_rate_vs_time'][1]],
         ['activity', 'rotation-phase sinusoid (P = %.1f d): amplitude (min)' % ro['Prot'], f6(ro['amp_min']), '', 'delta chi2 %.2f for 2 dof' % ro['dchi2']],
         ['activity', 'slope with the rotation term (min/yr)', f6(ro['slope_with_rot'][0]), f6(ro['slope_with_rot'][1]), ''],
         ['activity', 'residual periodogram 30-150 d: best period (d)', f6(ro['scan_best_P']), '', 'delta chi2 %.1f, permutation FAP %.3f' % (ro['scan_max_dchi2'], ro['scan_fap'])],
         ['activity', 'transits with an in-transit bump > 3 sigma (fraction)', f6(sp['frac_in_gt3']), '', 'out-of-transit controls %.4f; KS p = %.2f' % (sp['frac_out_gt3'], sp['ks_p'])]]
write('data/significance_307210830.csv', [
    'L 98-59 d (TIC 307210830): significance of the TESS duration trend and audits (scripts/transit_duration/joint/; t9_boot.py, t10_inject.py, t13_audit.py, t20_summary.py, t22_significance.py), %s' % DATE,
    'model comparison: evolving against constant geometry, same data and noise model; delta chi2 and delta BIC (separate MAP fits) positive favour the first model; likelihood-ratio z nominal, not a calibrated test',
    'look-elsewhere (illustrative): p_global = 1 - (1 - p_local)^N, N = 801 (planets with a screening statistic) or 3204 (801 x 4 flag statistics); two-sided normal p; assumes Gaussian tails and N independent, correctly calibrated trials (not established)',
    'injections: constant-geometry copy of d (b = 0.923, posterior k, a/R*, LD) at phases U(0.08, 0.92), inserted into quality-selected, normalised PDCSAP flux, or into calibrated pixels after the own background subtraction and before aperture summation (whole background-subtracted image multiplied); the subsequent transit-processing and fitting steps are rerun; mean dT14/dt uncertainty = sd / sqrt(n)',
    'telemetry audit: per-transit T14 (shape fixed at the posterior median) minus a linear trend; Pearson correlation of each covariate with the residuals divided by their individual uncertainties, 4000-permutation p-value; time covariates (dt_*) as log10 days',
    'uncertainty column: 1-sigma where given'],
    ['group', 'statistic', 'value', 'uncertainty', 'note'], rows)

# ---------------------------------------------------------------- JWST
DJ = {'ns_main_v2': 'NIRSpec pair, NRS1 + NRS2, main', 'ns_gp': 'NIRSpec pair, GP noise', 'ns_aRfree': 'NIRSpec pair, a/R* free',
      'ns_aR33.6': 'NIRSpec pair, a/R* = 33.6', 'ns_aR37.1': 'NIRSpec pair, a/R* = 37.1', 'ns_ldGressier': 'NIRSpec pair, LD fixed (Gressier et al. 2024)',
      'ns_sqrtlaw': 'NIRSpec pair, square-root law', 'ns_power2': 'NIRSpec pair, power-2 law', 'ns_noregs': 'NIRSpec pair, no detector regressors',
      'ns_poly1': 'NIRSpec pair, time polynomial order 1', 'ns_poly3': 'NIRSpec pair, time polynomial order 3', 'ns_ramp': 'NIRSpec pair, exponential ramp',
      'ns_matchcov': 'NIRSpec pair, matched pre/post coverage', 'ns_nrs1': 'NRS1 only', 'ns_gp_nrs1': 'NRS1 only, GP noise',
      'ex_hw5_nrs1': 'NRS1 only, box half-width 5 px', 'ns_matchcov_nrs1': 'NRS1 only, matched coverage', 'ns_nrs2': 'NRS2 only',
      'ns_gp_nrs2': 'NRS2 only, GP noise', 'ex_hw5_nrs2': 'NRS2 only, box half-width 5 px', 'ns_matchcov_nrs2': 'NRS2 only, matched coverage',
      'ex_blue_nrs2': 'NRS2 blue half (3.83-4.45 um)', 'ex_red_nrs2': 'NRS2 red half (4.45-5.15 um)', 'pipe_nrs2': 'NRS2 pipeline white light (stage 3)',
      'pipe_nrs2_ramp': 'NRS2 pipeline white light, exponential ramp', 'all_visits': 'NIRISS 1201 and 4098o1 (x1dints) + NIRSpec pair; 4098o1 excluded only from the duration-versus-time regression',
      'all_visits_gp': 'NIRISS 1201 and 4098o1 (x1dints) + NIRSpec pair, GP noise; 4098o1 excluded only from the duration-versus-time regression', 'all_visits_2d': 'NIRISS 1201 and 4098o1 (own 2D) + NIRSpec pair; both NIRISS visits excluded only from the duration-versus-time regression (T14 error > 1 min)',
      'all_visits_2d_nrs1': 'NIRISS 1201 and 4098o1 (own 2D) + NRS1; both NIRISS visits excluded only from the duration-versus-time regression (T14 error > 1 min)', 'all_visits_2d_nrs2': 'NIRISS 1201 and 4098o1 (own 2D) + NRS2; both NIRISS visits excluded only from the duration-versus-time regression (T14 error > 1 min)', 'ni2d_pair': 'NIRISS pair, own 2D',
      'ctrl_c_nrs1': 'planet c control, GO 2512 pair, NRS1', 'ni_1201': 'NIRISS 1201 alone (x1dints)'}
rows = []
for n, d in S['jwst'].items():
    if n not in DJ: continue
    fj = json.load(open(os.path.join(O, 'jwst', 'fits', n + '.json')))
    dc = fj.get('const_same_noise', {}).get('dchi2')
    for v, q in d['T14'].items():
        rows.append([n, DJ[n], v, 'T14_min', f6(q[0]), f6(q[1]), '', '%.4f +- %.4f' % (d['b'][v][0], d['b'][v][1])])
    rows.append([n, DJ[n], '', 'a/R*', f6(d['aR'][0]), f6(d['aR'][1]), '', ''])
    if d.get('dT14'): rows.append([n, DJ[n], 'second minus first', 'delta_T14_min', f6(d['dT14'][0]), f6(d['dT14'][1]), f6(dc), 'delta b %.4f +- %.4f' % (d['db'][0], d['db'][1])])
    if n in JS:
        s = JS[n]; rows.append([n, DJ[n], ','.join(s['visits']), 'slope_min_per_yr', f6(s['slope'][0]), f6(s['slope'][1]), '', 'JWST minus TESS %.2f sigma' % s['tension_sigma']])
    rows.append([n, DJ[n], '', 'fit:n_over_tau', f6(d['n_over_tau']), '', '', 'max R-hat %.3f' % d['rhat_max']])
write('data/jwst_joint_307210830.csv', [
    'L 98-59 d (TIC 307210830): joint fits of JWST white-light curves from own calints/x1dints extractions and stage-3 pipeline comparisons (scripts/transit_duration/joint/j1_extract.py, j1b_soss.py, j1c_soss2d.py, j2_fit.py, t21_jwst_slope.py), %s' % DATE,
    'visits: 1201 = GTO 1201 obs. 311 NIRISS SOSS (2022.99); 1224 = GTO 1224 obs. 1 NIRSpec G395H (2023.48); 4098 = GO 4098 obs. 2 NIRSpec G395H (2024.20); 4098o1 = GO 4098 obs. 1 NIRISS SOSS (2023.99); c2512a/b = GO 2512 obs. 15/16 (planet c)',
    'bands: NRS1 2.870-3.710 um, NRS2 3.830-5.150 um (same wavelength bounds applied independently to each visit; the selected detector columns are not identical), box half-width 3 px unless stated; NIRISS order 1',
    'model: shared a/R* (stellar-density prior ln rho ~ N(9.50, 0.09) unless stated), k and LD (quadratic, free) per detector, b and mid-time per visit; quadratic in time plus trace x/y shift and width regressors, profiled; white noise times beta, or celerite Matern-3/2 GP',
    'values: posterior median and standard deviation; delta_chi2_const = chi2 increase when all b are tied (same noise model, MAP); slope: weighted line through the visits per posterior sample (visits with T14 error > 1 min dropped)'],
    ['fit', 'description', 'visits', 'quantity', 'value', 'error', 'delta_chi2_const', 'note'], rows)
