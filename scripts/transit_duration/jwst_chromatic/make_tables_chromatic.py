"""Write the tables of section 6.6.11 (JWST per-detector and per-wavelength durations of L 98-59 d) from the outputs of
run_chromatic.sh: data/jwst_chromatic_307210830.csv (per-bin and per-detector fits), data/jwst_chromatic_summary_307210830.csv
(derived quantities) and data/jwst_inventory_307210830.csv (MAST inventory). Own code.
Usage (repository root): python scripts/transit_duration/jwst_chromatic/make_tables_chromatic.py [fits_dir]
Inputs in $TDC_WORK (default work/transit_duration_chromatic): fits/*.jsonl (c3_binfit.py; or fits_dir), out/chromatic.json
(c6_analyze.py), out/sysvariants.json (c9), out/white_diag.json (c7), out/jwst_summary.json (c11), out/predictions.json (m5),
out/mast_inventory.csv (m1)."""
import csv, json, os, sys, numpy as np
H = os.environ.get('TDC_WORK', 'work/transit_duration_chromatic')
FD = sys.argv[1] if len(sys.argv) > 1 else os.path.join(H, 'fits')
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA = os.path.join(ROOT, 'data')
R = json.load(open(H + '/out/chromatic.json')); S = json.load(open(H + '/out/sysvariants.json'))
W = json.load(open(H + '/out/white_diag.json')); M = json.load(open(H + '/out/jwst_summary.json')); P = json.load(open(H + '/out/predictions.json'))
VIS = {'d1224': ('d', 'GTO 1224 obs. 1', 2023.481), 'd4098': ('d', 'GO 4098 obs. 2', 2024.195), 'd1201': ('d', 'GTO 1201 obs. 311', 2022.992),
       'c2512a': ('c', 'GO 2512 obs. 15', 2023.531), 'c2512b': ('c', 'GO 2512 obs. 16', 2023.520)}
LDN = {'prior': 'Kurucz prior (sigma 0.1)', 'free': 'free', 'ramp': 'Kurucz prior + exponential ramp (15 min)'}
def r(x, n=4): return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)

# 1. per-bin and per-detector fits
head = ['fit', 'planet', 'visit', 'year', 'channel', 'ld', 'row', 'lam_lo_um', 'lam_hi_um', 'ncol', 'T14_min', 'T14_err_min', 'b', 'b_err', 'k', 'k_err',
        'u1', 'u1_err', 'u2', 'u2_err', 'sigma_ppm', 'beta', 'n_over_tau', 'used', 'T14_wmean_min', 'T14_wmean_err_min', 'dchi2_free_b_vs_shared_b', 'dof', 'p_shared_b',
        'slope_T14_min_per_um', 'slope_err_min_per_um']
rows = []
for key, v in R.items():
    if key.startswith('between'): continue
    vis, ldm = key.rsplit('_', 1); pl, vname, yr = VIS[vis]
    ncol = {}
    for det in ('nrs1', 'nrs2', 'soss1', 'soss2'):
        f = os.path.join(FD, '%s_%s_%s.jsonl' % (vis, det, ldm))
        if os.path.exists(f):
            for s in open(f): j = json.loads(s); ncol[(det.upper(), j['lo'])] = j['ncol']
    for x in v['bins']:
        rows.append([key, pl, vname, yr, x['det'], ldm, 'bin', x['lo'], x['hi'], ncol.get((x['det'], x['lo']), ''), r(x['T14'][0], 3), r(x['T14'][1], 3), r(x['b'][0], 5), r(x['b'][1], 5),
                     r(x['k'][0], 5), r(x['k'][1], 5), r(x['u1'][0], 3), r(x['u1'][1], 3), r(x['u2'][0], 3), r(x['u2'][1], 3), r(x['sig_ppm'], 1), r(x['beta'], 3), r(x['n_over_tau'], 1),
                     int(x['used'])] + [''] * 7)
    for det, d in v['per_det'].items():
        lo = min(x['lo'] for x in v['bins'] if x['det'] == det and x['used']); hi = max(x['hi'] for x in v['bins'] if x['det'] == det and x['used'])
        rows.append([key, pl, vname, yr, det, ldm, 'channel, shared b', lo, hi, '', r(d['T14_shared'][0], 3), r(d['T14_shared'][1], 3), r(d['b_shared'][0], 5), r(d['b_shared'][1], 5)] + [''] * 9 + [1,
                     r(d['T14_wmean'][0], 3), r(d['T14_wmean'][1], 3), r(d['dchi2_within'], 2), d['dof_within'], r(d['p_within'], 3), r(d['slope_T14_per_um'][0], 3), r(d['slope_T14_per_um'][1], 3)])
    s = v['shared_all']
    rows.append([key, pl, vname, yr, 'all used bins', ldm, 'visit, shared b', '', '', '', r(s['T14'], 3), '', r(s['b'][0], 5), r(s['b'][1], 5)] + [''] * 9 + [1, '', '', r(s['dchi2_free_vs_shared'], 2), s['dof'], r(s['p'], 3), '', ''])
hdr = ['# L 98-59 d (TIC 307210830), section 6.6.11: JWST transit fits per wavelength bin and per detector (own spectroscopic light curves; scripts/transit_duration/jwst_chromatic/c1_extract_cols.py, c4_soss_bins.py, c3_binfit.py, c6_analyze.py; run 2026-10-21)',
       '# fit: <visit>_<ld>; visits d1224 = GTO 1224 obs. 1 (NIRSpec G395H, 2023.48), d4098 = GO 4098 obs. 2 (NIRSpec G395H, 2024.20), d1201 = GTO 1201 obs. 311 (NIRISS SOSS, 2022.99), c2512a/b = GO 2512 obs. 15/16 (planet c, NIRSpec NRS1 only)',
       '# channel: NRS1 (2.87-3.71 um, 8 bins of 0.105 um), NRS2 (3.83-5.15 um, 8 bins of 0.165 um), SOSS1 (order 1, 0.9-2.8 um, 8 bins), SOSS2 (order 2, 2 bins, unconstrained, used = 0); bin edges identical for every visit',
       '# model per bin: batman quadratic LD, P, mid-time and a/R* fixed per visit at the MAP of the joint white-light fit of section 6.6.10 (d: a/R* 36.25; c: 22.14); k, b, u1, u2 free; quadratic in time + trace x/y shift and width regressors (SOSS: background and order-1 tilt), profiled; first 10 min removed, 5-integration bins; white noise sigma x beta',
       '# ld: prior = Gaussian sigma 0.1 about Kurucz ATLAS9 values ([M/H] = -0.5, extrapolated to 3415 K; c2_ld.py), free = uniform (-0.5 < u1, u2 < 1, 0 < u1 + u2 < 1), ramp = prior + exponential ramp (tau 15 min) from the exposure start (GTO 1224 only)',
       '# row = bin: posterior median and standard deviation (emcee, 20 walkers, 1500 steps, 3000 for d1224_nrs1_prior, first third discarded); T14 = first to fourth contact; n_over_tau = post-burn steps / max autocorrelation time',
       '# row = channel, shared b / visit, shared b: one b for all used bins of the channel / visit from the summed per-bin profile chi2(b) (k, u1, u2 re-optimised per bin; the per-bin objective includes the Gaussian limb-darkening penalty when ld = prior or ramp), b and its approximate error from a local quadratic approximation of the summed profile; T14 at that b with the mean of the per-bin MAP k; T14_err transforms the b error with that mean k held fixed and does not propagate the full k-b covariance; P, mid-time and a/R* held fixed (cancellation of their uncertainties in wavelength or epoch differences not demonstrated); absolute T14 on a different geometry scale from section 6.6.10',
       '# T14_wmean = weighted mean of the per-bin T14; dchi2_free_b_vs_shared_b, dof, p_shared_b = test of one b against one b per bin (profile-based p-values: nominal, conditional diagnostics); slope = weighted line of the per-bin T14 against wavelength. For planet c the b grid step (0.015) is coarse against sigma_b (0.007) and the shared-b chi2 values are not reliable']
with open(os.path.join(DATA, 'jwst_chromatic_307210830.csv'), 'w', newline='') as fo:
    fo.write('\n'.join(hdr) + '\n'); w = csv.writer(fo); w.writerow(head); w.writerows(rows)

# 2. derived quantities
out = []
def add(group, quantity, ld, value, error='', note=''): out.append([group, quantity, ld, r(value, 4) if value != '' else '', r(error, 4) if error != '' else '', note])
for vis in ('d1224', 'd4098'):
    for ldm in ('prior', 'free', 'ramp'):
        v = R.get('%s_%s' % (vis, ldm))
        if not v: continue
        g = 'per-detector, %s' % VIS[vis][1]
        add(g, 'NRS1 minus NRS2, shared-b T14 (min)', ldm, *v['det_diff_T14_shared']); add(g, 'NRS1 minus NRS2, weighted-mean T14 (min)', ldm, *v['det_diff_T14_wmean'])
        m = v['model_compare']
        add(g, 'delta chi2, step (one b per detector) vs constant b', ldm, m['dchi2_step_vs_const'], note='approximate shared-b profiles (local quadratic minima)')
        # c6_analyze.py evaluates the models with b linear in wavelength on linearly interpolated profiles and the constant
        # and step models with local quadratic minima; the mixed comparisons are withheld (values not exported)
        for lab in ('delta chi2, b linear in wavelength vs constant b', 'delta chi2, step added to the linear model', 'delta chi2, linear term added to the step model'):
            add(g, lab, ldm, '', note='Withheld: inconsistent profile approximations produce negative improvements for nested models.')
        add(g, 'adjacent bins across the detector gap, NRS1 3.605-3.710 minus NRS2 3.830-3.995 um, T14 (min)', ldm, *v['adjacent']['dT14'])
        add(g, 'within-detector lines extrapolated to %.2f um, NRS1 minus NRS2 (min)' % v['extrapolated_to_gap']['lam_gap'], ldm, *v['extrapolated_to_gap']['diff'])
for ldm in ('prior', 'free'):
    D = R['between_%s' % ldm]; g = 'between NIRSpec visits, GO 4098 obs. 2 minus GTO 1224 obs. 1'
    for det, d in D['per_det'].items():
        add(g, '%s, weighted mean of the per-bin delta T14 (min)' % det, ldm, *d['dT14_wmean'], note='chi2 %.1f for %d dof' % (d['chi2'], d['dof']))
        add(g, '%s, slope of the per-bin delta T14 (min/um)' % det, ldm, *d['slope_per_um']); add(g, '%s, delta of the shared-b T14 (min)' % det, ldm, *d['dT14_shared'])
    a = D['all']
    add(g, 'all 16 bins, weighted mean delta T14 (min)', ldm, *a['dT14_wmean'], note='chi2 about a constant %.1f for %d dof, p = %.3f' % (a['chi2_const'], a['dof'], a['p_const']))
    add(g, 'all 16 bins, step NRS2 minus NRS1 (min)', ldm, *a['step_NRS2_minus_NRS1'], note='chi2 %.1f' % a['chi2_step'])
    add(g, 'all 16 bins, linear in wavelength (min/um)', ldm, *a['slope_per_um'], note='chi2 %.1f' % a['chi2_linear'])
    add(g, 'all 16 bins, step given the linear term (min)', ldm, *a['step_given_linear'], note='chi2 %.1f' % a['chi2_step_plus_linear'])
    for q in ('wmean', 'shared'):
        A, B = R['c2512a_%s' % ldm]['per_det']['NRS1']['T14_%s' % q], R['c2512b_%s' % ldm]['per_det']['NRS1']['T14_%s' % q]
        add('planet c control, NRS1', 'GO 2512 obs. 16 (2023-07-09) minus obs. 15 (2023-07-13), %s T14 (min)' % ('weighted-mean' if q == 'wmean' else 'shared-b'), ldm, B[0] - A[0], float(np.hypot(A[1], B[1])))
for ldm in ('prior', 'free'):
    o = M[ldm]; g = 'detector systematics and averages'
    add(g, 'change of NRS1 minus NRS2 between the visits (min)', ldm, *o['interaction'], note='formal errors only')
    add(g, 'empirical systematic floor per detector and visit (min)', ldm, o['sigma_sys_per_detector'], note='sqrt(mean(d_v^2 - s_v^2) / 2) over the two visits; plug-in estimate from two differences, rough uncertainty about 50 % (not a calibrated interval); assumes zero-mean offsets independent between detectors and visits with a common variance (applied to NIRISS by assumption); errors that include the floor hold it fixed')
    for v, yr in (('1224', 2023.481), ('4098', 2024.195)):
        va = o['visit_avg'][v]; add(g, 'detector-averaged T14, %s (min)' % yr, ldm, va[0], va[2], note='statistical error %.2f; error column includes the floor / sqrt(2)' % va[1])
    pa = o['pair_dT14_avg']; add(g, 'detector-averaged delta T14, 2024.195 minus 2023.481 (min)', ldm, pa[0], pa[2], note='statistical error %.2f; TESS trend predicts -0.85 +- 0.08' % pa[1])
    for lab, d in o['slopes'].items():
        for k, s in d.items(): add('JWST-only slopes (%s)' % lab, k + ' (min/yr)', ldm, *s, note='NIRISS 1201 point: SOSS1 shared-b T14 with the full floor added' if 'NIRISS' in k and lab != 'formal' else '')
g = 'systematics-model variants (c9_sysvariants.py; white light, Kurucz-prior LD, least squares)'
for v in ('1224', '4098'):
    for det in ('NRS1', 'NRS2', 'NRS2pipe'):
        for name, x in S[v][det].items():
            add(g, '%s %s, %s: T14 (min)' % (VIS['d' + v][1], 'NRS2 pipeline (stage 3)' if det == 'NRS2pipe' else det, name), 'prior', *x['T14'])
for v in ('1224', '4098'):
    o = W[v]; g = 'white-light diagnostics, %s (c7_white_diag.py)' % VIS['d' + v][1]
    add(g, 'delta chi2, common b minus one b per detector', 'prior', o['dchi2_common_minus_separate'])
    add(g, 'one b per detector: NRS1 T14 (min)', 'prior', o['sep']['T14_1']); add(g, 'one b per detector: NRS2 T14 (min)', 'prior', o['sep']['T14_2'])
    for i, c in enumerate(o['contacts_min_from_start']): add(g, 'contact %d, minutes after the exposure start' % (i + 1), 'prior', c)
    q = o['ratio']; add(g, 'NRS1/NRS2 ratio light curve: part of the per-detector T14 difference supported (min)', 'prior', *q['T14diff_supported'], note='of %.2f min; red-noise factor %.2f' % (q['T14diff_sep'], q['red_beta']))
    add(g, 'NRS1/NRS2 ratio: delta chi2, common minus per-detector geometry', 'prior', q['dchi2_common_minus_sep'])
    for wn in ('pre', 'ingress', 'bottom', 'egress', 'post'):
        add(g, 'NRS1/NRS2 ratio minus common-geometry model, mean in window %s (ppm)' % wn, 'prior', q['window_mean_minus_common_ppm'][wn], q['window_mean_ppm'][wn][1])
    for det, s in o['spot_box'].items(): add(g, 'largest 10-min in-transit residual excess, %s (ppm)' % det.upper(), 'prior', s['max_box_ppm'], s['sigma_box_ppm'])
    n = sum(1 for d in o['resid_corr_5min'].values() for x in d.values() if x[1] is not None); k = sum(1 for d in o['resid_corr_5min'].values() for x in d.values() if x[1] is not None and x[1] < 0.05)
    lst = '; '.join('%s vs %s r = %.2f, p = %.3f' % (nm, kk, x[0], x[1]) for nm, d in o['resid_corr_5min'].items() for kk, x in d.items() if x[1] is not None and x[1] < 0.05)
    add(g, 'residual-regressor correlation tests (5-min bins, 2000 permutations) with p < 0.05', 'prior', k, note='of %d tests%s' % (n, ('; ' + lst) if lst else ''))
for v in ('1224', '4098'):
    a = S[v]['NRS2']['poly2+regs (section 6.6.10)']['T14'][0]; a0 = S[v]['NRS2']['poly2']['T14'][0]; b = S[v]['NRS2pipe']['poly2']['T14'][0]
    add('own minus pipeline NRS2 white light, %s' % VIS['d' + v][1], 'own (quadratic + regressors) minus pipeline (quadratic) T14 (min)', 'prior', a - b)
    add('own minus pipeline NRS2 white light, %s' % VIS['d' + v][1], 'own minus pipeline, both quadratic only, T14 (min)', 'prior', a0 - b)
for tr in P['transits']:
    for sc, d in tr['scales'].items():
        add('prediction, GO 10498 %s, mid-transit %s UTC (%.3f)' % (tr['label'], tr['tc_utc'], tr['year']), 'T14 on the scale %s, TESS trend continued (min)' % sc, '', *d['trend'])
        add('prediction, GO 10498 %s, mid-transit %s UTC (%.3f)' % (tr['label'], tr['tc_utc'], tr['year']), 'T14 on the scale %s, no change (min)' % sc, '', *d['none'])
hdr = ['# L 98-59 d (TIC 307210830), section 6.6.11: derived quantities of the JWST per-detector and per-wavelength duration analysis (scripts/transit_duration/jwst_chromatic/c6_analyze.py, c7_white_diag.py, c9_sysvariants.py, c11_summary.py, m5_predictions.py; run 2026-10-21)',
       '# ld: prior = Kurucz-prior limb darkening, free = free limb darkening, ramp = prior + exponential ramp; value and error in the unit given in the quantity; errors are formal (statistical) unless the note says otherwise; shared-b T14 values and their errors are approximate and conditional (header of data/jwst_chromatic_307210830.csv, section 6.6.11)',
       '# systematics-model variants of the pipeline NRS2 white light (NRS2pipe): the product has no regressors, so the variants named +regs use the polynomial (and ramp) terms only',
       '# the detector averages and the systematic floor were defined after the NRS1-NRS2 discrepancy had been seen (section 6.6.11); predictions: TESS trend -1.185 +- 0.112 min/yr from the GO 4098 visit (2024.195), Cadieux et al. (2025) linear ephemeris']
with open(os.path.join(DATA, 'jwst_chromatic_summary_307210830.csv'), 'w', newline='') as fo:
    fo.write('\n'.join(hdr) + '\n'); w = csv.writer(fo); w.writerow(['group', 'quantity', 'ld', 'value', 'error', 'note']); w.writerows(out)

# 3. MAST inventory
inv = os.path.join(H, 'out', 'mast_inventory.csv')
if os.path.exists(inv):
    hdr = ['# L 98-59 (TIC 307210830): JWST observations in MAST (CAOM) within 3 arcmin, by 12 target-name variants and 6 proposal ids (scripts/transit_duration/jwst_chromatic/m1_mast_inventory.py), queried 2026-10-21 about 16:50 UTC',
           '# t_min, t_max, t_obs_release: MJD (UTC); start_utc, end_utc, release_utc: the same in UTC; rows with calib_level -1 and no times are planned observations; dataRights as reported by MAST at the query time']
    with open(os.path.join(DATA, 'jwst_inventory_307210830.csv'), 'w') as fo: fo.write('\n'.join(hdr) + '\n' + open(inv).read())
print('tables written to', DATA, len(rows), len(out))
