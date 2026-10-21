"""Write the tables of README sections 6.6 and 6.7 from the outputs of run_l9859.sh and run_toi620.sh.
Usage: python scripts/transit_duration/make_tables_td.py [--work $TD_WORK] [--out data] [--date YYYY-MM-DD]
Writes data/durations_307210830.csv, duration_slopes_307210830.csv, groupfit_307210830.csv, transit_times_307210830.csv,
jwst_fits_307210830.csv, jwst_differences_307210830.csv, durations_296739893.csv, duration_slopes_296739893.csv,
transit_times_296739893.csv and diffimage_296739893.csv (formats in the comment headers)."""
import argparse, csv, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--work', default=os.environ.get('TD_WORK', os.path.join('work', 'transit_duration')))
ap.add_argument('--psf-work', default=None, help='directory holding tpf620/psf_decomposition.json (default: --work)')
ap.add_argument('--out', default='data')
ap.add_argument('--date', default='2026-10-21')
A = ap.parse_args()
W = A.work


def year(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451544.5) / 365.25


def J(name):
    return json.load(open(os.path.join(W, name)))


def write(name, comments, head, rows):
    with open(os.path.join(A.out, name), 'w', newline='') as f:
        for c in comments:
            f.write('# ' + c + '\n')
        w = csv.writer(f, lineterminator='\n'); w.writerow(head)
        for r in rows:
            w.writerow(['' if (isinstance(x, float) and not np.isfinite(x)) else x for x in r])
    print('wrote', name, len(rows), 'rows')


def f(x, n):
    return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else ('%.' + str(n) + 'f') % x


def regress_rows(target, planet, variant, jf):
    R = jf['regress']; n = sum(1 for u in jf['units'] if not u['edge_b'])
    out = []
    for q, nd in (('T14', 3), ('b', 5), ('w', 4)):
        r = R[q]
        out.append([target, planet, variant, q, n, f(r['slope'], nd), f(r['err_scaled'], nd), f(r['z_scaled'], 2), f(r['chi2_lin'], 2), f(r['chi2_const'], 2), r['dof'],
                    f(r['ols_slope'], nd), f(r['boot_err'], nd), f(r['boot_z'], 2)])
    return out


def unit_rows(planet, variant, jf):
    k, aR, b0 = jf['template'][:3]; rows = []
    for u in jf['units']:
        T14e = max(0.05, (u['T14_hi'] - u['T14_lo']) / 2)
        rows.append([planet, variant, u['unit'] if len(u['sectors']) == 1 and jf['config'][5] == 'sector' else '', ' '.join(str(s) for s in u['sectors']), f(u['year'], 4), u['ntr'],
                     f(u['T14'], 2), f(u['T14_lo'], 2), f(u['T14_hi'], 2), f(T14e, 2), f(u['b'], 4), f(u['b_lo'], 4), f(u['b_hi'], 4), f(u['w'], 3), f(u['w_lo'], 3), f(u['w_hi'], 3),
                     f(u['A'], 3), f(u['dchi2_transit'], 1), f(k, 4), f(aR, 3), f(b0, 4), 0 if u['edge_b'] else 1])
    return rows


UHEAD = ['planet', 'variant', 'sector', 'sectors', 'year', 'n_transits', 'T14_min', 'T14_lo_min', 'T14_hi_min', 'T14_err_min', 'b', 'b_lo', 'b_hi', 'w', 'w_lo', 'w_hi',
         'amplitude', 'dchi2_transit', 'k_template', 'aR_template', 'b_template', 'in_regression']
UCOM = ['per sector: T14 from a fit of b alone (k, a/R*, limb darkening fixed at the template); T14 = first to fourth contact, circular orbit; T14_lo, T14_hi and b_lo, b_hi = the',
        'range with chi2 - chi2_min <= s^2, s^2 = max(1, chi2_min/dof); T14_err_min = max(0.05, (T14_hi - T14_lo)/2), used in the regressions; w = time-stretch factor of the template',
        'with free amplitude (shape-agnostic width; 1 = template); amplitude = depth scale of the stretched template; dchi2_transit = chi2(no transit) - chi2(best b)',
        'in_regression = 0 where the best b is at the grid edge (excluded from the regressions)']

# ---------------------------------------------------------------- L 98-59
VAR_D = [('2m_pdc_ttv', 'ps_d_2m_pdc_ttv.json', 'baseline: 2-min PDCSAP, transits stacked on the sinusoidal TTV model'),
         ('2m_pdc_pt', 'ps_d_2m_pdc_pt.json', '2-min PDCSAP, per-transit mid-times'),
         ('2m_pdc_lin', 'ps_d_2m_pdc_lin.json', '2-min PDCSAP, linear ephemeris'),
         ('2m_pdc_sec', 'ps_d_2m_pdc_sec.json', '2-min PDCSAP, linear ephemeris plus a mean offset per sector'),
         ('2m_sap_ttv', 'ps_d_2m_sap.json', '2-min SAP, TTV model'),
         ('20s_pdc_ttv', 'ps_d_20s_pdc.json', '20-s PDCSAP (sectors 28-98), TTV model'),
         ('2m_pdc_ttv_nobc', 'ps_d_clean.json', '2-min PDCSAP, TTV model, only windows without points of b or c')]
CTRL = [('b', '2m_pdc_lin', 'ps_b_2m_pdc.json', 'planet b, 2-min PDCSAP, linear ephemeris (template k = 0.0243, a/R* = 15.15, b = 0.51)'),
        ('c', '2m_pdc_ttv', 'ps_c_2m_pdc.json', 'planet c, 2-min PDCSAP, TTV model (template k = 0.0393, a/R* = 21.05, b = 0.41)'),
        ('d_injected_phase0.30', '2m_pdc_ttv', 'ps_inj30.json', 'copy of d with constant b = 0.915 injected at phase 0.30 of d, full chain'),
        ('d_injected_phase0.62', '2m_pdc_ttv', 'ps_inj62.json', 'copy of d with constant b = 0.915 injected at phase 0.62 of d, full chain')]
rows = []; srows = []
for v, fn, _ in VAR_D:
    jf = J(fn); rows += unit_rows('d', v, jf); srows += regress_rows('L 98-59', 'd', v, jf)
for p, v, fn, _ in CTRL:
    jf = J(fn); rows += unit_rows(p, v, jf); srows += regress_rows('L 98-59', p, v, jf)
com = ['L 98-59 = TIC 307210830: per-sector transit duration T14, impact parameter b and width w of planet d, of the controls b and c and of injected copies of d',
       'scripts/transit_duration/run_l9859.sh (v6_persector.py), SPOC light curves (MAST), ' + A.date,
       'template of d (all variants and injections): k = 0.0453, a/R* = 34.80, b = 0.9257, u = (0.1895, 0.4604), the constant-geometry common fit of the 2-min PDCSAP transits (P = 7.450729 d)'] + UCOM + \
      ['variants: ' + '; '.join('%s = %s' % (v, d) for v, _, d in VAR_D)] + ['controls: ' + '; '.join('%s = %s' % (p, d) for p, _, _, d in CTRL)]
write('durations_307210830.csv', com, UHEAD, rows)
SHEAD = ['target', 'planet', 'variant', 'quantity', 'n_sectors', 'slope_per_yr', 'slope_err_per_yr', 'z', 'chi2_linear', 'chi2_constant', 'dof', 'ols_slope_per_yr', 'bootstrap_err_per_yr', 'bootstrap_z']
SCOM = ['weighted linear fit against year (reference 2022.0) of the per-sector values with in_regression = 1; errors T14_err_min, (b_hi - b_lo)/2, (w_hi - w_lo)/2 (floors 0.05 min, 1e-4, 1e-3);',
        'slope_err = formal error times sqrt(max(1, chi2_linear/dof)); z = slope/slope_err; chi2_constant = chi2 about the weighted mean; ols_slope, bootstrap_err = unweighted fit and its',
        'scatter over 2000 bootstrap resamples of the sectors (resamples spanning less than 1 yr skipped); units: T14 in min/yr, b and w per yr']
write('duration_slopes_307210830.csv', ['L 98-59 = TIC 307210830: slopes of the per-sector T14, b and w of data/durations_307210830.csv, ' + A.date] + SCOM, SHEAD, srows)

# group fits
rows = []
for lab, fn in (('free', 'grp_d_2m_pdc_ttv_prior_year.json'), ('33.6', 'grp_d_aR33.6.json')):
    G = J(fn); g = G['groups']; ng = len(g['year'])
    xA = G['A']['x']; xF = G['F']['x']
    for j in range(ng):
        rows.append([lab, f(g['year'][j], 3), ' '.join(str(s) for s in g['sectors'][j]), g['ntr'][j], f(xA[0], 5), f(xA[1], 2), f(xA[2 + j], 4), f(G['A']['T14_min'][j], 2),
                     f(xF[3 * j], 5), f(xF[3 * j + 1], 2), f(xF[3 * j + 2], 4), f(G['F']['T14_min'][j], 2), f(G['F']['T23_min'][j], 2),
                     f(G['dchi2_A_vs_C'], 1), f(G['dchi2_L_vs_C'], 1), f(G['dchi2_F_vs_A'], 1), f(G['L']['x'][3], 5)])
G0 = J('grp_d_2m_pdc_ttv_prior_year.json')
write('groupfit_307210830.csv',
      ['L 98-59 d: joint fits of the 2-min PDCSAP transits in four epoch groups (scripts/transit_duration/v5_groups.py), ' + A.date,
       'transits stacked on the sinusoidal TTV model (period %.0f d, semi-amplitude %.2f min; chi2 of the per-transit times %.1f linear, %.1f with the sinusoid, %d transits);' % (G0['ttv']['Pttv'], G0['ttv']['amp_min'], G0['ttv']['chi2_lin'], G0['ttv']['chi2_ttv'], G0['ttv']['n']),
       'quadratic limb darkening with a Gaussian prior (0.17, 0.44) +- 0.1; per-transit quadratic baseline profiled; aR_mode = free (common a/R* fitted) or 33.6 (a/R* fixed by a prior 33.6 +- 0.02)',
       'model A: common k and a/R*, b per group (k_A, aR_A, b_A, T14_A); model F: k, a/R* and b per group (k_F, aR_F, b_F, T14_F, T23_F); for aR_mode 33.6 only k and b are free in F',
       'dchi2_A_vs_C = chi2(one geometry) - chi2(A), 3 extra parameters; dchi2_L_vs_C: b linear in time (1 extra parameter, slope db_dt_L per yr); dchi2_F_vs_A: k and a/R* per group'],
      ['aR_mode', 'year', 'sectors', 'n_transits', 'k_A', 'aR_A', 'b_A', 'T14_A_min', 'k_F', 'aR_F', 'b_F', 'T14_F_min', 'T23_F_min', 'dchi2_A_vs_C', 'dchi2_L_vs_C', 'dchi2_F_vs_A', 'db_dt_L'], rows)


def ttv_times(T, Pmin=300, Pmax=500):
    rows_ = T['rows']; ep = rows_[:, 1]; tc = rows_[:, 2]; e = rows_[:, 3]; ok = np.isfinite(e) & (rows_[:, 6] == 0); wt = 1 / e[ok] ** 2
    Al = np.vstack([np.ones(ok.sum()), ep[ok]]).T; lin = np.linalg.solve(Al.T @ (Al * wt[:, None]), Al.T @ (wt * tc[ok])); tl = lin[0] + lin[1] * ep
    def fit(Pt):
        B = np.vstack([np.ones(ok.sum()), ep[ok], np.sin(2 * np.pi * tc[ok] / Pt), np.cos(2 * np.pi * tc[ok] / Pt)]).T
        c = np.linalg.solve(B.T @ (B * wt[:, None]), B.T @ (wt * tc[ok])); r = tc[ok] - B @ c; return np.sum(r * r * wt), c
    grid = np.linspace(Pmin, Pmax, 401); Pt = grid[int(np.argmin([fit(g)[0] for g in grid]))]; cs = fit(Pt)[1]
    tt = cs[0] + cs[1] * ep + cs[2] * np.sin(2 * np.pi * tc / Pt) + cs[3] * np.cos(2 * np.pi * tc / Pt)
    return rows_, ok, lin, tl, Pt, cs, tt


T = np.load(os.path.join(W, 'tim_L9859_d_2m_pdc.npz'), allow_pickle=True)
rr, ok, lin, tl, Pt, cs, tt = ttv_times(T)
rows = [[int(r[0]), int(r[1]), f(year(r[2]), 4), f(r[2], 6), f(r[3] * 1440, 2), f((r[2] - tl[i]) * 1440, 2), f((tt[i] - tl[i]) * 1440, 2), f(r[4], 1), int(ok[i])] for i, r in enumerate(rr)]
write('transit_times_307210830.csv',
      ['L 98-59 d: per-transit mid-times from the 2-min PDCSAP light curves with the common shape (scripts/transit_duration/v4_timing.py), ' + A.date,
       'linear ephemeris (weighted, transits with determined errors): T0 = BTJD %.6f, P = %.7f d (epoch 0 = BTJD 1362.74002); TTV model = linear ephemeris plus a sinusoid of period %.1f d' % (lin[0], lin[1], Pt),
       'oc_min = tc - linear ephemeris; ttv_model_min = TTV model - linear ephemeris; dchi2 = chi2(no transit) - chi2(transit) of the window; used = 1 if the error is determined and the best shift is not at the +-20-min edge'],
      ['sector', 'epoch', 'year', 'tc_btjd', 'tc_err_min', 'oc_min', 'ttv_model_min', 'dchi2', 'used'], rows)

# JWST
VIS = {'d_1201': ('d', '1201', 'GTO', 'D. Lafreniere', 311, 'NIRISS SOSS', 'jw01201-o311_t019_niriss_clear-gr700xd-substrip96_whtlt.ecsv', 'order 1'),
       'd_1224': ('d', '1224', 'GTO', 'S. Birkmann', 1, 'NIRSpec G395H', 'jw01224-o001_t009_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'd_4098b': ('d', '4098', 'GO', 'B. Benneke', 2, 'NIRSpec G395H', 'jw04098-o002_t011_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'c_1201': ('c', '1201', 'GTO', 'D. Lafreniere', 302, 'NIRISS SOSS', 'jw01201-o302_t019_niriss_clear-gr700xd-substrip96_whtlt.ecsv', 'order 1'),
       'c_2512a': ('c', '2512', 'GO', 'N. Batalha', 15, 'NIRSpec G395H', 'jw02512-o015_t003_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'c_2512b': ('c', '2512', 'GO', 'N. Batalha', 16, 'NIRSpec G395H', 'jw02512-o016_t003_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'b_3942a': ('b', '3942', 'GO', 'M. Damiano', 1, 'NIRSpec G395H', 'jw03942-o001_t001_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'b_3942b': ('b', '3942', 'GO', 'M. Damiano', 2, 'NIRSpec G395H', 'jw03942-o002_t001_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2'),
       'b_3942c': ('b', '3942', 'GO', 'M. Damiano', 3, 'NIRSpec G395H', 'jw03942-o003_t001_nirspec_f290lp-g395h-s1600a1-sub2048_whtlt.ecsv', 'NRS2')}
FITS = [('jfit_d_1201.json', 'single', 'free', 'nominal'), ('jfit_d_1224.json', 'single', 'free', 'nominal'), ('jfit_d_4098b.json', 'single', 'free', 'nominal'),
        ('jfit_d_1201_aR33.6.json', 'single', '33.6', 'nominal'), ('jfit_d_1224_aR33.6.json', 'single', '33.6', 'nominal'), ('jfit_d_4098b_aR33.6.json', 'single', '33.6', 'nominal'),
        ('jfit_d_nrs2_joint.json', 'joint', 'free', 'nominal'), ('jfit_d_nrs2_joint_t10.json', 'joint', 'free', 'first 10 min removed'),
        ('jfit_d_nrs2_joint_mim.json', 'joint', 'free', 'GO 4098 cut to the pre-ingress coverage of GTO 1224'), ('jfit_d_nrs2_joint_aR33.6.json', 'joint', '33.6', 'nominal'),
        ('jfit_c_1201.json', 'single', 'free', 'nominal'), ('jfit_c_2512_joint.json', 'joint', 'free', 'nominal'), ('jfit_b_3942_joint.json', 'joint', 'free', 'nominal')]
rows = []; drows = []
for fn, mode, arm, var in FITS:
    R = J(fn)
    for i, n in enumerate(R['names']):
        key = n.replace('t10', '').replace('mim', ''); p, prog, kind, pi, obs, inst, prod, det = VIS[key]; r = R[n]
        rows.append([p, prog, kind, pi, obs, inst, det, prod, fn.replace('.json', ''), mode, arm, var, f(r['year'], 3), f(r['tc'][0], 6), f(max(r['tc'][1], r['tc'][2]) * 1440, 2),
                     f(r['T14_min'][0], 2), f(r['T14_min'][1], 2), f(r['T14_min'][2], 2), f(r['T23_min'][0], 2), f(r['b'][0], 4), f(r['b'][1], 4), f(r['b'][2], 4),
                     f(r['k'][0], 4), f(max(r['k'][1], r['k'][2]), 4), f(r['aR'][0], 2), f(r['aR'][1], 2), f(r['aR'][2], 2), f(r['depth_ppm'][0], 0), f(r['u1'][0], 3), f(r['u2'][0], 3),
                     f(R['beta'][i], 2), f(R['s_white_ppm'][i], 0)])
    for k_ in R:
        if k_.startswith('dT14_'):
            pair = k_[5:]; db = R['db_' + pair]
            drows.append([fn.replace('.json', ''), R['names'][0][0], pair.replace('_minus_', ' minus '), arm, var, f(R[k_][0], 2), f(R[k_][1], 2), f(R[k_][2], 2), f(db[0], 4), f(db[1], 4), f(db[2], 4)])
JCOM = ['transit fits of the public JWST stage-3 white-light curves (MAST *_whtlt.ecsv) of L 98-59 b, c and d, refitted here (scripts/transit_duration/j1_prep.py, j2_fit.py), ' + A.date,
        'model: batman, circular orbit, quadratic limb darkening free (q1, q2 uniform), times (quadratic polynomial in time plus an exponential ramp); first 3 min removed; 5-sigma clip about a 31-point running median',
        'errors: white scatter of the best-fit residuals times the red-noise factor beta (median over bins of 5-80 points, at least 1); emcee 40 walkers, 3000 steps, second half kept; values = median, +/- = 84th - 50th and 50th - 16th percentiles',
        'fit = single (one visit, k, a/R*, b free) or joint (common k, a/R*, limb darkening; b, mid-time and ramp per visit); aR_mode = free or 33.6 (fixed); T14 first to fourth contact',
        'programme, PI and observation as listed in MAST; detector: NIRSpec NRS2 (the NRS1 white light of the products used is not usable for these visits), NIRISS SOSS order 1']
write('jwst_fits_307210830.csv', JCOM,
      ['planet', 'programme', 'programme_type', 'pi', 'observation', 'instrument', 'detector', 'product', 'fit', 'mode', 'aR_mode', 'variant', 'year', 'tc_btjd', 'tc_err_min',
       'T14_min', 'T14_plus', 'T14_minus', 'T23_min', 'b', 'b_plus', 'b_minus', 'k', 'k_err', 'aR', 'aR_plus', 'aR_minus', 'depth_ppm', 'u1', 'u2', 'beta', 'white_rms_ppm'], rows)
write('jwst_differences_307210830.csv',
      ['differences between visits in the joint JWST fits of data/jwst_fits_307210830.csv (posterior median and 84th - 50th, 50th - 16th percentiles), ' + A.date,
       'dT14_min = T14(second visit) - T14(first visit); db = b(second) - b(first); names: d_1224 = GTO 1224 (2023-06-25), d_4098b = GO 4098 observation 2 (2024-03-12), c_2512a/b = GO 2512 observations 15/16, b_3942a/b/c = GO 3942 observations 1/2/3'],
      ['fit', 'planet', 'pair', 'aR_mode', 'variant', 'dT14_min', 'dT14_plus', 'dT14_minus', 'db', 'db_plus', 'db_minus'], drows)

# ---------------------------------------------------------------- TOI-620
VAR_T = [('2m_pdc_sec', 'ps_T620_2m_pdc_sec.json', 'baseline: 2-min PDCSAP, linear ephemeris plus a mean offset per sector'),
         ('2m_pdc_pt', 'ps_T620_2m_pdc_pt.json', '2-min PDCSAP, per-transit mid-times'),
         ('2m_sap_pt', 'ps_T620_2m_sap_pt.json', '2-min SAP, per-transit mid-times'),
         ('own_r1.5_sec', 'ps_T620_r15_sec.json', 'own TPF photometry, r = 1.5 px, per-sector offsets'),
         ('own_r3.0_sec', 'ps_T620_r30_sec.json', 'own TPF photometry, r = 3.0 px, per-sector offsets (b-only T14 biased by the neighbour flux; w is not)')]
rows = []; srows = []
for v, fn, _ in VAR_T:
    jf = J(fn); rows += unit_rows('b', v, jf); srows += regress_rows('TOI-620', 'b', v, jf)
tmpl = J('ps_T620_2m_pdc_sec.json')['template']
write('durations_296739893.csv',
      ['TOI-620 = TIC 296739893: per-sector transit duration T14, impact parameter b and width w of TOI-620 b (P = 5.098784 d), scripts/transit_duration/run_toi620.sh (v6_persector.py), ' + A.date,
       'template = common fit of the 2-min PDCSAP transits: k = %.4f, a/R* = %.2f, b = %.4f, u = (%.3f, %.3f) (SAP: its own common fit; own apertures: k = 0.0602, a/R* = 16.64, b = 0.9326)' % tuple(tmpl)] + UCOM +
      ['variants: ' + '; '.join('%s = %s' % (v, d) for v, _, d in VAR_T)], UHEAD, rows)
write('duration_slopes_296739893.csv', ['TOI-620 = TIC 296739893: slopes of the per-sector T14, b and w of data/durations_296739893.csv, ' + A.date] + SCOM, SHEAD, srows)
T = np.load(os.path.join(W, 'tim_T620_2m_pdc.npz'), allow_pickle=True)
rr, ok, lin, tl, _, _, _ = ttv_times(T)
sec = rr[:, 0]; off = {}
for s in np.unique(sec):
    m = (sec == s) & ok; ww = 1 / rr[m, 3] ** 2
    off[s] = (np.sum((rr[m, 2] - tl[m]) * ww) / ww.sum() * 1440, 1440 / np.sqrt(ww.sum()))
pt = {(u['unit']): u for u in J('ps_T620_2m_pdc_sec_tr.json')['units']}
rows = []
for i, r in enumerate(rr):
    u = pt.get(i, None)
    rows.append([int(r[0]), int(r[1]), f(year(r[2]), 4), f(r[2], 6), f(r[3] * 1440, 2), f((r[2] - tl[i]) * 1440, 2), f(off[r[0]][0], 2), f(off[r[0]][1], 2),
                 f(u['T14'], 2) if u else '', f(u['T14_lo'], 2) if u else '', f(u['T14_hi'], 2) if u else '', f(u['b'], 4) if u else '', f(u['w'], 3) if u else '', f(u['w_lo'], 3) if u else '', f(u['w_hi'], 3) if u else '', int(ok[i])])
chi2r = np.sum(((rr[ok, 2] - tl[ok]) / rr[ok, 3]) ** 2) / (ok.sum() - 2)
write('transit_times_296739893.csv',
      ['TOI-620 b: per-transit mid-times (2-min PDCSAP, common shape; scripts/transit_duration/v4_timing.py) and per-transit T14, b, w (v6_persector.py, unit transit, stacked on per-sector offsets), ' + A.date,
       'linear ephemeris (weighted): T0 = BTJD %.6f, P = %.7f d (epoch 0 = BTJD 1518.004514); chi2_red of the O-C = %.2f for %d transits' % (lin[0], lin[1], chi2r, ok.sum()),
       'oc_min = tc - linear ephemeris; sector_offset_min, sector_offset_err_min = weighted mean O-C of the sector; T14, b, w as in data/durations_296739893.csv but per transit'],
      ['sector', 'epoch', 'year', 'tc_btjd', 'tc_err_min', 'oc_min', 'sector_offset_min', 'sector_offset_err_min', 'T14_min', 'T14_lo_min', 'T14_hi_min', 'b', 'w', 'w_lo', 'w_hi', 'used'], rows)
PSF = os.path.join(A.psf_work or W, 'tpf620', 'psf_decomposition.json')
P2 = json.load(open(PSF)) if os.path.exists(PSF) else {}
if not P2: print('no', PSF, '- diffimage_296739893.csv not written')
rows = []
for s in sorted(P2):
    d = P2[s]; m = json.load(open(os.path.join(A.psf_work or W, 'tpf620', s + '.json')))
    rows.append([int(s.lstrip('s')), m['ntransits'], f(d['sep_px'], 2), f(d['flux_ratio_n_t'], 3), f(d['amp_target'], 2), f(d['amp_target_err'], 2), f(d['amp_neighbour'], 2), f(d['amp_neighbour_err'], 2),
                 f(d['share_neighbour'], 3), f(d['free_offset_px'], 2)])
if P2: write('diffimage_296739893.csv',
      ['TOI-620 b: two-source decomposition of the mean out-of-transit minus in-transit difference image of each sector (scripts/transit_duration/p1_tpf.py, p2_psf.py; SPOC 2-min target pixel files), ' + A.date,
       'two Gaussian PSFs (common width) fitted to the median image: the target and the neighbour TIC 296739875; the difference image fitted with the two PSFs plus a constant (amplitudes in e-/s per pixel at the peak)',
       'share_neighbour = fraction of the fitted transit signal (amplitude times PSF sum) on the neighbour; free_offset_px = distance of a free-position single-PSF fit of the difference image from the target'],
      ['sector', 'n_transits', 'separation_px', 'flux_ratio_neighbour_target', 'amp_target', 'amp_target_err', 'amp_neighbour', 'amp_neighbour_err', 'share_neighbour', 'free_offset_px'], rows)
