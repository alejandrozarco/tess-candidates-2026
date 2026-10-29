#!/usr/bin/env python3
"""TOI-421 (TIC 94986319): write the tables of README section 6.12 from the outputs of fit_toi421.py, planet_b_toi421.py,
ffi_toi421.py, nbody_toi421.py (work/toi421/out/), cheops/summarize.py (work/toi421/cheops/summary_cheops.json) and the
JWST fits (work/toi421/jwst/fits/*_var*.json), and the curves of the figures from the recorded predictions
(data/predictions_94986319.json, data/predictions_cheops_jwst_94986319.json).
usage: python scripts/toi421/make_tables_toi421.py   (needs numpy, scipy and batman)
writes data/{sectors,pertransit,trend,injection,planetb,cheops,jwst,nbody,forecast}_94986319.csv
"""
import csv, glob, json, os
import numpy as np
from scipy.stats import norm
import batman

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.path.join(ROOT, 'work', 'toi421'); O = os.path.join(W, 'out'); D = os.path.join(ROOT, 'data')
J = lambda name: json.load(open(os.path.join(O, name + '.json')))
PC, KC, ARC = 16.06750764, 0.05459, 31.91
PB, KB, ARB = 5.19759807, 0.0279, 13.76
YEAR = lambda btjd: 2014.0 + (btjd + 2457000.0 - 2456658.5) / 365.25
LOCAL_Z = 4.5          # adopted local significance of the c slope (section 6.12.4)
NSCREEN = 271          # planets in the screen that selected TOI-421 c


def T14(b, aR, k, P):
    return P / np.pi * np.arcsin(np.sqrt(max((1 + k) ** 2 - b * b, 0)) / np.sqrt(aR * aR - b * b)) * 1440


def boxdepth(b, k=KB, aR=ARB, P=PB):
    t = np.linspace(-30, 30, 601) / 1440.; q = batman.TransitParams(); q.t0 = 0; q.per = P; q.rp = k; q.a = aR
    q.inc = np.degrees(np.arccos(min(b / aR, 1))); q.ecc = 0; q.w = 90; q.u = [0.401, 0.222]; q.limb_dark = 'quadratic'
    return float((1 - batman.TransitModel(q, t).light_curve(q)).mean() * 1e6)


def r(x, n=4):
    return '' if x is None or (isinstance(x, float) and not np.isfinite(x)) else (round(float(x), n) if isinstance(x, (float, np.floating)) else x)


def write(name, comments, head, rows):
    with open(os.path.join(D, name), 'w', newline='') as f:
        for c in comments:
            f.write('# ' + c + '\n')
        w = csv.writer(f, lineterminator='\n'); w.writerow(head)
        for row in rows:
            w.writerow([r(v) for v in row])
    print('wrote', name, len(rows), 'rows')


def T14err(b, eb, aR, k, P):
    h = 1e-4
    return abs(T14(b + h, aR, k, P) - T14(b - h, aR, k, P)) / (2 * h) * eb


def sectors():
    info = J('lc_header_info'); base = J('baseline'); dur = J('duration'); ffi = J('ffi'); pb = J('planet_b')
    tr = {x['sec']: x for x in dur['per_sector']}; Dk = {x['sec']: x for x in base['D']}
    rows = []
    for x in base['B']:
        s = x['sec']; h = info['s%03d_2m' % s]
        rows.append(['c', 'SPOC 2-min PDCSAP', s, x['year'], h['cam'], h['ccd'], 120, h['procver'], h['crowdsap'],
                     len([a for a in base['A'] if a['sec'] == s]), x['b'], x['eb'], x['T14'], T14err(x['b'], x['eb'], ARC, KC, PC), Dk[s]['b'], Dk[s]['eb'], Dk[s]['k'], Dk[s]['ek'],
                     tr[s]['T14'], tr[s]['eT14'], tr[s]['depth_ppm'], tr[s]['edepth']])
    for x in J('cadence')['20s_late_per_sector']:
        s = x['sec']; h = info['s%03d_20s' % s]
        rows.append(['c', 'SPOC 20-s PDCSAP', s, x['year'], h['cam'], h['ccd'], 20, h['procver'], h['crowdsap'], '', x['b'], x['eb'], x['T14'], T14err(x['b'], x['eb'], ARC, KC, PC)] + [''] * 8)
    for R in ('R2.0', 'R2.5', 'R3.0'):
        for x in ffi[R]['c_per_sector']:
            s = x['sec']; h = info['s%03d_2m' % s]
            rows.append(['c', 'own TESScut photometry r = %s px' % R[1:], s, x['year'], h['cam'], h['ccd'],
                         {5: 1800, 6: 1800, 32: 600}.get(s, 200), 'TESScut (no PDC)', '', '', x['b'], x['eb'], x['T14'], T14err(x['b'], x['eb'], ARC, KC, PC)] + [''] * 8)
    for tag, x in pb['per_sector'].items():
        s = int(tag[1:4]); h = info['s%03d_2m' % s]
        rows.append(['b', 'SPOC 2-min PDCSAP', s, x['year'], h['cam'], h['ccd'], 120, h['procver'], h['crowdsap'], x['ntr'],
                     x['x'][0], x['ex'][0], T14(x['x'][0], ARB, KB, PB), T14err(x['x'][0], x['ex'][0], ARB, KB, PB)] + [''] * 8)
    for tag, x in ffi['R2.5']['b_per_sector'].items():
        s = int(tag[1:4]); h = info['s%03d_2m' % s]
        rows.append(['b', 'own TESScut photometry r = 2.5 px', s, x['year'], h['cam'], h['ccd'], {5: 1800, 6: 1800, 32: 600}.get(s, 200),
                     'TESScut (no PDC)', '', x['ntr'], x['x'][0], x['ex'][0], T14(x['x'][0], ARB, KB, PB), T14err(x['x'][0], x['ex'][0], ARB, KB, PB)] + [''] * 8)
    write('sectors_94986319.csv',
          ['TOI-421 = TIC 94986319: impact parameter per TESS sector; scripts/toi421/{fit,planet_b,ffi}_toi421.py, make_tables_toi421.py',
           'planet c: k 0.05459, a/R* 31.91, u 0.401/0.222, P 16.06750764 d; planet b: k 0.0279, a/R* 13.76, same u, P 5.19759807 d; circular; errors from the Jacobian scaled by sqrt(chi2_red) if > 1',
           'T14_model_min = contact 1 to contact 4 of the fixed template at the fitted b; b_kfree/k: planet c with k free per sector; trap_*: trapezoid fit per sector (shared depth, T14, T12; convention-free)',
           'camera/CCD and crowding (CROWDSAP) from the SPOC 2-min headers; own photometry rows: same sectors, cadence of the full-frame images'],
          ['planet', 'product', 'sector', 'year', 'camera', 'ccd', 'cadence_s', 'pipeline', 'crowdsap', 'n_transits', 'b', 'b_err',
           'T14_model_min', 'T14_model_err_min', 'b_kfree', 'b_kfree_err', 'k', 'k_err', 'trap_T14_min', 'trap_T14_err_min', 'trap_depth_ppm', 'trap_depth_err_ppm'], rows)


def pertransit():
    base = J('baseline'); dur = J('duration'); ttv = J('ttv')['linear_ephem']; sp = J('spots')
    tr = {(x['sec'], x['n']): x for x in dur['per_transit']}; res = {(x['sec'], x['n']): x for x in sp['residual_bins']}
    rows = []
    for i, a in enumerate(base['A']):
        k = (a['sec'], a['n']); t = tr[k]; q = res[k]
        rows.append([a['sec'], a['n'], a['year'], a['b'], a['eb'], a['T14'], ttv['oc_min'][i], ttv['eoc_min'][i],
                     t['T14'], t['eT14'], t['T23'], t['eT23'], t['T12'], t['eT12'], t['depth_ppm'], t['edepth'],
                     q['in_max_pos_z'], q['in_mean_z2'], q['out_mean_z2']])
    write('pertransit_94986319.csv',
          ['TOI-421 c = TIC 94986319: the 12 complete TESS SPOC 2-min PDCSAP transits (sectors 5, 6, 32, 98, 107); scripts/toi421/fit_toi421.py',
           'b, b_err, T14_model_min: per-transit fit (b, mid-time, quadratic baseline free; k 0.05459, a/R* 31.91, u 0.401/0.222 fixed); epoch n from T0 BTJD 2195.30635',
           'OC_min: mid-time minus the weighted linear ephemeris of these 12 times (T0 BTJD %.5f, P %.6f d)' % (ttv['T0'], ttv['P']),
           'trap_*: trapezoid fit per transit (depth, T14, T12 free; T23 = T14 - 2 T12); resid_*: residuals of the per-transit fit in 10-min bins (z of the bin mean), in transit and out of transit'],
          ['sector', 'epoch', 'year', 'b', 'b_err', 'T14_model_min', 'OC_min', 'OC_err_min', 'trap_T14_min', 'trap_T14_err_min',
           'trap_T23_min', 'trap_T23_err_min', 'trap_T12_min', 'trap_T12_err_min', 'trap_depth_ppm', 'trap_depth_err_ppm',
           'resid_in_max_pos_z', 'resid_in_mean_z2', 'resid_out_mean_z2'], rows)


LABELS = {
    'SAP_crowd_corrected': 'SAP, own crowding correction with the SPOC CROWDSAP', 'SAP_no_crowd_corr': 'SAP, no crowding correction',
    'PDC_crowd_shift_early-0.01_late+0.01': 'PDCSAP, crowding shifted: early sectors -0.01, late +0.01',
    'PDC_crowd_shift_early--0.01_late+-0.01': 'PDCSAP, crowding shifted: early sectors +0.01, late -0.01',
    'PDC_crowd_shift_early-0.02_late+0.02': 'PDCSAP, crowding shifted: early sectors -0.02, late +0.02',
    'template_0.401_0.222': 'u = 0.401, 0.222 (template)', 'claret_like_0.44_0.20': 'u = 0.44, 0.20', 'low_0.30_0.25': 'u = 0.30, 0.25',
    'high_0.50_0.15': 'u = 0.50, 0.15', 'linear_0.60_0': 'linear law, u = 0.60', 'per_transit_at_global_free_LD': 'u fitted jointly to all transits, then fixed',
    'binned_600s_all': 'all sectors binned to 600 s', 'binned_1800s_all': 'all sectors binned to 1800 s',
    'late_binned_1800_early_2m': 'late sectors binned to 1800 s, early at 120 s', 'early_binned_1800_late_2m': 'early sectors binned to 1800 s, late at 120 s'}


def trend():
    rows = []

    def add(test, variant, s, note='', qty='db/dt (1/yr)'):
        rows.append([test, variant, qty, s['slope'], s['eslope_scaled'], s['z'], s['chi2'], s['dof'], s.get('early'), s.get('e_early'),
                     s.get('late'), s.get('e_late'), s.get('step'), s.get('estep'), s.get('zstep'), note])
    b = J('baseline')
    add('baseline', 'per transit, SPOC 2-min PDCSAP (adopted)', b['A_slope'], '12 transits')
    add('baseline', 'per sector, b free', b['B_slope'])
    add('baseline', 'per sector, b and k free', b['D_slope'])
    add('baseline', 'per sector, k with b free', b['D_k_slope'], qty='dk/dt (1/yr)')
    for k, v in J('loso').items():
        if k.startswith('without'):
            add('leave one sector out', k.replace('_', ' ').replace('S', 'sector '), v)
    a = J('aR')
    for k, v in a.items():
        if k.startswith('fixed_aR_'):
            add('a/R* fixed', 'a/R* = %s' % k[9:], v['slope'])
    g = a['global_free_aR_shared_b']
    rows.append(['a/R* free', 'one a/R* and one b for all 12 transits', 'a/R*', g['aR'], g['eaR'], '', '', '', '', '', '', '', '', '', '', 'b = %.4f +- %.4f' % (g['b'], g['eb'])])
    for p in a['profile_two_epoch']:
        rows.append(['a/R* profile', 'a/R* = %.0f, b per epoch group' % p['aR'], 'b late - b early', p['db'], abs(p['db'] / p['z']), p['z'],
                     '', '', p['b_early'], p['eb_early'], p['b_late'], p['eb_late'], p['db'], abs(p['db'] / p['z']), p['z'], ''])
    for k, v in a.items():
        if k.startswith('per_sector_free_aR'):
            add('a/R* free per sector', {'per_sector_free_aR_prior_29.2_1.5': 'prior a/R* 29.2 +- 1.5', 'per_sector_free_aR_prior_31.9_3.0': 'prior a/R* 31.9 +- 3.0', 'per_sector_free_aR_noprior': 'no prior'}[k], v['slope'],
                'a/R* slope %+.3f +- %.3f /yr' % (v['aR_slope']['slope'], v['aR_slope']['eslope_scaled']))
    l = J('ld')
    for k, v in l.items():
        if isinstance(v, dict) and 'slope' in v and 'z' in v:
            add('limb darkening', LABELS.get(k, k), v)
    e = l['epoch_free_LD_each']
    rows.append(['limb darkening', 'LD free in each epoch group', 'b late - b early', e['late']['b'] - e['early']['b'], np.hypot(e['late']['eb'], e['early']['eb']),
                 (e['late']['b'] - e['early']['b']) / np.hypot(e['late']['eb'], e['early']['eb']), '', '', e['early']['b'], e['early']['eb'], e['late']['b'], e['late']['eb'],
                 '', '', '', 'u early %.2f/%.2f, late %.2f/%.2f' % (e['early']['u1'], e['early']['u2'], e['late']['u1'], e['late']['u2'])])
    dl = J('dilution')
    for k, v in dl.items():
        if isinstance(v, dict) and 'z' in v:
            add('dilution', LABELS.get(k, k), v)
    add('dilution', 'per transit, k free: b', dl['per_transit_free_k']['b'])
    add('dilution', 'per transit, k free: k', dl['per_transit_free_k']['k'], qty='dk/dt (1/yr)')
    du = J('duration')
    for k in ('T14', 'T23', 'T12'):
        add('duration (trapezoid)', 'per transit', du['slope_' + k], qty='d%s/dt (min/yr)' % k)
    add('duration (trapezoid)', 'per transit', du['slope_depth_ppm'], qty='ddepth/dt (ppm/yr)')
    add('duration (trapezoid)', 'per transit, SAP', du['SAP_slope_T14'], qty='dT14/dt (min/yr)')
    ee, ll = du['epoch']['early'], du['epoch']['late']
    for k in ('T14', 'T23', 'T12'):
        d = ll[k] - ee[k]; s = np.hypot(ll['e' + k], ee['e' + k])
        rows.append(['duration (trapezoid)', 'epoch groups 2018-2020 / 2025-2026 (convention-free)', k + ' (min)', d, s, d / s, '', '',
                     ee[k], ee['e' + k], ll[k], ll['e' + k], d, s, d / s, 'joint trapezoid fit of each group'])
    d = ll['depth_ppm'] - ee['depth_ppm']; s = np.hypot(ll['edepth'], ee['edepth'])
    rows.append(['duration (trapezoid)', 'epoch groups 2018-2020 / 2025-2026', 'depth (ppm)', d, s, d / s, '', '', ee['depth_ppm'], ee['edepth'],
                 ll['depth_ppm'], ll['edepth'], d, s, d / s, ''])
    c = J('cadence')
    add('cadence', '20-s light curves for sectors 98 and 107', c['20s_late']['slope'])
    for k in ('binned_600s_all', 'binned_1800s_all', 'late_binned_1800_early_2m', 'early_binned_1800_late_2m'):
        add('cadence', LABELS.get(k, k), c[k])
    sp = J('spots')
    add('stellar activity', 'per transit, 4-sigma clipping', sp['clipped_4sigma'])
    zl = [x['z'] for x in sp['leave_one_transit_out']]
    rows.append(['stellar activity', 'leave one transit out (12 fits)', 'z range', '', '', '', '', '', '', '', '', '', '', '', '',
                 'z from %.2f to %.2f' % (max(zl), min(zl))])
    t = J('ttv')
    add('timing', 'per sector, one mid-time offset per sector on the linear ephemeris', t['per_sector_linear_ephem_fold']['slope'])
    rows.append(['timing', 'linear ephemeris of the 12 per-transit mid-times', 'chi2 of O-C', t['linear_ephem']['chi2'], '', '', t['linear_ephem']['chi2'],
                 t['linear_ephem']['n'] - 2, '', '', '', '', '', '', '', 'O-C from %+.1f to %+.1f min' % (min(t['linear_ephem']['oc_min']), max(t['linear_ephem']['oc_min']))])
    f = J('ffi')
    for R in ('R2.0', 'R2.5', 'R3.0'):
        add('own TESScut photometry (no PDC)', 'per transit, r = %s px' % R[1:], f[R]['c_slope'])
    for k in ('inject', 'inject_trend'):
        s = J(k)['summary']
        rows.append(['injection', 'constant b = 0.748' if k == 'inject' else 'injected trend -0.0051 /yr', 'per-transit slope z',
                     s['slope_mean'], s['slope_std'], s['z_mean'], '', '', '', '', '', '', '', '', '',
                     'n %d offsets; z std %.2f, max abs z %.2f; sector-fit z mean %.2f' % (s['n'], s['z_std'], s['z_absmax'], s['sector_z_mean'])])
    p_loc = 2 * norm.sf(LOCAL_Z); p_glob = 1 - (1 - p_loc) ** NSCREEN
    rows.append(['look-elsewhere', 'adopted local z %.1f over %d screened planets (independent trials, Gaussian tails)' % (LOCAL_Z, NSCREEN), 'global z',
                 norm.isf(p_glob / 2), '', '', '', '', '', '', '', '', '', '', '', 'two-sided p local %.2e, global %.2e' % (p_loc, p_glob)])
    write('trend_94986319.csv',
          ['TOI-421 c = TIC 94986319: rates, epoch-group values and variants; scripts/toi421/fit_toi421.py, ffi_toi421.py, make_tables_toi421.py',
           'slopes: weighted linear fits against year (errors scaled by sqrt(chi2/dof) if > 1); early = 2018-2020 (sectors 5, 6, 32), late = 2025-2026 (98, 107); step = late - early',
           'template unless stated: k 0.05459, a/R* 31.91, u 0.401/0.222, circular; b depends on the adopted a/R*, the trapezoid durations do not'],
          ['test', 'variant', 'quantity', 'value', 'err', 'z', 'chi2', 'dof', 'early', 'early_err', 'late', 'late_err', 'step', 'step_err', 'step_z', 'note'], rows)


def injection():
    rows = []
    for k in ('inject', 'inject_trend'):
        for x in J(k)['rows']:
            rows.append([J(k)['mode'], x['offset'], x['n_early'], x['n_late'], x['slope'], x['z'], x['zstep'], x['sector_slope'], x['sector_z']])
    write('injection_94986319.csv',
          ['TOI-421 c = TIC 94986319: synthetic transits of c injected into the real SPOC 2-min PDCSAP data at orbital-phase offsets (d) away from the real transits, fitted as the real ones; scripts/toi421/fit_toi421.py inject, inject_trend',
           'mode const: b = 0.748 at every epoch; mode trend: b = 0.748 - 0.0051 (year - 2022.4)'],
          ['mode', 'offset_d', 'n_early', 'n_late', 'slope_per_yr', 'z', 'step_z', 'sector_slope_per_yr', 'sector_z'], rows)


def planetb():
    pb = J('planet_b'); f = J('ffi')['R2.5']['b_epochs']; rows = []
    for src, d in (('SPOC 2-min PDCSAP', pb), ('own TESScut photometry r = 2.5 px', f)):
        for kind in ('batman', 'box', 'trap'):
            if kind not in d:
                continue
            for ep in ('early', 'late'):
                x = d[kind][ep]
                if kind == 'batman':
                    rows.append([src, ep, x['year'], x['ntr'], 'b (k 0.0279, a/R* 13.76)', x['x'][0], x['ex'][0], x['dt_min'], ''])
                elif kind == 'box':
                    rows.append([src, ep, x['year'], x['ntr'], 'depth in a 60-min box (ppm)', x['x'][0] * 1e6, x['ex'][0] * 1e6, '', ''])
                else:
                    rows.append([src, ep, x['year'], x['ntr'], 'trapezoid T14 (min)', x['x'][1] * 1440, x['ex'][1] * 1440, x['dt_min'],
                                 'depth %.0f +- %.0f ppm' % (x['x'][0] * 1e6, x['ex'][0] * 1e6)])
    be, bl = pb['batman']['early'], pb['batman']['late']
    sl = (bl['x'][0] - be['x'][0]) / (bl['year'] - be['year']); esl = np.hypot(bl['ex'][0], be['ex'][0]) / (bl['year'] - be['year'])
    rows.append(['SPOC 2-min PDCSAP', 'late - early', '', '', 'db/dt (1/yr)', sl, esl, '', 'z %.1f' % (sl / esl)])
    for kind, lab in (('box', 'depth in a 60-min box (ppm)'),):
        e, l = pb[kind]['early'], pb[kind]['late']
        d = (l['x'][0] - e['x'][0]) * 1e6; s = np.hypot(l['ex'][0], e['ex'][0]) * 1e6
        rows.append(['SPOC 2-min PDCSAP', 'late - early', '', '', lab, d, s, '', 'z %.1f' % (d / s)])
    write('planetb_94986319.csv',
          ['TOI-421 b = TIC 94986319 b: early (sectors 5, 6, 32) and late (98, 107) folds; one mid-time offset per epoch group, quadratic baseline per window, c transits masked; scripts/toi421/planet_b_toi421.py, ffi_toi421.py',
           'b at the fixed template k 0.0279, a/R* 13.76, u 0.401/0.222 (circular); near grazing b mainly measures T14 and depth at this template; dt_min = mid-time offset of the group'],
          ['product', 'group', 'year', 'n_transits', 'quantity', 'value', 'err', 'dt_min', 'note'], rows)


def cheops():
    S = json.load(open(os.path.join(W, 'cheops', 'summary_cheops.json')))['runs']
    pred = json.load(open(os.path.join(D, 'predictions_cheops_jwst_94986319.json')))['predictions']
    PK = {1: 'b: CHEOPS 2021-11-19', 2: 'b: CHEOPS 2021-12-10', 3: 'b: CHEOPS 2022-01-10', 4: 'b: CHEOPS 2022-11-28',
          5: 'c: CHEOPS 2022-12-04', 6: 'c: CHEOPS 2022-12-20'}
    G = {'v1': [1], 'v2': [2], 'v3': [3], 'v4': [4], 'jb2021': [1, 2, 3], 'jball': [1, 2, 3, 4], 'v5': [5], 'v6': [6], 'jcall': [5, 6]}
    CONST = {'b': (0.9350, 0.0070), 'c': (0.7654, 0.0046)}
    rows = []
    for name in sorted(S, key=lambda n: (n[0], n.split('_')[1], list(G).index(n.split('_')[2]), n)):
        x = S[name]; pl, prod, grp = name.split('_')[:3]; var = '_'.join(name.split('_')[3:])
        vs = G[grp]; ps = [pred[PK[v]] for v in vs]; w = np.array([1 / p['sig'] ** 2 for p in ps])
        pb_ = float(np.sum(w * [p['b'] for p in ps]) / w.sum()); psig = float(np.mean([p['sig'] for p in ps]))
        yr = float(np.mean([p['year'] for p in ps])); b, bm, bp = x['b']; sb = 0.5 * (bm + bp)
        k = x.get('k') or [None, None, None]
        rows.append([pl, prod, grp, ' '.join(str(v) for v in vs), var, yr, b, bm, bp, k[0], k[1], k[2], pb_, psig, (b - pb_) / np.hypot(sb, psig),
                     CONST[pl][0], (b - CONST[pl][0]) / np.hypot(sb, CONST[pl][1]), x.get('tau_max')])
    write('cheops_94986319.csv',
          ['TOI-421 b and c = TIC 94986319: impact parameters from the six public CHEOPS visits of programme CH_PR100024 (2021-11 to 2022-12); scripts/toi421/cheops/ (emcee; median and 16th/84th percentiles)',
           'products: PIPE = PSF photometry published by Krenn et al. (2024), VizieR J/A+A/686/A301/cheops1; DRP = CHEOPS DRP 14.1.2 DEFAULT aperture (DACE, V0300); variants: k fixed (b: 0.0279; c: 0.05459) or free; LD atlas/phoenix (Claret 2021), ldcu (Krenn et al. 2024 prior centre), tess (TESS-band values, diagnostic only), free',
           'template: b a/R* 13.76, c a/R* 31.91, circular, as in the TESS fits; groups: v1-v6 single visits, jb2021 = visits 1-3, jball = 1-4, jcall = 5-6',
           'pred_2026_10_28: linear-trend value recorded before these fits (data/predictions_cheops_jwst_94986319.json); const: TESS 2018-2020 value of this record (b: 0.935 +- 0.007, SPOC early fold; c: 0.7654 +- 0.0046)'],
          ['planet', 'product', 'group', 'visits', 'variant', 'year', 'b', 'b_minus', 'b_plus', 'k', 'k_minus', 'k_plus',
           'pred_2026_10_28', 'pred_err', 'z_vs_pred', 'const', 'z_vs_const', 'max_autocorr_steps'], rows)


def jwst():
    rows = []
    for fn in sorted(glob.glob(os.path.join(W, 'jwst', 'fits', '*_var*.json'))):
        x = json.load(open(fn)); p, d = x['params'], x['derived']; beta = max(x['diag']['beta_time_avg'].values())
        b = p['b']; sbb = 0.5 * (b['plus'] + b['minus']) * beta
        yr = 2023.692 if x['label'].startswith('V1') else 2023.849
        rows.append([x['label'], x['variant'], yr, b['median'], b['minus'], b['plus'], beta, sbb, p['k']['median'], 0.5 * (p['k']['minus'] + p['k']['plus']),
                     d['T14_h']['median'] * 60, 0.5 * (d['T14_h']['minus'] + d['T14_h']['plus']) * 60, d['aR']['median'],
                     d['tc_minus_pred_min']['median'], x['diag']['redchi2_raw_err'], x['diag']['jitter_ppm'],
                     (b['median'] - 0.963) / np.hypot(sbb, 0.005), (b['median'] - 0.935) / np.hypot(sbb, 0.007)])
    write('jwst_94986319.csv',
          ['TOI-421 b = TIC 94986319 b: JWST GO 1935 white-light transits (V1 2023-09-10 NIRSpec G395M NRS1; V2 2023-11-06 NIRISS SOSS order 1); scripts/toi421/jwst/ (emcee, 48 walkers)',
           'light curves: eureka = Stage-4 white-light curve of Davenport et al. (2025), Zenodo 10.5281/zenodo.15176438; mast = own sum of the MAST stage-2 x1dints',
           'variants: a = a/R* 13.76, LD free; b = a/R* 13.76, LD ExoTiC-LD MPS1; bs = Stagger LD; bl = linear baseline, no ramp; c = a/R* free; tc_OC_min relative to P 5.19759807 d, T0 BTJD 1441.283952',
           'sigma_x_beta = half-width of the 68 % interval times the time-averaging red-noise factor beta; z_vs_pred: prediction of 2026-10-28 (0.963 +- 0.005); z_vs_const: SPOC 2018-2020 early fold (0.935 +- 0.007)'],
          ['label', 'variant', 'year', 'b', 'b_minus', 'b_plus', 'beta', 'b_sigma_x_beta', 'k', 'k_err', 'T14_min', 'T14_err_min', 'aR',
           'tc_OC_min', 'redchi2_raw', 'jitter_ppm', 'z_vs_pred', 'z_vs_const'], rows)


def nbody():
    n = J('nodal_nbody'); rows = []
    for x in n['runs'] + [dict(m, scan='mass') for m in n['mass_scan']]:
        rows.append(['mass scan' if x.get('scan') else 'geometry scan', x['i_bc'], x['mb'], x['mc'], 'same' if x['same_hemi'] else 'opposite', x['sgn'],
                     x['dOmega_deg'], n['periods'][str(x['i_bc'])][0] if not x.get('scan') else '', x['db_c_dt_conv'], x['db_b_dt_conv'], x['ratio'], x['Lc_over_Lb']])
    write('nbody_94986319.csv',
          ['TOI-421 b and c = TIC 94986319: REBOUND WHFast runs (40 yr; node period from 400-yr runs); scripts/toi421/nbody_toi421.py',
           'M* 0.852 Msun, R* 0.871 Rsun, masses 7.17 and 16.42 Mearth (Carleo et al. 2020) unless scanned; sky inclinations from the 2018-2020 TESS b values (b 0.944 at a/R* 13.76, c 0.766 at a/R* 31.91); dOmega = sky-plane node offset that gives the mutual inclination i_bc',
           'hemisphere: both transit chords on the same / opposite side of the stellar centre; rates converted to b at the fixed templates (a/R* 31.91 for c, 13.76 for b); ratio = -d(cos i_b)/d(cos i_c)'],
          ['scan', 'i_bc_deg', 'm_b_Mearth', 'm_c_Mearth', 'hemisphere', 'node_sign', 'dOmega_deg', 'node_period_yr', 'db_c_dt_per_yr',
           'db_b_dt_per_yr', 'ratio', 'Lc_over_Lb'], rows)


def forecast():
    P = json.load(open(os.path.join(D, 'predictions_94986319.json')))
    c = P['c_fit']; C = np.array(c['cov']); bf = P['b_fit']
    rows = []
    for yr in np.round(np.arange(2018.5, 2031.01, 0.1), 2):
        v = np.array([1, yr - 2022]); bc = float(v @ [c['b2022'], c['slope']]); ebc = float(np.sqrt(v @ C @ v))
        bb = bf['b_early'] + bf['slope'] * (yr - bf['year_early'])
        rows.append(['curve', yr, bc, ebc, T14(bc, ARC, KC, PC), T14(0.7654, ARC, KC, PC), bb, '', T14(bb, ARB, KB, PB) if bb < 1 + KB else 0,
                     boxdepth(bb) if bb < 1 + KB else 0, boxdepth(bf['b_early'])])
    for x in P['rows']:
        rows.append(['recorded prediction', x['year'], x['c_b'], x['c_eb'], x['c_T14'], x['c_T14_const_early'], x['b_b'], x['b_eb'], x['b_T14'], x['b_box60_depth'], x['b_box60_depth_const']])
    write('forecast_94986319.csv',
          ['TOI-421 = TIC 94986319: linear trends of the recorded predictions (data/predictions_94986319.json, written 2026-10-29 22:54 UTC) evaluated on a grid for the figures; scripts/toi421/make_tables_toi421.py',
           'c: b = b2022 + slope (year - 2022) with the fit covariance; T14 of the template (a/R* 31.91, k 0.05459); const: b = 0.7654 (2018-2020); planet b: linear between the epoch folds, T14 and 60-min box depth of the template (a/R* 13.76, k 0.0279); 0 = no transit (b > 1 + k)',
           'these are extrapolations of linear fits, not predictions of a dynamical model'],
          ['kind', 'year', 'c_b', 'c_b_err', 'c_T14_min', 'c_T14_const_min', 'b_b', 'b_b_err', 'b_T14_min', 'b_box60_depth_ppm', 'b_box60_depth_const_ppm'], rows)


if __name__ == '__main__':
    sectors(); pertransit(); trend(); injection(); planetb(); cheops(); jwst(); nbody(); forecast()
