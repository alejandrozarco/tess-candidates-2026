#!/usr/bin/env python
"""TOI-421 (TIC 94986319) CHEOPS impact-parameter fits, planets b and c: collect runs/fit_*.json into tables
(markdown) and a machine-readable summary (summary_cheops.json). Conventions: see cheops_fit.py docstring.
Comparison values: predictions recorded on 2026-10-28 in data/predictions_cheops_jwst_94986319.json (not used in any fit).
usage: python scripts/toi421/cheops/summarize.py   (reads work/toi421/cheops/runs/fit_*.json; writes
work/toi421/cheops/summary_cheops.json and tables_cheops.md)"""
import json, glob, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
HERE = os.path.join(ROOT, 'work', 'toi421', 'cheops')
RUNS = os.path.join(HERE, 'runs')
PRED = json.load(open(os.path.join(ROOT, 'data', 'predictions_cheops_jwst_94986319.json')))['predictions']
PKEY = {1: 'b: CHEOPS 2021-11-19', 2: 'b: CHEOPS 2021-12-10', 3: 'b: CHEOPS 2022-01-10', 4: 'b: CHEOPS 2022-11-28',
        5: 'c: CHEOPS 2022-12-04', 6: 'c: CHEOPS 2022-12-20'}
# constant alternative: TESS 2018-2020 values (b: SPOC 2-min early-epoch fold, work/toi421/out/planet_b.json;
# c: early group mean of the per-transit fits, work/toi421/out/baseline.json)
CONST = {'b': (0.9350, 0.0070), 'c': (0.7654, 0.0046)}
GROUPS = {'v1': [1], 'v2': [2], 'v3': [3], 'v4': [4], 'jb2021': [1, 2, 3], 'jball': [1, 2, 3, 4],
          'v5': [5], 'v6': [6], 'jcall': [5, 6]}
DATES = {1: '2021-11-19', 2: '2021-12-10', 3: '2022-01-10', 4: '2022-11-28', 5: '2022-12-04', 6: '2022-12-20'}


def pred_for(group):
    vs = GROUPS[group]
    p = [PRED[PKEY[v]] for v in vs]
    w = np.array([1 / x['sig'] ** 2 for x in p])
    return float(np.sum(w * [x['b'] for x in p]) / w.sum()), float(np.mean([x['sig'] for x in p]))


def load_all():
    F = {}
    for f in sorted(glob.glob(os.path.join(RUNS, 'fit_*.json'))):
        o = json.load(open(f)); F[o['spec']['name']] = o
    return F


def fmt(r, nd=3):
    return f"{r[0]:.{nd}f} -{r[1]:.{nd}f} +{r[2]:.{nd}f}"


def zscore(bm, sm, ref, sref):
    return (bm - ref) / np.hypot(sm, sref)


def table(F, planet, product, variants, groups):
    rows = []
    hdr = ("| fit | visits (date) | variant | b | k | tc offset (min) | jitter (ppm) | chi2r nom. | beta | "
           "z vs trend | z vs const |")
    rows.append(hdr); rows.append('|' + '---|' * 11)
    for g in groups:
        for vn in variants:
            nm = f"{planet}_{product}_{g}_{vn}"
            if nm not in F: continue
            o = F[nm]; r = o['results']; vs = GROUPS[g]
            b = r['b']; sb = 0.5 * (b[1] + b[2])
            pt, ps = pred_for(g); c0, cs = CONST[planet]
            k = fmt(r['k'], 4) if 'k' in r else 'fixed'
            dtc = ', '.join(f"{r[f'dtc_{v}_min'][0]:+.1f}+-{0.5*(r[f'dtc_{v}_min'][1]+r[f'dtc_{v}_min'][2]):.1f}" for v in vs)
            jit = ', '.join(f"{r[f'jitter_{v}_ppm'][0]:.0f}" for v in vs)
            chi = ', '.join(f"{o['diagnostics'][str(v)]['chi2_red_nominal']:.2f}" for v in vs)
            bet = ', '.join(f"{o['diagnostics'][str(v)]['beta_median']:.2f}" for v in vs)
            rows.append(f"| {g} | {', '.join(f'{v} ({DATES[v]})' for v in vs)} | {vn} | {fmt(b)} | {k} | {dtc} | {jit} | "
                        f"{chi} | {bet} | {zscore(b[0], sb, pt, ps):+.1f} | {zscore(b[0], sb, c0, cs):+.1f} |")
    return '\n'.join(rows)


def ld_spread(F, planet, product, group, variants):
    v = [F[f"{planet}_{product}_{group}_{vn}"]['results']['b'][0] for vn in variants
         if f"{planet}_{product}_{group}_{vn}" in F]
    return (max(v) - min(v), v) if v else (np.nan, v)


if __name__ == '__main__':
    F = load_all()
    VB = ['kfix_LDatlas', 'kfix_LDfree', 'kfree_LDatlas', 'kfree_LDfree', 'kfix_LDphoenix', 'kfix_LDldcu', 'kfix_LDtess']
    VC = ['kfree_LDatlas', 'kfree_LDfree', 'kfix_LDatlas', 'kfix_LDfree', 'kfree_LDphoenix', 'kfree_LDldcu', 'kfree_LDtess']
    out = []
    for prod in ('PIPE', 'DRP'):
        out.append(f"\n### Planet b, {prod}\n")
        out.append(table(F, 'b', prod, VB, ['v1', 'v2', 'v3', 'v4', 'jb2021', 'jball']))
        out.append(f"\n### Planet c, {prod}\n")
        out.append(table(F, 'c', prod, VC, ['v5', 'v6', 'jcall']))
    out.append("\n### LD spread (max - min of b over CHEOPS LD variants: Claret ATLAS, Claret PHOENIX, Krenn LDCU, free)\n")
    for prod in ('PIPE', 'DRP'):
        for pl, grp, vv in (('b', 'jball', ['kfix_LDatlas', 'kfix_LDphoenix', 'kfix_LDldcu', 'kfix_LDfree']),
                            ('b', 'jb2021', ['kfix_LDatlas', 'kfix_LDphoenix', 'kfix_LDldcu', 'kfix_LDfree']),
                            ('b', 'v4', ['kfix_LDatlas', 'kfix_LDphoenix', 'kfix_LDldcu', 'kfix_LDfree']),
                            ('c', 'jcall', ['kfree_LDatlas', 'kfree_LDphoenix', 'kfree_LDldcu', 'kfree_LDfree'])):
            s, v = ld_spread(F, pl, prod, grp, vv)
            tv = F.get(f"{pl}_{prod}_{grp}_{vv[0].split('_')[0]}_LDtess")
            tb = tv['results']['b'][0] if tv else np.nan
            out.append(f"- {prod} planet {pl} {grp}: spread {s:.4f} (values {', '.join(f'{x:.4f}' for x in v)}); "
                       f"TESS-band template LD gives {tb:.4f}")
    txt = '\n'.join(out)
    print(txt)
    summ = {nm: dict(b=o['results']['b'], k=o['results'].get('k'), header=o['header'],
                     tau_max=max(o['mcmc']['tau']) if o['mcmc']['tau'] else None, acc=o['mcmc']['acceptance'])
            for nm, o in F.items()}
    json.dump(dict(header="TOI-421 (TIC 94986319) planets b and c, CHEOPS visits 1-6 (CH_PR100024), impact-parameter "
                          "fits; per-run conventions in each 'header'", runs=summ),
              open(os.path.join(HERE, 'summary_cheops.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'tables_cheops.md'), 'w').write(
        "TOI-421 (TIC 94986319) planets b and c, CHEOPS visits 1-6 (CH_PR100024): fit tables generated by summarize.py; "
        "conventions in cheops_fit.py.\n" + txt + '\n')
