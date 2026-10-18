"""Summarise the TRICERATOPS runs: FPP/NFPP statistics over seeds per configuration and mean scenario probabilities.

usage: python scripts/triceratops_summary.py [--dir data/triceratops]
input:  <dir>/res_*.json and <dir>/probs_*.csv (triceratops_run.py); or, as shipped, <dir>/runs.csv + probs_*.csv
output: printed summary (the numbers quoted in README.md) and <dir>/summary.csv
"""
import argparse, collections, glob, json, os
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--dir', default='data/triceratops'); args = p.parse_args()
R = collections.defaultdict(list)
res = sorted(glob.glob(os.path.join(args.dir, 'res_*.json')))
if res:
    for f in res:
        d = json.load(open(f)); R[(d['tic'], d['label'], d['cc_hypothetical'])].append((d['seed'], d['FPP'], d['NFPP'], os.path.join(args.dir, 'probs_' + os.path.basename(f)[4:-5] + '.csv')))
else:
    for _, d in pd.read_csv(os.path.join(args.dir, 'runs.csv'), comment='#').iterrows():
        R[(int(d.tic), d.label, bool(d.cc_hypothetical))].append((int(d.seed), float(d.FPP), float(d.NFPP), os.path.join(args.dir, d.probs_file)))
rows = []
for key in sorted(R, key=lambda k: (k[0], k[1], k[2])):
    v = sorted(R[key]); F = np.array([x[1] for x in v]); N = np.array([x[2] for x in v])
    print(f"\n=== TIC {key[0]}  {key[1]}  {'ASSUMED contrast curve' if key[2] else 'no imaging constraint'}  seeds={len(v)}")
    print(f"  FPP  mean {F.mean():.3g} +- {F.std(ddof=1) if len(F) > 1 else 0:.2g}; median {np.median(F):.3g}; range [{F.min():.2g}, {F.max():.2g}]")
    print(f"  NFPP mean {N.mean():.3g}; range [{N.min():.2g}, {N.max():.2g}]")
    P = pd.concat([pd.read_csv(x[3]).assign(seed=x[0]) for x in v]); P['who'] = np.where(P.ID == key[0], 'target', P.ID.astype(str))
    g = P.groupby(['who', 'scenario'], sort=False).prob.agg(['mean', 'std', 'min', 'max']); print(g[g['max'] > 1e-6].to_string(float_format=lambda x: f'{x:.3g}'))
    print('  best-fit TP planet radius (Earth radii) per seed:', np.round(P[P.scenario == 'TP'].R_p.values, 2))
    rows.append([key[0], key[1], key[2], len(F), f'{F.mean():.3e}', f'{F.std(ddof=1) if len(F) > 1 else 0:.2e}', f'{np.median(F):.3e}', f'{F.min():.2e}', f'{F.max():.2e}', f'{N.mean():.2e}', f'{N.max():.2e}'])
pd.DataFrame(rows, columns=['tic', 'label', 'cc_hypothetical', 'n_seeds', 'FPP_mean', 'FPP_std', 'FPP_median', 'FPP_min', 'FPP_max', 'NFPP_mean', 'NFPP_max']).to_csv(os.path.join(args.dir, 'summary_recomputed.csv'), index=False)
