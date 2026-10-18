"""One TRICERATOPS calc_probs run (N = 1e6 draws) on the prepared light curve; writes the scenario table and a result record.

usage: python scripts/triceratops_run.py <tic> <seed> [--tag _J1.0] [--label nominal] [--sigx 1.0]
                                        [--nb ID:M:R:Teff] [--cc file [--cc-measured]] [--work ./work] [--out data/triceratops]
input:  <work>/<tic>/triceratops_target.pkl (triceratops_target.py), <work>/<tic>/lc_<tic><tag>.npz (triceratops_prep.py)
output: <out>/probs_<tic>_<label>[_HYPOcc]_seed<seed>.csv (scenario table) and <out>/res_<...>.json (FPP, NFPP, settings)
Options: --sigx scales the flux error; --nb overrides mass, radius and Teff of one TIC star (used for the 2.3"
neighbour of TIC 402710280, which has no TIC parameters); --cc passes a contrast-curve file: without --cc-measured it
is treated as ASSUMED (data/triceratops/hypothetical_speckle_832nm.txt, reference runs only; file tag _HYPOcc); with
--cc-measured it is an observed curve (TIC 455017950: data/triceratops/cc_HIP74344_SOAR2019_I.txt, SOAR HRCam I band,
Tokovinin et al. 2020; run with --label SOARcc; the runs without it have --label nocc).
Settings fixed for all runs: parallel=True, exptime=10 min, nsamples=20, dilution_corrected=False. ~2 GB RSS.
"""
import argparse, json, os, pickle, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('seed', type=int); p.add_argument('--tag', default=''); p.add_argument('--label', default='nominal')
p.add_argument('--sigx', type=float, default=1.0); p.add_argument('--nb', default=None); p.add_argument('--cc', default=None); p.add_argument('--out', default='data/triceratops')
p.add_argument('--cc-measured', action='store_true', help='the --cc file is an observed contrast curve')
args = p.parse_args(); wd = workdir(args, args.tic); tic = args.tic; os.makedirs(args.out, exist_ok=True)
t = pickle.load(open(os.path.join(wd, 'triceratops_target.pkl'), 'rb'))
if args.nb:
    i, M, R, Te = args.nb.split(':'); j = t.stars.index[t.stars.ID.astype(str) == str(i)][0]
    t.stars.loc[j, 'mass'] = float(M); t.stars.loc[j, 'rad'] = float(R); t.stars.loc[j, 'Teff'] = float(Te)
z = np.load(os.path.join(wd, f'lc_{tic}{args.tag}.npz')); sig = float(z['sig']) * args.sigx; P = float(z['P'])
np.random.seed(args.seed); t0 = time.time()
kw = dict(contrast_curve_file=args.cc, filt='TESS') if args.cc else {}
t.calc_probs(time=z['t'], flux_0=z['f'], flux_err_0=sig, P_orb=P, N=1000000, parallel=True, exptime=10 / 1440, nsamples=20, dilution_corrected=False, verbose=0, **kw)
HYPO = bool(args.cc) and not args.cc_measured
tag = f'{tic}_{args.label}{"_HYPOcc" if HYPO else ""}_seed{args.seed}'
t.probs.to_csv(os.path.join(args.out, f'probs_{tag}.csv'), index=False)
json.dump(dict(tic=tic, seed=args.seed, label=args.label, lc=args.tag, cc_hypothetical=HYPO, contrast_curve=os.path.basename(args.cc) if args.cc else None, sigx=args.sigx, nb=args.nb, FPP=float(t.FPP), NFPP=float(t.NFPP),
               degenerate=bool(t.FPP_degenerate), sigma=sig, npts=len(z['t']), P=P, runtime_s=round(time.time() - t0), lc_note=str(z['note']), bg_file=os.path.basename(t.trilegal_fname)),
          open(os.path.join(args.out, f'res_{tag}.json'), 'w'), indent=1)
print(tag, 'FPP=%.4g NFPP=%.4g' % (t.FPP, t.NFPP), 'degenerate', t.FPP_degenerate, '%ds' % (time.time() - t0), flush=True)
