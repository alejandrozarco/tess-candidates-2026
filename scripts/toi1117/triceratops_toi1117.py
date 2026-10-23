"""TOI-1117 c (TIC 295541511), README section 6.10: TRICERATOPS false-positive probabilities.

usage: python scripts/toi1117/triceratops_toi1117.py [--work ./work] [--seeds 1,2] [--N 300000]   (TRICERATOPS environment)
input:  <work>/295541511/lc_fit_pdc.npz (fit_toi1117.py), data/ephemeris_295541511.csv, data/model_295541511.csv,
        data/triceratops/cc_295541511_approx_Zorro832.txt (approximate contrast curve, see below)
output: data/triceratops_input_295541511.csv   the folded, binned light curve given to TRICERATOPS
        data/triceratops/probs_295541511_<label>_seed<seed>.csv  scenario tables
        data/triceratops/gaia_background_295541511.csv  the Gaia DR3 background population drawn by TRICERATOPS
        data/triceratops_295541511.csv         FPP and NFPP per run, and the largest scenario probabilities
Light curve: SPOC 2-min PDCSAP sectors 39, 100-103 (detrended as in fit_toi1117.py), folded on the adopted ephemeris, 10-min bins
within +-0.22 d; flux error = standard deviation of the bins more than 0.06 d from mid-transit. Target: tr.target(ID, sectors 39,
100-103, the SPOC 2-min sectors of the light curve), apertures = the SPOC optimal apertures (tg.get_spoc_apertures()); transit
depth = depth of the model curve rho_2min in data/model_295541511.csv (1 - its minimum: the depth in the PDCSAP light curve);
dilution_corrected = True (PDCSAP flux is corrected for crowding); N draws, parallel = False, exptime = 10 min, nsamples = 20. Background population:
TRICERATOPS Gaia option. Labels: nocc = no contrast curve; ccapprox = the approximate curve.
The approximate contrast curve is not a measured curve: it was written by hand from the text of Lockley et al. (2025): no companion
within 5-8 mag of the target from 0.02" to 1.2" in Gemini South/Zorro speckle imaging and 4.2 mag at 1" from SOAR (sections 2.3-2.4);
the measured curves were not used. Requires triceratops 1.1.0 and network access to MAST and Gaia.
"""
import argparse, os, sys, time, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import TIC, SPOC2 as SECTORS, phase, read_csv, write_csv, workdir, DATE
import triceratops.triceratops as tr
from triceratops.funcs import query_gaia_background

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--seeds', default='1,2'); p.add_argument('--N', type=int, default=300000)
args = p.parse_args(); ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); DATA = os.path.join(ROOT, 'data')
CC = os.path.join(DATA, 'triceratops', 'cc_295541511_approx_Zorro832.txt')
E = read_csv(os.path.join(DATA, 'ephemeris_295541511.csv')); ia = list(E['role']).index('adopted'); Tref, P = float(E['Tref_btjd'][ia]), float(E['P_d'][ia])
dep = float(1.0 - np.min(read_csv(os.path.join(DATA, 'model_295541511.csv'))['rho_2min']))

z = np.load(os.path.join(workdir(args.work), 'lc_fit_pdc.npz')); k = z['sec'] != 13
ph = phase(z['t'][k], P, Tref); f = z['f'][k]; w = np.abs(ph) < 0.22; ph, f = ph[w], f[w]
b = np.floor((ph + 0.22) / (10 / 1440)).astype(int); ub = np.unique(b)
tb = np.array([ph[b == u].mean() for u in ub]); fb = np.array([f[b == u].mean() for u in ub]); nb = np.array([(b == u).sum() for u in ub])
sig = float(np.std(fb[np.abs(tb) > 0.06]))
write_csv(os.path.join(DATA, 'triceratops_input_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): light curve given to TRICERATOPS; scripts/toi1117/triceratops_toi1117.py, {DATE}',
           f'SPOC 2-min PDCSAP sectors 39, 100-103 detrended as in fit_toi1117.py (b masked), folded on Tref = BTJD {Tref:.5f}, P = {P:.7f} d; 10-min phase bins within +-0.22 d',
           f'flux error used for every bin: {sig:.6f} (standard deviation of the bins more than 0.06 d from mid-transit); n = cadences per bin; depth given to TRICERATOPS: {dep * 1e6:.0f} ppm'],
          ['dt_days', 'flux', 'n'], [(f'{a:.6f}', f'{c:.7f}', n_) for a, c, n_ in zip(tb, fb, nb)])

tg = tr.target(ID=TIC, sectors=np.array(SECTORS))
for a in range(5):
    if isinstance(tg.trilegal_fname, str):
        break
    time.sleep(20); tg.trilegal_fname = query_gaia_background(tg.stars.ra.values[0], tg.stars.dec.values[0], TIC, verbose=0)
if not isinstance(tg.trilegal_fname, str):
    sys.exit('no background population')
dst = os.path.join(DATA, 'triceratops', f'gaia_background_{TIC}.csv'); os.replace(tg.trilegal_fname, dst); tg.trilegal_fname = dst
aps = tg.get_spoc_apertures()
if any(ap is None for ap in aps):
    raise RuntimeError("Missing SPOC aperture")
tg.calc_depths(tdepth=dep, all_ap_pixels=aps, dilution_corrected=True)
rows = []
for label, ccf in (('nocc', None), ('ccapprox', CC)):
    for seed in [int(x) for x in args.seeds.split(',')]:
        np.random.seed(seed); t0 = time.time()
        tg.calc_probs(time=tb, flux_0=fb, flux_err_0=sig, P_orb=P, N=args.N, parallel=False, exptime=10 / 1440, nsamples=20, dilution_corrected=True, verbose=0,
                      contrast_curve_file=ccf)
        pr = tg.probs.copy(); pr.to_csv(os.path.join(DATA, 'triceratops', f'probs_295541511_{label}_seed{seed}.csv'), index=False)
        top = pr.sort_values('prob', ascending=False)
        nonhost = top[top.ID.astype(str) != str(TIC)].head(3)
        rows.append((label, seed, args.N, f'{tg.FPP:.4f}', f'{tg.NFPP:.4f}', f'{float(pr.prob[pr.scenario == "TP"].sum()):.4f}',
                     '; '.join(f'{r.scenario} on TIC {r.ID}: {r.prob:.4f}' for r in nonhost.itertuples()), int(time.time() - t0)))
        print(rows[-1], flush=True)
st = tg.stars[['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth']].head(8)
write_csv(os.path.join(DATA, 'triceratops_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): TRICERATOPS 1.1.0 runs; scripts/toi1117/triceratops_toi1117.py, {DATE}; input light curve triceratops_input_295541511.csv',
           f'P = {P:.7f} d, depth {dep * 1e6:.0f} ppm (PDCSAP model depth, dilution_corrected = True), SPOC optimal apertures of sectors {", ".join(map(str, SECTORS))}; label nocc = no contrast curve, ccapprox = triceratops/cc_295541511_approx_Zorro832.txt (approximate, written from the text of Lockley et al. 2025; not the measured curve)',
           'P_TP = probability of the transiting-planet scenario on the target; largest_off_target = the three largest scenario probabilities on other stars; runtime_s = wall-clock time of the run',
           'nearest TRICERATOPS stars (TIC, Tmag, separation ", flux ratio, depth required on that star): ' + '; '.join(f'{int(r.ID)} {r.Tmag:.2f} {r._3:.1f} {r.fluxratio:.4f} {r.tdepth:.4f}' for r in st.itertuples())],
          ['label', 'seed', 'N', 'FPP', 'NFPP', 'P_TP', 'largest_off_target', 'runtime_s'], rows)
