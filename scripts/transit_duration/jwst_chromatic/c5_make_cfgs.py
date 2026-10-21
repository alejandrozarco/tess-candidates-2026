"""C5: write the per-visit bin-fit configurations for c3_binfit.py, with mid-time and a/R* from the MAP of the joint
white-light fits of section 6.6.10 (ns_main for d, ctrl_c_nrs1 for c, ni_1201 for the NIRISS visit of d). Own code.
Inputs: $TDJ_WORK/out/jwst/fits/{ns_main,ctrl_c_nrs1,ni_1201}.json (../joint/run_joint.sh step 10) and the column
extractions $TDC_WORK/out/cols_*.npz / soss_1201_bins.npz (c1_extract_cols.py, c4_soss_bins.py).
Outputs: $TDC_WORK/cfg/<name>.json; the fits are written to $TDC_WORK/fits/<name>.jsonl.
Settings of the published run: 1500 emcee steps per bin, except d1224_nrs1_prior (3000 steps, the first run); the
exponential-ramp variant (tau 15 min) for the GTO 1224 visit only. Usage: c5_make_cfgs.py"""
import json, os, numpy as np
H = os.path.abspath(os.environ.get('TDC_WORK', 'work/transit_duration_chromatic'))
F = os.path.join(os.path.abspath(os.environ.get('TDJ_WORK', 'work/transit_duration_joint')), 'out', 'jwst', 'fits')
PL = {'b': (2.2531140, 1366.1706), 'c': (3.6906764, 1367.2730), 'd': (7.450729, 1362.7400)}   # P (d), T0 (BTJD): Cadieux et al. (2025)
G = 6.674e-11
def aR(lnrho, P): return float((G * np.exp(lnrho) * (P * 86400.) ** 2 / (3 * np.pi)) ** (1 / 3.))
def tc(pl, tmean, dtc): P, T0 = PL[pl]; n = np.round((tmean - T0) / P); return float(T0 + n * P + dtc / 1440.)
os.makedirs(H + '/cfg', exist_ok=True); os.makedirs(H + '/fits', exist_ok=True)
ns = json.load(open(F + '/ns_main.json'))['map']; cc = json.load(open(F + '/ctrl_c_nrs1.json'))['map']; ni = json.load(open(F + '/ni_1201.json'))['map']
spec = []
for v in ('1224', '4098'):
    for det in ('nrs1', 'nrs2'):
        spec.append(dict(name='d%s_%s' % (v, det), planet='d', file=H + '/out/cols_%s_%s.npz' % (v, det), ld_key=det.upper(), lnrho=ns['lnrho'], dtc=ns['dtc_' + v], b_guess=0.93, k_guess=0.045, bgrid=[0.890, 0.9705, 0.0015]))
for v in ('c2512a', 'c2512b'):
    spec.append(dict(name='%s_nrs1' % v, planet='c', file=H + '/out/cols_%s_nrs1.npz' % v, ld_key='NRS1', lnrho=cc['lnrho'], dtc=cc['dtc_' + v], b_guess=0.38, k_guess=0.040, bgrid=[0.0, 0.6005, 0.015], b_starts=[0.2, 0.38, 0.5]))
for order in ('SOSS1', 'SOSS2'):
    spec.append(dict(name='d1201_%s' % order.lower(), planet='d', file=H + '/out/soss_1201_bins.npz', ld_key=order, lnrho=ns['lnrho'], dtc=ni['dtc_1201'], b_guess=0.93, k_guess=0.043, bgrid=[0.890, 0.9705, 0.0015], regs=['bkg', 'tilt']))
for s in spec:
    if not os.path.exists(s['file']): print('missing', s['file']); continue
    z = np.load(s['file']); P = PL[s['planet']][0]
    modes = ['prior', 'free'] + (['ramp'] if s['name'].startswith('d1224') else [])
    for ldm in modes:
        name = s['name'] + '_' + ldm
        c = dict(name=name, file=s['file'], P=P, tc_btjd=tc(s['planet'], z['t'].mean(), s['dtc']), aR=aR(s['lnrho'], P), ld='prior' if ldm == 'ramp' else ldm, ld_sigma=0.1,
                 ld_file=H + '/out/ld_kurucz.json', ld_key=s['ld_key'], regs=s.get('regs', ['xshift', 'yshift', 'width']), bin=5, trim_min=10, poly=2,
                 nsteps=3000 if name == 'd1224_nrs1_prior' else 1500, bgrid=s['bgrid'], b_guess=s['b_guess'], k_guess=s['k_guess'], out=H + '/fits/%s.jsonl' % name)
        if 'b_starts' in s: c['b_starts'] = s['b_starts']
        if ldm == 'ramp': c['ramp_min'] = 15
        json.dump(c, open(H + '/cfg/%s.json' % name, 'w'), indent=1); print('cfg', name, 'tc %.5f aR %.2f' % (c['tc_btjd'], c['aR']))
