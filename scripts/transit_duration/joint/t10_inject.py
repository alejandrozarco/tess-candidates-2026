"""T10: null injections. A constant-geometry copy of d (posterior-median k, a/R*, LD and b at 2022.0) is inserted at
a random orbital phase into quality-selected, normalised PDCSAP flux (mode 'lc'), or into calibrated pixels after the
own background subtraction and before aperture summation (mode 'pix'; the whole background-subtracted image is
multiplied, including contaminating sources); the subsequent transit-processing and fitting steps (windows, masks,
flare and outlier clipping, noise model, const and lin_b least-squares fits) are rerun. Upstream calibration, PDC and
background estimation are not tested. Records the fitted db/dt, dT14/dt, their
errors and delta-chi2. Mode 'lc' injects into SPOC PDCSAP; mode 'pix' injects into the target pixels and uses own
aperture photometry (t11). Checkpointed jsonl. Own code.
Usage: t10_inject.py <lc|pix> <n> <seed> <out.jsonl> <main_posterior.json>"""
import sys, os, json, numpy as np, time, tempfile, shutil
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from t5lib import *; from t3_windows import make_windows
mode, n, seed, outp, postf = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
post = json.load(open(postf))['post']
lnrho = post['lnrho'][0]; aR = float(aR_from_lnrho(lnrho, PL['d'][0])); u1, u2 = ld_q2u(post['q1'][0], post['q2'][0])
base = dict(k=post['k'][0], aR=aR, b=post['b0'][0], u1=float(u1), u2=float(u2))
done = sum(1 for _ in open(outp)) if os.path.exists(outp) else 0
rng = np.random.default_rng(seed); phases = rng.uniform(0.08, 0.92, n)
lcs = load_lc('2m') if mode == 'lc' else None
if mode == 'pix':
    import t11_tpf
    tmpd = tempfile.mkdtemp(dir=ROOT + '/out')
t0 = time.time()
for i in range(done, n):
    inj = dict(base, phase=float(phases[i]))
    if mode == 'lc':
        segs, fl = make_windows(lcs, 'pdc', 'mask', inj)
    else:
        t11_tpf.phot('spoc', tmpd, inj); segs, fl = make_windows(load_lc('2m', tmpd), 'pdc', 'mask', inj, apply_inj=False)   # injected at pixel level only
    D = prepare(segs); r = fit_pair(D)
    rec = dict(i=i, phase=inj['phase'], ntr=len(segs), **{k: v for k, v in r.items() if k not in ('x', 'names')})
    with open(outp, 'a') as fo: fo.write(json.dumps(rec) + '\n')
    print(i, 'phase %.3f ntr %d bdot %.5f +- %.5f dT14 %.3f +- %.3f dchi2 %.1f  %.0f s' % (inj['phase'], len(segs), r['bdot'], r['bdot_err'], r['dT14'], r['dT14_err'], r['dchi2'], time.time() - t0), flush=True)
if mode == 'pix': shutil.rmtree(tmpd, ignore_errors=True)
