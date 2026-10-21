"""T9: block bootstrap of the evolving-geometry fit: resample whole transits, or whole sectors, with replacement
(multiplicity weights), refit lin_b by least squares, record db/dt and dT14/dt. Also leave-one-campaign-out and
leave-one-sector-out jackknives. Checkpointed (resume). Own code.
Usage: t9_boot.py <win.npy> <mode transit|sector|jack> <nrep> <out.npz>"""
import sys, os, numpy as np, time
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from t5lib import *
win, mode, nrep, outp = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
segs = list(np.load(win, allow_pickle=True))
D0 = prepare(segs)            # clipping + noise model from the full data, then held fixed
keep, sig = D0.keep, D0.sig
S0 = Spec('lin_b', D=D0, ld0=LD0); x_full, c2, lp, cov = lsq(S0, default_x0(S0))
def fit_with(mult):
    D = TData(segs, mult=mult); D.build(keep); D.set_sig(sig)
    S = Spec('lin_b', D=D, ld0=LD0); x, c2, lp, cov = lsq(S, x_full)
    T = S.t14_at(x, np.array([TREF - 0.5, TREF + 0.5]))
    return x[S.ix['bdot']], T[1] - T[0], x
res = dict(np.load(outp, allow_pickle=True)) if os.path.exists(outp) else {}
bd, dt, lab = list(res.get('bdot', [])), list(res.get('dT14', [])), list(res.get('label', []))
rng = np.random.default_rng({'transit': 11, 'sector': 22, 'jack': 33}[mode] + len(bd))
secs = np.array([s['sector'] for s in segs]); usec = np.unique(secs); camp = np.array([campaign(s) for s in secs])
jobs = []
if mode == 'jack':
    jobs = [('camp%d' % c, (camp != c).astype(float)) for c in range(4)] + [('sec%d' % s, (secs != s).astype(float)) for s in usec]
t0 = time.time()
i = len(bd)
while i < (len(jobs) if mode == 'jack' else nrep):
    if mode == 'transit':
        mult = np.bincount(rng.integers(0, len(segs), len(segs)), minlength=len(segs)).astype(float); label = 'b%d' % i
    elif mode == 'sector':
        pick = rng.choice(usec, len(usec), replace=True); cnt = {s: np.sum(pick == s) for s in usec}
        mult = np.array([cnt[s] for s in secs], float); label = 'b%d' % i
    else:
        label, mult = jobs[i]
    b, d, x = fit_with(mult); bd.append(b); dt.append(d); lab.append(label); i += 1
    if i % 10 == 0 or mode == 'jack':
        np.savez(outp, bdot=np.array(bd), dT14=np.array(dt), label=np.array(lab), full_bdot=x_full[S0.ix['bdot']])
        print(mode, i, 'bdot %.5f dT14 %.3f' % (b, d), '%.0f s' % (time.time() - t0), flush=True)
np.savez(outp, bdot=np.array(bd), dT14=np.array(dt), label=np.array(lab), full_bdot=x_full[S0.ix['bdot']])
bd, dt = np.array(bd), np.array(dt)
if mode != 'jack':
    print('%s bootstrap n=%d: bdot %.5f +- %.5f (16-84: %.5f..%.5f); dT14/dt %.3f +- %.3f min/yr; fraction dT14>=0: %.4f' % (
        mode, len(bd), np.median(bd), np.std(bd), *np.percentile(bd, [16, 84]), np.median(dt), np.std(dt), np.mean(dt >= 0)))
else:
    for l, b, d in zip(lab, bd, dt): print(l, '%.5f %.3f' % (b, d))
