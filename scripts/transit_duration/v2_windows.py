"""Cut per-transit windows for a target planet, masking other planets; upward flare clipping later.
Usage: v2_windows.py <system> <planet> <cad> <flux> <out.npy> [half-width in d, default 0.16]"""
import sys, numpy as np, json
sys.path.insert(0, __file__.rsplit('/',1)[0]); from lib import load
SYS = {
 'L9859': dict(tic=307210830, pl={  # Cadieux+2025 Table (TBJD), T14 h
   'b': (2.2531140, 1366.17056, 1.01), 'c': (3.6906764, 1367.27303, 1.28), 'd': (7.450729, 1362.74002, 0.84)}),
 'L9859inj30': dict(tic=307210830, pl={'b': (2.2531140, 1366.17056, 1.01), 'c': (3.6906764, 1367.27303, 1.28), 'd': (7.450729, 1362.74002, 0.84), 'x': (7.450729, 1362.74002 + 0.30 * 7.450729, 0.84)}),
 'L9859inj62': dict(tic=307210830, pl={'b': (2.2531140, 1366.17056, 1.01), 'c': (3.6906764, 1367.27303, 1.28), 'd': (7.450729, 1362.74002, 0.84), 'x': (7.450729, 1362.74002 + 0.62 * 7.450729, 0.84)}),
 'TOI620': dict(tic=296739893, pl={'b': (5.098784, 1518.004514, 1.45)}),
}
sysn, pln, cad, flux, out = sys.argv[1:6]
S = SYS[sysn]; P, T0, D = S['pl'][pln]; D /= 24.
half = float(sys.argv[6]) if len(sys.argv) > 6 else 0.16
segs = []
for lc in load(S['tic'], cad, flux):
    t, f, e = lc['t'], lc['f'], lc['e']
    other = np.zeros(len(t), bool)
    for q, (Pq, Tq, Dq) in S['pl'].items():
        if q == pln: continue
        ph = (t - Tq + 0.5 * Pq) % Pq - 0.5 * Pq
        other |= np.abs(ph) < (0.5 * Dq / 24. + 0.025)   # +36 min allowance for TTVs/ingress
    n0, n1 = np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P)
    for n in np.arange(n0, n1 + 1):
        tc = T0 + n * P
        sel = (np.abs(t - tc) < half)
        if sel.sum() < 10: continue
        tt, ff, ee, oo = t[sel], f[sel], e[sel], other[sel]
        intr = np.abs(tt - tc) < 0.5 * D
        pre = (tt < tc - 0.5 * D - 0.01) & ~oo; post = (tt > tc + 0.5 * D + 0.01) & ~oo
        cad_d = np.median(np.diff(tt))
        n_exp_in = D / cad_d
        cov_in = (intr & ~oo).sum() / n_exp_in
        ok = (cov_in > 0.8) and (pre.sum() * cad_d > 0.04) and (post.sum() * cad_d > 0.04)
        keep = ~oo
        segs.append(dict(sector=int(lc['sector']), epoch=int(n), tc_lin=tc, t=tt[keep], f=ff[keep], e=ee[keep], ok=bool(ok), cov_in=float(cov_in),
                         other_in_window=bool(oo.any()), crowd=float(lc['meta']['crowd'])))
good = [s for s in segs if s['ok']]
print(sysn, pln, cad, flux, 'windows', len(segs), 'usable', len(good), 'with other-planet points masked', sum(s['other_in_window'] for s in good))
np.save(out, np.array(good, dtype=object), allow_pickle=True)
from collections import Counter
print(sorted(Counter(s['sector'] for s in good).items()))
