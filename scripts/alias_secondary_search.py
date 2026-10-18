"""Blind box search in all sectors, period-alias check (P/k, k = 1..12, and 2P) and secondary-eclipse limit.

usage: python scripts/alias_secondary_search.py <tic> [--harmonics] [--work ./work]
input:  <work>/<tic>/cut/lc_s<sector>.npz (aperture_photometry.py), r = 2.0 px aperture, quality 0
Light curves are flattened with a biweight filter (wotan, window 1.0 d). A box of width 0.8 x T14 is slid over every
sector (step 0.01 d); events with depth/error > 5 are printed. For each alias P/k every predicted epoch inside the
data is classified TRANSIT (depth > 0.5 x nominal and > 4 sigma), ABSENT (depth < 0.3 x nominal with error < 0.25 x
nominal) or ambiguous; a gap is reported when less than 70 % of the window has data. The secondary limit is the box
depth at phase 0.5 plus 3 x the empirical scatter of box depths over all phases (primaries masked); the box depth of
every covered individual phase-0.5 epoch and their inverse-variance weighted mean (white-noise limit) are also printed.
--harmonics: fit cos/sin of the orbital phase and of twice the phase (reflection/beaming, ellipsoidal) plus one offset
per sector to the out-of-transit light curve (|dt| > 0.2 d), errors from 200 bootstrap resamplings of 1-d chunks,
after 1.0-d and 5.0-d biweight detrending, and repeat the fit after injecting -200 ppm cos(2 phi) and +100 ppm
sin(phi) into the undetrended flux to measure how much of such a signal survives the detrending.
"""
import argparse, os
import numpy as np
from wotan import flatten
from targets import TARGETS, workdir, add_common_args

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--harmonics', action='store_true')
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic)
T0, P = c['T0'], c['P']; HW = c['T14h'] / 48; DEP = (c['fix'][0] if c['fix'] else c['k0']) ** 2
L = {}
for s in c['sectors']:
    z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); m = (z['q'] == 0) & np.isfinite(z['f2.0']); t = z['tm'][m]; f = z['f2.0'][m] / np.nanmedian(z['f2.0'][m])
    fl = flatten(t, f, method='biweight', window_length=1.0, break_tolerance=0.3, edge_cutoff=0.0); ok = np.isfinite(fl); L[s] = (t[ok], fl[ok], np.median(np.diff(t)))


def box(t, f, tc, hw=HW * 0.8):
    win = np.abs(t - tc) < hw; oot = (np.abs(t - tc) > HW + 0.05) & (np.abs(t - tc) < 0.8)
    if win.sum() < 3 or oot.sum() < 10:
        return np.nan, np.nan, win.sum()
    s = 1.4826 * np.median(np.abs(f[oot] - np.median(f[oot]))); return np.median(f[oot]) - np.mean(f[win]), s / np.sqrt(win.sum()), win.sum()


print('=== blind box search (events with depth/error > 5) ===')
for s, (t, f, cad) in L.items():
    best = sorted([(d / e, tc, d, e) for tc in np.arange(t.min() + 0.2, t.max() - 0.2, 0.01) for d, e, n in [box(t, f, tc)] if np.isfinite(d) and n * cad > 0.6 * 2 * HW * 0.8], reverse=True); seen = []
    for snr, tc, d, e in best:
        if snr < 5:
            break
        if any(abs(tc - x) < 0.4 for x in seen):
            continue
        seen.append(tc); print(f'  S{s} tc={tc:.3f} depth={d * 1e6:.0f}+-{e * 1e6:.0f} ppm SNR={snr:.1f} (phase {((tc - T0) / P) % 1:.3f})')
    print(f'  S{s}: span {t.min():.2f}-{t.max():.2f}, rms {np.std(f) * 1e6:.0f} ppm, n={len(t)}')
print('=== alias coverage P/k (timing slop +-0.04 d) ===')
for k in [0.5] + list(range(1, 13)):
    Pk = P / k; res = []
    for s, (t, f, cad) in L.items():
        for n in range(int(np.floor((t.min() - T0) / Pk)), int(np.floor((t.min() - T0) / Pk)) + int(30 / Pk) + 3):
            tp = T0 + n * Pk
            if not (t.min() + 0.1 < tp < t.max() - 0.1):
                continue
            cov = (np.abs(t - tp) < HW).sum() * cad / (2 * HW)
            if cov < 0.7:
                res.append((s, round(tp, 3), 'gap')); continue
            d, e, _ = max([box(t, f, tp + o) for o in np.linspace(-0.04, 0.04, 9)], key=lambda x: x[0] if np.isfinite(x[0]) else -1)
            res.append((s, round(tp, 3), 'TRANSIT' if d > 0.5 * DEP and d / e > 4 else ('ABSENT' if d < 0.3 * DEP and e < 0.25 * DEP else 'ambiguous'), round(d * 1e6), round(e * 1e6)))
    print(f'k={k:>4} P={Pk:.4f}: TRANSIT {sum(r[2] == "TRANSIT" for r in res)}, ABSENT {sum(r[2] == "ABSENT" for r in res)}, gaps {sum(r[2] == "gap" for r in res)}, ambiguous {sum(r[2] == "ambiguous" for r in res)}' + (f'  {res}' if k <= 4 else ''))
print('=== secondary eclipse / any-phase box (primaries masked) ===')
T = np.concatenate([L[s][0] for s in L]); F = np.concatenate([L[s][1] for s in L]); ph = ((T - T0) / P + 0.5) % 1 - 0.5; m = np.abs(ph * P) > 0.25; T, F, ph = T[m], F[m], ph[m]
res = []
for p0 in np.arange(-0.4995, 0.5, 0.001):
    d = np.abs(((ph - p0 + 0.5) % 1) - 0.5) * P; win = d < HW * 0.85
    if win.sum() >= 20:
        res.append((p0, (1 - np.mean(F[win])) * 1e6, win.sum(), len(np.unique(np.round((T[win] - T0) / P - p0)))))
res = np.array(res); d = res[:, 1]; e_emp = 1.4826 * np.median(np.abs(d - np.median(d))); j5 = np.argmin(np.abs(np.abs(res[:, 0]) - 0.5))
print(f'phase 0.5: depth {res[j5, 1]:.0f} ppm (npts {int(res[j5, 2])}, epochs {int(res[j5, 3])}); empirical 1-sigma over phases {e_emp:.0f} ppm; 3-sigma upper limit {max(res[j5, 1], 0) + 3 * e_emp:.0f} ppm')
j = np.argmax(d); print(f'deepest box anywhere off-transit: phase {res[j, 0]:+.4f} depth {d[j]:.0f} ppm ({d[j] / e_emp:.1f} sigma), epochs {int(res[j, 3])}')
print('=== individual phase-0.5 epochs (box of width 0.8 x T14) ===')
pe = []
for s, (t, f, cad) in L.items():
    for n in range(int(np.floor((t.min() - T0) / P)) - 1, int(np.floor((t.max() - T0) / P)) + 1):
        tp = T0 + (n + 0.5) * P
        if t.min() < tp < t.max():
            dd, ee, nn = box(t, f, tp); cov = (np.abs(t - tp) < HW).sum() * cad / (2 * HW)
            print(f'  S{s} cycle {n + 0.5:+.1f} BTJD {tp:.3f}: coverage {cov:.2f}, depth {dd * 1e6:.0f} +- {ee * 1e6:.0f} ppm')
            if np.isfinite(dd) and cov >= 0.7:
                pe.append((dd, ee))
if pe:
    wv = np.array([1 / e ** 2 for _, e in pe]); m_ = np.sum([d * w_ for (d, _), w_ in zip(pe, wv)]) / wv.sum(); e_ = wv.sum() ** -0.5
    print(f'weighted mean of {len(pe)} phase-0.5 epochs: {m_ * 1e6:.0f} +- {e_ * 1e6:.0f} ppm; 3-sigma white-noise limit {(max(m_, 0) + 3 * e_) * 1e6:.0f} ppm')
if args.harmonics:
    print('=== phase-curve harmonics (out-of-transit) ===')

    def lcs(wl, inj=False):
        T_, F_, S_ = [], [], []
        for s in c['sectors']:
            z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); m = (z['q'] == 0) & np.isfinite(z['f2.0']); t = z['tm'][m]; f = z['f2.0'][m] / np.nanmedian(z['f2.0'][m])
            if inj:
                phx = ((t - T0) / P) % 1; f = f * (1 - 200e-6 * np.cos(4 * np.pi * phx) + 100e-6 * np.sin(2 * np.pi * phx))
            fl = flatten(t, f, method='biweight', window_length=wl, break_tolerance=0.3, edge_cutoff=0.0); ok = np.isfinite(fl)
            T_.append(t[ok]); F_.append(fl[ok]); S_.append(np.full(ok.sum(), s))
        T_, F_, S_ = np.concatenate(T_), np.concatenate(F_), np.concatenate(S_); ph_ = ((T_ - T0) / P + 0.5) % 1 - 0.5; m = np.abs(ph_ * P) > 0.2
        return T_[m], F_[m], S_[m]

    def harm(Tx, Fx, Sx, label):
        phx = ((Tx - T0) / P) % 1; A = np.array([np.cos(2 * np.pi * phx), np.sin(2 * np.pi * phx), np.cos(4 * np.pi * phx), np.sin(4 * np.pi * phx)] + [(Sx == s).astype(float) for s in np.unique(Sx)]).T
        cf = np.linalg.lstsq(A, Fx, rcond=None)[0]; rng = np.random.default_rng(0); ch = np.floor(Tx).astype(int); uc = np.unique(ch); bs = []
        for _ in range(200):
            idx = np.concatenate([np.nonzero(ch == u_)[0] for u_ in rng.choice(uc, len(uc))]); bs.append(np.linalg.lstsq(A[idx], Fx[idx], rcond=None)[0][:4])
        e = np.array(bs).std(0)
        print(f'  {label}: cos(phi) {cf[0] * 1e6:+.0f}+-{e[0] * 1e6:.0f}  sin(phi) {cf[1] * 1e6:+.0f}+-{e[1] * 1e6:.0f}  cos(2phi) {cf[2] * 1e6:+.0f}+-{e[2] * 1e6:.0f}  sin(2phi) {cf[3] * 1e6:+.0f}+-{e[3] * 1e6:.0f} ppm')
        return cf[:4] * 1e6

    for wl in (1.0, 5.0):
        h0 = harm(*lcs(wl), f'biweight {wl} d'); h1 = harm(*lcs(wl, inj=True), f'biweight {wl} d, injected -200 ppm cos(2phi), +100 ppm sin(phi)')
        print(f'  recovered after {wl}-d detrending: cos(2phi) {h1[2] - h0[2]:+.0f} ppm of -200, sin(phi) {h1[1] - h0[1]:+.0f} ppm of +100')
