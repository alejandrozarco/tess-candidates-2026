"""T2: ephemeris, shape and per-sector depth of the ~1.0492-d contaminating eclipses (Cadieux+2025 Sect. 2).
Masks b, c, d and flares, detrends with a running median, then grid-searches P, T0, duration with a box. Own code."""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
from scipy.ndimage import median_filter
T, F, S = [], [], []
for lc in load_lc('2m'):
    g = lc['quality'] == 0
    t = lc['time'][g]; f = lc['pdcsap_flux'][g]; ok = np.isfinite(f); t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    m = np.ones(len(t), bool)
    for q, (P, T0, D) in PL.items():
        ph = (t - T0 + 0.5 * P) % P - 0.5 * P; m &= np.abs(ph) > (0.5 * D / 24 + 0.04)
    # running median over 0.5 d (360 cadences), computed on unmasked points
    tr = np.interp(t, t[m], median_filter(f[m], size=361, mode='nearest'))
    r = f / tr - 1
    s = 1.4826 * np.median(np.abs(r[m]))
    bad = r > 4 * s; bad |= np.roll(bad, 1) | np.roll(bad, -1) | np.roll(bad, 2)
    m &= ~bad
    T.append(t[m]); F.append(r[m]); S.append(np.full(m.sum(), lc['meta']['sector']))
T = np.concatenate(T); F = np.concatenate(F); S = np.concatenate(S)
print('points', len(T))
def boxdepth(P, T0, D):
    ph = (T - T0 + 0.5 * P) % P - 0.5 * P; inn = np.abs(ph) < D / 2
    return -np.mean(F[inn]), inn.sum()
from astropy.timeseries import BoxLeastSquares
bls = BoxLeastSquares(T, F, dy=np.full(len(F), np.std(F)))
Pg = np.arange(1.0488, 1.0497, 0.000002)
res = bls.power(Pg, np.arange(0.8, 2.41, 0.1) / 24, objective='snr')
j = np.argmax(res.power); P, T0, D, dep, snr = res.period[j], res.transit_time[j], res.duration[j], res.depth[j], res.power[j]
print('BLS: P %.6f T0 %.4f D %.2f h depth %.0f ppm snr %.1f' % (P, T0, D * 24, dep * 1e6, snr))
# local refinement on a fine grid with the box statistic
grid = []
for Pq in P + np.arange(-6e-6, 6.01e-6, 5e-7):
    for T0q in T0 + np.arange(-0.006, 0.0061, 0.0005):
        for Dq in D + np.arange(-0.2, 0.21, 0.05) / 24:
            dq, n = boxdepth(Pq, T0q, Dq); grid.append((dq / (np.std(F) / np.sqrt(n)), Pq, T0q, Dq, dq))
snr, P, T0, D, dep = max(grid)
print('refined snr %.1f P %.6f T0 %.4f D %.2f h depth %.0f ppm' % (snr, P, T0, D * 24, dep * 1e6))
# move T0 to middle of baseline
n_mid = np.round((np.median(T) - T0) / P); T0m = T0 + n_mid * P
per = {}
for s in np.unique(S):
    sel = S == s
    ph = (T[sel] - T0 + 0.5 * P) % P - 0.5 * P; inn = np.abs(ph) < D / 2
    e = np.std(F[sel][~inn]) / np.sqrt(max(inn.sum(), 1))
    per[int(s)] = (float(-np.mean(F[sel][inn]) * 1e6), float(e * 1e6), int(inn.sum()))
    print('S%03d depth %6.0f +- %4.0f ppm (n=%d)' % (s, *per[int(s)]))
# fold profile in 10-min bins
ph = ((T - T0 + 0.5 * P) % P - 0.5 * P) * 24
bins = np.arange(-3, 3.01, 1 / 6.); prof = [(0.5 * (a + b), 1e6 * np.mean(F[(ph >= a) & (ph < b)])) for a, b in zip(bins[:-1], bins[1:])]
json.dump(dict(P=P, T0=T0m, D_h=D * 24, depth_ppm=dep * 1e6, snr=snr, per_sector=per, profile=prof), open(ROOT + '/out/eb_ephem.json', 'w'), indent=1)
print('profile (h, ppm):', ' '.join('%.2f:%.0f' % p for p in prof))
