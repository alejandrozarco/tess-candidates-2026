"""C2: model limb darkening per wavelength bin from Kurucz ATLAS9 specific intensities (Castelli & Kurucz 2003 grid,
[M/H] = -0.5, vturb 2 km/s, new ODF; file im05k2new.pck, 17 mu values, plane-parallel;
http://kurucz.harvard.edu/grids/gridm05/im05k2new.pck, 31659975 bytes, retrieved 2026-10-21). Own code.
L 98-59: Teff 3415 K, log g 4.91, [Fe/H] -0.46 (Cadieux et al. 2025). The grid starts at 3500 K, so coefficients are
linearly extrapolated in Teff from 3500 and 3750 K (log g 5.0); the 3500 K values and the log g 4.5 values are kept
as a systematic spread. Within a bin, intensities are photon-weighted (I_nu / lambda, flat throughput). Quadratic law
fitted by least squares to I(mu)/I(1) interpolated on 100 points uniform in mu.
Usage: c2_ld.py <im05k2new.pck> <out.json> ; bins are defined in binsdef.py"""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from binsdef import BINS_UM
MU = np.array([1.000, .900, .800, .700, .600, .500, .400, .300, .250, .200, .150, .125, .100, .075, .050, .025, .010])
def read_block(fn, teff, logg):
    L = open(fn).read().split('\n'); out = []
    for i, s in enumerate(L):
        if s.startswith('EFF') or s.startswith('TEFF'):
            t = float(s.split()[1]); g = float(s.split()[3])
            if abs(t - teff) < 1 and abs(g - logg) < 1e-3:
                j = i + 1
                while j < len(L) and not (L[j].startswith('EFF') or L[j].startswith('TEFF')):
                    p = L[j].split()
                    if len(p) == 18:
                        try: out.append([float(p[0]), float(p[1])] + [float(x) for x in p[2:]])
                        except ValueError: pass
                    j += 1
                break
    a = np.array(out); lam_um = a[:, 0] / 1000.; I1 = a[:, 1]
    ratio = np.hstack([np.ones((len(a), 1)), a[:, 2:] / 1e5])          # I(mu)/I(1) at the 17 mu values
    return lam_um, I1, ratio
def fit_quad(prof):
    g = np.linspace(0.01, 1, 100); Ig = np.interp(g, MU[::-1], prof[::-1])
    A = np.vstack([-(1 - g), -(1 - g) ** 2]).T; u = np.linalg.lstsq(A, Ig - 1, rcond=None)[0]
    return [float(u[0]), float(u[1])]
def band_ld(model, lo, hi):
    lam, I1, R = model; s = (lam >= lo) & (lam <= hi)
    if s.sum() < 1:   # bin narrower than the grid spacing: nearest point
        s = np.zeros(len(lam), bool); s[np.argmin(np.abs(lam - 0.5 * (lo + hi)))] = True
    w = I1[s] / lam[s]; prof = np.sum(R[s] * w[:, None], 0) / np.sum(w)
    return fit_quad(prof)
fn, outp = sys.argv[1], sys.argv[2]
models = {(T, g): read_block(fn, T, g) for T in (3500, 3750) for g in (4.5, 5.0)}
for k, m in models.items(): print('model', k, 'n lambda', len(m[0]), 'range %.3f-%.1f um' % (m[0].min(), m[0].max()))
res = {}
bands = dict(BINS_UM); bands['TESS_0.6-1.0'] = [(0.6, 1.0)]; bands['K_2.0-2.4'] = [(2.0, 2.4)]
bands['NRS1_white'] = [(2.87, 3.71)]; bands['NRS2_white'] = [(3.83, 5.15)]
for name, edges in bands.items():
    rows = []
    for lo, hi in edges:
        u = {k: band_ld(m, lo, hi) for k, m in models.items()}
        ext = [u[(3500, 5.0)][j] + (3415 - 3500) / 250. * (u[(3750, 5.0)][j] - u[(3500, 5.0)][j]) for j in (0, 1)]
        sp = [max(abs(ext[j] - u[(3500, 5.0)][j]), abs(u[(3500, 4.5)][j] - u[(3500, 5.0)][j])) for j in (0, 1)]
        rows.append(dict(lo=lo, hi=hi, u1=ext[0], u2=ext[1], u1_3500=u[(3500, 5.0)][0], u2_3500=u[(3500, 5.0)][1], spread=sp))
        print('%-14s %.3f-%.3f  u1 %.3f u2 %.3f  (3500 K: %.3f %.3f; spread %.3f %.3f)' % (name, lo, hi, ext[0], ext[1], u[(3500, 5.0)][0], u[(3500, 5.0)][1], sp[0], sp[1]))
    res[name] = rows
json.dump(dict(source='Kurucz ATLAS9 im05k2new.pck ([M/H]=-0.5), extrapolated to 3415 K, log g 5.0', bins=res), open(outp, 'w'), indent=1)
