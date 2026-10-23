"""TOI-1117 (TIC 295541511), README section 6.10: refit of the published HARPS radial velocities (Lockley et al. 2025) and
comparison of the RV orbit of planet c with the transit ephemeris.

usage: python scripts/toi1117/rv_toi1117.py [--work ./work] [--steps 30000]
input:  <work>/295541511/dl/arxiv/2506.05521.tar.gz (download_toi1117.py; the arXiv source of Lockley et al. 2025, whose
        appendix table lists the HARPS velocities), data/ephemeris_295541511.csv, <work>/295541511/chain_*.npz
output: data/rv_harps_295541511.csv     the HARPS table as published (removed = 1: the two rows in bold, removed by the authors)
        data/rvfit_295541511.csv        posteriors of the two RV fits, and the mass and density of c
        data/rv_ephemeris_295541511.csv transit times against the RV-predicted conjunctions of c
Model: offset + jitter + three circular Keplerians (b, c, d), emcee (44 walkers, seeded); planet b's period and transit time
carry Gaussian priors from its transits (Lockley et al. 2025: P = 2.22816 +- 0.00001 d, T0 = BTJD 2386.964 +- 0.004).
Mode 'free': P_c in (4, 5) d and T_c in (2490, 2505) BTJD uniform (RV-only orbit of c). Mode 'transit': P_c and T_c with
Gaussian priors from the adopted transit ephemeris (effectively fixed) -> K_c. Time: RJD - 57000 = BTJD (RJD = BJD - 2400000).
Masses: M sin i = K (P / 2 pi G)^(1/3) M*^(2/3) (M_p << M*), M* = 0.97 +- 0.02 Msun; true mass with the inclination samples of the
transit fit; density with the planet radius samples of the transit fit (R* = 1.06 +- 0.03 Rsun).
"""
import argparse, os, re, sys, tarfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_csv, write_csv, workdir, dldir, a_rs, MS, MS_E, RS, RS_E, DATE
import emcee

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--steps', type=int, default=30000)
p.add_argument('--adopt', default='rho_2min_pdc'); p.add_argument('--free-fit', default='free_2min_pdc'); args = p.parse_args()
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')

# 1. the HARPS table from the arXiv source (read as text only)
with tarfile.open(os.path.join(dldir(args.work, 'arxiv'), '2506.05521.tar.gz')) as tf:
    tex = tf.extractfile(tf.getmember('main.tex')).read().decode('utf-8', 'replace')
rows = []
for l in tex.split('\n'):
    bold = '\\textbf' in l; s = re.sub(r'\\textbf\{([^}]*)\}', r'\1', l)
    m = re.match(r'\s*([56]\d{4}\.\d+)\s*&\s*(-?[\d.]+)\s*&\s*([\d.]+)\s*&\s*([\d.]+)\s*&\s*(-?[\d.]+)\s*&\s*([\d.]+)', s)
    if m:
        rows.append(list(m.groups()) + [int(bold)])
write_csv(os.path.join(DATA, 'rv_harps_295541511.csv'),
          ['TOI-1117 = TIC 295541511: HARPS radial velocities as published by Lockley et al. (2025, MNRAS 541, 919; table "tab:RVs" of the arXiv source 2506.05521)',
           f'transcribed by scripts/toi1117/rv_toi1117.py, {DATE}; rjd = BJD - 2400000 (as in the table); rv, rv_err, fwhm, bis in m/s; sindex dimensionless',
           'removed = 1: the two rows printed in bold, which the authors removed because of their large uncertainties; not used here either'],
          ['rjd', 'rv_ms', 'rv_err_ms', 'fwhm_ms', 'bis_ms', 'sindex', 'removed'], rows)
R = np.array([[float(x) for x in r[:3]] + [r[6]] for r in rows]); R = R[R[:, 3] == 0]
t, y, e = R[:, 0] - 57000.0, R[:, 1], R[:, 2]
print(len(rows), 'rows,', len(t), 'used')

E = read_csv(os.path.join(DATA, 'ephemeris_295541511.csv')); ia = list(E['role']).index('adopted')
Tph, Pph, sT, sP, cTP = E['Tref_btjd'][ia], E['P_d'][ia], E['Tref_sd_d'][ia], E['P_sd_d'][ia], E['cov_Tref_P_d2'][ia]
def sig_ph(n):
    return np.sqrt(sT ** 2 + n ** 2 * sP ** 2 + 2 * n * cTP)


def model(th, tt):
    g, lnj, Pb, Tb, Kb, Pc, Tc, Kc, Pd, Td, Kd = th
    return g - Kb * np.sin(2 * np.pi * (tt - Tb) / Pb) - Kc * np.sin(2 * np.pi * (tt - Tc) / Pc) - Kd * np.sin(2 * np.pi * (tt - Td) / Pd)


def lnp(th, mode):
    g, lnj, Pb, Tb, Kb, Pc, Tc, Kc, Pd, Td, Kd = th
    if not (-24900 < g < -24890 and -3 < lnj < 3 and 0 < Kb < 50 and 0 < Kc < 50 and 0 < Kd < 50 and 4 < Pc < 5 and 2490 < Tc < 2505 and 8 < Pd < 9 and 2485 < Td < 2505):
        return -np.inf
    lp = -0.5 * ((Pb - 2.22816) / 1e-5) ** 2 - 0.5 * ((Tb - 2386.964) / 0.004) ** 2
    if mode == 'transit':
        n = np.round((Tc - Tph) / Pph); lp += -0.5 * ((Pc - Pph) / sP) ** 2 - 0.5 * ((Tc - (Tph + n * Pph)) / sig_ph(n)) ** 2
    s2 = e ** 2 + np.exp(2 * lnj); r = y - model(th, t)
    return lp - 0.5 * np.sum(r ** 2 / s2 + np.log(2 * np.pi * s2))


G, MSUN, MEARTH, REARTH = 6.674e-11, 1.989e30, 5.972e24, 6.371e6
names = ['offset_ms', 'ln_jitter', 'P_b_d', 'T_b_btjd', 'K_b_ms', 'P_c_d', 'T_c_btjd', 'K_c_ms', 'P_d_d', 'T_d_btjd', 'K_d_ms']
fmt = ['.3f', '.2f', '.6f', '.4f', '.2f', '.5f', '.3f', '.2f', '.4f', '.3f', '.2f']
out, ch = [], {}
for mode in ('free', 'transit'):
    p0 = np.array([-24894.4, 1.5, 2.22816, 2386.964, 4.45, 4.579, 2498.9, 3.45, 8.665, 2495.6, 3.4]); nd, nw = len(p0), 44
    if mode == 'transit':
        p0[5] = Pph; p0[6] = Tph + np.round((2498.9 - Tph) / Pph) * Pph
    rng = np.random.default_rng(1); sc = np.array([0.2, 0.1, 1e-5, 0.003, 0.3, 0.002 if mode == 'free' else 1e-6, 0.3 if mode == 'free' else 0.001, 0.3, 0.01, 0.4, 0.3])
    S = emcee.EnsembleSampler(nw, nd, lnp, args=(mode,)); S.random_state = np.random.RandomState(1).get_state()
    S.run_mcmc(p0 + sc * rng.standard_normal((nw, nd)), args.steps, progress=False)
    c = S.get_chain(discard=args.steps // 3, thin=20, flat=True); ch[mode] = c
    try:
        tau = float(np.nanmax(S.get_autocorr_time(discard=args.steps // 3, quiet=True)))
    except Exception:
        tau = np.nan
    for i, (nm, f_) in enumerate(zip(names, fmt)):
        lo, me, hi = np.percentile(c[:, i], [16, 50, 84]); out.append((mode, nm, format(me, f_), format(me - lo, f_), format(hi - me, f_)))
    rs = np.random.default_rng(2); Ms = (MS + MS_E * rs.standard_normal(len(c))) * MSUN
    for pl, iK, iP in (('b', 4, 2), ('c', 7, 5), ('d', 10, 8)):
        m = c[:, iK] * (c[:, iP] * 86400 / (2 * np.pi * G)) ** (1 / 3) * Ms ** (2 / 3) / MEARTH; lo, me, hi = np.percentile(m, [16, 50, 84])
        out.append((mode, f'Msini_{pl}_earth', f'{me:.2f}', f'{me - lo:.2f}', f'{hi - me:.2f}'))
        if mode == 'transit' and pl == 'c':
            msini_c = m
    out.append((mode, 'acceptance', f'{np.mean(S.acceptance_fraction):.3f}', '', '')); out.append((mode, 'tau_steps_max', f'{tau:.0f}', '', ''))
    out.append((mode, 'n_rv', str(len(t)), '', '')); out.append((mode, 'K_c_over_err', f'{np.median(c[:, 7]) / np.std(c[:, 7]):.1f}', '', ''))
    print(mode, 'K_c', np.percentile(c[:, 7], [16, 50, 84]), 'P_c', np.percentile(c[:, 5], [16, 50, 84]), 'tau', tau, flush=True)

# 2. mass and density of c with the inclination and radius of the transit fits
for cfg, lab in ((args.adopt, 'adopted transit fit'), (args.free_fit, 'free-density transit fit')):
    z = np.load(os.path.join(workdir(args.work), f'chain_{cfg}.npz'))['chain']; rs = np.random.default_rng(4)
    j = rs.integers(0, len(z), len(msini_c)); T0, P, lnr, b, lnrho, q1, q2 = z[j].T
    inc = np.arccos(b / a_rs(np.exp(lnrho), P)); Rp = np.exp(lnr) * (RS + RS_E * rs.standard_normal(len(j))) * 6.957e8 / REARTH
    M = msini_c / np.sin(inc); rho = M * MEARTH / (4 / 3 * np.pi * (Rp * REARTH) ** 3) / 1000
    for nm, x in (('inc_deg', np.degrees(inc)), ('M_c_earth', M), ('R_c_earth', Rp), ('rho_c_gcc', rho)):
        lo, me, hi = np.percentile(x, [16, 50, 84]); out.append((f'transit + {lab} ({cfg})', nm, f'{me:.2f}', f'{me - lo:.2f}', f'{hi - me:.2f}'))
write_csv(os.path.join(DATA, 'rvfit_295541511.csv'),
          [f'TOI-1117 (TIC 295541511): circular three-planet fits to the published HARPS velocities (rv_harps_295541511.csv, removed = 0; {len(t)} points); scripts/toi1117/rv_toi1117.py, {DATE}',
           'mode free: orbit of c from the RVs alone; mode transit: P_c and T_c with Gaussian priors from the adopted transit ephemeris (ephemeris_295541511.csv); planet b: Gaussian priors on P and T0 from Lockley et al. (2025)',
           'value = posterior median, err_lo/err_hi = distances to the 16th/84th percentiles; T_* = time of inferior conjunction (BTJD = BJD_TDB - 2457000); ln_jitter in ln(m/s)',
           'Msini = minimum mass (M* = 0.97 +- 0.02 Msun); rows "transit + ...": inclination and radius samples of the named transit fit combined with the Msini_c samples of mode transit (true mass and density of c, conditional on that fit)',
           'tau_steps_max = longest autocorrelation time (steps); K_c_over_err = median / standard deviation of K_c'],
          ['mode', 'quantity', 'value', 'err_lo', 'err_hi'], out)

# 3. transit times against the RV-predicted conjunctions of c (free mode, full P-T covariance, matched epoch numbers)
c = ch['free']; Pc, Tc = c[:, 5], c[:, 6]
cmp_rows = []
nn = np.arange(-400, 300); sd_rv = np.array([np.std(Tc + np.round((Tph + n * Pph - Tc) / Pc) * Pc) for n in nn]); n_rv = int(nn[np.argmin(sd_rv)])
T = read_csv(os.path.join(DATA, 'transit_times_295541511.csv'))
for lab, n in [(f'sector {int(s)}', int(np.round((tc - Tph) / Pph))) for s, tc in zip(T['sector'], T['tc_btjd'])] + [('epoch of the smallest RV-predicted uncertainty', n_rv)]:
    tph = Tph + n * Pph; tt = Tc + np.round((tph - Tc) / Pc) * Pc; mu, sd = float(np.median(tt)), float(np.std(tt))
    cmp_rows.append((lab, n, f'{tph:.4f}', f'{sig_ph(n) * 1440:.1f}', f'{mu:.3f}', f'{sd:.3f}', f'{(tph - mu) / sd:+.2f}'))
cmp_rows.append(('period', '', f'{Pph:.7f}', f'{sP * 86400:.2f} s', f'{np.median(Pc):.5f}', f'{np.std(Pc):.5f}', f'{(Pph - np.median(Pc)) / np.std(Pc):+.2f}'))
write_csv(os.path.join(DATA, 'rv_ephemeris_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): adopted transit ephemeris (Tref = BTJD {Tph:.5f}, P = {Pph:.7f} d) against the RV-only orbit of c (rvfit_295541511.csv, mode free); scripts/toi1117/rv_toi1117.py, {DATE}',
           'epoch counted from Tref; transit_btjd = linear transit ephemeris at that epoch (transit_err_min its 1-sigma); rv_conj_btjd, rv_conj_sd_d = median and standard deviation of the RV-predicted conjunction at the same epoch (P_c - T_c samples)',
           'z = (transit - RV) / RV standard deviation; last row: periods (transit_err in seconds); the epochs listed are those of the per-sector transit times and the epoch where the RV prediction is most precise'],
          ['epoch_label', 'epoch', 'transit_btjd', 'transit_err_min', 'rv_conj_btjd', 'rv_conj_sd_d', 'z'], cmp_rows)
for r in cmp_rows:
    print(r)
