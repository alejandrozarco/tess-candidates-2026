"""TOI-1117 c (TIC 295541511), README section 6.10: adopted ephemeris, per-sector mid-transit times, predicted transits,
and the folded light curve and models for the figure.

usage: python scripts/toi1117/ephemeris_toi1117.py [--work ./work] [--adopt rho_2min_pdc] [--alt rho_all_pdc] [--free free_2min_pdc]
input:  <work>/295541511/chain_<config>.npz (mcmc_toi1117.py), <work>/295541511/lc_fit_pdc.npz (fit_toi1117.py)
output: data/ephemeris_295541511.csv         adopted and alternative linear ephemerides with covariances
        data/transit_times_295541511.csv     per-sector mid-transit times (shape fixed), O - C, depth scale
        data/predicted_transits_295541511.csv predicted transits 2026-11-01 to 2027-12-31 (UTC, geocentric)
        data/folded_295541511.csv            detrended windows (PDCSAP) with phase and the adopted model at each cadence
        data/model_295541511.csv             model curves (adopted fit at 2-min and 30-min exposure; free-density fit)
Per-sector mid-times: the shape is fixed at the maximum-likelihood sample of the adopted fit; the mid-time offset of all
transits of a sector is scanned in 0.5-min steps (chi2 + 1 interval), then the depth scale at that offset. The model is averaged
over each exposure at the centres of nsub equal parts.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import nsub, import_batman, a_rs, phase, write_csv, workdir, RA, DEC, PB, TB, DB, DATE
batman = import_batman()
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--adopt', default='rho_2min_pdc')
p.add_argument('--alt', default='rho_all_pdc'); p.add_argument('--free', default='free_2min_pdc'); args = p.parse_args()
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
WD = workdir(args.work)
CO = SkyCoord(RA, DEC, unit='deg'); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)


def chain(cfg):
    z = np.load(os.path.join(WD, f'chain_{cfg}.npz')); return z['chain'], z['best']


def params(th):
    T0, P, lnr, b, lnrho, q1, q2 = th; a = a_rs(np.exp(lnrho), P); sq = np.sqrt(q1)
    pm = batman.TransitParams(); pm.t0 = 0.; pm.per = P; pm.rp = np.exp(lnr); pm.a = a; pm.inc = np.degrees(np.arccos(b / a)); pm.ecc = 0.; pm.w = 90.
    pm.limb_dark = 'quadratic'; pm.u = [2 * sq * q2, sq * (1 - 2 * q2)]; return pm, T0, P


def model(pm, dt, cad):
    """exposure-averaged model at offsets dt (d) from mid-transit"""
    n = nsub(cad); off = ((np.arange(n) + 0.5) / n - 0.5) * cad; tg = (np.asarray(dt)[:, None] + off[None, :]).ravel()
    return batman.TransitModel(pm, tg).light_curve(pm).reshape(-1, n).mean(1)


def t14(pm):
    inc = np.radians(pm.inc); b = pm.a * np.cos(inc)
    return pm.per / np.pi * np.arcsin(min(1, np.sqrt(max((1 + pm.rp) ** 2 - b ** 2, 0)) / (pm.a * np.sin(inc))))


def to_utc(bjd_tdb):
    tb = Time(bjd_tdb, format='jd', scale='tdb', location=GEO); tg = tb - tb.light_travel_time(CO, 'barycentric')
    tg = Time(tg.jd, format='jd', scale='tdb', location=GEO); tg = tb - tg.light_travel_time(CO, 'barycentric')
    return tg.utc


# ephemerides
eph = []
for cfg in (args.adopt, args.alt):
    ch, best = chain(cfg); T0s, Ps = ch[:, 0], ch[:, 1]; cov = np.cov(T0s, Ps)
    n101 = int(np.round((4115.17 - np.median(T0s)) / np.median(Ps))); t101 = T0s + n101 * Ps
    eph.append((cfg, 'adopted' if cfg == args.adopt else 'alternative', f'{np.median(T0s):.5f}', f'{np.median(T0s) - np.percentile(T0s, 16):.5f}', f'{np.percentile(T0s, 84) - np.median(T0s):.5f}',
                f'{np.std(T0s):.5f}', f'{np.median(Ps):.7f}', f'{np.median(Ps) - np.percentile(Ps, 16):.7f}', f'{np.percentile(Ps, 84) - np.median(Ps):.7f}', f'{np.std(Ps):.7f}',
                f'{cov[0, 1]:.4e}', f'{cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]):.3f}', n101, f'{np.median(t101):.5f}', f'{np.std(t101):.5f}'))
write_csv(os.path.join(DATA, 'ephemeris_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): linear ephemerides from the transit fits of mcmc_toi1117.py (circular, stellar-density prior); scripts/toi1117/ephemeris_toi1117.py, {DATE}',
           f'adopted: {args.adopt} (SPOC 2-min sectors 39, 100-103, PDCSAP); alternative: {args.alt} (also TESS-SPOC sector 13 at 1800 s)',
           'Tref_btjd = posterior median of the mid-transit time at the reference epoch (near the data-weighted centre); _lo/_hi = distances to the 16th/84th percentiles; _sd = posterior standard deviation',
           'cov_Tref_P_d2 and corr: posterior covariance and correlation of Tref and P; t(n) = Tref + n P, sigma(n)^2 = Tref_sd^2 + n^2 P_sd^2 + 2 n cov',
           'n_s101 = epoch (from Tref) of the transit near BTJD 4115.17 (sector 101); tc_s101 = its posterior median time and standard deviation; BTJD = BJD_TDB - 2457000'],
          ['config', 'role', 'Tref_btjd', 'Tref_lo_d', 'Tref_hi_d', 'Tref_sd_d', 'P_d', 'P_lo_d', 'P_hi_d', 'P_sd_d', 'cov_Tref_P_d2', 'corr', 'n_s101', 'tc_s101_btjd', 'tc_s101_sd_d'], eph)

# per-sector mid-times with the adopted shape
ch, best = chain(args.adopt); pm, T0b, Pb = params(best)
z = np.load(os.path.join(WD, 'lc_fit_pdc.npz')); t, f, e, sec, cad = z['t'], z['f'], z['e'], z['sec'], z['cad']
rows = []
for s in np.unique(sec):
    k = sec == s; c = float(np.median(cad[k])); nm = int(np.round((np.median(t[k]) - T0b) / Pb)); tc_lin = T0b + nm * Pb
    sig_lin = np.std(ch[:, 0] + nm * ch[:, 1]) * 1440
    ph = phase(t[k], Pb, T0b)
    dts = np.arange(-0.06, 0.06001, 0.0005 / 1.44)                      # 0.5-min steps
    chi = np.array([np.sum(((f[k] - model(pm, ph - d, c)) / e[k]) ** 2) for d in dts]); i = int(np.argmin(chi)); ok = chi - chi[i] <= 1.0
    m0 = model(pm, ph - dts[i], c); sc = np.linspace(0, 2.5, 501)
    cs = np.array([np.sum(((f[k] - (1 - x * (1 - m0))) / e[k]) ** 2) for x in sc]); j = int(np.argmin(cs)); okk = cs - cs[j] <= 1
    rows.append((int(s), f'{c * 1440:.1f}', nm, f'{tc_lin + dts[i]:.5f}', f'{dts[i] * 1440:+.1f}', f'{(dts[i] - dts[ok].min()) * 1440:.1f}', f'{(dts[ok].max() - dts[i]) * 1440:.1f}',
                 f'{sig_lin:.1f}', f'{sc[j]:.2f}', f'{(sc[okk].max() - sc[okk].min()) / 2:.2f}', f'{np.sum(((f[k] - 1) / e[k]) ** 2) - cs[j]:.1f}', 0 if s == 13 else 1))
    print(rows[-1], flush=True)
write_csv(os.path.join(DATA, 'transit_times_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): one mid-transit time per sector (all transits of the sector together), shape fixed at the maximum-likelihood sample of {args.adopt}; scripts/toi1117/ephemeris_toi1117.py, {DATE}',
           f'epoch counted from Tref = BTJD {T0b:.5f} (maximum-likelihood sample); o_c_min = offset of the best mid-time from that sample\'s linear ephemeris (P = {Pb:.7f} d); err_lo/hi = chi2 + 1 interval (formal, white errors)',
           'eph_err_min = posterior standard deviation of the adopted linear ephemeris at that epoch; depth_scale = depth relative to the adopted model at the best mid-time (chi2 + 1 half-width)',
           'dchi2_null = chi2(no transit) - chi2(scaled model); in_fit = 1 if the sector is in the adopted fit (sector 13, 1800 s, is not); model averaged over each exposure (midpoint sampling)',
           'injection tests (checks_295541511.csv, timing rows) show that formal timing errors of single sectors underestimate the scatter of recovered mid-times by factors of about 2-3'],
          ['sector', 'cadence_min', 'epoch', 'tc_btjd', 'o_c_min', 'err_lo_min', 'err_hi_min', 'eph_err_min', 'depth_scale', 'depth_scale_err', 'dchi2_null', 'in_fit'], rows)

# predictions
T0s, Ps = ch[:, 0], ch[:, 1]; T0m, Pm = float(np.median(T0s)), float(np.median(Ps)); T14 = t14(pm)
start = Time('2026-11-01T00:00:00', scale='utc').tdb.jd - 2457000; stop = Time('2028-01-01T00:00:00', scale='utc').tdb.jd - 2457000
pred = []
for n in range(int(np.ceil((start - T0m) / Pm)), int(np.floor((stop - T0m) / Pm)) + 1):
    tt = T0s + n * Ps; mu, sd = float(np.median(tt)), float(np.std(tt))
    ut = to_utc(np.array([mu, mu - T14 / 2, mu + T14 / 2]) + 2457000)
    nb = np.round((mu - TB) / PB); dtb = (mu - (TB + nb * PB)) * 24
    pred.append((n, f'{mu:.5f}', f'{mu + 2457000:.5f}', ut[0].iso[:16], f'{sd * 1440:.1f}', ut[1].iso[:16], ut[2].iso[:16], f'{dtb:+.2f}',
                 int(abs(dtb) < (T14 * 24 + DB * 24) / 2)))
write_csv(os.path.join(DATA, 'predicted_transits_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): predicted mid-transit times 2026-11-01 to 2027-12-31 from the adopted ephemeris ({args.adopt}): Tref = BTJD {T0m:.5f}, P = {Pm:.7f} d; scripts/toi1117/ephemeris_toi1117.py, {DATE}',
           'err_min = formal 1-sigma from the posterior (Tref, P covariance); the formal error does not include timing variations or the underestimate of TESS timing errors (factor 2-3 in injections); an approximate planning allowance is +-10 min; its statistical coverage is not calibrated',
           f'*_utc = geocentric UTC (barycentric light-travel time toward the target removed; a topocentric observer differs by < 0.1 s); ingress/egress = mid-time -/+ T14/2 with T14 = {T14 * 24:.2f} h (maximum-likelihood sample of the adopted fit; the posterior-median T14 is in fit_295541511.csv)',
           f'dt_b_h = time from the nearest predicted transit of planet b (P = {PB:.7f} d, T0 = BTJD {TB:.5f}, T14 = {DB * 24:.2f} h); b_overlap = 1 if the two transits overlap in time',
           'TESS coverage of these dates is not checked'],
          ['epoch', 'tc_btjd', 'tc_bjd_tdb', 'tc_utc', 'err_min', 'ingress_utc', 'egress_utc', 'dt_b_h', 'b_overlap'], pred)
print(len(pred), pred[0], pred[-1])

# folded light curve and model curves for the figure
ph = phase(t, Pb, T0b); mod = np.ones(len(t))
for c in np.unique(np.round(cad, 6)):
    k = np.round(cad, 6) == c; mod[k] = model(pm, ph[k], float(np.median(cad[k])))
write_csv(os.path.join(DATA, 'folded_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): detrended PDCSAP windows (+-0.45 d) used in the transit fits; scripts/toi1117/ephemeris_toi1117.py from the output of fit_toi1117.py, {DATE}',
           'sectors 39, 100-103: SPOC 2-min; sector 13: TESS-SPOC FFI (1800 s); planet b masked; wotan biweight 1.0 d with the windows of c masked; 5-sigma clipped; quality flag 0 cadences only',
           f'dt_days = phase relative to the maximum-likelihood sample of {args.adopt} (Tref = BTJD {T0b:.5f}, P = {Pb:.7f} d); flux_err = robust out-of-transit scatter of the sector; model = that sample averaged over the exposure'],
          ['sector', 'cadence_s', 'time_btjd', 'dt_days', 'flux', 'flux_err', 'model'],
          [(int(s), int(round(c * 86400)), f'{tt:.6f}', f'{x:.6f}', f'{y:.6f}', f'{ee:.6f}', f'{m:.7f}') for s, c, tt, x, y, ee, m in zip(sec, cad, t, ph, f, e, mod)])
dt = np.arange(-0.45, 0.45001, 0.0005)
chf, bestf = chain(args.free); pmf, _, _ = params(bestf)
pma, _, _ = params(chain(args.alt)[1])
write_csv(os.path.join(DATA, 'model_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): model curves for the figure; scripts/toi1117/ephemeris_toi1117.py, {DATE}',
           f'rho_2min, rho_30min: maximum-likelihood sample of {args.adopt} (k = {pm.rp:.4f}, b = {pm.a * np.cos(np.radians(pm.inc)):.3f}, a/R* = {pm.a:.2f}, u = [{pm.u[0]:.2f}, {pm.u[1]:.2f}]) averaged over a 120-s and an 1800-s exposure',
           f'free_2min: maximum-likelihood sample of {args.free} (k = {pmf.rp:.4f}, b = {pmf.a * np.cos(np.radians(pmf.inc)):.3f}, a/R* = {pmf.a:.2f}, u = [{pmf.u[0]:.2f}, {pmf.u[1]:.2f}]) at 120 s',
           f'all_30min: maximum-likelihood sample of {args.alt} (k = {pma.rp:.4f}, b = {pma.a * np.cos(np.radians(pma.inc)):.3f}, a/R* = {pma.a:.2f}) at 1800 s; dt_days = time from mid-transit'],
          ['dt_days', 'rho_2min', 'rho_30min', 'free_2min', 'all_30min'],
          [(f'{a:.4f}', f'{b:.7f}', f'{c:.7f}', f'{d:.7f}', f'{g:.7f}') for a, b, c, d, g in zip(dt, model(pm, dt, 120 / 86400), model(pm, dt, 1800 / 86400), model(pmf, dt, 120 / 86400), model(pma, dt, 1800 / 86400))])
