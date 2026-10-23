"""TOI-1117 c (TIC 295541511), README section 6.10: circular transit fits (batman + emcee) to the detrended windows.

usage: python scripts/toi1117/mcmc_toi1117.py [--work ./work] [--steps 8000] [--configs free_2min_pdc,rho_2min_pdc,rho_all_pdc,rho_all_sap]
input:  <work>/295541511/lc_fit_{pdc,sap}.npz (fit_toi1117.py)
output: data/fit_295541511.csv  posterior medians and 16th/84th percentiles of every configuration
        <work>/295541511/chain_<config>.npz  thinned posterior samples (for ephemeris_toi1117.py, rv_toi1117.py)
Configurations: mode 'rho' = Gaussian prior on the stellar density from M = 0.97 +- 0.02 Msun, R = 1.06 +- 0.03 Rsun
(Lockley et al. 2025) and Gaussian priors on the quadratic limb-darkening coefficients (u1 = 0.40 +- 0.10, u2 = 0.21 +- 0.10);
mode 'free' = log-uniform density (the density of a circular orbit implied by the transit shape) and uniform q1, q2 (Kipping
2013). Data 2min = SPOC 2-min sectors 39, 100-103; all = also TESS-SPOC sector 13 (1800 s). Parameters: T_ref (near the
data-weighted centre, epoch 366 from BTJD 1670.1655), P, ln(Rp/R*), b, ln(rho*), q1, q2; circular orbit.
The model is averaged over each exposure at the centres of nsub equal parts (3 at 120 s, 45 at 1800 s); it is evaluated only
within the transit duration of each sample plus one cadence and is 1 elsewhere. 32 walkers; the first third of the steps
is discarded. Seeded (numpy default_rng(3)).
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import nsub, import_batman, a_rs, rho_prior, phase, write_csv, workdir, RS, RS_E, DATE, PC
batman = import_batman()
import emcee

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--steps', type=int, default=8000)
p.add_argument('--configs', default='free_2min_pdc,rho_2min_pdc,rho_all_pdc,rho_all_sap'); args = p.parse_args()
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data')
P0 = PC; TREF0 = 1670.1655 + 366 * PC
RHO_MU, RHO_SD = rho_prior()


class Fit:
    def __init__(self, mode, which, kind):
        z = np.load(os.path.join(workdir(args.work), f'lc_fit_{kind}.npz')); t, f, e, sec, cad = z['t'], z['f'], z['e'], z['sec'], z['cad']
        if which == '2min':
            k = sec != 13; t, f, e, sec, cad = t[k], f[k], e[k], sec[k], cad[k]
        self.mode, self.t, self.f, self.e, self.sec, self.cad = mode, t, f, e, sec, cad
        self.groups = []                                    # (index array, n, offsets) per cadence
        for c in np.unique(np.round(cad, 6)):
            idx = np.where(np.round(cad, 6) == c)[0]; c = float(np.median(cad[idx])); n = nsub(c); off = ((np.arange(n) + 0.5) / n - 0.5) * c
            self.groups.append((idx, n, off, c))
        pm = batman.TransitParams(); pm.t0 = 0.; pm.per = P0; pm.rp = 0.02; pm.a = 10.; pm.inc = 89.; pm.ecc = 0.; pm.w = 90.
        pm.limb_dark = 'quadratic'; pm.u = [0.4, 0.2]; self.pm = pm

    def unpack(self, th):
        T0, P, lnr, b, lnrho, q1, q2 = th
        a = a_rs(np.exp(lnrho), P); sq = np.sqrt(q1)
        return T0, P, np.exp(lnr), b, a, [2 * sq * q2, sq * (1 - 2 * q2)]

    def model(self, th):
        T0, P, k, b, a, u = self.unpack(th)
        pm = self.pm; pm.per = P; pm.rp = k; pm.a = a; pm.inc = np.degrees(np.arccos(min(b / a, 1.0))); pm.u = u
        half = P / np.pi * np.arcsin(min(1.0, (1 + k) / a)) / 2 * 1.05
        m = np.ones(len(self.t)); ph = phase(self.t, P, T0)
        for idx, n, off, c in self.groups:
            w = idx[np.abs(ph[idx]) < half + c]
            if len(w) == 0:
                continue
            tg = (ph[w][:, None] + off[None, :]).ravel()
            m[w] = batman.TransitModel(pm, tg).light_curve(pm).reshape(-1, n).mean(1)
        return m

    def chi2(self, th):
        return float(np.sum(((self.f - self.model(th)) / self.e) ** 2))

    def lnp(self, th):
        T0, P, lnr, b, lnrho, q1, q2 = th
        if not (abs(P - P0) < 0.01 and abs(T0 - TREF0) < 0.1 and -6 < lnr < -2 and 0 <= b < 1 + np.exp(lnr) and -3 < lnrho < 3 and 0 < q1 < 1 and 0 < q2 < 1):
            return -np.inf
        lp = 0.0
        if self.mode == 'rho':
            sq = np.sqrt(q1); u1, u2 = 2 * sq * q2, sq * (1 - 2 * q2)
            lp += -0.5 * ((np.exp(lnrho) - RHO_MU) / RHO_SD) ** 2 - 0.5 * (((u1 - 0.40) / 0.1) ** 2 + ((u2 - 0.21) / 0.1) ** 2)
        return lp - 0.5 * self.chi2(th)


rows = []
for cfg in args.configs.split(','):
    mode, which, kind = cfg.split('_'); F = Fit(mode, which, kind)
    p0 = np.array([TREF0, P0, np.log(0.021), 0.5, np.log(RHO_MU), 0.36, 0.3]); nw, nd = 32, 7; rng = np.random.default_rng(3)
    start = p0 + np.array([0.002, 3e-5, 0.05, 0.2, 0.2, 0.05, 0.05]) * rng.standard_normal((nw, nd)); start[:, 3] = np.abs(start[:, 3])
    S = emcee.EnsembleSampler(nw, nd, F.lnp); S.random_state = np.random.RandomState(3).get_state()
    S.run_mcmc(start, args.steps, progress=False)
    ch = S.get_chain(discard=args.steps // 3, thin=10, flat=True); lpost = S.get_log_prob(discard=args.steps // 3, thin=10, flat=True)
    try:
        tau = S.get_autocorr_time(discard=args.steps // 3, quiet=True)
    except Exception:
        tau = np.full(nd, np.nan)
    T0, P, lnr, b, lnrho, q1, q2 = ch.T; k = np.exp(lnr); rho = np.exp(lnrho); a = a_rs(rho, P); inc = np.degrees(np.arccos(b / a))
    T14 = P / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (a * np.sin(np.radians(inc))), 0, 1)) * 24
    Rp = k * RS * 109.076 * (1 + RS_E / RS * np.random.default_rng(5).standard_normal(len(k)))
    sq = np.sqrt(q1); u1, u2 = 2 * sq * q2, sq * (1 - 2 * q2)
    # maximum-likelihood sample (chi2 without priors) and the flat model
    i_ml = int(np.argmin([F.chi2(x) for x in ch[np.argsort(lpost)[-400:]]])); best = ch[np.argsort(lpost)[-400:]][i_ml]
    chi2_best = F.chi2(best); chi2_flat = float(np.sum(((F.f - 1) / F.e) ** 2))
    np.savez(os.path.join(workdir(args.work), f'chain_{cfg}.npz'), chain=ch, lnpost=lpost, best=best, names=np.array(['Tref', 'P', 'ln_k', 'b', 'ln_rho', 'q1', 'q2']))
    pc = lambda x: np.percentile(x, [16, 50, 84])
    cov = np.cov(T0, P)
    for name, x, fmt in (('Tref_btjd', T0, '.5f'), ('P_d', P, '.7f'), ('k_rp_rs', k, '.4f'), ('b', b, '.3f'), ('rho_star_gcc', rho, '.2f'), ('a_rs', a, '.2f'),
                         ('inc_deg', inc, '.2f'), ('T14_h', T14, '.3f'), ('depth_k2_ppm', k ** 2 * 1e6, '.0f'), ('Rp_earth', Rp, '.2f'), ('u1', u1, '.2f'), ('u2', u2, '.2f')):
        lo, me, hi = pc(x); rows.append((cfg, name, format(me, fmt), format(me - lo, fmt), format(hi - me, fmt), format(best[0], '.5f') if name == 'Tref_btjd' else ''))
    rows.append((cfg, 'sigma_Tref_min', f'{np.std(T0) * 1440:.2f}', '', '', ''))
    rows.append((cfg, 'sigma_P_s', f'{np.std(P) * 86400:.3f}', '', '', ''))
    rows.append((cfg, 'cov_Tref_P_d2', f'{cov[0, 1]:.4e}', '', '', ''))
    rows.append((cfg, 'P_rho_lt_prior_mean', f'{np.mean(rho < RHO_MU):.3f}', '', '', ''))
    rows.append((cfg, 'chi2_maxlike', f'{chi2_best:.2f}', '', '', ''))
    rows.append((cfg, 'dchi2_flat_minus_maxlike', f'{chi2_flat - chi2_best:.2f}', '', '', ''))
    rows.append((cfg, 'n_points', str(len(F.t)), '', '', ''))
    rows.append((cfg, 'acceptance', f'{np.mean(S.acceptance_fraction):.3f}', '', '', ''))
    rows.append((cfg, 'tau_steps_max', f'{np.nanmax(tau):.0f}', '', '', ''))
    rows.append((cfg, 'steps', str(args.steps), '', '', ''))
    print(cfg, 'P', pc(P), 'b', pc(b), 'T14', pc(T14), 'k', pc(k), 'rho', pc(rho), 'chi2', chi2_best, 'dchi2 flat', chi2_flat - chi2_best, 'tau', tau, flush=True)

mu, sd = RHO_MU, RHO_SD
write_csv(os.path.join(DATA, 'fit_295541511.csv'),
          [f'TOI-1117 c (TIC 295541511): circular transit fits (batman, quadratic limb darkening, emcee); scripts/toi1117/mcmc_toi1117.py, {DATE}',
           'config = mode_data_flux: mode rho = stellar-density prior N(%.3f, %.3f) g/cm^3 from M = 0.97 +- 0.02 Msun and R = 1.06 +- 0.03 Rsun (Lockley et al. 2025) and limb-darkening priors u1 = 0.40 +- 0.10, u2 = 0.21 +- 0.10; mode free = density log-uniform, q1 and q2 uniform' % (mu, sd),
           'data 2min = SPOC 2-min sectors 39, 100-103; all = also TESS-SPOC sector 13 (1800 s); flux pdc = PDCSAP_FLUX, sap = SAP_FLUX (SPOC and TESS-SPOC), detrended as in fit_toi1117.py, windows +-0.45 d, planet b masked',
           'model averaged over each exposure at the centres of nsub equal parts (midpoint sampling); circular orbit; per-point error = robust out-of-transit scatter of each sector',
           'value = posterior median; err_lo, err_hi = distances to the 16th and 84th percentiles; Rp_earth includes the 0.03 Rsun error of R*; depth_k2_ppm = (Rp/R*)^2 (not the observed depth)',
           'T14_h = first-to-fourth-contact duration; chi2_maxlike = chi2 (without priors) of the maximum-likelihood sample among the 400 highest-posterior samples; dchi2_flat_minus_maxlike = chi2(no transit) - chi2_maxlike',
           'maxlike_Tref_btjd = T_ref of that sample; tau_steps_max = longest emcee autocorrelation time (steps) of the seven parameters; BTJD = BJD_TDB - 2457000'],
          ['config', 'quantity', 'value', 'err_lo', 'err_hi', 'maxlike_Tref_btjd'], rows)
