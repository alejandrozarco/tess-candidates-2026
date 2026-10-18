"""Joint transit fit with the target flux fraction as a free parameter (MCMC), and the photoeccentric constraint.

usage: python scripts/fit_free_dilution.py <tic> --fit <fit json> [--aperture 2.0] [--window 0.8] [--steps 3000]
                                              [--walkers 40] [--rho-prior] [--work ./work] [--out data]
input:  <work>/<tic>/cut/lc_s<sector>.npz (aperture_photometry.py) and <work>/<tic>/<fit json> (fit_transits.py; the
        per-epoch mid-transit times centre the windows and start the chains, the per-epoch residual rms is the fixed
        white noise), `star_sed` in targets.py (adopted R*, M*), `free_dilution` in targets.py (f_t prior, start values,
        stellar-density prior)
output: <work>/<tic>/free_dilution[_rho].json, <work>/<tic>/free_dilution[_rho]_chain.npy, <out>/fit_free_dilution[_rho]_<tic>.csv
        (suffix _rho with --rho-prior)
Model (emcee + batman): raw aperture flux of each epoch (quality 0, within +-window d of the fitted mid-time)
F = f_t T(t) + (1 - f_t), multiplied by a quadratic baseline per epoch that is solved linearly at every step. Free:
k = Rp/R* (undiluted), b, log10 rho_circ (uniform in -0.5 to 1.5; circular orbit, a/R* from rho_circ and P), Kipping
q1, q2 with Gaussian priors on u1, u2 (the limb darkening of targets.py, +- 0.1), f_t with a Gaussian prior from the
dilution_correction.py fractions (TIC 89247735: N(0.735, 0.03), 0.70-0.78 by epoch; TIC 455017950: N(0.966, 0.010)),
and one mid-time per epoch (+-0.03 d). With --rho-prior the circular-orbit density gets the Gaussian prior
`rho_prior` of targets.py (TIC 455017950: 2.40 +- 0.30 g cm^-3, applied to rho_circ with uniform-in-rho_circ weighting)
and the photoeccentric step is skipped. Exposure time = median cadence of each epoch (200 s or 600 s), supersampling
max(3, exposure / 120 s). The first third of each chain is discarded, the rest thinned by 10.
Rp = k R* with R* drawn from N(R, eR). Photoeccentric constraint (Dawson & Johnson 2012): g = (rho_circ/rho*)^(1/3)
with rho* = M/R^3 drawn from N(M, eM), N(R, eR); e and omega drawn uniformly in [0, 0.95) and [0, 2 pi) and weighted
by the kernel density of g at g(e, omega) = (1 + e sin omega)/sqrt(1 - e^2). Savage-Dickey ratio = posterior density
of e at e = 0 (first two bins of width 0.01) divided by the uniform prior density.
"""
import argparse, csv, json, os, time
import numpy as np
from scipy.stats import gaussian_kde
from targets import TARGETS, workdir, add_common_args, import_batman, sector_of

batman = import_batman()
import emcee

p = add_common_args(argparse.ArgumentParser(description=__doc__))
p.add_argument('--fit', required=True); p.add_argument('--aperture', default='2.0'); p.add_argument('--window', type=float, default=0.8)
p.add_argument('--steps', type=int, default=3000); p.add_argument('--walkers', type=int, default=40); p.add_argument('--out', default='data')
p.add_argument('--date', default=''); p.add_argument('--rho-prior', action='store_true', help='Gaussian stellar-density prior of targets.py')
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); st = c['star_sed']; FD = c['free_dilution']
P = c['P']; G = 6.674e-8; FT_MU, FT_SD = FD['ft']; U_MU, U_SD = c['u'], 0.1
RHO = FD['rho_prior'] if args.rho_prior else None; SUF = '_rho' if RHO else ''
FJ = json.load(open(os.path.join(wd, args.fit))); USE = list(FJ)
DATE = args.date or __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
D = {}
for s in USE:
    z = np.load(os.path.join(wd, 'cut', f'lc_s{sector_of(s)}.npz')); tc = FJ[s]['tc']
    t = np.asarray(z['tm'], dtype=float); f = np.asarray(z[f'f{args.aperture}'], dtype=float); q = z['q']   # native byte order for batman
    m = (np.abs(t - tc) < args.window) & np.isfinite(f) & (q == 0); t, f = t[m], f[m] / np.median(f[m])
    pp = batman.TransitParams(); pp.t0 = tc; pp.per = P; pp.rp = 0.1; pp.a = 33; pp.inc = 89.9; pp.ecc = 0; pp.w = 90; pp.u = [0.5, 0.2]; pp.limb_dark = 'quadratic'
    x = t - tc; exp_s = round(float(np.median(np.diff(t))) * 86400)
    D[s] = dict(t=t, f=f, p=pp, mod=batman.TransitModel(pp, t, supersample_factor=max(3, round(exp_s / 120)), exp_time=exp_s / 86400.), A=np.c_[np.ones(len(t)), x, x * x],
                tc0=tc, sig=FJ[s]['rms'] * 1e-6)
NAMES = ['k', 'b', 'log_rho_circ', 'q1', 'q2', 'f_t'] + [f'dtc_{s}' for s in USE]
print(f'{sum(len(d["t"]) for d in D.values())} points in epochs {USE}; noise per epoch (ppm)', {s: round(d['sig'] * 1e6) for s, d in D.items()})


def lnprior(th):
    k, b, lr, q1, q2, ft = th[:6]; dtc = th[6:]
    if not (0.02 < k < 0.3 and 0 <= b < 1 + k and -0.5 < lr < 1.5 and 0 < q1 < 1 and 0 < q2 < 1 and 0.5 < ft < 1.0 and np.all(np.abs(dtc) < 0.03)):
        return -np.inf
    u1 = 2 * np.sqrt(q1) * q2; u2 = np.sqrt(q1) * (1 - 2 * q2)
    lp = -0.5 * ((u1 - U_MU[0]) / U_SD) ** 2 - 0.5 * ((u2 - U_MU[1]) / U_SD) ** 2 - 0.5 * ((ft - FT_MU) / FT_SD) ** 2
    if RHO:
        lp += -0.5 * ((10 ** lr - RHO[0]) / RHO[1]) ** 2 + np.log(10 ** lr)
    return lp


def lnlike(th):
    k, b, lr, q1, q2, ft = th[:6]; dtc = th[6:]
    aR = (G * 10 ** lr * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3.); cosi = b / aR
    if cosi >= 1:
        return -np.inf
    u1 = 2 * np.sqrt(q1) * q2; u2 = np.sqrt(q1) * (1 - 2 * q2); L = 0.0
    for i, d in enumerate(D.values()):
        pp = d['p']; pp.t0 = d['tc0'] + dtc[i]; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(cosi)); pp.u = [u1, u2]
        mm = ft * d['mod'].light_curve(pp) + (1 - ft); A = d['A'] * mm[:, None]; cc = np.linalg.lstsq(A, d['f'], rcond=None)[0]
        L += -0.5 * np.sum(((d['f'] - A @ cc) / d['sig']) ** 2)
    return L


def lnprob(th):
    lp = lnprior(th)
    if not np.isfinite(lp):
        return -np.inf
    lk = lnlike(th)
    return lp + lk if np.isfinite(lk) else -np.inf


p0 = np.array([FD['p0'][0], FD['p0'][1], np.log10(FD['p0'][2]), FD['p0'][3], FD['p0'][4], FT_MU] + [0.0] * len(USE)); nd = len(p0); nw = args.walkers
sc = np.array([0.002, 0.1, 0.03, 0.05, 0.05, 0.01] + [0.0005] * len(USE)); rng = np.random.default_rng(7)
pos = p0 + sc * rng.normal(size=(nw, nd)); pos[:, 1] = np.abs(pos[:, 1])
for i in range(nw):
    while not np.isfinite(lnprob(pos[i])):
        pos[i] = p0 + sc * rng.normal(size=nd); pos[i, 1] = abs(pos[i, 1])
t0 = time.time(); S = emcee.EnsembleSampler(nw, nd, lnprob); S._random.seed(7); S.run_mcmc(pos, args.steps, progress=False)
ch = S.get_chain(discard=args.steps // 3, thin=10, flat=True); lp = S.get_log_prob(discard=args.steps // 3, thin=10, flat=True)
try:
    tau = [float(x) for x in S.get_autocorr_time(tol=0)]
except Exception:
    tau = None
PCT = [2.5, 16, 50, 84, 97.5]
q = lambda v: [float(x) for x in np.percentile(v, PCT)]
best = ch[np.argmax(lp)]
res = dict(nwalkers=nw, nsteps=args.steps, runtime_s=time.time() - t0, acceptance=float(np.mean(S.acceptance_fraction)), tau=tau,
           npts=int(sum(len(d['t']) for d in D.values())), best=dict(zip(NAMES, best.tolist())), best_chi2=float(-2 * lnlike(best)))
for i, n in enumerate(NAMES):
    res[n] = q(ch[:, i])
k, b, rc = ch[:, 0], ch[:, 1], 10 ** ch[:, 2]; aR = (G * rc * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3.); inc = np.arccos(b / aR)
res['depth_true_ppm'] = q(k * k * 1e6); res['rho_circ'] = q(rc); res['aR'] = q(aR); res['inc_deg'] = q(np.degrees(inc))
res['T14_h'] = q(P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b ** 2) / (aR * np.sin(inc))) * 24)
res['T23_h'] = q(P / np.pi * np.arcsin(np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))) * 24)
Rs = np.random.default_rng(1).normal(st['R'], st['eR'], len(k)); Rp = k * Rs * 109.076; res['Rp_Rearth'] = q(Rp); res['Rp_RJ'] = q(Rp / 11.209)
# photoeccentric (Dawson & Johnson 2012); not with the density prior, which sets rho_circ
pe = None
if not RHO:
    r2 = np.random.default_rng(2); M = r2.normal(st['M'], st['eM'], len(rc)); R = r2.normal(st['R'], st['eR'], len(rc)); rs = M / R ** 3 * 1.4100
    g = (rc / rs) ** (1 / 3.); kde = gaussian_kde(g)
    E = r2.uniform(0, 0.95, 400000); W = r2.uniform(0, 2 * np.pi, 400000); wt = kde((1 + E * np.sin(W)) / np.sqrt(1 - E ** 2)); wt /= wt.sum()
    e = E[r2.choice(len(E), 200000, p=wt)]
    h, _ = np.histogram(e, bins=np.linspace(0, 0.95, 96), density=True)
    pe = res['photoeccentric'] = dict(rho_star=q(rs), ratio_rho_circ_over_rho_star=q(rc / rs), g=q(g), e=q(e), e_upper_68=float(np.percentile(e, 68.3)),
                                      e_upper_95=float(np.percentile(e, 95)), P_e_lt_0p1=float(np.mean(e < 0.1)), P_e_gt_0p3=float(np.mean(e > 0.3)),
                                      savage_dickey_e0=float(h[:2].mean() * 0.95), g_equals_1_sigma_above_median=float((1 - np.median(g)) / np.std(g)))
np.save(os.path.join(wd, f'free_dilution{SUF}_chain.npy'), ch); json.dump(res, open(os.path.join(wd, f'free_dilution{SUF}.json'), 'w'), indent=1)
for n in ['k', 'depth_true_ppm', 'f_t', 'b', 'aR', 'rho_circ', 'inc_deg', 'T14_h', 'T23_h', 'Rp_Rearth', 'Rp_RJ']:
    v = res[n]; print(f'{n}: {v[2]:.4g} +{v[3] - v[2]:.3g} -{v[2] - v[1]:.3g} (95 %: {v[0]:.4g}-{v[4]:.4g})')
if pe:
    print(f'rho_circ/rho* = {pe["ratio_rho_circ_over_rho_star"][2]:.2f}; e = {pe["e"][2]:.2f} +{pe["e"][3] - pe["e"][2]:.2f} -{pe["e"][2] - pe["e"][1]:.2f}, '
          f'e < {pe["e_upper_68"]:.2f} (68 %), < {pe["e_upper_95"]:.2f} (95 %); P(e < 0.1) = {pe["P_e_lt_0p1"]:.2f}; Savage-Dickey {pe["savage_dickey_e0"]:.2f}')
print(f'best chi2 {res["best_chi2"]:.1f} for {res["npts"]} points; acceptance {res["acceptance"]:.2f}; tau {tau}; {res["runtime_s"]:.0f} s')
os.makedirs(args.out, exist_ok=True); fn = os.path.join(args.out, f'fit_free_dilution{SUF}_{args.tic}.csv')
rows = [(n, res[n], '') for n in ['k', 'depth_true_ppm', 'f_t', 'b', 'aR', 'rho_circ', 'inc_deg', 'T14_h', 'T23_h', 'q1', 'q2']]
rows += [(n, res[n], 'BTJD offset from the fit_transits.py time (d)') for n in NAMES[6:]]
rows += [('Rp_Rearth', res['Rp_Rearth'], f'R* = {st["R"]} +- {st["eR"]} Rsun'), ('Rp_RJ', res['Rp_RJ'], f'R* = {st["R"]} +- {st["eR"]} Rsun')]
if pe:
    rows += [('rho_star', pe['rho_star'], f'M/R^3, M = {st["M"]} +- {st["eM"]}, R = {st["R"]} +- {st["eR"]}'),
             ('rho_circ_over_rho_star', pe['ratio_rho_circ_over_rho_star'], ''), ('g', pe['g'], '(rho_circ/rho*)^(1/3)'),
             ('e', pe['e'], f'photoeccentric; e < {pe["e_upper_68"]:.2f} (68 %), < {pe["e_upper_95"]:.2f} (95 %); P(e < 0.1) = {pe["P_e_lt_0p1"]:.2f}; '
                            f'P(e > 0.3) = {pe["P_e_gt_0p3"]:.2f}; Savage-Dickey ratio at e = 0: {pe["savage_dickey_e0"]:.2f}')]
RHOTXT = f'rho_circ ~ N({RHO[0]}, {RHO[1]}) g cm^-3 (stellar-density prior)' if RHO else 'log10 rho_circ uniform'
with open(fn, 'w', newline='') as fo:
    fo.write(f'# TIC {args.tic} = Gaia DR3 {c["gaia_dr3"]}: joint transit fit with free target flux fraction f_t (scripts/fit_free_dilution.py), {DATE}\n'
             f'# data: TESScut r = {args.aperture} px raw aperture flux, epochs {",".join(USE)}, quality 0, +-{args.window} d, {res["npts"]} points; '
             f'per-epoch noise fixed at the fit_transits.py residual rms\n'
             f'# priors: f_t ~ N({FT_MU}, {FT_SD}); u1 ~ N({U_MU[0]}, {U_SD}), u2 ~ N({U_MU[1]}, {U_SD}); {RHOTXT}; circular orbit; P = {P} d fixed\n'
             f'# emcee: {nw} walkers x {args.steps} steps, first {args.steps // 3} discarded, thinned by 10; acceptance {res["acceptance"]:.2f}; '
             f'best chi2 {res["best_chi2"]:.1f} for {res["npts"]} points\n'
             f'# columns p2.5 ... p97.5 = percentiles of the posterior; k = Rp/R* corrected for dilution; depth_true_ppm = k^2; rho_circ in g cm^-3; '
             f'T14_h, T23_h in hours\n'
             + ('# photoeccentric rows: Dawson & Johnson (2012) with uniform priors e in [0, 0.95), omega in [0, 2 pi)\n' if pe else
                '# no photoeccentric rows: the density prior sets rho_circ (see the fit without the prior)\n'))
    w = csv.writer(fo); w.writerow(['quantity', 'p2.5', 'p16', 'p50', 'p84', 'p97.5', 'note'])
    for n, v, note in rows:
        w.writerow([n] + [f'{x:.6g}' for x in v] + [note])
print('written', fn)
