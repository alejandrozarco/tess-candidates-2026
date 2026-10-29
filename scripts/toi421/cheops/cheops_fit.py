#!/usr/bin/env python
"""TOI-421 (TIC 94986319, Gaia DR3 2984582227215748864): CHEOPS impact-parameter fits for planets b and c.

Conventions (fixed, identical to the TESS series):
  planet b: P = 5.19759807 d, T0 = BTJD 1441.283952, a/R* fixed 13.76, k fixed 0.0279 in the primary variant
  planet c: P = 16.0675076 d, T0 = BTJD 2195.30635, a/R* fixed 31.910, k free (primary) or fixed 0.05459
  BJD_TDB = BTJD + 2457000; circular orbits; batman quadratic LD; model integrated over the 60 s exposure
  (supersample 7). Detrending: multiplicative linear basis (constant, time, time^2, centroid x, x^2, y, y^2,
  background, sin/cos of n*roll for n = 1..3), terms chosen per visit by forward BIC selection; the linear
  coefficients are free and marginalised analytically (flat priors). White-noise jitter per visit free.
  Free: b, tc, k (unless fixed), LD (only in LD-free variants: Kipping q1, q2 ~ U(0,1)), jitter.

Products:
  PIPE : Krenn et al. 2024 PIPE PSF photometry (VizieR J/A+A/686/A301/cheops1, raw, author-cleaned)
  DRP  : CHEOPS DRP 14.1.2 DEFAULT aperture (file keys CH_PR100024_TG0xxxxx_V0300 via DACE), own cleaning,
         contamination column subtracted.

Usage: python scripts/toi421/cheops/cheops_fit.py RUNSPEC_JSON OUTDIR   (run specs are built by run_all.py; inputs from
download_cheops.py in work/toi421/cheops/data)
"""
import sys, json, os, time
import numpy as np
import pandas as pd
import batman
import emcee
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
DATA = os.path.join(ROOT, 'work', 'toi421', 'cheops', 'data')

PLANETS = {
    'b': dict(P=5.19759807, T0=1441.283952, aR=13.76, kfix=0.0279),
    'c': dict(P=16.0675076, T0=2195.30635, aR=31.910, kfix=0.05459),
}
VISITS = {  # visit number (Krenn+2024 Table 1), file key, planet
    1: ('CH_PR100024_TG013601_V0300', 'b'),
    2: ('CH_PR100024_TG014701_V0300', 'b'),
    3: ('CH_PR100024_TG014702_V0300', 'b'),
    4: ('CH_PR100024_TG014703_V0300', 'b'),
    5: ('CH_PR100024_TG015301_V0300', 'c'),
    6: ('CH_PR100024_TG015302_V0300', 'c'),
}
# Theoretical quadratic LD, CHEOPS band, Teff 5291 K, logg 4.48, [Fe/H] -0.02
LD = {
    'atlas':   (0.5153, 0.1864),  # Claret 2021 RNAAS 5,13 (VizieR J/other/RNAAS/5.13 table8), ATLAS, LSM, vturb 2, interpolated
    'phoenix': (0.5700, 0.1438),  # Claret 2021 table2, PHOENIX-COND (solar Z), LSM, interpolated in Teff at logg 4.5
    'ldcu':    (0.5013, 0.1968),  # Krenn+2024 LDCU prior centre q1 0.487, q2 0.359 converted to u1, u2
    'tess':    (0.401, 0.222),    # TESS-band template values of our TESS series (diagnostic only, wrong band)
}
EXPTIME = 60.0 / 86400.
SUPERSAMPLE = 7
GAP = 5.0 / 1440.        # a gap is > 5 min without data
EDGE = 2.0 / 1440.       # drop points within 2 min of a gap edge (SAA / Earth occultation)
UNITS = ['t', 't2', 'x', 'x2', 'y', 'y2', 'bg', 'r1', 'r2', 'r3']


def visit_window():
    idx = pd.read_csv(os.path.join(DATA, 'dace_cheops_full_index.csv'))
    w = {}
    for v, (fk, _) in VISITS.items():
        r = idx[idx.file_key == fk].iloc[0]
        w[v] = (r.date_mjd_start + 2400000.5 - 2457000 - 0.01, r.date_mjd_end + 2400000.5 - 2457000 + 0.01)
    return w


def drop_gap_edges(t):
    keep = np.ones(len(t), bool)
    for i in np.where(np.diff(t) > GAP)[0]:
        keep[(t > t[i] - EDGE) & (t <= t[i])] = False
        keep[(t >= t[i + 1]) & (t < t[i + 1] + EDGE)] = False
    keep[(t < t[0] + EDGE) | (t > t[-1] - EDGE)] = False
    return keep


def read_vizier_tsv(path):
    lines = [l.rstrip('\n') for l in open(path) if not l.startswith('#') and l.strip()]
    hdr = lines[0].split('\t'); rows = []
    for l in lines[3:]:
        try: rows.append([float(z) for z in l.split('\t')])
        except ValueError: pass
    return pd.DataFrame(rows, columns=[h.strip() for h in hdr])


def load(product, visit):
    win = visit_window()[visit]
    if product == 'PIPE':
        d = read_vizier_tsv(os.path.join(DATA, 'VizieR_J_A+A_686_A301_cheops1.tsv'))
        t = d.Time.values - 2457000.
        m = (t > win[0]) & (t < win[1])
        d = d[m]; t = t[m]
        out = dict(t=t, f=d.NFlux.values, e=d.e_NFlux.values, roll=d.Roll.values,
                   x=d.CenX.values, y=d.CenY.values, bg=d.Bg.values)
        nraw = len(t)
    else:
        d = pd.read_csv(os.path.join(DATA, f'{VISITS[visit][0]}_DRP_default.csv'))
        nraw = len(d)
        d = d[d.photom_event_vect == 0]
        f = d.photom_flux_vect.values - d.photom_conta_lc_vect.values
        e = d.photom_flux_vect_err.values
        bg = d.photom_background_vect.values
        x, y = d.photom_centroid_x_vect.values, d.photom_centroid_y_vect.values
        ok = (bg < 4 * np.median(bg)) & (np.hypot(x - np.median(x), y - np.median(y)) < 1.0)
        med = np.median(f[ok])
        out = dict(t=d.obj_date_bjd_vect.values[ok] + 2400000. - 2457000., f=f[ok] / med, e=e[ok] / med,
                   roll=d.photom_roll_angle_vect.values[ok], x=x[ok], y=y[ok], bg=bg[ok])
    o = np.argsort(out['t'])
    out = {k: v[o] for k, v in out.items()}
    keep = drop_gap_edges(out['t'])
    out = {k: v[keep] for k, v in out.items()}
    out['nraw'] = nraw
    return out


def basis(d, units):
    t = d['t']; tm = 0.5 * (t.min() + t.max()); ts = (t - tm) / (0.5 * (t.max() - t.min()))
    dx = d['x'] - np.median(d['x']); dy = d['y'] - np.median(d['y'])
    bg = (d['bg'] - np.median(d['bg'])) / (np.std(d['bg']) + 1e-30)
    phi = np.radians(d['roll'])
    cols = [np.ones_like(t)]
    for u in units:
        if u == 't': cols.append(ts)
        elif u == 't2': cols.append(ts ** 2)
        elif u == 'x': cols.append(dx)
        elif u == 'x2': cols.append(dx ** 2)
        elif u == 'y': cols.append(dy)
        elif u == 'y2': cols.append(dy ** 2)
        elif u == 'bg': cols.append(bg)
        elif u[0] == 'r':
            n = int(u[1]); cols += [np.sin(n * phi), np.cos(n * phi)]
    A = np.vstack(cols).T
    # normalise non-constant columns for conditioning
    s = A.std(0); s[0] = 1.; s[s == 0] = 1.
    return A / s


class Visit:
    def __init__(self, d, planet):
        self.d = d; self.pl = PLANETS[planet]
        P, T0 = self.pl['P'], self.pl['T0']
        tmid = 0.5 * (d['t'].min() + d['t'].max())
        n = np.round((tmid - T0) / P)
        self.epoch = int(n); self.tpred = T0 + n * P
        self.params = batman.TransitParams()
        p = self.params
        p.t0 = self.tpred; p.per = P; p.rp = self.pl['kfix']; p.a = self.pl['aR']; p.inc = 89.
        p.ecc = 0.; p.w = 90.; p.u = [0.5, 0.2]; p.limb_dark = 'quadratic'
        self.m = batman.TransitModel(p, d['t'], supersample_factor=SUPERSAMPLE, exp_time=EXPTIME)
        self.units = []
        self.A = basis(d, [])

    def set_units(self, units):
        self.units = list(units); self.A = basis(self.d, self.units)

    def transit(self, b, tc, k, u):
        p = self.params
        cosi = b / self.pl['aR']
        if cosi >= 1: return None
        p.inc = np.degrees(np.arccos(cosi)); p.t0 = tc; p.rp = k; p.u = list(u)
        return self.m.light_curve(p)

    def lnlike(self, b, tc, k, u, lns, return_model=False):
        T = self.transit(b, tc, k, u)
        if T is None: return -np.inf
        f, e = self.d['f'], self.d['e']
        W = 1. / (e ** 2 + np.exp(2 * lns))
        D = T[:, None] * self.A
        DW = D * W[:, None]
        M = D.T @ DW
        try:
            c = np.linalg.solve(M, DW.T @ f)
        except np.linalg.LinAlgError:
            return -np.inf
        mod = D @ c
        r = f - mod
        sign, logdet = np.linalg.slogdet(M)
        ll = -0.5 * np.sum(W * r ** 2) + 0.5 * np.sum(np.log(W / (2 * np.pi))) - 0.5 * logdet
        if return_model:
            return ll, dict(model=mod, transit=T, sys=self.A @ c, coef=c, resid=r, W=W)
        return ll


def chi2_profile(v, b, tc, k, u):
    T = v.transit(b, tc, k, u)
    if T is None: return 1e30
    D = T[:, None] * v.A; W = 1 / v.d['e'] ** 2
    c = np.linalg.lstsq(D * np.sqrt(W)[:, None], v.d['f'] * np.sqrt(W), rcond=None)[0]
    return np.sum(W * (v.d['f'] - D @ c) ** 2)


def best_transit(v, kfree, u, x0=None):
    """optimise (b, tc[, k]) for the current basis, chi2 with nominal errors."""
    k0 = v.pl['kfix']
    def fun(x):
        b, dtc = x[0], x[1]; k = x[2] if kfree else k0
        if b < 0 or b > 1 + k or abs(dtc) > 0.04 or k <= 0.005 or k > 0.2: return 1e30
        return chi2_profile(v, b, v.tpred + dtc, k, u)
    best = None
    starts = [x0] if x0 is not None else []
    for b0 in (0.3, 0.6, 0.75, 0.85, 0.95):
        for dt0 in (-0.006, 0.0, 0.006):
            starts.append([b0, dt0] + ([k0] if kfree else []))
    for s in starts:
        r = minimize(fun, s, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-4, maxiter=4000))
        if best is None or r.fun < best.fun: best = r
    return best


def select_basis(v, kfree, u):
    """forward BIC selection of detrending units with the transit re-optimised at every trial."""
    n = len(v.d['t']); chosen = []
    v.set_units(chosen)
    r = best_transit(v, kfree, u)
    npar = v.A.shape[1] + len(r.x)
    bic = r.fun + npar * np.log(n)
    log = [dict(units=[], chi2=float(r.fun), bic=float(bic))]
    while True:
        trials = []
        for un in UNITS:
            if un in chosen: continue
            v.set_units(chosen + [un])
            rr = best_transit(v, kfree, u, x0=list(r.x))
            npar = v.A.shape[1] + len(rr.x)
            trials.append((rr.fun + npar * np.log(n), un, rr))
        if not trials: break
        trials.sort(key=lambda z: z[0])
        if trials[0][0] < bic:
            bic, un, r = trials[0]; chosen.append(un)
            log.append(dict(units=list(chosen), chi2=float(r.fun), bic=float(bic)))
        else:
            break
    v.set_units(chosen)
    return chosen, log, r


def q_to_u(q1, q2):
    s = np.sqrt(q1); return 2 * s * q2, s * (1 - 2 * q2)


def run(spec, outdir):
    t0 = time.time()
    planet = spec['planet']; pl = PLANETS[planet]
    visits = spec['visits']; product = spec['product']
    kfree = spec['kfree']; ldmode = spec['ld']
    rng = np.random.default_rng(spec.get('seed', 1))
    bas = json.load(open(os.path.join(outdir, f"basis_{product}.json")))
    V = []
    for vi in visits:
        d = load(product, vi)
        key = str(vi)
        keep = np.array(bas[key]['clip_mask'])
        d2 = {k: (val[keep] if isinstance(val, np.ndarray) else val) for k, val in d.items()}
        v = Visit(d2, planet)
        v.set_units(bas[key]['units'])
        V.append(v)
    # parameter vector: b, [k], [q1, q2], then per visit (dtc, lns)
    names = ['b'] + (['k'] if kfree else []) + (['q1', 'q2'] if ldmode == 'free' else [])
    for vi in visits: names += [f'dtc_{vi}', f'lns_{vi}']
    ng = len(names) - 2 * len(visits)

    def unpack(x):
        i = 0; b = x[i]; i += 1
        k = x[i] if kfree else pl['kfix']; i += kfree
        if ldmode == 'free':
            q1, q2 = x[i], x[i + 1]; i += 2; u = q_to_u(q1, q2)
        else:
            q1 = q2 = None; u = LD[ldmode]
        return b, k, q1, q2, u, x[i:]

    def lnprob(x):
        b, k, q1, q2, u, rest = unpack(x)
        if not (0 <= b <= 1 + k and 0.005 < k < 0.2): return -np.inf
        if q1 is not None and not (0 < q1 < 1 and 0 < q2 < 1): return -np.inf
        ll = 0.
        for j, v in enumerate(V):
            dtc, lns = rest[2 * j], rest[2 * j + 1]
            if abs(dtc) > 0.04 or not (-14 < lns < -4.5): return -np.inf
            ll += v.lnlike(b, v.tpred + dtc, k, u, lns)
            if not np.isfinite(ll): return -np.inf
        return ll

    # starting point: optimise lnprob
    b0 = 0.94 if planet == 'b' else 0.75
    x0 = [b0] + ([pl['kfix']] if kfree else []) + ([0.35, 0.36] if ldmode == 'free' else [])
    for v in V: x0 += [0.0, np.log(1e-4)]
    x0 = np.array(x0)
    best = None
    for bstart in ([0.85, 0.92, 0.96, 1.0] if planet == 'b' else [0.6, 0.7, 0.75, 0.8]):
        xs = x0.copy(); xs[0] = bstart
        for j, v in enumerate(V):
            r = best_transit(v, False, LD['atlas'] if ldmode == 'free' else LD[ldmode],
                             x0=[bstart, 0.0]) if planet == 'b' else None
            if r is not None: xs[ng + 2 * j] = r.x[1]
        if planet == 'c':
            for j, v in enumerate(V):
                r = best_transit(v, True, LD['atlas']); xs[ng + 2 * j] = r.x[1]
        rr = minimize(lambda x: -lnprob(x), xs, method='Nelder-Mead',
                      options=dict(maxiter=20000, maxfev=20000, xatol=1e-7, fatol=1e-3))
        if best is None or rr.fun < best.fun: best = rr
    xb = best.x
    ndim = len(xb); nw = max(32, 4 * ndim)
    scale = np.array([0.01] + ([0.001] if kfree else []) + ([0.05, 0.05] if ldmode == 'free' else [])
                     + [0.0005, 0.1] * len(V))
    p0 = []
    while len(p0) < nw:
        p = xb + scale * rng.standard_normal(ndim)
        if planet == 'b' and not kfree: p[0] = np.clip(p[0], 0.0, 1.0 + pl['kfix'] - 1e-4)
        if np.isfinite(lnprob(p)): p0.append(p)
    nsteps = spec.get('nsteps', 6000); burn = spec.get('burn', 2000)
    sampler = emcee.EnsembleSampler(nw, ndim, lnprob,
                                    moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    sampler.run_mcmc(np.array(p0), nsteps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)
    lp = sampler.get_log_prob(discard=burn, flat=True)
    try:
        tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    except Exception:
        tau = np.full(ndim, np.nan)
    pct = lambda a: [float(np.percentile(a, 50)), float(np.percentile(a, 50) - np.percentile(a, 16)),
                     float(np.percentile(a, 84) - np.percentile(a, 50))]
    res = {}
    for i, nm in enumerate(names): res[nm] = pct(chain[:, i])
    if ldmode == 'free':
        iq = names.index('q1'); u1, u2 = q_to_u(chain[:, iq], chain[:, iq + 1])
        res['u1'] = pct(u1); res['u2'] = pct(u2)
    for j, vi in enumerate(visits):
        res[f'tc_{vi}_BTJD'] = pct(V[j].tpred + chain[:, names.index(f'dtc_{vi}')])
        res[f'dtc_{vi}_min'] = pct(1440 * chain[:, names.index(f'dtc_{vi}')])
        res[f'jitter_{vi}_ppm'] = pct(1e6 * np.exp(chain[:, names.index(f'lns_{vi}')]))
    ib = names.index('b')
    res['b_95pct_interval'] = [float(np.percentile(chain[:, ib], 2.5)), float(np.percentile(chain[:, ib], 97.5))]
    # diagnostics at MAP
    xm = chain[np.argmax(lp)]
    b, k, q1, q2, u, rest = unpack(xm)
    diag = {}
    for j, (vi, v) in enumerate(zip(visits, V)):
        ll, mm = v.lnlike(b, v.tpred + rest[2 * j], k, u, rest[2 * j + 1], return_model=True)
        r = mm['resid']; e = v.d['e']
        chi2n = float(np.sum((r / e) ** 2)); dof = len(r) - v.A.shape[1] - 2
        # time-averaging beta factor (10-30 min bins)
        betas = []
        for nb in (10, 15, 20, 30):
            m = len(r) // nb
            if m < 4: continue
            rb = r[:m * nb].reshape(m, nb).mean(1)
            betas.append(np.std(rb) / (np.std(r) / np.sqrt(nb)) * np.sqrt(m / (m - 1.)))
        intr = (mm['transit'] < 1)
        diag[str(vi)] = dict(file_key=VISITS[vi][0], n_raw=int(v.d.get('nraw', -1)), n_used=int(len(r)),
                             n_in_transit=int(intr.sum()), units=v.units, n_basis=int(v.A.shape[1]),
                             rms_ppm=float(1e6 * np.std(r)), median_err_ppm=float(1e6 * np.median(e)),
                             chi2_red_nominal=chi2n / dof, beta_median=float(np.median(betas)) if betas else None,
                             epoch=v.epoch, tpred_BTJD=float(v.tpred))
        np.savez_compressed(os.path.join(outdir, f"lc_{spec['name']}_v{vi}.npz"),
                            info=np.array(f"TOI-421 planet {planet}, visit {vi} ({VISITS[vi][0]}), product {product}, "
                                          f"run {spec['name']}; columns: t BTJD (BJD_TDB-2457000), flux, err, "
                                          f"systematics model, transit model (MAP)"),
                            t=v.d['t'], f=v.d['f'], e=v.d['e'], sys=mm['sys'], transit=mm['transit'])
    out = dict(
        header=(f"TOI-421 (TIC 94986319) planet {planet}; CHEOPS visits {visits} "
                f"({', '.join(VISITS[x][0] for x in visits)}); product {product}; "
                f"k {'free' if kfree else 'fixed ' + str(pl['kfix'])}; LD {ldmode}"
                f"{'' if ldmode == 'free' else ' u=' + str(LD[ldmode])}; P {pl['P']} d, T0 BTJD {pl['T0']}, "
                f"a/R* fixed {pl['aR']}; circular; batman quadratic LD, 60 s exposure supersample {SUPERSAMPLE}; "
                f"multiplicative linear detrending marginalised, BIC-selected units; jitter free; emcee"),
        spec=spec, names=names, results=res, diagnostics=diag,
        map=dict(zip(names, map(float, xm))), lnp_map=float(lp.max()),
        mcmc=dict(nwalkers=nw, nsteps=nsteps, burn=burn, tau=[float(x) for x in tau],
                  acceptance=float(np.mean(sampler.acceptance_fraction))),
        runtime_s=time.time() - t0)
    json.dump(out, open(os.path.join(outdir, f"fit_{spec['name']}.json"), 'w'), indent=1)
    np.savez_compressed(os.path.join(outdir, f"chain_{spec['name']}.npz"), chain=chain[::5], names=np.array(names),
                        info=np.array(out['header']))
    return out


def select_all(product, outdir):
    bas = dict(header=f"TOI-421 (TIC 94986319) CHEOPS detrending-basis selection, product {product}; forward BIC with "
                      f"transit (b, tc; k free for c, fixed 0.0279 for b; LD Claret 2021 ATLAS) re-optimised per trial; "
                      f"then one 5x MAD residual clip; keys = visit numbers (Krenn+2024 Table 1)")
    for vi, (fk, planet) in VISITS.items():
        pl = PLANETS[planet]
        d = load(product, vi)
        v = Visit(d, planet)
        ch, log, r = select_basis(v, planet == 'c', LD['atlas'])
        b_, dtc_ = r.x[0], r.x[1]; k_ = r.x[2] if planet == 'c' else pl['kfix']
        T = v.transit(b_, v.tpred + dtc_, k_, LD['atlas'])
        D = T[:, None] * v.A
        c = np.linalg.lstsq(D, d['f'], rcond=None)[0]
        res = d['f'] - D @ c
        mad = 1.4826 * np.median(np.abs(res - np.median(res)))
        keep = np.abs(res - np.median(res)) < 5 * mad
        bas[str(vi)] = dict(file_key=fk, planet=planet, n_loaded=int(len(d['t'])), n_raw=int(d['nraw']),
                            units=ch, bic_log=log, best_b=float(b_), best_dtc_min=float(1440 * dtc_),
                            best_k=float(k_), n_clip=int((~keep).sum()), clip_mask=keep.tolist())
        print(product, vi, fk, planet, 'n', len(d['t']), 'units', ch, 'b %.3f k %.4f dtc %.1f min clip %d'
              % (b_, k_, 1440 * dtc_, (~keep).sum()), flush=True)
    json.dump(bas, open(os.path.join(outdir, f"basis_{product}.json"), 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1] == 'select':
        select_all(sys.argv[2], sys.argv[3]); sys.exit()
    spec = json.load(open(sys.argv[1])) if sys.argv[1].endswith('.json') else json.loads(sys.argv[1])
    o = run(spec, sys.argv[2])
    print(spec['name'], 'b =', o['results']['b'], 'k' in o['results'] and o['results'].get('k'), f"{o['runtime_s']:.0f}s")
