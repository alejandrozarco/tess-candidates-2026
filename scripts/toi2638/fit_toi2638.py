"""TOI-2638.01 (TIC 458342769): ephemeris, per-sector transit depths in every TESS product, per-transit depths,
per-sector trapezoid durations, joint models of the change, trends and their red-noise tests (README section 6.9).

usage: python scripts/toi2638/fit_toi2638.py [--work ./work] [--out data] [--nboot 300] [--nnull 1000] [--nsecboot 2000]
input:  the light curves of download_toi2638.py in <work>/458342769/dl/lc/
output: <out>/ephemeris_458342769.csv, depths_458342769.csv, pertransit_458342769.csv, shape_458342769.csv,
        jointfit_458342769.csv, trends_458342769.csv, predictions_458342769.csv; <work>/458342769/fit.json
        (formats in the file headers)
Steps (windows and model as in common.py):
 1. ephemeris: common shape fitted to the SPOC PDCSAP windows plus QLP sector 10 at the TOI-table ephemeris; one
    mid-time offset per sector (shape fixed, amplitude free); weighted linear fit of the sector mid-times, errors
    scaled by sqrt(chi2/dof) when it exceeds 1; the windows are then rebuilt on this ephemeris.
 2. common shape (k, b, a/R*) on SPOC PDCSAP + QLP sector 10 ('null' model); per-sector depth = amplitude of this
    template times its instantaneous peak depth, by weighted least squares, in every product; error = max(white
    error x beta, bootstrap over transits); beta from the residuals binned to 15, 30 and 60 min.
 3. per-transit depths (SPOC PDCSAP) with the quality flags and background of each window (SPOC files).
 4. per-sector trapezoid (depth, T14, ingress fraction, centre) on SPOC PDCSAP, exposure-averaged, errors from a
    bootstrap over transits (150 draws).
 5. joint models over all sectors of a product (least squares, several starting points): null (one shape), M0 (one
    shape, free depth scale per sector, sector 64 fixed at 1), M1 (common k and a/R*, one b per sector), M2 (common k
    and a/R*, b linear in time), M5 (one shape, depth scale linear in time); white errors per window.
 6. trends: weighted linear fits of the per-sector depths and trapezoid durations against time; sector-block
    bootstrap; leave-one-sector-out; red-noise null test (constant-depth model plus the observed residuals of each
    sector shifted cyclically by a random amount, per-sector depths refitted, slope recorded).
 7. predictions from M2 and from the linear depth and duration trends.
Requires numpy, scipy, astropy, batman-package.
"""
import argparse, json, os, sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (TIC, P_TOI, T0_TOI, REF_YEAR, DATE, SPOC_SECTORS, load, crowd_correct, windows, Model, peak_depth,
                    model_t14, trapezoid, wlinfit, hdr, year, dldir)

ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--work', default='./work'); ap.add_argument('--out', default='data')
ap.add_argument('--nboot', type=int, default=300); ap.add_argument('--nnull', type=int, default=1000); ap.add_argument('--nsecboot', type=int, default=2000)
args = ap.parse_args(); os.makedirs(args.out, exist_ok=True)
rng = np.random.default_rng(2638)
OUT = {}

# ---------------------------------------------------------------- data
pdc = load(args.work, 'SPOC'); sapr = load(args.work, 'SPOC', 'SAP_FLUX'); sapc = crowd_correct(load(args.work, 'SPOC', 'SAP_FLUX'))
tsp = load(args.work, 'TESSSPOC'); qlp = load(args.work, 'QLP'); tglc = load(args.work, 'TGLC')
print('sectors: SPOC', sorted(pdc), 'TESS-SPOC', sorted(tsp), 'QLP', sorted(qlp), 'TGLC', sorted(tglc))


def with_s10(d):
    x = dict(d); x[10] = qlp[10]; return x


# ---------------------------------------------------------------- generic fitting helpers
class Set:
    """windows of one product with one Model per sector"""
    def __init__(self, d, P, T0):
        self.W = windows(d, P, T0); self.secs = sorted(self.W); self.P, self.T0 = P, T0
        self.M = {s: Model(self.W[s]['t'], self.W[s]['cad'], P, T0) for s in self.secs}
        self.N = sum(len(self.W[s]['t']) for s in self.secs)

    def resid(self, fn):
        return np.concatenate([(self.W[s]['f'] - (1 - fn(s))) / self.W[s]['e'] for s in self.secs])


def fit(fun, starts, lo, hi):
    best = None
    for x0 in starts:
        r = least_squares(fun, np.clip(x0, np.array(lo) + 1e-9, np.array(hi) - 1e-9), bounds=(lo, hi), x_scale='jac', max_nfev=5000)
        if best is None or r.cost < best.cost:
            best = r
    try:
        err = np.sqrt(np.diag(np.linalg.inv(best.jac.T @ best.jac)))
    except np.linalg.LinAlgError:
        err = np.full(len(best.x), np.nan)
    return best.x, err, 2 * best.cost


def null_fit(S):
    f = lambda p: S.resid(lambda s: S.M[s].depth(*p))
    return fit(f, [[0.055, 0.92, 11], [0.07, 1.0, 9], [0.05, 0.85, 13]], [0.01, 0, 2], [0.5, 1.5, 40])


def amp(y, x, e, eps=None):
    """weighted least-squares amplitude of template x in y = 1 - f"""
    w = 1 / e ** 2; A = np.sum(w * x * y) / np.sum(w * x * x); return A, 1 / np.sqrt(np.sum(w * x * x))


def beta_factor(r, t, cad, ep):
    """median over 15, 30, 60 min of rms(binned residuals) / white expectation, within windows; at least 1"""
    out = []
    for bmin in (15, 30, 60):
        nb = max(2, int(round(bmin / (cad * 1440))))
        rb = []
        for n in np.unique(ep):
            rr = r[ep == n]; M = len(rr) // nb
            if M >= 2:
                rb.append(rr[:M * nb].reshape(M, nb).mean(1))
        if rb:
            rb = np.concatenate(rb); out.append(np.std(rb) / (np.std(r) / np.sqrt(nb)))
    return max(1.0, float(np.median(out))) if out else 1.0


# ---------------------------------------------------------------- 1. ephemeris
S0 = Set(with_s10(pdc), P_TOI, T0_TOI)
x0, _, _ = null_fit(S0)
rows_t = []
for s in S0.secs:
    W, M = S0.W[s], S0.M[s]
    f = lambda p: (W['f'] - (1 - p[0] * M.depth(*x0, dt0=p[1]))) / W['e']
    r = least_squares(f, [1.0, 0.0], bounds=([0, -0.03], [5, 0.03]), x_scale=[0.1, 0.001])
    ef = np.sqrt(np.linalg.inv(r.jac.T @ r.jac)[1, 1])
    eps = np.unique(W['ep']); bt = []
    for _ in range(100):
        pick = rng.choice(eps, len(eps)); idx = np.concatenate([np.where(W['ep'] == e)[0] for e in pick])
        g = lambda p: (W['f'][idx] - (1 - p[0] * M.depth(*x0, dt0=p[1])[idx])) / W['e'][idx]
        bt.append(least_squares(g, r.x, bounds=([0, -0.03], [5, 0.03]), x_scale=[0.1, 0.001]).x[1])
    nmid = int(np.round(np.median(W['ep'])))
    rows_t.append(dict(sector=s, product='QLP' if s == 10 else 'SPOC PDCSAP', epoch_toi=nmid, tc=T0_TOI + nmid * P_TOI + r.x[1],
                       err=max(ef, np.std(bt)), formal=ef, boot=float(np.std(bt)), ntr=W['ntr']))
nref = int(np.round(np.average([r_['epoch_toi'] for r_ in rows_t], weights=[1 / r_['err'] ** 2 for r_ in rows_t])))
E = np.array([r_['epoch_toi'] - nref for r_ in rows_t]); T = np.array([r_['tc'] for r_ in rows_t]); eT = np.array([r_['err'] for r_ in rows_t])
(c0, c1), (e0, e1), chi = wlinfit(E, T, eT); dof = len(E) - 2; sc = max(1.0, np.sqrt(chi / dof))
T0, P = float(c0), float(c1); eT0, eP = float(e0 * sc), float(e1 * sc)
for r_, e_ in zip(rows_t, E):
    r_['epoch'] = int(e_); r_['oc_min'] = (r_['tc'] - T0 - e_ * P) * 1440
print(f'ephemeris: T0 = {T0:.5f} +- {eT0:.5f}, P = {P:.7f} +- {eP:.7f} d, chi2 = {chi:.2f} / {dof}; O-C (min)', [round(r_['oc_min'], 1) for r_ in rows_t])
OUT['ephemeris'] = dict(T0=T0, eT0=eT0, P=P, eP=eP, chi2=chi, dof=dof, nref=nref)
with open(os.path.join(args.out, 'ephemeris_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: per-sector mid-transit times and linear ephemeris; scripts/toi2638/fit_toi2638.py, {DATE}',
            'one mid-time per sector: common transit shape fixed (null model of jointfit_458342769.csv at the TOI-table ephemeris), amplitude and time offset free, '
            'on the detrended +-0.30 d windows; SPOC PDCSAP (2 min) and QLP KSPSAP (sector 10, 1800 s); tc_btjd = BJD_TDB - 2457000 at the window-median epoch',
            'tc_err_min = max(formal, bootstrap over transits (100 draws)); o_c_min relative to the linear ephemeris of the rows with row_type = ephemeris',
            f'linear ephemeris: weighted fit of the sector mid-times, errors scaled by sqrt(chi2/dof) = {sc:.2f} when above 1; chi2 = {chi:.2f} for {dof} degrees of freedom; '
            f'epoch 0 = TOI-table epoch {nref}; TOI table (ExoFOP, TESS project, 2025-07-09): P = {P_TOI:.7f} d, epoch BTJD {T0_TOI:.6f}',
            'row_type,sector,product,epoch,n_transits,tc_btjd,tc_err_min,formal_err_min,bootstrap_err_min,o_c_min,T0_btjd,T0_err_d,P_d,P_err_d'])
    for r_ in rows_t:
        f.write(f"sector,{r_['sector']},{r_['product']},{r_['epoch']},{r_['ntr']},{r_['tc']:.5f},{r_['err'] * 1440:.2f},{r_['formal'] * 1440:.2f},{r_['boot'] * 1440:.2f},{r_['oc_min']:+.2f},,,,\n")
    f.write(f'ephemeris,,,0,,,,,,,{T0:.5f},{eT0:.5f},{P:.7f},{eP:.7f}\n')

# ---------------------------------------------------------------- 2. common shape and per-sector depths
VAR = {'SPOC PDCSAP': pdc, 'SPOC SAP crowding-corrected': sapc, 'SPOC SAP uncorrected': sapr, 'TESS-SPOC PDCSAP': tsp, 'QLP KSPSAP': qlp,
       'TGLC cal_psf_flux': tglc}
SETS = {k: Set(v, P, T0) for k, v in VAR.items()}
SNULL = Set(with_s10(pdc), P, T0)
xn, en, chin = null_fit(SNULL); k0, b0, a0 = xn; D0 = peak_depth(k0, b0, a0, P)
print(f'common shape (PDCSAP + QLP S10): k = {k0:.4f}, b = {b0:.3f}, a/R* = {a0:.2f}, peak depth {D0 * 1e6:.0f} ppm, chi2 {chin:.1f}')
OUT['null_shape'] = dict(k=k0, b=b0, aR=a0, peak_depth=D0)
DEP = []
resid_store = {}
for name, S in SETS.items():
    for s in S.secs:
        W = S.W[s]; x = S.M[s].depth(k0, b0, a0) / D0; y = 1 - W['f']; A, eA = amp(y, x, W['e'])
        r = y - A * x; be = beta_factor(r, W['t'], W['cad'], W['ep'])
        eps = np.unique(W['ep']); bt = []
        for _ in range(args.nboot):
            pick = rng.choice(eps, len(eps)); idx = np.concatenate([np.where(W['ep'] == e)[0] for e in pick]); bt.append(amp(y[idx], x[idx], W['e'][idx])[0])
        eb = float(np.std(bt)) if len(eps) > 1 else np.nan
        DEP.append(dict(product=name, sector=s, year=W['year'], cad=W['cad'] * 86400, ntr=W['ntr'], depth=A * 1e6, ew=eA * 1e6, beta=be,
                        eboot=eb * 1e6, err=max(eA * be, eb if np.isfinite(eb) else 0) * 1e6))
        resid_store[(name, s)] = (x, r, A)
        print(f'{name:28s} S{s:3d} {W["year"]:.2f} ntr {W["ntr"]} depth {A * 1e6:6.0f} +- {DEP[-1]["err"]:4.0f} ppm (white {eA * 1e6:.0f}, beta {be:.2f}, boot {eb * 1e6:.0f})')
with open(os.path.join(args.out, 'depths_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: transit depth per sector and light-curve product; scripts/toi2638/fit_toi2638.py, {DATE}',
            f'depth_ppm = amplitude of the common transit shape (k = {k0:.4f}, b = {b0:.3f}, a/R* = {a0:.2f}; null model of jointfit_458342769.csv, SPOC PDCSAP + QLP sector 10) '
            'times its instantaneous peak depth, exposure-averaged model, weighted least squares on the detrended +-0.30 d windows',
            'err_ppm = max(white_err_ppm x beta, bootstrap_err_ppm); beta = red-noise factor of the residuals binned to 15, 30 and 60 min; bootstrap over transits '
            f'({args.nboot} draws); year = median time of the windows',
            'products: SPOC 2-min PDCSAP_FLUX; SPOC SAP_FLUX divided by its median and corrected with the header CROWDSAP, (f - (1 - CROWDSAP)) / CROWDSAP; '
            'SPOC SAP_FLUX without that correction (depth diluted by the other stars in the SPOC aperture); TESS-SPOC FFI PDCSAP_FLUX; QLP KSPSAP_FLUX (no deblending '
            'applied here); TGLC cal_psf_flux (sector 10); QUALITY = 0',
            'product,sector,year,cadence_s,n_transits,depth_ppm,err_ppm,white_err_ppm,beta,bootstrap_err_ppm'])
    for d in DEP:
        f.write(f"{d['product']},{d['sector']},{d['year']:.3f},{d['cad']:.0f},{d['ntr']},{d['depth']:.0f},{d['err']:.0f},{d['ew']:.0f},{d['beta']:.2f},"
                f"{'' if not np.isfinite(d['eboot']) else round(d['eboot'])}\n")

# ---------------------------------------------------------------- 3. per-transit depths, SPOC PDCSAP, with flags and background
from astropy.io import fits
import glob
PT = []
S = SETS['SPOC PDCSAP']
for s in S.secs:
    fn = glob.glob(os.path.join(dldir(args.work), 'lc', f'tess*-s{s:04d}-0000000458342769-*-s_lc.fits'))[0]
    d = fits.open(fn)[1].data; traw = np.array(d['TIME'], float); q = np.array(d['QUALITY'], int); bk = np.array(d['SAP_BKG'], float)
    W = S.W[s]; x = S.M[s].depth(k0, b0, a0) / D0; y = 1 - W['f']
    for n in np.unique(W['ep']):
        m = W['ep'] == n; A, eA = amp(y[m], x[m], W['e'][m]); tc = T0 + n * P
        win = np.isfinite(traw) & (np.abs(traw - tc) < 0.035); fl = [int(v) for v in sorted(set(q[win])) if v != 0]
        g = np.isfinite(traw) & (q == 0) & np.isfinite(bk); dbk = np.nanmedian(bk[g & (np.abs(traw - tc) < 0.03)]) - np.nanmedian(bk[g & (np.abs(traw - tc) > 0.05) & (np.abs(traw - tc) < 0.3)])
        PT.append(dict(sector=s, epoch=int(n - 0), tc=tc, depth=A * 1e6, err=eA * 1e6, nin=int(win.sum()), nflag=int((q[win] != 0).sum()), flags=fl, dbk=float(dbk)))
chi_sec = {}
for s in S.secs:
    a = np.array([(r_['depth'], r_['err']) for r_ in PT if r_['sector'] == s]); w = 1 / a[:, 1] ** 2; mu = np.sum(a[:, 0] * w) / w.sum()
    chi_sec[s] = (float(np.sum(((a[:, 0] - mu) / a[:, 1]) ** 2)), len(a) - 1, mu)
    print(f'per-transit S{s}: weighted mean {mu:.0f} ppm, chi2 {chi_sec[s][0]:.1f} / {len(a) - 1}; range {a[:, 0].min():.0f}-{a[:, 0].max():.0f}')
with open(os.path.join(args.out, 'pertransit_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: depth of each transit in the SPOC 2-min PDCSAP light curves; scripts/toi2638/fit_toi2638.py, {DATE}',
            'depth_ppm = amplitude of the common shape (as depths_458342769.csv) in one +-0.30 d window, white error from the out-of-transit scatter of that window',
            'n_cadences_in = cadences (any quality) within +-0.035 d of the predicted mid-time; n_flagged_in = of these, QUALITY != 0 (excluded from the fit); flags = the non-zero QUALITY values',
            'bkg_in_minus_out = median SAP_BKG (e-/s, whole aperture) within +-0.03 d minus 0.05-0.30 d from mid-transit, QUALITY = 0; epoch relative to the ephemeris of ephemeris_458342769.csv',
            'chi2_sector, dof_sector = chi2 of the per-transit depths about their weighted mean in that sector',
            'sector,epoch,tc_btjd,depth_ppm,err_ppm,n_cadences_in,n_flagged_in,flags,bkg_in_minus_out,chi2_sector,dof_sector'])
    for r_ in PT:
        cs = chi_sec[r_['sector']]
        f.write(f"{r_['sector']},{r_['epoch']},{r_['tc']:.5f},{r_['depth']:.0f},{r_['err']:.0f},{r_['nin']},{r_['nflag']},{' '.join(map(str, r_['flags']))},{r_['dbk']:.2f},{cs[0]:.1f},{cs[1]}\n")

# ---------------------------------------------------------------- 4. per-sector trapezoid (dilution-free duration)
TR = {}
for s in S.secs:
    W = S.W[s]; ph = (W['t'] - T0 + 0.5 * P) % P - 0.5 * P; y = 1 - W['f']; e = W['e']

    def trfit(idx):
        best = None
        for x_ in ([0.003, 0.058, 0.9, 0], [0.003, 0.075, 0.6, 0]):
            r = least_squares(lambda p: (y[idx] - trapezoid(p, ph[idx], W['cad'])) / e[idx], x_, bounds=([0, 0.02, 0.05, -0.01], [0.02, 0.12, 1.0, 0.01]))
            if best is None or r.cost < best.cost:
                best = r
        return best.x
    allidx = np.arange(len(y)); p0 = trfit(allidx); eps = np.unique(W['ep']); bt = []
    for _ in range(150):
        pick = rng.choice(eps, len(eps)); bt.append(trfit(np.concatenate([np.where(W['ep'] == e_)[0] for e_ in pick])))
    bt = np.array(bt)
    TR[s] = dict(year=W['year'], ntr=W['ntr'], depth=p0[0] * 1e6, edepth=np.std(bt[:, 0]) * 1e6, T14=p0[1] * 1440, eT14=np.std(bt[:, 1]) * 1440,
                 f12=p0[2], ef12=np.std(bt[:, 2]), dt0=p0[3] * 1440)
    print(f'trapezoid S{s}: depth {TR[s]["depth"]:.0f} +- {TR[s]["edepth"]:.0f}, T14 {TR[s]["T14"]:.1f} +- {TR[s]["eT14"]:.1f} min, f12 {p0[2]:.2f} +- {TR[s]["ef12"]:.2f}')

# ---------------------------------------------------------------- 5. joint models
def joint(S, ref=64):
    secs = S.secs; ns = len(secs); yr = {s: S.W[s]['year'] - REF_YEAR for s in secs}; ref = ref if ref in secs else secs[len(secs) // 2]
    others = [s for s in secs if s != ref]; res = {}
    lc = lambda s, k, b, aR, D=1.0: D * S.M[s].depth(k, b, aR)
    xn_, en_, cn_ = fit(lambda p: S.resid(lambda s: lc(s, *p)), [[0.055, 0.92, 11], [0.07, 1.0, 9], [0.05, 0.85, 13]], [0.01, 0, 2], [0.5, 1.5, 40])
    res['null'] = (xn_, en_, cn_, ['k', 'b', 'aR']); k_, b_, a_ = xn_

    def m0(p):
        D = dict(zip(others, p[3:])); D[ref] = 1.0; return S.resid(lambda s: lc(s, *p[:3], D[s]))
    res['M0'] = fit(m0, [[k_, b_, a_] + [1] * (ns - 1), [0.07, 1.0, 8] + [1] * (ns - 1)], [0.01, 0, 2] + [0.05] * (ns - 1), [0.5, 1.5, 40] + [5] * (ns - 1)) + (['k', 'b', 'aR'] + [f'D_s{s}' for s in others],)

    def m1(p):
        bs = dict(zip(secs, p[2:])); return S.resid(lambda s: lc(s, p[0], bs[s], p[1]))
    res['M1'] = fit(m1, [[k_, a_] + [b_] * ns, [0.06, 12] + [0.9] * ns, [0.08, 8] + [1.0] * ns], [0.01, 2] + [0] * ns, [0.5, 40] + [1.5] * ns) + (['k', 'aR'] + [f'b_s{s}' for s in secs],)

    def m2(p):
        return S.resid(lambda s: lc(s, p[0], p[2] + p[3] * yr[s], p[1]))
    res['M2'] = fit(m2, [[k_, a_, b_, -0.02], [0.06, 12, 0.9, -0.02], [0.08, 8, 1.0, -0.03]], [0.01, 2, 0, -0.5], [0.5, 40, 1.5, 0.5]) + (['k', 'aR', 'b_ref', 'db_dt'],)

    def m5(p):
        return S.resid(lambda s: lc(s, *p[:3], 1 + p[3] * yr[s]))
    res['M5'] = fit(m5, [[k_, b_, a_, 0.1], [0.07, 1.0, 9, 0.1], [0.05, 0.85, 13, 0.1]], [0.01, 0, 2, -0.3], [0.5, 1.5, 40, 0.3]) + (['k', 'b', 'aR', 'dD_dt'],)
    return res


def prewhiten(d):
    """remove the 7.31 c/d signal of the 17-arcsec neighbour (first and second harmonic, fitted out of transit per sector)"""
    out = {}
    for s, x in d.items():
        t, fl = x['t'], x['f'] - 1; m = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) > 0.06; best = None
        for fr in np.linspace(7.301, 7.321, 201):
            X = np.c_[np.ones(len(t)), np.sin(2 * np.pi * fr * t), np.cos(2 * np.pi * fr * t), np.sin(4 * np.pi * fr * t), np.cos(4 * np.pi * fr * t)]
            c, *_ = np.linalg.lstsq(X[m], fl[m], rcond=None); a = np.hypot(c[1], c[2])
            if best is None or a > best[0]:
                best = (a, fr, c, X)
        y = dict(x); y['f'] = x['f'] - best[3][:, 1:] @ best[2][1:]; y['pw'] = (best[1], best[0]); out[s] = y
    return out


JV = {'SPOC PDCSAP + QLP sector 10': Set(with_s10(pdc), P, T0), 'SPOC PDCSAP': SETS['SPOC PDCSAP'], 'SPOC SAP crowding-corrected': SETS['SPOC SAP crowding-corrected'],
      'TESS-SPOC PDCSAP': SETS['TESS-SPOC PDCSAP'], 'QLP KSPSAP': SETS['QLP KSPSAP'],
      'SPOC PDCSAP prewhitened + QLP sector 10': Set(with_s10(prewhiten(pdc)), P, T0)}
JR = {}
for name, S_ in JV.items():
    JR[name] = joint(S_)
    print(name, {m: round(v[2], 1) for m, v in JR[name].items()}, 'M2', np.round(JR[name]['M2'][0], 4), 'M5', np.round(JR[name]['M5'][0], 4))
with open(os.path.join(args.out, 'jointfit_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: joint transit models over all sectors of a product; scripts/toi2638/fit_toi2638.py, {DATE}',
            'models: null = one shape (k, b, aR); M0 = one shape, free depth scale D_s per sector (sector 64 fixed at 1); M1 = common k and aR, one b per sector; '
            f'M2 = common k and aR, b = b_ref + db_dt (year - {REF_YEAR:.0f}); M5 = one shape, depth scaled by 1 + dD_dt (year - {REF_YEAR:.0f})',
            'batman, circular orbit, quadratic limb darkening (0.17, 0.27) fixed, exposure-averaged (midpoint sampling); detrended +-0.30 d windows; white errors (out-of-transit scatter per window); '
            'chi2 = sum of squared normalised residuals; npts = number of points; err = formal (inverse of J^T J), not scaled',
            f'ephemeris of ephemeris_458342769.csv; QLP KSPSAP is not deblended here; prewhitened = the 7.31 c/d signal of the 17" neighbour removed from each SPOC sector (two harmonics, fitted out of transit)',
            'variant,model,npts,npar,chi2,parameter,value,err'])
    for name, R_ in JR.items():
        N = JV[name].N
        for mname, (x, e, c, names) in R_.items():
            for pn, v, ev in zip(names, x, e):
                f.write(f'{name},{mname},{N},{len(x)},{c:.2f},{pn},{v:.5f},{ev:.5f}\n')

# ---------------------------------------------------------------- 6. trends
TRD = []


def slope_rows(name, rows, key='depth', ekey='err', unit='ppm/yr'):
    yrs = np.array([r_['year'] for r_ in rows]) - REF_YEAR; v = np.array([r_[key] for r_ in rows]); e = np.array([r_[ekey] for r_ in rows])
    (a, b), (ea, eb), chi = wlinfit(yrs, v, e); dof = len(v) - 2; sc = max(1.0, np.sqrt(chi / dof)) if dof > 0 else 1.0
    return dict(name=name, a=a, b=b, ea=ea, eb=eb, chi=chi, dof=dof, sc=sc, yrs=yrs, v=v, e=e)


def depth_rows(prod, secs=None):
    return [d for d in DEP if d['product'] == prod and (secs is None or d['sector'] in secs)]


main_rows = depth_rows('SPOC PDCSAP') + [d for d in DEP if d['product'] == 'QLP KSPSAP' and d['sector'] == 10]
for label, rows in (('SPOC PDCSAP + QLP sector 10', main_rows), ('SPOC PDCSAP', depth_rows('SPOC PDCSAP')),
                    ('SPOC SAP crowding-corrected', depth_rows('SPOC SAP crowding-corrected')), ('SPOC SAP uncorrected', depth_rows('SPOC SAP uncorrected')),
                    ('TESS-SPOC PDCSAP', depth_rows('TESS-SPOC PDCSAP')), ('QLP KSPSAP', depth_rows('QLP KSPSAP'))):
    L = slope_rows(label, rows)
    TRD.append(dict(quantity='depth slope', variant=label, value=L['b'], err=L['eb'] * L['sc'], z=L['b'] / (L['eb'] * L['sc']), unit='ppm/yr',
                    note=f"weighted linear fit of {len(rows)} sectors; chi2 = {L['chi']:.1f} for {L['dof']} dof; error scaled by {L['sc']:.2f}; depth at {REF_YEAR:.0f}: {L['a']:.0f} +- {L['ea'] * L['sc']:.0f} ppm; "
                         f"relative slope {L['b'] / L['a']:+.3f} per yr"))
    print(TRD[-1])
L = slope_rows('main', main_rows); OUT['depth_trend'] = dict(a=L['a'], b=L['b'], ea=L['ea'] * L['sc'], eb=L['eb'] * L['sc'])
# sector-block bootstrap (resample sectors with replacement; at least 3 distinct sectors)
sl = []
n_ = len(main_rows)
while len(sl) < args.nsecboot:
    pick = rng.integers(0, n_, n_)
    if len(set(pick)) < 3:
        continue
    sl.append(wlinfit(L['yrs'][pick], L['v'][pick], L['e'][pick])[0][1])
TRD.append(dict(quantity='depth slope (sector-block bootstrap)', variant='SPOC PDCSAP + QLP sector 10', value=L['b'], err=float(np.std(sl)), z=L['b'] / np.std(sl), unit='ppm/yr',
                note=f'{args.nsecboot} resamples of the {n_} sectors with replacement (at least 3 distinct); err = standard deviation of the resampled slopes; fraction of resampled slopes <= 0: {np.mean(np.array(sl) <= 0):.4f}'))
loo = [wlinfit(np.delete(L['yrs'], i), np.delete(L['v'], i), np.delete(L['e'], i))[0][1] for i in range(n_)]
TRD.append(dict(quantity='depth slope (leave one sector out)', variant='SPOC PDCSAP + QLP sector 10', value=float(np.min(loo)), err=float('nan'), z=float('nan'), unit='ppm/yr',
                note=f'minimum and maximum over the {n_} omitted sectors: {np.min(loo):.0f} to {np.max(loo):.0f} ppm/yr (omitted sector of the extremes: {main_rows[int(np.argmin(loo))]["sector"]}, {main_rows[int(np.argmax(loo))]["sector"]})'))
# red-noise null: constant depth (weighted mean amplitude) + cyclically shifted residuals per sector
Amean = np.average(L['v'], weights=1 / L['e'] ** 2) / 1e6
keys = [('SPOC PDCSAP', d['sector']) if d['product'] == 'SPOC PDCSAP' else ('QLP KSPSAP', 10) for d in main_rows]
nulls = []
for _ in range(args.nnull):
    vv = []
    for key in keys:
        x, r, A = resid_store[key]; W = SETS[key[0]].W[key[1]]; y = Amean * x + np.roll(r, rng.integers(len(r)))
        vv.append(amp(y, x, W['e'])[0] * 1e6)
    nulls.append(wlinfit(L['yrs'], np.array(vv), L['e'])[0][1])
nulls = np.array(nulls)
TRD.append(dict(quantity='depth slope (red-noise null test)', variant='SPOC PDCSAP + QLP sector 10', value=L['b'], err=float(np.std(nulls)), z=L['b'] / np.std(nulls), unit='ppm/yr',
                note=f'{args.nnull} simulations: constant depth ({Amean * 1e6:.0f} ppm) plus the observed residuals of each sector shifted cyclically by a random amount; per-sector depths refitted; '
                     f'err = standard deviation of the simulated slopes (mean {nulls.mean():.1f}); fraction >= observed: {np.mean(nulls >= L["b"]):.4f}'))
for name in JR:
    c_null, c_m5, c_m2 = JR[name]['null'][2], JR[name]['M5'][2], JR[name]['M2'][2]
    TRD.append(dict(quantity='delta chi2 null - M5 (depth linear in time)', variant=name, value=c_null - c_m5, err=float('nan'), z=float('nan'), unit='',
                    note=f"white errors; 1 extra parameter; M5 dD_dt = {JR[name]['M5'][0][3]:+.3f} +- {JR[name]['M5'][1][3]:.3f} per yr (relative to the {REF_YEAR:.0f} depth)"))
    TRD.append(dict(quantity='delta chi2 M5 - M2 (shape test: depth scaling vs b drift)', variant=name, value=c_m5 - c_m2, err=float('nan'), z=float('nan'), unit='',
                    note=f"same number of parameters; M2 db_dt = {JR[name]['M2'][0][3]:+.4f} +- {JR[name]['M2'][1][3]:.4f} per yr, k = {JR[name]['M2'][0][0]:.4f}, aR = {JR[name]['M2'][0][1]:.2f}, "
                         f"b_ref = {JR[name]['M2'][0][2]:.3f}; calibration in simulations_458342769.csv"))
Ltr = slope_rows('trap', [dict(year=TR[s]['year'], T14=TR[s]['T14'], eT14=TR[s]['eT14']) for s in sorted(TR)], key='T14', ekey='eT14')
OUT['t14_trend'] = dict(a=Ltr['a'], b=Ltr['b'], ea=Ltr['ea'], eb=Ltr['eb'], chi=Ltr['chi'], dof=Ltr['dof'])
TRD.append(dict(quantity='trapezoid T14 slope', variant='SPOC PDCSAP', value=Ltr['b'], err=Ltr['eb'], z=Ltr['b'] / Ltr['eb'], unit='min/yr',
                note=f"weighted fit of {len(TR)} sectors with the bootstrap errors (not scaled); chi2 = {Ltr['chi']:.1f} for {Ltr['dof']} dof; T14 at {REF_YEAR:.0f}: {Ltr['a']:.1f} min; "
                     'the calibrated error is in simulations_458342769.csv'))
with open(os.path.join(args.out, 'trends_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: trends of the transit depth and duration, and model comparisons; scripts/toi2638/fit_toi2638.py, {DATE}',
            'per-sector depths of depths_458342769.csv; trapezoid T14 of shape_458342769.csv; delta chi2 from jointfit_458342769.csv; z = value / err',
            'quantity,variant,value,err,z,unit,note'])
    for r_ in TRD:
        ev = '' if not np.isfinite(r_['err']) else format(r_['err'], '.4g'); zv = '' if not np.isfinite(r_['z']) else format(r_['z'], '.2f')
        f.write(f"{r_['quantity']},{r_['variant']},{r_['value']:.4g},{ev},{zv},{r_['unit']},\"{r_['note']}\"\n")

# shape table: trapezoid + M1 b per sector + M2 prediction
M1x = JR['SPOC PDCSAP + QLP sector 10']['M1']; M1b = dict(zip(JV['SPOC PDCSAP + QLP sector 10'].secs, M1x[0][2:])); M1eb = dict(zip(JV['SPOC PDCSAP + QLP sector 10'].secs, M1x[1][2:]))
k2, a2, bref, bd = JR['SPOC PDCSAP + QLP sector 10']['M2'][0]
with open(os.path.join(args.out, 'shape_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: transit shape per sector; scripts/toi2638/fit_toi2638.py, {DATE}',
            'trapezoid_*: trapezoid (depth, T14, ingress fraction f12 = 2 T12 / T14, centre) fitted to the SPOC PDCSAP windows of one sector, exposure-averaged; errors = standard deviation '
            'of 150 bootstrap resamples of the transits (not calibrated; see simulations_458342769.csv); sector 10 (QLP, 1800 s) is not fitted',
            f'M1_b: impact parameter per sector of model M1 (common k = {M1x[0][0]:.4f}, a/R* = {M1x[0][1]:.2f}; SPOC PDCSAP + QLP sector 10) with formal errors; '
            f'M2_b, M2_depth_ppm, M2_T14_min: model M2 (k = {k2:.4f}, a/R* = {a2:.2f}, b = {bref:.4f} {bd:+.5f} (year - {REF_YEAR:.0f})) at the sector year; depth = instantaneous peak depth',
            'sector,year,n_transits,trapezoid_depth_ppm,trapezoid_depth_err_ppm,trapezoid_T14_min,trapezoid_T14_err_min,trapezoid_f12,trapezoid_f12_err,trapezoid_dt0_min,M1_b,M1_b_err,M2_b,M2_depth_ppm,M2_T14_min'])
    for s in JV['SPOC PDCSAP + QLP sector 10'].secs:
        yr_ = JV['SPOC PDCSAP + QLP sector 10'].W[s]['year']; b_ = bref + bd * (yr_ - REF_YEAR)
        t = TR.get(s); tr = (f"{t['depth']:.0f},{t['edepth']:.0f},{t['T14']:.1f},{t['eT14']:.1f},{t['f12']:.2f},{t['ef12']:.2f},{t['dt0']:+.1f}" if t else ',,,,,,')
        f.write(f"{s},{yr_:.3f},{JV['SPOC PDCSAP + QLP sector 10'].W[s]['ntr']},{tr},{M1b[s]:.4f},{M1eb[s]:.4f},{b_:.4f},{peak_depth(k2, b_, a2, P) * 1e6:.0f},{model_t14(k2, b_, a2, P):.1f}\n")

# ---------------------------------------------------------------- 7. predictions
PR = []
for name in ('SPOC PDCSAP + QLP sector 10', 'SPOC SAP crowding-corrected'):
    k_, a_, b_r, bd_ = JR[name]['M2'][0]
    for y_ in (2019.0, 2021.0, 2023.0, 2025.0, 2026.0, 2027.0, 2027.25, 2028.0, 2030.0):
        b_ = b_r + bd_ * (y_ - REF_YEAR)
        PR.append(dict(source=f'M2 ({name})', year=y_, b=b_, depth=peak_depth(k_, b_, a_, P) * 1e6, T14=model_t14(k_, b_, a_, P), grazing=int(b_ > 1 - k_)))
    for lab, bb in (('b = 1 - k (end of grazing)', 1 - k_), ('b = 1 + k (first transit)', 1 + k_)):
        PR.append(dict(source=f'M2 ({name}): {lab}', year=REF_YEAR + (bb - b_r) / bd_, b=bb, depth=float('nan'), T14=float('nan'), grazing=-1))
Ld, Lt = OUT['depth_trend'], OUT['t14_trend']
for y_ in (2027.0, 2027.25, 2027.5):
    dy = y_ - REF_YEAR
    PR.append(dict(source='linear depth trend (SPOC PDCSAP + QLP sector 10)', year=y_, b=float('nan'), depth=Ld['a'] + Ld['b'] * dy,
                   edepth=float(np.hypot(Ld['ea'], Ld['eb'] * dy)), T14=float('nan'), grazing=-1))
    PR.append(dict(source='linear trapezoid T14 trend (SPOC PDCSAP; bootstrap errors)', year=y_, b=float('nan'), depth=float('nan'), T14=Lt['a'] + Lt['b'] * dy,
                   eT14=float(np.hypot(Lt['ea'], Lt['eb'] * dy)), grazing=-1))
with open(os.path.join(args.out, 'predictions_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: extrapolations of the fitted models; scripts/toi2638/fit_toi2638.py, {DATE}',
            'M2 rows: model M2 of jointfit_458342769.csv (b linear in time), instantaneous peak depth and first-to-fourth-contact T14 of the circular model; grazing = 1 if b > 1 - k; '
            'the rows "b = 1 - k" and "b = 1 + k" give the year at which the linear b(t) reaches that value',
            'linear-trend rows: weighted linear fits of the per-sector depths (errors scaled by sqrt(chi2/dof)) and trapezoid durations (bootstrap errors, not calibrated); '
            'errors propagate the intercept and slope errors only (correlation neglected); these are extrapolations, not forecasts of a model',
            'source,year,b,depth_ppm,depth_err_ppm,T14_min,T14_err_min,grazing'])
    for r_ in PR:
        g = lambda k_, fmt: '' if (k_ not in r_ or not np.isfinite(r_[k_])) else format(r_[k_], fmt)
        f.write(f"{r_['source']},{r_['year']:.2f},{g('b', '.4f')},{g('depth', '.0f')},{g('edepth', '.0f')},{g('T14', '.1f')},{g('eT14', '.1f')},{'' if r_['grazing'] < 0 else r_['grazing']}\n")

OUT['M2'] = {n: list(map(float, JR[n]['M2'][0])) for n in JR}; OUT['M5'] = {n: list(map(float, JR[n]['M5'][0])) for n in JR}
OUT['trap'] = {str(s): v for s, v in TR.items()}
json.dump(OUT, open(os.path.join(args.work, str(TIC), 'fit.json'), 'w'), indent=1, default=float)
print('written', args.out)
