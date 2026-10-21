"""TOI-6482.01 (TIC 189300053): common-shape fit of six TESS QLP transits, mid-transit times, linear ephemeris,
per-transit depths, odd/even and phase-window tests (README section 6.8).

usage: python scripts/toi6482/fit_toi6482.py [--work ./work] [--out data] [--no-download] [--trappist FILE]
input:  QLP light curves hlsp_qlp_tess_ffi_s00NN-0000000189300053_tess_v01_llc.fits (NN = 08, 35, 62, 89, 99), downloaded
        from MAST (HLSP QLP) to <work>/189300053/dl/ unless already present (--no-download: use the files there).
        --trappist FILE: the TRAPPIST-South light curve of 2025-04-21, reduction v2 provided by the observer
        (K. Barkaoui; not distributed in this repository); without it the ground-based refit is skipped and the
        TRAPPIST row of the transit-time table is copied from the existing table if present.
output: <out>/transit_times_189300053.csv, ephemeris_189300053.csv, phase_windows_189300053.csv,
        depth_tests_189300053.csv, transits_189300053.csv (formats in the file headers).
Model: trapezoid (depth, total duration T14, ingress duration shared by the six transits), averaged over each exposure
(10 sub-steps; 1800 s in sector 8, 600 s in sector 35, 200 s in sectors 62, 89, 99), times a quadratic baseline per
transit, over +-0.7 d around the preliminary ephemeris T0 = BTJD 1522.18710, P = 20.0479129 d; SAP_FLUX with QUALITY = 0,
divided by its sector median; unweighted least squares. Mid-times: refit per transit with the shape fixed (mid-time and
baseline free). Timing error = larger of (formal error x beta) and the prayer-bead scatter; beta = max(1, mean over bins
of 3, 5 and 10 cadences of rms(binned residuals) / (rms / sqrt(n))); prayer bead = 50 cyclic shifts of the residuals by
multiples of N//50 cadences. Ephemeris: weighted linear fit, errors scaled by sqrt(chi2/dof) when it exceeds 1.
Depths at fixed shape and mid-time (per transit) or at fixed shape and predicted time (phase windows 0.5, 1/3, 2/3),
with the quadratic baseline solved simultaneously. TRAPPIST: same shape, 60-s exposures, baseline = constant + time +
airmass + a free flux offset after the meridian flip (largest jump of dx_MOVE/dy_MOVE); mid-time from a 0.0005-d grid
search refined by least squares, prayer bead over all cyclic shifts with each shifted light curve fitted the same way.
Requires numpy, scipy, astropy, astroquery (download only).
"""
import argparse, csv, glob, os, sys
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2dist
from astropy.io import fits
from astropy.time import Time

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
p.add_argument('--no-download', action='store_true'); p.add_argument('--trappist', default=None)
args = p.parse_args()
TIC = 189300053; SECTORS = [8, 35, 62, 89, 99]; CAD = {8: 1800, 35: 600, 62: 200, 89: 200, 99: 200}
T0_PRE, P_PRE = 1522.18710, 20.0479129          # preliminary ephemeris (window centres only)
EPOCHS = [(8, 0), (35, 37), (62, 74), (89, 110), (99, 126), (99, 127)]
HW = 0.7; NSUB = 10; NREF = 100; DATE = '2026-10-21'
dl = os.path.join(args.work, str(TIC), 'dl'); os.makedirs(dl, exist_ok=True); os.makedirs(args.out, exist_ok=True)


def qlp_files():
    fs = {s: glob.glob(os.path.join(dl, '**', f'hlsp_qlp_tess_ffi_s{s:04d}-{TIC:016d}_tess_v01_llc.fits'), recursive=True) for s in SECTORS}
    return {s: v[0] for s, v in fs.items() if v}


fs = qlp_files()
if len(fs) < len(SECTORS) and not args.no_download:
    from astroquery.mast import Observations
    obs = Observations.query_criteria(provenance_name='QLP', target_name=str(TIC))
    pl = Observations.get_product_list(obs)
    pl = pl[[str(x).endswith('_llc.fits') for x in pl['productFilename']]]
    Observations.download_products(pl, download_dir=dl)
    fs = qlp_files()
missing = [s for s in SECTORS if s not in fs]
if missing:
    sys.exit(f'missing QLP light curves for sectors {missing} in {dl}')

LC = {}
for s in SECTORS:
    d = fits.getdata(fs[s], 1); t = np.array(d['TIME'], float); f = np.array(d['SAP_FLUX'], float); q = np.array(d['QUALITY'])
    ok = (q == 0) & np.isfinite(t) & np.isfinite(f); LC[s] = (t[ok], f[ok] / np.median(f[ok]))


def trap(x, T, ri):
    return np.clip((T / 2 - np.abs(x)) / ri, 0, 1)


def shape(x, T, ri, dur, n=NSUB):
    off = (np.arange(n) + 0.5) / n - 0.5
    return trap(x[:, None] + off * dur, T, ri).mean(1)


def beta_factor(r):
    sig = np.std(r); b = []
    for nb in (3, 5, 10):
        M = len(r) // nb; bb = r[:M * nb].reshape(M, nb).mean(1); b.append(np.std(bb) / (sig / np.sqrt(nb)))
    return max(1.0, float(np.mean(b)))


def depth_fit(t, f, tc, dur, T, ri, x0):
    """depth (and white error) at fixed shape and mid-time, multiplicative quadratic baseline solved simultaneously"""
    x = t - x0; s = shape(t - tc, T, ri, dur); A = np.column_stack([np.ones_like(x), x, x**2])
    z = np.linalg.lstsq(np.column_stack([A, -s]), f, rcond=None)[0]
    for _ in range(3):
        B = A @ z[:3]; A2 = np.column_stack([A, -s * B]); z = np.linalg.lstsq(A2, f, rcond=None)[0]
    r = f - A2 @ z; C = np.linalg.inv(A2.T @ A2) * (r @ r) / (len(f) - 4)
    return z[3], np.sqrt(C[3, 3]), beta_factor(r), B


# --- segments and joint fit
SEG = []
for s, e in EPOCHS:
    g = T0_PRE + e * P_PRE; t, f = LC[s]; m = np.abs(t - g) < HW; SEG.append((s, e, t[m], f[m], g))


def seg_model(i, tc, c, dep, T, ri):
    s, e, t, f, g = SEG[i]; x = t - g
    return (c[0] + c[1] * x + c[2] * x**2) * (1 - dep * shape(t - tc, T, ri, CAD[s] / 86400))


def res_joint(p):
    return np.concatenate([SEG[i][3] - seg_model(i, p[3 + i], p[9 + 3 * i:12 + 3 * i], *p[:3]) for i in range(6)])


p0 = [0.0094, 0.213, 0.045] + [sg[4] for sg in SEG] + [1, 0, 0] * 6
J = least_squares(res_joint, p0, x_scale='jac', xtol=1e-12, ftol=1e-12, gtol=1e-12)
dep, T14, Ti = J.x[:3]; s2 = (J.fun @ J.fun) / (len(J.fun) - len(J.x)); CJ = np.linalg.inv(J.jac.T @ J.jac) * s2
edep, eT14, eTi = np.sqrt(np.diag(CJ))[:3]

# --- per-transit mid-times at fixed shape
rows = []; LCROWS = []
for i, (s, e, t, f, g) in enumerate(SEG):
    dur = CAD[s] / 86400
    rr = lambda q, y: y - (q[1] + q[2] * (t - g) + q[3] * (t - g)**2) * (1 - dep * shape(t - q[0], T14, Ti, dur))
    q0 = [J.x[3 + i]] + list(J.x[9 + 3 * i:12 + 3 * i])
    S = least_squares(rr, q0, args=(f,), x_scale='jac', xtol=1e-12, ftol=1e-12, gtol=1e-12)
    C = np.linalg.inv(S.jac.T @ S.jac) * (S.fun @ S.fun) / (len(f) - 4); fe = np.sqrt(C[0, 0]) * 1440
    r = S.fun; b = beta_factor(r); mod = f - r; step = max(1, len(f) // 50)
    pb = [least_squares(rr, S.x, args=(mod + np.roll(r, k * step),), x_scale='jac').x[0] for k in range(50)]
    pbe = np.std(pb) * 1440; inn = np.abs(t - S.x[0]) < T14 / 2
    rows.append(dict(sector=s, epoch=e, cadence_s=CAD[s], tc=S.x[0], err=max(fe * b, pbe), fe=fe, beta=b, pb=pbe, n=len(t),
                     cov=inn.sum() * dur / T14, base=S.x[1:]))

# --- linear ephemeris
E = np.array([r['epoch'] for r in rows], float); tc = np.array([r['tc'] for r in rows]); et = np.array([r['err'] for r in rows]) / 1440
A = np.column_stack([np.ones_like(E), E]); W = 1 / et**2; Cs = np.linalg.inv(A.T @ (A * W[:, None])); cc = Cs @ (A.T @ (W * tc))
chi2 = float(np.sum((tc - A @ cc)**2 * W)); dof = len(E) - 2; sc = max(1.0, np.sqrt(chi2 / dof)); Cs = Cs * sc**2
oc = (tc - A @ cc) * 1440
pred = lambda n: (cc[0] + n * cc[1], np.sqrt(np.array([1, n]) @ Cs @ np.array([1, n])))
Tref, eTref = pred(NREF)

# --- per-transit depths, odd/even
for r, (s, e, t, f, g) in zip(rows, SEG):
    d_, de_, b_, B = depth_fit(t, f, r['tc'], CAD[s] / 86400, T14, Ti, g); r.update(dep=d_, edep=de_, dbeta=b_)
    model = 1 - dep * shape(t - r['tc'], T14, Ti, CAD[s] / 86400)
    x = t - g; base = r['base'][0] + r['base'][1] * x + r['base'][2] * x**2
    lin = cc[0] + e * cc[1]
    for ti, fi, mi in zip(t, f / base, model):
        LCROWS.append((s, e, ti, (ti - lin) * 24, fi, mi))


def wmean(v, e):
    w = 1 / np.asarray(e)**2; return float(np.sum(w * v) / w.sum()), float(np.sqrt(1 / w.sum()))


dv = np.array([r['dep'] for r in rows]) * 1e3; de = np.array([r['edep'] for r in rows]) * 1e3; par = E % 2
even = wmean(dv[par == 0], de[par == 0]); odd = wmean(dv[par == 1], de[par == 1]); allm = wmean(dv, de)
chi2_dep = float(np.sum((dv - allm[0])**2 / de**2))
p_dep = float(chi2dist.sf(chi2_dep, len(dv) - 1))
ty = (tc - 3000.0) / 365.25; Ad = np.column_stack([np.ones_like(ty), ty]); Wd = 1 / de**2   # weighted linear depth trend (ppt/yr)
Cd = np.linalg.inv(Ad.T @ (Ad * Wd[:, None])); cd = Cd @ (Ad.T @ (Wd * dv)); chi2_lin = float(np.sum((dv - Ad @ cd)**2 * Wd))

# --- phase windows (0 = the transits; 0.5 and 1/3, 2/3: aliases P/2, P/3 and secondary eclipse)
PW = []
for ph, alias in ((0.5, 'P/2'), (1 / 3, 'P/3'), (2 / 3, 'P/3')):
    for s in SECTORS:
        t0, f0 = LC[s]
        for n in np.arange(np.floor((t0.min() - cc[0]) / cc[1]) - 1, np.ceil((t0.max() - cc[0]) / cc[1]) + 1):
            tw = cc[0] + (n + ph) * cc[1]; m = np.abs(t0 - tw) < HW
            if m.sum() < 10:
                continue
            d_, de_, b_, _ = depth_fit(t0[m], f0[m], tw, CAD[s] / 86400, T14, Ti, tw)
            cov = (np.abs(t0[m] - tw) < T14 / 2).sum() * CAD[s] / 86400 / T14
            PW.append(dict(alias=alias, phase=ph, sector=s, t=tw, n=int(m.sum()), cov=cov, dep=d_ * 1e3, ew=de_ * 1e3, beta=b_, e=de_ * b_ * 1e3))
sel = lambda a: [w for w in PW if w['alias'] == a]
p2 = wmean([w['dep'] for w in sel('P/2')], [w['e'] for w in sel('P/2')]); p3 = wmean([w['dep'] for w in sel('P/3')], [w['e'] for w in sel('P/3')])
sec = wmean([w['dep'] for w in sel('P/2')], [w['ew'] for w in sel('P/2')])
dm_shape = [2.5 * np.log10((ti / (T14 - ti))**2 / dep) for ti in (Ti, Ti + eTi, Ti + 3 * eTi)]

# --- TRAPPIST-South refit (optional; data file provided by the observer)
TR = None
if args.trappist:
    d = np.genfromtxt(args.trappist, names=True); tt = d['BJDTDB'] - 2457000; ff = d['DIFF_FLUX'] / np.median(d['DIFF_FLUX']); am = d['AIRMASS']
    jump = np.hypot(np.diff(d['dx_MOVE']), np.diff(d['dy_MOVE'])); kf = int(np.argmax(jump)); flip = (np.arange(len(tt)) > kf).astype(float)
    ntr = int(round((np.median(tt) - cc[0]) / cc[1])); tp, etp = pred(ntr)

    def tfit(m, cols):
        def design(dt):
            s = shape(tt[m] - tp - dt, T14, Ti, 60 / 86400); X = [np.ones(m.sum())]
            if 'time' in cols: X.append(tt[m] - tp)
            if 'air' in cols: X.append(am[m] - am.mean())
            if 'flip' in cols: X.append(flip[m])
            return np.column_stack(X), s

        def rsd(q, y):
            X, s = design(q[0]); z = np.linalg.lstsq(np.column_stack([X, -s]), y, rcond=None)[0]
            for _ in range(3):
                B = X @ z[:-1]; A2 = np.column_stack([X, -s * B]); z = np.linalg.lstsq(A2, y, rcond=None)[0]
            return y - A2 @ z, z, A2
        def best(y):
            """grid search over the mid-time (0.0005-d steps), then least-squares refinement from the best grid point"""
            grid = np.arange(-0.02, 0.02, 0.0005); g0 = grid[np.argmin([np.sum(rsd([g], y)[0]**2) for g in grid])]
            return least_squares(lambda q: rsd(q, y)[0], [g0], diff_step=1e-3, x_scale=[0.003], xtol=1e-12, ftol=1e-12)
        y = ff[m]; S = best(y)
        r, z, A2 = rsd(S.x, y); n = m.sum(); k = A2.shape[1] + 1; s2 = (r @ r) / (n - k)
        fe = np.sqrt(s2 / (S.jac.T @ S.jac)[0, 0]) * 1440; b = beta_factor(r); mod = y - r
        pb = [best(mod + np.roll(r, sh)).x[0] for sh in range(1, n)]   # every cyclic shift, each fitted like the data
        ed = np.sqrt((np.linalg.inv(A2.T @ A2) * s2)[-1, -1])
        return dict(tc=tp + S.x[0], fe=fe, beta=b, pb=np.std(pb) * 1440, err=max(fe * b, np.std(pb) * 1440), dep=z[-1] * 1e3, edep=ed * 1e3, n=int(n))
    full = tfit(np.ones(len(tt), bool), ('time', 'air', 'flip')); post = tfit(flip == 1, ('time', 'air'))
    TR = dict(epoch=ntr, pred=tp, epred=etp * 1440, full=full, post=post, tflip=0.5 * (tt[kf] + tt[kf + 1]))

# --- write tables
def hdr(f, lines):
    """comment lines, then the column line (last element) without '#'"""
    for l in lines[:-1]:
        f.write('# ' + l + '\n')
    f.write(lines[-1] + '\n')
ttfile = os.path.join(args.out, 'transit_times_189300053.csv')
old_tr = None
if TR is None and os.path.exists(ttfile):
    with open(ttfile) as f:
        old_tr = [l for l in f if l.startswith('TRAPPIST')]
with open(ttfile, 'w') as f:
    hdr(f, [f'TOI-6482.01 = TIC 189300053: mid-transit times; scripts/toi6482/fit_toi6482.py, {DATE}',
            'TESS: QLP SAP_FLUX (MAST HLSP, v01), QUALITY = 0; common trapezoid shape, mid-time and quadratic baseline per transit (+-0.7 d)',
            'tc_err_min = max(formal_err_min x beta, prayer_bead_min); depth_ppt, depth_err_ppt = depth at fixed shape and mid-time (white error)',
            f'o_c_min relative to the linear ephemeris of data/ephemeris_189300053.csv (TESS rows only, in_fit = 1); tc_btjd = BJD_TDB - 2457000',
            'TRAPPIST-South row (in_fit = 0): Rc, 2025-04-21, reduction v2 provided by the observer (K. Barkaoui); same shape, baseline constant + time',
            '  + airmass + flux offset after the meridian flip; o_c_err_min includes the ephemeris uncertainty at that epoch; the flip falls during ingress',
            'instrument,sector,epoch,cadence_s,tc_btjd,tc_err_min,formal_err_min,beta,prayer_bead_min,o_c_min,o_c_err_min,in_transit_coverage,depth_ppt,depth_err_ppt,in_fit'])
    for r, o in zip(rows, oc):
        f.write(f"TESS QLP,{r['sector']},{r['epoch']},{r['cadence_s']},{r['tc']:.5f},{r['err']:.1f},{r['fe']:.1f},{r['beta']:.2f},{r['pb']:.1f},{o:+.1f},{r['err']:.1f},{min(r['cov'],1):.2f},{r['dep']*1e3:.2f},{r['edep']*1e3:.2f},1\n")
    if TR:
        fu = TR['full']; o = (fu['tc'] - TR['pred']) * 1440
        f.write(f"TRAPPIST-South,,{TR['epoch']},60,{fu['tc']:.5f},{fu['err']:.1f},{fu['fe']:.1f},{fu['beta']:.2f},{fu['pb']:.1f},{o:+.1f},{np.hypot(fu['err'],TR['epred']):.1f},1.00,{fu['dep']:.2f},{fu['edep']:.2f},0\n")
    elif old_tr:
        f.writelines(old_tr)

with open(os.path.join(args.out, 'ephemeris_189300053.csv'), 'w') as f:
    hdr(f, [f'TOI-6482.01 = TIC 189300053: linear ephemeris and common transit shape from six TESS QLP transits (sectors 8-99); scripts/toi6482/fit_toi6482.py, {DATE}',
            'T0_btjd at epoch 0 (sector-8 transit); errors and cov_T0_P (d^2) scaled by sqrt(chi2/dof) (scale); chi2 with the unscaled timing errors',
            f'Tref_btjd at epoch n_ref, cov_Tref_P (d^2) its covariance with P (correlation {(Cs[0,1]+NREF*Cs[1,1])/(eTref*np.sqrt(Cs[1,1])):+.3f}); T0 and P are correlated ({Cs[0,1]/np.sqrt(Cs[0,0]*Cs[1,1]):+.3f}), so propagate from Tref:',
            '  t(n) = Tref_btjd + (n - n_ref) P_d, sigma(n)^2 = Tref_err_d^2 + (n - n_ref)^2 P_err_d^2 + 2 (n - n_ref) cov_Tref_P (equivalently from T0 with cov_T0_P)',
            'depth (no dilution correction), T14 and ingress with errors from the unweighted joint fit; P_prelim, T0_prelim = window centres',
            'T0_btjd,T0_err_d,P_d,P_err_d,cov_T0_P,chi2,dof,scale,n_ref,Tref_btjd,Tref_err_d,cov_Tref_P,depth_ppt,depth_err_ppt,T14_h,T14_err_h,ingress_h,ingress_err_h,T0_prelim_btjd,P_prelim_d'])
    covRP = Cs[0, 1] + NREF * Cs[1, 1]
    f.write(f"{cc[0]:.7f},{np.sqrt(Cs[0,0]):.7f},{cc[1]:.8f},{np.sqrt(Cs[1,1]):.9f},{Cs[0,1]:.6e},{chi2:.2f},{dof},{sc:.3f},{NREF},{Tref:.7f},{eTref:.7f},{covRP:.4e},"
            f"{dep*1e3:.2f},{edep*1e3:.2f},{T14*24:.3f},{eT14*24:.3f},{Ti*24:.3f},{eTi*24:.3f},{T0_PRE},{P_PRE}\n")

with open(os.path.join(args.out, 'phase_windows_189300053.csv'), 'w') as f:
    hdr(f, [f'TOI-6482.01 = TIC 189300053: depth at fixed transit shape in every +-0.7-d window at orbital phase 0.5, 1/3 and 2/3 with >= 10 QUALITY = 0 cadences; {DATE}',
            'alias = period tested (P/2: phase 0.5, also the secondary-eclipse window; P/3: phases 1/3 and 2/3); t_btjd = window centre on the linear ephemeris',
            'in_window_coverage = cadences within T14/2 of the centre x cadence / T14; depth_err_ppt includes beta, depth_err_white_ppt does not; transit depth 9.4 ppt',
            'alias,phase,sector,t_btjd,n_cadences,in_window_coverage,depth_ppt,depth_err_ppt,depth_err_white_ppt,beta'])
    for w in PW:
        f.write(f"{w['alias']},{w['phase']:.3f},{w['sector']},{w['t']:.2f},{w['n']},{min(w['cov'],1):.2f},{w['dep']:+.2f},{w['e']:.2f},{w['ew']:.2f},{w['beta']:.2f}\n")

with open(os.path.join(args.out, 'depth_tests_189300053.csv'), 'w') as f:
    hdr(f, [f'TOI-6482.01 = TIC 189300053: depth tests on the TESS QLP data (shape fixed at the common fit); {DATE}',
            'per-transit depths: data/transit_times_189300053.csv; means are inverse-variance weighted with white errors except where noted',
            'shape_limit: largest magnitude difference of a blended eclipsing source that can give the flat-bottomed shape, 2.5 log10((T12/T13)^2/depth) (Vanderburg et al. 2019)',
            'test,value,error,unit,note'])
    rowsD = [('depth_even', *even, 'ppt', 'epochs 0 74 110 126'), ('depth_odd', *odd, 'ppt', 'epochs 37 127'),
             ('odd_minus_even', odd[0] - even[0], np.hypot(odd[1], even[1]), 'ppt', f'{abs(odd[0]-even[0])/np.hypot(odd[1],even[1]):.1f} sigma'),
             ('depth_weighted_mean_per_transit', *allm, 'ppt', f'chi2 = {chi2_dep:.1f} for 5 dof about the mean (p = {p_dep:.2f})'),
             ('depth_linear_trend', cd[1], np.sqrt(Cd[1, 1]), 'ppt/yr', f'weighted linear fit of the per-transit depths against time (white errors); {abs(cd[1])/np.sqrt(Cd[1,1]):.2f} sigma; chi2 = {chi2_lin:.1f} for 4 dof'),
             ('secondary_phase_0.5', *sec, 'ppt', f'3 windows; 3-sigma upper limit {sec[0]+3*sec[1]:.2f} ppt; with beta: error {wmean([w["dep"] for w in sel("P/2")],[w["e"] for w in sel("P/2")])[1]:.2f}'),
             ('alias_P/2', *p2, 'ppt', '3 windows; errors with beta'), ('alias_P/3', *p3, 'ppt', '10 windows; errors with beta'),
             ('shape_limit_dm_nominal', dm_shape[0], '', 'mag', 'ingress as fitted'), ('shape_limit_dm_ingress_+1sigma', dm_shape[1], '', 'mag', ''),
             ('shape_limit_dm_ingress_+3sigma', dm_shape[2], '', 'mag', ''), ('blend_dm_100pct_eclipse', 2.5 * np.log10(1 / dep), '', 'mag', ''),
             ('blend_dm_50pct_eclipse', 2.5 * np.log10(0.5 / dep), '', 'mag', '')]
    for name, v, e, u, n in rowsD:
        f.write(f"{name},{v:.2f},{'' if e == '' else f'{e:.2f}'},{u},{n}\n")

with open(os.path.join(args.out, 'transits_189300053.csv'), 'w') as f:
    hdr(f, [f'TOI-6482.01 = TIC 189300053: the six TESS QLP transit windows (QUALITY = 0), flux divided by the fitted per-transit baseline; {DATE}',
            'dt_lin_h = hours from the linear ephemeris; model = common trapezoid at the fitted mid-time, exposure-averaged',
            'sector,epoch,t_btjd,dt_lin_h,flux,model'])
    for s, e, ti, dh, fi, mi in LCROWS:
        f.write(f'{s},{e},{ti:.6f},{dh:.4f},{fi:.6f},{mi:.6f}\n')

# --- report
print(f'shape: depth {dep*1e3:.2f} +- {edep*1e3:.2f} ppt, T14 {T14*24:.3f} +- {eT14*24:.3f} h, ingress {Ti*24:.3f} +- {eTi*24:.3f} h')
for r, o in zip(rows, oc):
    print(f"S{r['sector']:<3d} E{r['epoch']:<4d} tc {r['tc']:.5f} err {r['err']:.1f} (formal {r['fe']:.2f} beta {r['beta']:.3f} pb {r['pb']:.2f}) O-C {o:+.1f} cov {r['cov']:.2f} depth {r['dep']*1e3:.2f} +- {r['edep']*1e3:.2f}")
print(f'T0 {cc[0]:.5f} +- {np.sqrt(Cs[0,0]):.5f}  P {cc[1]:.6f} +- {np.sqrt(Cs[1,1]):.6f}  cov {Cs[0,1]:.4g}  chi2 {chi2:.2f}/{dof}  scale {sc:.3f}  Tref({NREF}) {Tref:.5f} +- {eTref:.5f}')
print(f'depth: chi2 {chi2_dep:.2f}/5 about the mean (p {p_dep:.3f}); linear trend {cd[1]:+.3f} +- {np.sqrt(Cd[1,1]):.3f} ppt/yr (chi2 {chi2_lin:.2f}/4)')
print(f'even {even[0]:.2f} +- {even[1]:.2f}  odd {odd[0]:.2f} +- {odd[1]:.2f}  secondary {sec[0]:+.2f} +- {sec[1]:.2f}  P/2 {p2[0]:+.2f} +- {p2[1]:.2f}  P/3 {p3[0]:+.2f} +- {p3[1]:.2f}  shape-limit dm {np.round(dm_shape,2)}')
for n in [113, 114] + list(range(141, 147)):
    tp, ep = pred(n); print(f'predicted epoch {n}: BTJD {tp:.3f} +- {ep*1440:.1f} min ({Time(tp + 2457000, format="jd", scale="tdb").utc.iso[:16]} UTC, barycentric)')
if TR:
    fu, po = TR['full'], TR['post']
    print(f"TRAPPIST epoch {TR['epoch']} (flip at BTJD {TR['tflip']:.4f}): tc {fu['tc']:.5f} err {fu['err']:.1f} min (formal {fu['fe']:.2f} x beta {fu['beta']:.2f} = {fu['fe']*fu['beta']:.2f}; prayer bead {fu['pb']:.2f}), "
          f"O-C {(fu['tc']-TR['pred'])*1440:+.1f} +- {np.hypot(fu['err'],TR['epred']):.1f} min, depth {fu['dep']:.2f} +- {fu['edep']:.2f} ppt; "
          f"after the flip only: O-C {(po['tc']-TR['pred'])*1440:+.1f} +- {po['err']:.1f} min (timing error; formal {po['fe']:.2f} x beta {po['beta']:.2f}, prayer bead {po['pb']:.2f}), "
          f"+- {np.hypot(po['err'],TR['epred']):.1f} min with the ephemeris uncertainty ({TR['epred']:.2f} min), depth {po['dep']:.2f} +- {po['edep']:.2f} ppt, N = {po['n']}; N = {fu['n']} in all")
