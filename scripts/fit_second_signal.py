"""Joint fit of the three events of the candidate second transit signal of TIC 311179599 (README section 5.10).

usage: python scripts/fit_second_signal.py [--aperture 2.0] [--work ./work] [--out data] [--steps 600000]
input:  <work>/311179599/cut/lc_s<sector>.npz of sectors 12, 100 and 103 (tesscut_download.py --sectors 12,100,103,
        then aperture_photometry.py)
output: data/transit_times_311179599c.csv, fit_311179599c.csv, folded_311179599c.csv and model_311179599c.csv for the
        default aperture (formats in the file headers); for another aperture only <work>/311179599/fit_c_r<ap>.json.
Model: batman transit (circular orbit at the reference period 73.87 d, quadratic limb darkening u = [0.33, 0.27]) with a
common k = Rp/R*, a/R* and b and one mid-time per event, times one quadratic baseline per event over +-0.6 d (solved
linearly); quality-0 cadences; the transits of the 24.8-d signal are masked (+-0.1 d); exposure-time supersampling.
Per-event noise = white rms x max(1, beta at 30 and 60 min). Posterior: Metropolis sampler (600000 steps, the first
10000 discarded, every 5th kept; random seed 2). Linear ephemerides through the three median mid-times for n = 1
(cycles -33, 0, 1) and n = 2 (cycles -66, 0, 2). Per-event depth: the common shape scaled per event (linear fit with the
event's own baseline). The durations and densities are those at the reference period; for n = 2 the circular-orbit
density scales as P^-2 at fixed a/R* shape (factor 4 higher).
"""
import argparse, csv, json, os
import numpy as np
from scipy.optimize import least_squares
from targets import TARGETS, import_batman

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--aperture', default='2.0'); p.add_argument('--work', default='./work')
p.add_argument('--out', default='data'); p.add_argument('--steps', type=int, default=600000)
args = p.parse_args(); TIC = 311179599; c = TARGETS[TIC]; S = c['second']; AP = args.aperture; batman = import_batman()
wd = os.path.join(args.work, str(TIC)); PREF = S['P_ref']; U = c['u']; G = 6.674e-8
EV = S['events']; L = list(EV); BTC = np.array(S['mask_b'])
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
D = {}
for lab_, (s, tc0) in EV.items():
    z = np.load(os.path.join(wd, 'cut', f'lc_s{s}.npz')); f = z[f'f{AP}']; t = z['tm']
    inb = np.min(np.abs(t[:, None] - BTC[None, :]), 1) < 0.1
    m = (np.abs(t - tc0) < 0.6) & np.isfinite(f) & (z['q'] == 0) & ~inb; t, y = t[m], f[m] / np.median(f[m]); cad = np.median(np.diff(t))
    pp = batman.TransitParams(); pp.t0 = tc0; pp.per = PREF; pp.rp = 0.031; pp.a = 100; pp.inc = 89.9; pp.ecc = 0; pp.w = 90; pp.u = U; pp.limb_dark = 'quadratic'
    ns = max(3, int(round(cad * 1440 / 2)))   # midpoint sub-exposures (batman spans exp_time end to end)
    D[lab_] = dict(t=t, y=y, x=t - tc0, p=pp, mod=batman.TransitModel(pp, t, supersample_factor=ns, exp_time=cad * (ns - 1) / ns), cad=cad, tc0=tc0, sector=s)


def lcm(d, tc, k, aR, b):
    pp = d['p']; pp.t0 = tc; pp.rp = k; pp.a = aR; pp.inc = np.degrees(np.arccos(min(b / aR, 1))); return d['mod'].light_curve(pp)


def base(d, mm):
    x = d['x']; A = np.c_[mm, mm * x, mm * x * x]; cc = np.linalg.lstsq(A, d['y'], rcond=None)[0]; return cc, A


def rsd(d, mm):
    cc, A = base(d, mm); return d['y'] - A @ cc


SIG = {l: 1.0 for l in L}
def R(pp):
    k, aR, b = pp[:3]; return np.concatenate([rsd(D[l], lcm(D[l], pp[3 + i], k, aR, b)) / SIG[l] for i, l in enumerate(L)])
lo = np.array([0.005, 10, 0] + [EV[l][1] - 0.1 for l in L]); hi = np.array([0.2, 500, 1.1] + [EV[l][1] + 0.1 for l in L])
best = None
for a0 in [60, 100, 150]:
    for b0 in [0.2, 0.6]:
        F = least_squares(R, [0.031, a0, b0] + [EV[l][1] for l in L], bounds=(lo, hi), x_scale=[0.003, 10, 0.1] + [0.003] * 3)
        if best is None or F.cost < best.cost: best = F
BETA = {}
for i, l in enumerate(L):
    d = D[l]; r = rsd(d, lcm(d, best.x[3 + i], *best.x[:3])); sd = np.std(r); bet = []
    for hrs in [0.5, 1.0]:
        nb = max(1, int(round(hrs / 24 / d['cad'])))
        if nb < 2: continue
        nn = len(r) // nb; bet.append(np.std(r[:nn * nb].reshape(nn, nb).mean(1)) / (sd / np.sqrt(nb)))
    SIG[l] = sd * max([1.0] + bet); BETA[l] = max([1.0] + bet)
best = least_squares(R, best.x, bounds=(lo, hi), x_scale=[0.003, 10, 0.1] + [0.003] * 3)


def lnp(pp):
    if np.any(pp < lo) or np.any(pp > hi) or pp[2] > 1 + pp[0]: return -np.inf
    r = R(pp); return -0.5 * r @ r


rng = np.random.default_rng(2); cur = best.x.copy(); cur[2] = min(cur[2], 0.9); lc_ = lnp(cur); step = np.array([0.0015, 6, 0.06, 0.004, 0.004, 0.006]); ch = []; acc = 0
for i in range(args.steps):
    pr = cur + step * rng.normal(size=6); lp = lnp(pr)
    if np.log(rng.random()) < lp - lc_: cur, lc_ = pr, lp; acc += 1
    if i > 10000 and i % 5 == 0: ch.append(cur.copy())
ch = np.array(ch); k, aR, b = ch[:, 0], ch[:, 1], ch[:, 2]
T14 = PREF / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k)**2 - b**2, 0, None)) / (aR * np.sqrt(1 - (b / aR)**2)), 0, 1)) * 24
rho = 3 * np.pi * aR**3 / (G * (PREF * 86400)**2)
pc = lambda a: [float(v) for v in np.percentile(a, [2.5, 16, 50, 84, 97.5])]
out = dict(aperture=AP, acceptance=acc / args.steps, sig_ppm={l: float(SIG[l] * 1e6) for l in L}, beta=BETA, best=dict(zip(['k', 'aR', 'b'] + ['tc' + l for l in L], map(float, best.x))),
           k2_ppm=pc(k**2 * 1e6), b=pc(b), T14_h=pc(T14), aR=pc(aR), rho_circ_Pref=pc(rho), tc={l: pc(ch[:, 3 + i]) for i, l in enumerate(L)})
tcm = {l: float(np.median(ch[:, 3 + i])) for i, l in enumerate(L)}; tce = {l: float(0.5 * (np.percentile(ch[:, 3 + i], 84) - np.percentile(ch[:, 3 + i], 16))) for i, l in enumerate(L)}
CYC = {1: {'C': -33, 'A': 0, 'B': 1}, 2: {'C': -66, 'A': 0, 'B': 2}}
for n, cyc in CYC.items():
    E_ = np.array([cyc[l] for l in L], float); T = np.array([tcm[l] for l in L]); W = 1 / np.array([tce[l] for l in L])**2
    A = np.c_[np.ones(3), E_]; cv = np.linalg.inv(A.T @ (W[:, None] * A)); cc = cv @ A.T @ (W * T); oc = (T - A @ cc) * 1440
    out[f'eph_n{n}'] = dict(T0=float(cc[0]), P=float(cc[1]), eT0=float(np.sqrt(cv[0, 0])), eP=float(np.sqrt(cv[1, 1])), oc_min={l: float(v) for l, v in zip(L, oc)},
                            chi2=float(np.sum(oc**2 / (np.array([tce[l] for l in L]) * 1440)**2)))
# per-event depth with the best common shape and each event's own mid-time
kb, ab, bb = best.x[:3]; dep = {}
for i, l in enumerate(L):
    d = D[l]; mm = lcm(d, best.x[3 + i], kb, ab, bb); sh = 1 - mm; x = d['x']
    A = np.c_[np.ones_like(x), x, x * x, -sh]; cc = np.linalg.lstsq(A, d['y'], rcond=None)[0]; cov = np.linalg.inv(A.T @ A) * SIG[l]**2
    dep[l] = (float(cc[3] * kb**2 * 1e6), float(np.sqrt(cov[3, 3]) * kb**2 * 1e6))
out['depth_per_event_ppm'] = dep
os.makedirs(wd, exist_ok=True); json.dump(out, open(os.path.join(wd, f'fit_c_r{AP}.json'), 'w'), indent=1)
print(json.dumps({kk: v for kk, v in out.items() if kk not in ('tc',)}, indent=None))
if AP != '2.0':
    raise SystemExit
hdr = lambda lines: ''.join('# ' + s_ + '\n' for s_ in lines)
e1, e2 = out['eph_n1'], out['eph_n2']
CAM = S['camera_ccd']
with open(os.path.join(args.out, 'transit_times_311179599c.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 311179599 = Gaia DR3 5801882036668103680, candidate second transit signal c: mid-transit times of the three events (scripts/fit_second_signal.py, TESScut r = {AP} px aperture photometry), {DATE}',
                  'tc_btjd = median of the joint-fit posterior (BJD_TDB - 2457000); tc_lo_d, tc_hi_d = 16th and 84th percentiles; tc_err_min = half the 16-84 % interval',
                  f'epoch_n1, epoch_n2 = cycle numbers for P = (t_B - t_A)/n; O-C against the three-event linear ephemerides n = 1: T0 = BTJD {e1["T0"]:.5f} +- {e1["eT0"]:.5f}, P = {e1["P"]:.5f} +- {e1["eP"]:.5f} d (chi2 = {e1["chi2"]:.2f}, 1 degree of freedom); n = 2: P = {e2["P"]:.5f} +- {e2["eP"]:.5f} d (same O-C)',
                  'depth_ppm = per-event depth: the best common shape (k^2 of the joint fit) scaled by a linear fit with the event baseline; raw aperture flux, not dilution corrected',
                  'noise_ppm = per-event white rms x beta (time-averaging factor at 30 and 60 min) used in the fit']))
    w = csv.writer(fo); w.writerow(['event', 'sector', 'camera_ccd', 'cadence_min', 'epoch_n1', 'epoch_n2', 'tc_btjd', 'tc_bjd_tdb', 'tc_lo_d', 'tc_hi_d', 'tc_err_min', 'oc_min', 'depth_ppm', 'depth_err_ppm', 'noise_ppm', 'beta'])
    for l in sorted(L, key=lambda x: EV[x][1]):
        q_ = out['tc'][l]; w.writerow([l, EV[l][0], CAM[l], f'{D[l]["cad"] * 1440:.1f}', CYC[1][l], CYC[2][l], f'{q_[2]:.6f}', f'{q_[2] + 2457000:.6f}', f'{q_[1]:.6f}', f'{q_[3]:.6f}',
                   f'{tce[l] * 1440:.1f}', f'{e1["oc_min"][l]:+.1f}', f'{dep[l][0]:.0f}', f'{dep[l][1]:.0f}', f'{SIG[l] * 1e6:.0f}', f'{BETA[l]:.2f}'])
with open(os.path.join(args.out, 'fit_311179599c.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 311179599, candidate second transit signal c: joint fit of the three events (scripts/fit_second_signal.py, r = {AP} px), {DATE}',
                  f'percentiles 2.5/16/50/84/97.5 of the Metropolis posterior; T14 and rho_circ (g cm^-3, circular orbit) at the reference period {PREF} d (rho_circ x 4 for P/2); acceptance {out["acceptance"]:.3f}; best fit k = {kb:.5f}, a/R* = {ab:.2f}, b = {bb:.3f}']))
    w = csv.writer(fo); w.writerow(['quantity', 'p2.5', 'p16', 'p50', 'p84', 'p97.5'])
    for key, lab_ in (('k2_ppm', 'depth_k2_ppm'), ('b', 'b'), ('T14_h', 'T14_h'), ('aR', 'a_over_Rstar'), ('rho_circ_Pref', 'rho_circ_gcc')):
        w.writerow([lab_] + [f'{v:.4g}' for v in out[key]])
    for l in L:
        w.writerow([f'tc_{l}_btjd'] + [f'{v:.6f}' for v in out['tc'][l]])
mod_lines = [f'TIC 311179599, candidate second transit signal c: per-event light curves within +-0.6 d of the fitted mid-times (scripts/fit_second_signal.py, r = {AP} px), {DATE}',
             'flux = aperture flux divided by the event quadratic baseline fitted jointly with the best common model; not dilution corrected; quality-0 cadences; the 24.8-d transits masked',
             f'model = best common model (k = {kb:.5f}, a/R* = {ab:.2f}, b = {bb:.3f}, u = {U}) integrated over the cadence; dt_days = time - tc (best-fit mid-time of the event)']
with open(os.path.join(args.out, 'folded_311179599c.csv'), 'w') as fo:
    fo.write(hdr(mod_lines)); w = csv.writer(fo); w.writerow(['event', 'sector', 'time_btjd', 'dt_days', 'flux', 'model'])
    for i, l in enumerate(sorted(L, key=lambda x: EV[x][1])):
        j = L.index(l); d = D[l]; mm = lcm(d, best.x[3 + j], kb, ab, bb); cc, A = base(d, mm); bl = np.c_[np.ones_like(d['x']), d['x'], d['x']**2] @ cc
        for tt, ff, m_ in zip(d['t'], d['y'] / bl, mm):
            w.writerow([l, d['sector'], f'{tt:.6f}', f'{tt - best.x[3 + j]:.6f}', f'{ff:.6f}', f'{m_:.6f}'])
xx_ = np.arange(-0.5, 0.5 + 1e-9, 0.0005); pm_ = batman.TransitParams(); pm_.t0 = 0; pm_.per = PREF; pm_.rp = kb; pm_.a = ab; pm_.inc = np.degrees(np.arccos(bb / ab)); pm_.ecc = 0; pm_.w = 90; pm_.u = U; pm_.limb_dark = 'quadratic'
fl = batman.TransitModel(pm_, xx_).light_curve(pm_)
with open(os.path.join(args.out, 'model_311179599c.csv'), 'w') as fo:
    fo.write(hdr([f'TIC 311179599, candidate second transit signal c: best common transit model (batman, u = {U}), k = {kb:.5f}, b = {bb:.3f}, a/R* = {ab:.2f}, P_ref = {PREF} d; no exposure-time smearing', 'dt_days = time from mid-transit']))
    w = csv.writer(fo); w.writerow(['dt_days', 'flux'])
    for a_, f_ in zip(xx_, fl): w.writerow([f'{a_:.5f}', f'{f_:.7f}'])
print('tables written')
