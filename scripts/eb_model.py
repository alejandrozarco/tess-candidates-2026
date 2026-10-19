"""TIC 11613676 (CD-30 448): geometric eclipse models fitted to all eclipse windows of sectors 3-106.

usage: python scripts/eb_model.py [--models M0,MA,N4] [--work ./work] [--out data]
input:  <work>/11613676/lc_s<sector>.npz (eb_eclipses.py), SPOC 2-min light curves of sectors 3 and 30 (other_download.py)
output: data/eclipse_model_11613676.csv (parameters, errors, chi2 of each model) and data/model_history_11613676.csv
        (inclination, impact parameters, model eclipse depths and durations of each model on a time grid 2018-2029).
Model: two spheres with linear limb darkening, Keplerian orbit, third light 0, exposure-time
integration (limb darkening u = 0.6 for star 1, 0.7 for star 2); r1 = R1/a, k = R2/R1, J = surface-brightness ratio. Every window (+-0.3 d around each primary and
secondary eclipse; SPOC 2-min SAP flux in sectors 3 and 30, TESScut r = 2.5 px otherwise) is divided by a linear
out-of-eclipse fit and multiplied by a linear baseline solved at every evaluation; errors per window from the
out-of-eclipse scatter. Models: M0 static geometry; MA apsidal motion only (d omega/dt free, inclination constant);
N4 e cos(omega), e sin(omega) linear in time and inclination i(t) = i_L + A cos(2 pi (t - t_max)/P_node). Errors are
from the Jacobian, scaled by chi2/dof. The fits are not formally good (chi2/N > 2.5); see README.
"""
import argparse, csv, glob, os
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
from astropy.time import Time
from targets import OTHER

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--models', default='M0,MA,N4'); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
args = p.parse_args(); TIC = 11613676; c = OTHER[TIC]; wd = os.path.join(args.work, str(TIC))
T0G, P0 = c['eph0']; PH2 = c['ph2']; TREF = 3979.6; HW = 0.3; U = (0.6, 0.7)
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
RR = (np.arange(400) + 0.5) / 400.0


def blocked(z, k, u):
    """fraction of the light of a linearly limb-darkened unit disc hidden by a disc of radius k at projected distance z"""
    z = np.atleast_1d(z)[:, None]; r = RR[None, :]; I = 1 - u * (1 - np.sqrt(1 - r**2)); wgt = 2 * r * I
    frac = np.where(r <= k - z, 1.0, np.where((r >= z + k) | (r <= z - k), 0.0, np.arccos(np.clip((r**2 + z**2 - k**2) / (2 * r * z + 1e-12), -1, 1)) / np.pi))
    return (wgt * frac).sum(1) / wgt.sum(1)


def flux(t, P, Tc, e, om, inc, r1, k, J):
    nu_c = np.pi / 2 - om; Ec = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu_c / 2), np.sqrt(1 + e) * np.cos(nu_c / 2))
    M = 2 * np.pi / P * (t - Tc) + Ec - e * np.sin(Ec); E = np.mod(M, 2 * np.pi)
    for _ in range(30): E = E - (E - e * np.sin(E) - np.mod(M, 2 * np.pi)) / (1 - e * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    rs = (1 - e**2) / (1 + e * np.cos(nu)); front = np.sin(nu + om) > 0
    z = rs * np.sqrt(1 - np.sin(inc)**2 * np.sin(nu + om)**2) / r1; L2 = J * k**2; f = np.full_like(t, 1 + L2)
    m1 = front & (z < 1 + k); m2 = ~front & (z < 1 + k)
    if m1.any(): f[m1] -= blocked(z[m1], k, U[0])
    if m2.any(): f[m2] -= L2 * blocked(z[m2] / k, 1 / k, U[1])
    return f / (1 + L2)


NAMES = ['P', 'Tc', 'ecosw0', 'decosw_dt', 'esinw0', 'desinw_dt', 'omdot_rad_yr', 'i_deg', 'A_deg', 'Pnode_yr', 'tmax_yr', 'r1', 'k', 'J']
LO = np.array([P0 - 0.001, T0G - 0.05, -0.4, -0.05, -0.4, -0.05, -2.0, 60.0, 0.1, 8.0, -3.0, 0.015, 0.05, 0.01])
HI = np.array([P0 + 0.001, T0G + 0.05, 0.4, 0.05, 0.4, 0.05, 2.0, 95.0, 15.0, 120.0, 2.0, 0.08, 1.0, 1.5])


def geom(pf, dt_yr):
    P, Tc, h0, h1, g0, g1, omd, iL, A, Pn, tmax, r1, k, J = pf
    h, g = h0 + h1 * dt_yr, g0 + g1 * dt_yr; e = np.hypot(h, g); om = np.arctan2(g, h) + omd * dt_yr
    return P, Tc, e, om, np.radians(iL + A * np.cos(2 * np.pi * (dt_yr - tmax) / Pn)), r1, k, J


def model_win(pf, w):
    nsub = max(1, int(round(w['texp'] / 0.0015))); off = ((np.arange(nsub) + 0.5) / nsub - 0.5) * w['texp']
    tt = (w['t'][:, None] + off[None, :]).ravel()
    return flux(tt, *geom(pf, (w['tp'] - TREF) / 365.25)).reshape(-1, nsub).mean(1)


# ---- data windows
LCS = {}
for s in c['sectors']:
    if s in c['spoc_2min']:
        fn = glob.glob(os.path.join(wd, 'dl', 'mastDownload', 'TESS', f'*s{s:04d}-{TIC:016d}*', '*-s_lc.fits'))[0]; d = fits.open(fn)[1].data
        LCS[s] = (np.array(d['TIME'], float), np.array(d['SAP_FLUX'], float), np.array(d['QUALITY']))
    else:
        z = np.load(os.path.join(wd, f'lc_s{s}.npz')); LCS[s] = (z['t'], z['r2.5'], z['q'])
WINS = []
for s, (t, f, q) in LCS.items():
    g = (q == 0) & np.isfinite(f) & np.isfinite(t); t, f = t[g], f[g]; texp = float(np.median(np.diff(t)))
    for lab, ph in (('pri', 0.0), ('sec', PH2)):
        for n in range(int(np.floor((t.min() - T0G) / P0)) - 1, int(np.ceil((t.max() - T0G) / P0)) + 2):
            tp = T0G + (n + ph) * P0; m = np.abs(t - tp) < HW
            if m.sum() < 0.5 * 2 * HW / texp or not np.sum(np.abs(t[m] - tp) < 0.08) > 0.5 * 0.16 / texp: continue
            tt, ff = t[m], f[m]; oot = np.abs(tt - tp) > 0.1; ffn = ff / np.polyval(np.polyfit(tt[oot] - tp, ff[oot], 1), tt - tp)
            WINS.append(dict(sec=s, lab=lab, tp=tp, t=tt, f=ffn, sig=1.4826 * np.median(np.abs(ffn[oot] - np.median(ffn[oot]))), texp=texp))
print(len(WINS), 'windows,', sum(len(w['t']) for w in WINS), 'points', flush=True)


def resid(x, idx, pf0):
    pf = pf0.copy(); pf[idx] = x; out = []
    for w in WINS:
        m = model_win(pf, w); A = np.c_[m, m * (w['t'] - w['tp'])]; cc = np.linalg.lstsq(A, w['f'], rcond=None)[0]; out.append((w['f'] - A @ cc) / w['sig'])
    return np.concatenate(out)


def start(**kw):
    s = dict(P=P0, Tc=T0G, ecosw0=0.085, decosw_dt=0.0, esinw0=0.0, desinw_dt=0.0, omdot_rad_yr=0.0, i_deg=89.0, A_deg=0.0, Pnode_yr=1e6, tmax_yr=0.0, r1=0.035, k=0.3, J=0.25); s.update(kw)
    return np.array([s[n] for n in NAMES], float)


MODELS = {
    'M0': (['P', 'Tc', 'ecosw0', 'esinw0', 'i_deg', 'r1', 'k', 'J'], [start(i_deg=i) for i in (88.3, 89.0, 89.6)], 'static geometry'),
    'MA': (['P', 'Tc', 'ecosw0', 'esinw0', 'omdot_rad_yr', 'i_deg', 'r1', 'k', 'J'], [start(i_deg=i, omdot_rad_yr=o) for i in (88.3, 89.0, 89.6) for o in (0.0, 0.3)], 'apsidal motion only (inclination constant)'),
    'N4': (['P', 'Tc', 'ecosw0', 'decosw_dt', 'esinw0', 'desinw_dt', 'i_deg', 'A_deg', 'Pnode_yr', 'tmax_yr', 'r1', 'k', 'J'],
           [start(ecosw0=0.0935, decosw_dt=0.0061, esinw0=0.0489, i_deg=87.0 - A, A_deg=A, Pnode_yr=Pn, r1=0.05, k=0.3, J=0.45, tmax_yr=0.1) for Pn in (15.0, 25.0, 40.0, 60.0) for A in (1.0, 2.5)],
           'e-vector linear in time, sinusoidal inclination'),
}
LO_A = dict(zip(NAMES, LO)); HI_A = dict(zip(NAMES, HI)); HI_A['Pnode_yr'] = 2e6
RES = {}
for name in args.models.split(','):
    free, starts, label = MODELS[name]; idx = [NAMES.index(n) for n in free]; best = None
    for s in starts:
        r = least_squares(resid, s[idx], args=(idx, s), bounds=([LO_A[n] for n in free], [HI_A[n] for n in free]), x_scale='jac', max_nfev=300)
        print(f'{name} start -> chi2 {2 * r.cost:.1f}', flush=True)
        if best is None or r.cost < best[0].cost: best = (r, s)
    r, s = best; pf = s.copy(); pf[idx] = r.x; chi2 = 2 * r.cost; nd = len(r.fun)
    err = np.zeros(len(NAMES)); err[idx] = np.sqrt(np.diag(np.linalg.inv(r.jac.T @ r.jac) * max(chi2 / (nd - len(idx)), 1)))
    RES[name] = (pf, err, chi2, nd, len(idx), label, free)
    print(f'== {name} ({label}): chi2 = {chi2:.1f} for {nd} points, {len(idx)} free;', {n: (round(pf[NAMES.index(n)], 6), round(err[NAMES.index(n)], 6)) for n in free}, flush=True)

with open(os.path.join(args.out, 'eclipse_model_11613676.csv'), 'w') as fo:
    fo.write(''.join('# ' + l + '\n' for l in [
        f'TIC 11613676 = CD-30 448: geometric eclipse models fitted to {len(WINS)} eclipse windows (scripts/eb_model.py), {DATE}',
        'Parameters: P (d), Tc (BTJD of primary mid-eclipse), ecosw0/esinw0 = e cos(omega), e sin(omega) at BTJD 3979.6 (2025.78) and their rates (per yr), omdot (rad/yr),',
        'i_deg = inclination (constant, or the mean i_L of the sinusoid), A_deg, Pnode_yr, tmax_yr = amplitude, period and time of maximum (yr from BTJD 3979.6) of the inclination sinusoid,',
        'r1 = R1/a, k = R2/R1, J = surface-brightness ratio; fixed: third light 0, linear limb darkening 0.6 (star 1) and 0.7 (star 2); error 0 = parameter not free; chi2 for n_points with n_free parameters']))
    w = csv.writer(fo); w.writerow(['model', 'description', 'chi2', 'n_points', 'n_free', 'parameter', 'value', 'error'])
    for name, (pf, err, chi2, nd, nf, label, free) in RES.items():
        for n in NAMES:
            if n in free or n in ('omdot_rad_yr',): w.writerow([name, label, f'{chi2:.1f}', nd, nf, n, f'{pf[NAMES.index(n)]:.6g}', f'{err[NAMES.index(n)]:.3g}'])

with open(os.path.join(args.out, 'model_history_11613676.csv'), 'w') as fo:
    fo.write(''.join('# ' + l + '\n' for l in [
        f'TIC 11613676 = CD-30 448: model eclipse history from the fits of eclipse_model_11613676.csv (scripts/eb_model.py), {DATE}',
        'year (Julian), i_deg, b_pri, b_sec (projected separation at the conjunctions in units of R1), model depth (%) and T14 (h) of the primary and the secondary eclipse (0 = no eclipse);',
        'the models are extrapolations outside 2018.7-2026.6 and are not predictions with errors']))
    w = csv.writer(fo); w.writerow(['model', 'year', 'i_deg', 'b_pri', 'b_sec', 'depth_pri_pct', 'T14_pri_h', 'depth_sec_pct', 'T14_sec_h'])
    for name, (pf, *_rest) in RES.items():
        for yr in np.arange(2018.0, 2029.01, 0.05):
            t = Time(yr, format='jyear').jd - 2457000; P, Tc, e, om, inc, r1, k, J = geom(pf, (t - TREF) / 365.25); n = np.round((t - Tc) / P); row = [name, f'{yr:.2f}', f'{np.degrees(inc):.3f}']
            row += [f'{(1 - e**2) / (1 + e * np.sin(om)) * np.cos(inc) / r1:.3f}', f'{(1 - e**2) / (1 - e * np.sin(om)) * np.cos(inc) / r1:.3f}']
            for ph in (0.0, 0.5 + 2 * e * np.cos(om) / np.pi):
                tq = Tc + (n + ph) * P; tt = np.linspace(tq - 0.35, tq + 0.35, 3001); f = flux(tt, P, Tc, e, om, inc, r1, k, J); ine = tt[f < 1 - 1e-5]
                row += [f'{(1 - f.min()) * 100:.3f}', f'{(ine.max() - ine.min()) * 24 if len(ine) else 0:.2f}']
            w.writerow(row)
print('written', os.path.join(args.out, 'eclipse_model_11613676.csv'), os.path.join(args.out, 'model_history_11613676.csv'))
