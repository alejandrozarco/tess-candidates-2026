"""TOI-2638.01 (TIC 458342769): checks of the transit-depth change (README section 6.9): SPOC header keywords and
target flux per sector, Gaia DR3 neighbours, eclipsing-binary tests, pulsation search, the 7.31 c/d signal of the
17" neighbour, the stellar-oblateness (J2) precession estimate and the future TESS sectors.

usage: python scripts/toi2638/checks_toi2638.py [--work ./work] [--out data] [--no-gaia]
input:  the light curves of download_toi2638.py; data/ephemeris_458342769.csv and data/jointfit_458342769.csv of
        fit_toi2638.py (read from --out); Gaia DR3 (ESA archive TAP) unless --no-gaia, in which case
        <out>/neighbours_458342769.csv is read
output: <out>/sectors_458342769.csv, neighbours_458342769.csv, checks_458342769.csv, precession_458342769.csv
        (formats in the file headers)
Methods:
 - flux level: median SAP_FLUX x CROWDSAP / FLFRCSAP (QUALITY = 0) = the target flux implied by the SPOC crowding
   model; a source that faded or brightened would change it at fixed keywords.
 - odd/even: depth of the common shape in odd and even transits per observing campaign (SPOC PDCSAP windows).
 - secondary eclipse: PDCSAP divided by a running median of the out-of-transit flux in 0.5-d bins, transits removed;
   box depths of 0.040, 0.056 and 0.075 d at phase 0.5 and over phases 0.15-0.85 (step 0.0025) against the mean of
   the neighbouring flux within two box widths.
 - phase curve: per sector, PDCSAP out of transit, quadratic in time plus -A_ell cos(4 pi phase) - A_ref cos(2 pi phase)
   + A_beam sin(2 pi phase); errors scaled by the red-noise factor of 2-h bins; weighted mean over sectors.
 - pulsations: Lomb-Scargle amplitude spectrum of the out-of-transit PDCSAP flux (0.05-360 c/d) per sector.
 - 7.31 c/d: amplitude in the target's PDCSAP flux per sector (best frequency within +-0.01 c/d of 7.311) and the
   strongest peak in the QLP light curves of the neighbour TIC 458342741 (1 c/d to 98 % of the Nyquist frequency).
 - J2 precession: q = (v_eq / v_crit)^2 with v_eq = vsini / sin i_star, v_crit = sqrt(G M / R); J2 = (2 k_ap / 3) q
   (k_ap = apsidal-motion constant); omega0 = (3/2) n J2 (R/a)^2; mutual precession frequency =
   omega0 cos(psi) sqrt(1 + (L/S)^2 + 2 (L/S) cos(psi)), L = orbital angular momentum of the companion,
   S = 0.06 M R^2 Omega_star; |db/dt| <= (a/R) sin(i_star) sin(psi) |omega0 cos(psi)|; the rate envelopes do not impose
   the observed transiting geometry; Gaia DR3 FLAME mass and radius and ESP-HS vsini; order-of-magnitude estimate for
   stated assumptions.
Requires numpy, scipy, astropy, batman-package, pyvo (Gaia), tess-point.
"""
import argparse, csv, os, sys
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (TIC, GAIA_TARGET, GAIA_NB17, RA, DEC, load, windows, Model, peak_depth, hdr, dldir, year, DATE, REF_YEAR)
import glob

ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('--work', default='./work'); ap.add_argument('--out', default='data')
ap.add_argument('--no-gaia', action='store_true')
args = ap.parse_args()


def rows(fn):
    return [r for r in csv.reader(l for l in open(os.path.join(args.out, fn)) if not l.startswith('#'))]


eph = [r for r in rows('ephemeris_458342769.csv')[1:] if r[0] == 'ephemeris'][0]; T0, P = float(eph[10]), float(eph[12])
J = {}
for r in rows('jointfit_458342769.csv')[1:]:
    J.setdefault((r[0], r[1]), []).append(float(r[6]))
k0, b0, a0 = J[('SPOC PDCSAP + QLP sector 10', 'null')]; D0 = peak_depth(k0, b0, a0, P)
M2 = J[('SPOC PDCSAP + QLP sector 10', 'M2')]; M2sap = J[('SPOC SAP crowding-corrected', 'M2')]
CH = []


def add(test, scope, value, err='', unit='', note=''):
    CH.append([test, scope, value, err, unit, note])


# ---------------------------------------------------------------- sectors and keywords
SEC = []
for prov, pat in (('SPOC', 'tess*-s00??-0000000458342769-*-s_lc.fits'), ('TESS-SPOC', 'hlsp_tess-spoc_*458342769*_lc.fits'), ('QLP', 'hlsp_qlp_*458342769*_llc.fits')):
    for fn in sorted(glob.glob(os.path.join(dldir(args.work), 'lc', pat))):
        with fits.open(fn) as h:
            h0, h1, d = h[0].header, h[1].header, h[1].data; q = d['QUALITY'] == 0; s = int(h0['SECTOR'])
            row = dict(product=prov, sector=s, cam=h0.get('CAMERA'), ccd=h0.get('CCD'), cad=float(np.nanmedian(np.diff(d['TIME']))) * 86400,
                       t0=float(np.nanmin(d['TIME'])), t1=float(np.nanmax(d['TIME'])), file=os.path.basename(fn), npix='', crowd='', flfr='', sap='', pdc='', bkg='', target='')
            if prov != 'QLP':
                ap_ = int(((np.array(h[2].data) & 2) > 0).sum()); sap = float(np.nanmedian(d['SAP_FLUX'][q])); cr, ff = h1['CROWDSAP'], h1['FLFRCSAP']
                row.update(npix=ap_, crowd=cr, flfr=ff, sap=sap, pdc=float(np.nanmedian(d['PDCSAP_FLUX'][q])), bkg=float(np.nanmedian(d['SAP_BKG'][q])), target=sap * cr / ff)
            SEC.append(row)
spoc_t = np.array([r['target'] for r in SEC if r['product'] == 'SPOC'])
add('target flux implied by the SPOC keywords', 'SPOC sectors 36-99', f'{spoc_t.min():.0f} to {spoc_t.max():.0f}', '', 'e-/s',
    f'median SAP_FLUX x CROWDSAP / FLFRCSAP; range {100 * (spoc_t.max() - spoc_t.min()) / np.mean(spoc_t):.1f} % of the mean ({np.mean(spoc_t):.0f} e-/s)')
with open(os.path.join(args.out, 'sectors_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: TESS light-curve products used, with SPOC header keywords; scripts/toi2638/checks_toi2638.py, {DATE}',
            'cadence_s = median time step; aperture_pixels = pixels with bit 2 of the SPOC aperture image; crowdsap, flfrcsap = SPOC header keywords (fraction of the flux in the '
            'aperture from the target; fraction of the target flux in the aperture); sap_e_s, pdcsap_e_s, bkg_e_s = medians of SAP_FLUX, PDCSAP_FLUX, SAP_BKG (QUALITY = 0)',
            'target_flux_e_s = sap_e_s x crowdsap / flfrcsap (the target flux implied by the SPOC crowding model); times BTJD',
            'product,sector,camera,ccd,cadence_s,start_btjd,end_btjd,aperture_pixels,crowdsap,flfrcsap,sap_e_s,pdcsap_e_s,bkg_e_s,target_flux_e_s,file'])
    for r in SEC:
        fm = lambda v, fmt: '' if v == '' else format(v, fmt)
        f.write(f"{r['product']},{r['sector']},{r['cam']},{r['ccd']},{r['cad']:.0f},{r['t0']:.3f},{r['t1']:.3f},{r['npix']},{fm(r['crowd'], '.4f')},{fm(r['flfr'], '.4f')},"
                f"{fm(r['sap'], '.0f')},{fm(r['pdc'], '.0f')},{fm(r['bkg'], '.1f')},{fm(r['target'], '.0f')},{r['file']}\n")

# ---------------------------------------------------------------- Gaia DR3 neighbours
NBF = os.path.join(args.out, 'neighbours_458342769.csv')
if not args.no_gaia:
    import pyvo
    tap = pyvo.dal.TAPService('https://gea.esac.esa.int/tap-server/tap')
    q = f"""SELECT s.source_id, s.ra, s.dec, s.pmra, s.pmdec, s.parallax, s.phot_g_mean_mag, s.bp_rp, s.phot_variable_flag, s.ruwe, s.non_single_star,
     DISTANCE(POINT(s.ra, s.dec), POINT({RA}, {DEC})) * 3600 AS sep FROM gaiadr3.gaia_source s
     WHERE 1 = CONTAINS(POINT(s.ra, s.dec), CIRCLE({RA}, {DEC}, 0.0333334)) AND s.phot_g_mean_mag < 19 ORDER BY sep"""
    t = tap.search(q).to_table(); G0 = float(t['phot_g_mean_mag'][0]); assert int(t['source_id'][0]) == GAIA_TARGET
    with open(NBF, 'w') as f:
        hdr(f, [f'TOI-2638.01 = TIC 458342769: Gaia DR3 sources within 120" with G < 19 (gaia_source, ESA archive TAP, epoch 2016.0); scripts/toi2638/checks_toi2638.py, {DATE}',
                'sep_arcsec from the target (first row); dG = G - G(target); pmra, pmdec in mas/yr; empty = not in Gaia DR3',
                'source_id,ra,dec,pmra,pmdec,parallax,G,bp_rp,dG,sep_arcsec,phot_variable_flag,ruwe,non_single_star'])
        for r in t:
            fm = lambda v, fmt: '' if np.ma.is_masked(v) or not np.isfinite(v) else format(float(v), fmt)
            f.write(f"{r['source_id']},{r['ra']:.8f},{r['dec']:.8f},{fm(r['pmra'], '.3f')},{fm(r['pmdec'], '.3f')},{fm(r['parallax'], '.4f')},{r['phot_g_mean_mag']:.4f},"
                    f"{fm(r['bp_rp'], '.3f')},{r['phot_g_mean_mag'] - G0:.3f},{r['sep']:.2f},{r['phot_variable_flag']},{fm(r['ruwe'], '.3f')},{r['non_single_star']}\n")
NB = rows('neighbours_458342769.csv'); hd = NB[0]; NB = [dict(zip(hd, r)) for r in NB[1:]]
for r in NB:
    r['dG'] = float(r['dG']); r['sep'] = float(r['sep_arcsec'])
fl42 = sum(10 ** (-0.4 * r['dG']) for r in NB[1:] if r['sep'] < 42)
add('Gaia DR3 flux within 42" relative to the target', 'G < 19, excluding the target', f'{fl42:.2f}', '', '', f'{sum(1 for r in NB[1:] if r["sep"] < 42)} sources; sum of 10^(-0.4 dG)')
nb = [r for r in NB if int(r['source_id']) == GAIA_NB17][0]
add('nearest bright neighbour', f'Gaia DR3 {GAIA_NB17}', f'{nb["sep"]:.1f}', '', 'arcsec', f'dG = {nb["dG"]:.2f}, BP-RP = {nb["bp_rp"]}, phot_variable_flag {nb["phot_variable_flag"]}')
for r in NB[1:]:
    if r['phot_variable_flag'] == 'VARIABLE':
        add('Gaia DR3 VARIABLE flag within 120"', f'Gaia DR3 {r["source_id"]}', f'{r["sep"]:.1f}', '', 'arcsec', f'dG = {r["dG"]:.2f}; flux {100 * 10 ** (-0.4 * r["dG"]):.2f} % of the target')

# ---------------------------------------------------------------- odd/even
pdc = load(args.work, 'SPOC'); W = windows(pdc, P, T0)
camp = {'2021': [36, 37], '2023': [63, 64], '2025': [90], '2026': [99]}
for c, ss in camp.items():
    res = {}
    for par, lab in ((0, 'even'), (1, 'odd')):
        x, y, e = [], [], []
        for s in ss:
            m = (W[s]['ep'] % 2) == par; x.append(Model(W[s]['t'][m], W[s]['cad'], P, T0).depth(k0, b0, a0) / D0); y.append(1 - W[s]['f'][m]); e.append(W[s]['e'][m])
        x, y, e = map(np.concatenate, (x, y, e)); w = 1 / e ** 2; A = np.sum(w * x * y) / np.sum(w * x * x); eA = 1 / np.sqrt(np.sum(w * x * x)); res[lab] = (A * 1e6, eA * 1e6)
    dz = (res['odd'][0] - res['even'][0]) / np.hypot(res['odd'][1], res['even'][1])
    add('odd minus even depth', f'campaign {c} (sectors {"+".join(map(str, ss))})', f"{res['odd'][0] - res['even'][0]:.0f}", f"{np.hypot(res['odd'][1], res['even'][1]):.0f}", 'ppm',
        f"odd {res['odd'][0]:.0f} +- {res['odd'][1]:.0f}, even {res['even'][0]:.0f} +- {res['even'][1]:.0f} ppm; {dz:+.1f} sigma; white errors")

# ---------------------------------------------------------------- secondary eclipse and phase curve
allt, allr = [], []
for s, x in pdc.items():
    t, f_ = x['t'], x['f']; intr = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) < 0.07
    cen, med = [], []
    for a in np.arange(t.min(), t.max() + 0.5, 0.5):
        m = (t >= a) & (t < a + 0.5) & ~intr
        if m.sum() > 50:
            cen.append(np.median(t[m])); med.append(np.median(f_[m]))
    r = f_ / np.interp(t, cen, med) - 1; allt.append(t[~intr]); allr.append(r[~intr])
t, r = np.concatenate(allt), np.concatenate(allr); ph = ((t - T0) / P) % 1; sig = np.std(r); scan = []
for dur in (0.040, 0.056, 0.075):
    for p0 in np.arange(0.15, 0.85 + 1e-9, 0.0025):
        m = np.abs(ph - p0) * P < dur / 2; nbm = ~m & (np.abs(ph - p0) * P < dur * 2)
        if m.sum() < 50:
            continue
        dep = np.mean(r[nbm]) - np.mean(r[m]); e = sig / np.sqrt(m.sum()); scan.append((p0, dur, dep * 1e6, e * 1e6, dep / e))
scan = np.array(scan)
for dur in (0.040, 0.056, 0.075):
    z = scan[(np.abs(scan[:, 0] - 0.5) < 1e-6) & (scan[:, 1] == dur)][0]
    add('secondary eclipse depth at phase 0.5', f'box {dur:.3f} d, SPOC PDCSAP all sectors', f'{z[2]:.0f}', f'{z[3]:.0f}', 'ppm', 'white error; transits removed; 0.5-d running-median detrending')
j = np.argmax(scan[:, 4])
add('largest box S/N over phases 0.15-0.85', 'SPOC PDCSAP all sectors', f'{scan[j, 4]:.2f}', '', '', f'phase {scan[j, 0]:.4f}, box {scan[j, 1]:.3f} d, depth {scan[j, 2]:.0f} ppm; {len(scan)} trials; '
    f'99.9th percentile of |S/N| over the scan {np.percentile(np.abs(scan[:, 4]), 99.9):.2f}; rms of S/N {np.std(scan[:, 4]):.2f}')
PC = []
for s, x in sorted(pdc.items()):
    t, f_ = x['t'], x['f'] - 1; m = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) > 0.084; t, f_ = t[m], f_[m]; phi = 2 * np.pi * (t - T0) / P
    X = np.c_[np.ones_like(t), t - t.mean(), (t - t.mean()) ** 2, -np.cos(2 * phi), -np.cos(phi), np.sin(phi)]
    c, *_ = np.linalg.lstsq(X, f_, rcond=None); rr = f_ - X @ c; nb_ = 60; mm = len(rr) // nb_
    beta = max(1.0, np.std(rr[:mm * nb_].reshape(mm, nb_).mean(1)) / (np.std(rr) / np.sqrt(nb_))); cov = np.linalg.inv(X.T @ X) * np.var(rr) * beta ** 2
    PC.append((s, c[3] * 1e6, np.sqrt(cov[3, 3]) * 1e6, c[4] * 1e6, np.sqrt(cov[4, 4]) * 1e6, c[5] * 1e6, np.sqrt(cov[5, 5]) * 1e6, beta))
    add('phase curve', f'sector {s}', f'{c[3] * 1e6:.1f}', f'{np.sqrt(cov[3, 3]) * 1e6:.1f}', 'ppm', f'ellipsoidal amplitude; reflection {c[4] * 1e6:.1f} +- {np.sqrt(cov[4, 4]) * 1e6:.1f}, '
        f'beaming {c[5] * 1e6:.1f} +- {np.sqrt(cov[5, 5]) * 1e6:.1f} ppm; red-noise factor {beta:.2f}')
PCa = np.array(PC)
for col, lab in ((1, 'ellipsoidal'), (3, 'reflection')):
    w = 1 / PCa[:, col + 1] ** 2; mu = np.sum(PCa[:, col] * w) / w.sum(); emu = 1 / np.sqrt(w.sum())
    add(f'phase curve, weighted mean', f'{lab}, SPOC sectors 36-99', f'{mu:.1f}', f'{emu:.1f}', 'ppm', f'sector-to-sector standard deviation {np.std(PCa[:, col]):.0f} ppm; PDCSAP (PDC may suppress a 3.93-d signal)')
G, Ms, Rs = 6.674e-11, 1.989e30, 6.957e8
MFL, RFL = 2.0891666, 1.6712719                      # Gaia DR3 FLAME
a_m = (G * MFL * Ms * (P * 86400) ** 2 / (4 * np.pi ** 2)) ** (1 / 3); aR_dyn = a_m / (RFL * Rs)
add('a/R* from the Gaia DR3 FLAME mass and radius', f'M = {MFL:.2f} Msun, R = {RFL:.2f} Rsun', f'{aR_dyn:.2f}', '', '', f'a = {a_m / 1.496e11:.4f} au; fitted a/R*: null {a0:.2f}, M2 {M2[1]:.2f} (PDCSAP + QLP 10), {M2sap[1]:.2f} (SAP)')
Alim = 3 * max(np.std(PCa[:, 1], ddof=1) / np.sqrt(len(PCa)), 1 / np.sqrt(np.sum(1 / PCa[:, 2] ** 2)))
add('companion mass limit from the ellipsoidal amplitude', '3 x error of the mean (sector scatter / sqrt(6) or formal, the larger)', f'{Alim / 1e6 / (1.4 * aR_dyn ** -3) * MFL:.3f}', '', 'Msun',
    f'A_ell < {Alim:.0f} ppm; M2 / M1 = A / (alpha (R/a)^3 sin^2 i), alpha = 1.4, a/R* = {aR_dyn:.2f}, M1 = {MFL:.2f} Msun; order of magnitude; PDCSAP')

# ---------------------------------------------------------------- pulsations and the 7.31 c/d signal
for s, x in sorted(pdc.items()):
    t, f_ = x['t'], x['f'] - 1; m = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) > 0.084; t, f_ = t[m], f_[m]
    fr = np.linspace(0.05, 360, 400000); pw = LombScargle(t, f_, normalization='psd').power(fr); amp = np.sqrt(4 * pw / len(t)) * 1e6
    idx = np.argsort(amp)[::-1]; peaks = []
    for i in idx:
        if all(abs(fr[i] - fr[j_]) > 0.1 for j_ in peaks):
            peaks.append(i)
        if len(peaks) == 5:
            break
    best = None
    for f0 in np.linspace(7.301, 7.321, 201):
        X = np.c_[np.ones(len(t)), np.sin(2 * np.pi * f0 * t), np.cos(2 * np.pi * f0 * t)]; c, *_ = np.linalg.lstsq(X, f_, rcond=None); a_ = np.hypot(c[1], c[2])
        if best is None or a_ > best[0]:
            best = (a_, f0)
    hi1 = fr > 1.0; j1 = np.argmax(np.where(hi1, amp, 0))
    add('amplitude spectrum, out of transit', f'SPOC sector {s}', f'{amp[j1]:.0f}', '', 'ppm',
        f'highest peak above 1 c/d at {fr[j1]:.3f} c/d; five highest peaks overall (c/d: ppm) ' + '; '.join(f'{fr[i]:.3f}: {amp[i]:.0f}' for i in peaks) + f'; median amplitude 0.5-20 c/d {np.median(amp[(fr > 0.5) & (fr < 20)]):.0f} ppm, '
        f'above 20 c/d {np.median(amp[fr > 20]):.0f} ppm; amplitude at {best[1]:.4f} c/d (7.311 +- 0.01) {best[0] * 1e6:.0f} ppm')
for fn in sorted(glob.glob(os.path.join(dldir(args.work), 'nb17', 'hlsp_qlp_*.fits')), key=lambda x: int(x.split('_s00')[1][:2])):
    h = fits.open(fn); d = h[1].data; s = h[0].header['SECTOR']; col = 'KSPSAP_FLUX' if 'KSPSAP_FLUX' in d.names else 'DET_FLUX'
    t = np.array(d['TIME'], float); f_ = np.array(d[col], float); ok = np.isfinite(t) & np.isfinite(f_) & (d['QUALITY'] == 0); t, f_ = t[ok], f_[ok] / np.nanmedian(f_[ok])
    fr = np.linspace(1.0, 0.5 / np.median(np.diff(t)) * 0.98, 60000); pw = LombScargle(t, f_ - 1).power(fr, normalization='psd'); amp = np.sqrt(4 * pw / len(t)) * 1e6; j = np.argmax(amp)
    phs = ((t - T0) / P + 0.5) % 1 - 0.5; intr = np.abs(phs * P) < 0.03; near = (np.abs(phs * P) > 0.05) & (np.abs(phs * P) < 0.25)
    dep = (np.mean(f_[near]) - np.mean(f_[intr])) * 1e6; edep = np.std(f_[near]) / np.sqrt(intr.sum()) * 1e6
    k7 = np.abs(fr - 7.311) < 0.02; j7 = np.where(k7)[0][np.argmax(amp[k7])]
    add('neighbour TIC 458342741, QLP: strongest peak', f'sector {s}', f'{fr[j]:.3f}', '', 'c/d', f'amplitude {amp[j]:.0f} ppm (median amplitude {np.median(amp):.0f} ppm); '
        f'highest peak within 7.311 +- 0.02 c/d: {fr[j7]:.3f} c/d, {amp[j7]:.0f} ppm; '
        f'mean flux within +-0.03 d of the TOI-2638.01 transits lower than at 0.05-0.25 d by {dep:.0f} +- {edep:.0f} ppm')

# ---------------------------------------------------------------- J2 precession estimate
VSINI, EVSINI = 109.730064, 6.545823                 # Gaia DR3 ESP-HS vsini (km/s)
MJ = 1.898e27; yr_s = 3.15576e7; n = 2 * np.pi / (P * 86400); R = RFL * Rs; M = MFL * Ms
vcrit = np.sqrt(G * M / R) / 1e3
obs = abs(M2[3]); eobs = J.get(('SPOC PDCSAP + QLP sector 10', 'M2_err'), None)
PREC = []
for sinis in (1.0, 0.8, 0.6, 0.5, 0.4, 0.3):
    veq = VSINI / sinis; q = (veq / vcrit) ** 2
    if q > 0.45:
        continue
    Om = veq * 1e3 / R
    for kap in (0.003, 0.005, 0.007):
        J2 = 2 * kap / 3 * q; S_ = 0.06 * M * R ** 2 * Om
        for mp in (1, 10, 50):
            L_ = mp * MJ * np.sqrt(G * M * a_m); Od0 = 1.5 * n * J2 / aR_dyn ** 2
            dbmax = aR_dyn * sinis * 0.5 * Od0 * yr_s
            psi = np.radians(np.linspace(0.1, 89.9, 1799)); ok = aR_dyn * sinis * np.sin(psi) * np.cos(psi) * Od0 * yr_s >= obs
            rng_ = f'{np.degrees(psi[ok]).min():.0f}-{np.degrees(psi[ok]).max():.0f}' if ok.any() else 'none'
            PREC.append([sinis, veq, q, kap, J2, mp, L_ / S_, np.degrees(Od0 * yr_s), 2 * np.pi / (Od0 * np.cos(np.pi / 4) * np.sqrt(1 + (L_ / S_) ** 2 + 2 * (L_ / S_) * np.cos(np.pi / 4))) / yr_s, dbmax, rng_])
with open(os.path.join(args.out, 'precession_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: nodal precession of the orbit by the stellar oblateness (J2), order-of-magnitude grid; scripts/toi2638/checks_toi2638.py, {DATE}',
            f'star: Gaia DR3 FLAME M = {MFL:.2f} Msun, R = {RFL:.2f} Rsun; ESP-HS vsini = {VSINI:.0f} +- {EVSINI:.0f} km/s; v_crit = sqrt(GM/R) = {vcrit:.0f} km/s; a/R* = {aR_dyn:.2f} (Kepler, P = {P:.5f} d)',
            'q = (v_eq/v_crit)^2, v_eq = vsini / sin_istar (rows with q > 0.45 omitted); J2 = (2 k_ap / 3) q, k_ap = apsidal-motion constant (assumed 0.003-0.007); '
            'L_over_S = orbital angular momentum of a companion of mass m_c / stellar spin angular momentum (moment of inertia 0.06 M R^2, assumed)',
            'omega0_deg_yr = omega0 = (3/2) n J2 (R/a)^2; mutual precession frequency = omega0 cos(psi) sqrt(1 + (L/S)^2 + 2 (L/S) cos(psi)); P_prec_yr_psi45 = mutual precession period at psi = 45 deg',
            'dbdt_max_yr = (a/R) sin_istar (1/2) omega0 = maximum over psi and node of |db/dt| <= (a/R) sin(i_star) sin(psi) |omega0 cos(psi)|; the rate envelopes do not impose the observed transiting geometry',
            f'psi_range_deg = obliquities at which the maximum |db/dt| reaches the M2 value {obs:.4f} per yr (SPOC PDCSAP + QLP sector 10, jointfit_458342769.csv); none = not reached',
            'sin_istar,v_eq_km_s,q,k_ap,J2,m_c_MJ,L_over_S,omega0_deg_yr,P_prec_yr_psi45,dbdt_max_yr,psi_range_deg'])
    for r in PREC:
        f.write(f'{r[0]:.2f},{r[1]:.0f},{r[2]:.4f},{r[3]:.3f},{r[4]:.2e},{r[5]},{r[6]:.3f},{r[7]:.3f},{r[8]:.0f},{r[9]:.4f},{r[10]}\n')
print('J2 grid:', len(PREC), 'rows; observed |db/dt|', obs)

# ---------------------------------------------------------------- future TESS sectors
try:
    import tess_stars2px as TP
    res = TP.tess_stars2px_function_entry(TIC, RA, DEC); sc = TP.TESS_Spacecraft_Pointing_Data(); mid = {int(s): m for s, m in zip(sc.sectors, sc.midtimes)}
    from astropy.time import Time
    for s, cam, ccd in zip(res[3], res[4], res[5]):
        if s > 99:
            m = mid[int(s)] - 2457000; nt = np.arange(np.ceil((m - 13.7 - T0) / P), np.floor((m + 13.7 - T0) / P) + 1)
            add('future TESS sector (tess-point)', f'sector {s}', f'camera {cam} CCD {ccd}', '', '', f'nominal window {Time(m - 13.7 + 2457000, format="jd").iso[:10]} to '
                f'{Time(m + 13.7 + 2457000, format="jd").iso[:10]} (midtime +- 13.7 d); {len(nt)} transits of the linear ephemeris in the window')
except Exception as e:
    add('future TESS sector (tess-point)', 'not computed', '', '', '', repr(e))

with open(os.path.join(args.out, 'checks_458342769.csv'), 'w') as f:
    hdr(f, [f'TOI-2638.01 = TIC 458342769: checks of the transit signal and of its change; scripts/toi2638/checks_toi2638.py, {DATE}',
            f'ephemeris of ephemeris_458342769.csv; common shape of jointfit_458342769.csv (null, SPOC PDCSAP + QLP sector 10): k = {k0:.4f}, b = {b0:.3f}, a/R* = {a0:.2f}; methods in the script header',
            'test,scope,value,err,unit,note'])
    w = csv.writer(f, lineterminator='\n')
    for r in CH:
        w.writerow([str(v) for v in r])
for r in CH:
    print(r)
