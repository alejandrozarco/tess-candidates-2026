"""Spectral-energy-distribution fit of TIC 89247735: Teff, radius and luminosity from catalogue photometry.

usage: python scripts/sed_fit.py [--work ./work] [--out data]
input:  the catalogue photometry listed in PHOT below (Gaia DR3, SkyMapper DR4 PSF, APASS DR9, 2MASS, AllWISE; VizieR
        values at the Gaia DR3 position), the Gaia DR3 parallax; Castelli & Kurucz (2003) ODFNEW model spectra
        (Teff 3750-5000 K in 250 K steps, log g 4.5 and 5.0, [M/H] 0.0, +0.2, +0.5) and the filter curves and zero
        points from the SVO Theoretical Spectra and Filter Profile services, downloaded once into <work>/89247735/sed/
        (bounded timeouts, three attempts per file)
output: printed fits, <work>/89247735/sed.json and <out>/sed_89247735.csv (photometry, model magnitudes, residuals)
Model: log-linear interpolation of the model spectra in Teff (10 K grid), photon-counting synthetic magnitudes with the
SVO transmission curves and zero points (AB for SkyMapper), a free scale (R/d)^2 fitted as a weighted magnitude
offset, A_V = 0 (Cardelli, Clayton & Mathis 1989 law with R_V = 3.1 for
the A_V checks), error floor 0.02 mag; d = 1/parallax.
R = sqrt(scale) d, L = 4 pi d^2 sigma Teff^4 scale. Fits: all bands; without SkyMapper u and v (model ultraviolet
unreliable for K dwarfs; the reference fit); [M/H] and log g dependence; A_V = 0.05, 0.10, 0.20; Monte Carlo over
the photometric and parallax errors (60 draws). WISE W4 is an upper limit and is not fitted. The adopted parameters
of README section 3.1 combine this fit with the other Teff estimates listed there.
"""
import argparse, csv, json, os, re, time, urllib.request
import numpy as np

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
p.add_argument('--date', default='')
args = p.parse_args(); TIC = 89247735; SD = os.path.join(args.work, str(TIC), 'sed'); os.makedirs(SD, exist_ok=True)
DATE = args.date or __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
# (SVO filter id, magnitude, error, system, catalogue); W4 is an upper limit (no error)
PHOT = [('GAIA/GAIA3.G', 11.8509, 0.0028, 'Vega', 'Gaia DR3'), ('GAIA/GAIA3.Gbp', 12.5209, 0.0030, 'Vega', 'Gaia DR3'),
        ('GAIA/GAIA3.Grp', 11.0690, 0.0039, 'Vega', 'Gaia DR3'),
        ('SkyMapper/SkyMapper.u', 15.3792, 0.0142, 'AB', 'SkyMapper DR4 PSF'), ('SkyMapper/SkyMapper.v', 15.0495, 0.0114, 'AB', 'SkyMapper DR4 PSF'),
        ('SkyMapper/SkyMapper.g', 12.5626, 0.0107, 'AB', 'SkyMapper DR4 PSF'), ('SkyMapper/SkyMapper.r', 11.8553, 0.0111, 'AB', 'SkyMapper DR4 PSF'),
        ('SkyMapper/SkyMapper.i', 11.4054, 0.0097, 'AB', 'SkyMapper DR4 PSF'), ('SkyMapper/SkyMapper.z', 11.2245, 0.0120, 'AB', 'SkyMapper DR4 PSF'),
        ('Generic/Johnson.B', 13.451, 0.029, 'Vega', 'APASS DR9'), ('Generic/Johnson.V', 12.225, 0.018, 'Vega', 'APASS DR9'),
        ('2MASS/2MASS.J', 10.154, 0.024, 'Vega', '2MASS'), ('2MASS/2MASS.H', 9.606, 0.026, 'Vega', '2MASS'), ('2MASS/2MASS.Ks', 9.444, 0.019, 'Vega', '2MASS'),
        ('WISE/WISE.W1', 9.352, 0.024, 'Vega', 'AllWISE'), ('WISE/WISE.W2', 9.441, 0.020, 'Vega', 'AllWISE'), ('WISE/WISE.W3', 9.335, 0.037, 'Vega', 'AllWISE'),
        ('WISE/WISE.W4', 8.480, None, 'Vega', 'AllWISE (upper limit)')]
FITTED = [r[0] for r in PHOT if r[2] is not None]; NOUV = [f for f in FITTED if not f.endswith(('.u', '.v'))]
PLX, EPLX = 9.5397e-3, 0.0202e-3   # Gaia DR3 parallax (arcsec)
TEFFS = [3750, 4000, 4250, 4500, 4750, 5000]; LOGGS = [4.5, 5.0]; MHS = [0.0, 0.2, 0.5]


def get(url, fn):
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        return open(fn).read()
    for i in range(3):
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
            if '<TR>' in txt or '<TABLE' in txt:
                open(fn, 'w').write(txt); return txt
        except Exception as e:
            print('retry', url[:90], str(e)[:80]); time.sleep(5)
    raise RuntimeError(f'download failed: {url}')


def votab(txt):
    rows = re.findall(r'<TR>(.*?)</TR>', txt, re.S)
    return np.array([[float(v) for v in re.findall(r'<TD>(.*?)</TD>', r)] for r in rows])


lst = get('http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=Kurucz2003&REQUEST=queryData&teff_min=3750&teff_max=5000&logg_min=4.5&logg_max=5.0',
          os.path.join(SD, 'kurucz_list.xml'))
MOD = {}
for row in re.findall(r'<TR>(.*?)</TR>', lst, re.S):
    td = re.findall(r'<TD>(.*?)</TD>', row, re.S); te, lg, mh = int(float(td[0])), float(td[1]), float(td[3])
    if te in TEFFS and lg in LOGGS and mh in MHS:
        url = re.search(r'CDATA\[(.*?)\]', row).group(1)
        MOD[(te, lg, mh)] = votab(get(url, os.path.join(SD, f'kurucz_{te}_{lg:.1f}_{mh:+.1f}.xml')))   # wavelength (A), F_lambda at the surface
assert len(MOD) == len(TEFFS) * len(LOGGS) * len(MHS), f'{len(MOD)} model spectra'
FIL = {}
for f, m, e, sy, cat in PHOT:
    txt = get(f'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={f}', os.path.join(SD, 'filt_' + f.replace('/', '_') + '.xml')); a = votab(txt)
    zp = float(re.search(r'name="ZeroPoint"[^>]*value="([^"]+)"', txt).group(1))
    FIL[f] = (a[:, 0], a[:, 1], 3631.0 if sy == 'AB' else zp)
C = 2.99792458e18; PC = 3.0857e18; RSUN = 6.957e10; SIG = 5.6704e-5; LSUN = 3.828e33; d_pc = 1 / PLX


def synmag(lam, flam, f):
    wl, T, zp = FIL[f]; Fnu = np.interp(wl, lam, flam) * wl ** 2 / C
    return -2.5 * np.log10(np.trapezoid(Fnu * T / wl, wl) / np.trapezoid(T / wl, wl) / (zp * 1e-23))


def ccm(lam, Av, Rv=3.1):
    x = np.minimum(1e4 / lam, 3.3); y = x - 1.82   # optical/IR form of the law; held constant below 3030 A
    a = np.where(x < 1.1, 0.574 * x ** 1.61, 1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7)
    b = np.where(x < 1.1, -0.527 * x ** 1.61, 1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7)
    return Av * (a + b / Rv)


def model_flux(T, lg, mh):
    lo = max(t for t in TEFFS if t <= T); hi = min(t for t in TEFFS if t >= T); a = MOD[(lo, lg, mh)]; lam = a[:, 0]
    if hi == lo:
        return lam, a[:, 1]
    b = MOD[(hi, lg, mh)]; w = (T - lo) / (hi - lo)
    return lam, np.exp((1 - w) * np.log(np.clip(a[:, 1], 1e-30, None)) + w * np.log(np.clip(np.interp(lam, b[:, 0], b[:, 1]), 1e-30, None)))


def fit(use, floor=0.02, Av=0.0, mags=None):
    mags = mags or {f: m for f, m, e, sy, cat in PHOT}; err = {f: e for f, m, e, sy, cat in PHOT}
    ob = np.array([mags[f] for f in use]); w = 1 / np.array([max(err[f], floor) for f in use]) ** 2; best = None; grid = []
    for lg in LOGGS:
        for mh in MHS:
            for T in np.arange(3800, 4951, 10):
                lam, fl = model_flux(T, lg, mh)
                if Av:
                    fl = fl * 10 ** (-0.4 * ccm(lam, Av))
                sm = np.array([synmag(lam, fl, f) for f in use]); off = np.sum(w * (ob - sm)) / np.sum(w); chi2 = np.sum(w * (ob - sm - off) ** 2)
                grid.append((chi2, T, lg, mh, off))
                if best is None or chi2 < best[0]:
                    best = (chi2, T, lg, mh, off)
    return best, grid


def derived(T, off, d=d_pc):
    scale = 10 ** (-0.4 * off); return np.sqrt(scale) * d * PC / RSUN, 4 * np.pi * (d * PC) ** 2 * SIG * T ** 4 * scale / LSUN


out = {}
for lab_, use in (('all_bands', FITTED), ('no_uv', NOUV)):
    (chi2, T, lg, mh, off), grid = fit(use); R, L = derived(T, off); nd = len(use) - 3
    g = np.array([x[0] for x in grid if x[2] == lg and x[3] == mh]); Tg = np.array([x[1] for x in grid if x[2] == lg and x[3] == mh])
    ok = g <= chi2 + max(1, chi2 / nd)
    out[lab_] = dict(Teff=float(T), Teff_1sig_range=[float(Tg[ok].min()), float(Tg[ok].max())], logg=lg, mh=mh, chi2=float(chi2), ndof=nd, R=float(R), L=float(L),
                     off=float(off), per_grid={})
    print(f'{lab_}: Teff {T:.0f} K (1-sigma, chi2 scaled: {Tg[ok].min():.0f}-{Tg[ok].max():.0f}), log g {lg}, [M/H] {mh:+.1f}, R {R:.3f} Rsun, L {L:.4f} Lsun, chi2 {chi2:.1f} / {nd}')
    for lg2 in LOGGS:
        for mh2 in MHS:
            bb = min(x for x in grid if x[2] == lg2 and x[3] == mh2); R2, L2 = derived(bb[1], bb[4])
            out[lab_]['per_grid'][f'logg{lg2}_mh{mh2:+.1f}'] = dict(Teff=float(bb[1]), chi2=float(bb[0]), R=float(R2), L=float(L2))
            print(f'    log g {lg2} [M/H] {mh2:+.1f}: Teff {bb[1]:.0f} K, chi2 {bb[0]:.1f}, R {R2:.3f}, L {L2:.4f}')
for Av in (0.05, 0.10, 0.20):
    (chi2, T, lg, mh, off), _ = fit(NOUV, Av=Av); R, L = derived(T, off)
    out[f'no_uv_Av{Av:.2f}'] = dict(Teff=float(T), chi2=float(chi2), R=float(R), L=float(L)); print(f'A_V = {Av:.2f}: Teff {T:.0f} K, chi2 {chi2:.1f}, R {R:.3f}, L {L:.4f}')
rng = np.random.default_rng(3); mc = []
for i in range(60):
    mags = {f: m + (rng.normal(0, max(e, 0.02)) if e is not None else 0) for f, m, e, sy, cat in PHOT}
    b, _ = fit(NOUV, mags=mags); R, L = derived(b[1], b[4], 1 / (PLX + rng.normal(0, EPLX))); mc.append((b[1], R, L))
mc = np.array(mc); out['monte_carlo'] = {k: [float(np.mean(mc[:, i])), float(np.std(mc[:, i]))] for i, k in enumerate(('Teff', 'R', 'L'))}
print('Monte Carlo (60 draws, no u, v):', out['monte_carlo'])
json.dump(out, open(os.path.join(args.work, str(TIC), 'sed.json'), 'w'), indent=1)
ref = out['no_uv']; lam, fl = model_flux(ref['Teff'], ref['logg'], ref['mh'])
os.makedirs(args.out, exist_ok=True); fn = os.path.join(args.out, f'sed_{TIC}.csv')
with open(fn, 'w', newline='') as fo:
    a = out['all_bands']; mcv = out['monte_carlo']
    fo.write(f'# TIC {TIC} = Gaia DR3 6708071581514074368: SED fit (scripts/sed_fit.py), {DATE}\n'
             f'# model: Castelli & Kurucz (2003) ODFNEW spectra (SVO), photon-counting synthetic photometry, A_V = 0, error floor 0.02 mag, d = 1/parallax = {d_pc:.2f} pc\n'
             f'# reference fit (no SkyMapper u, v): Teff {ref["Teff"]:.0f} K (1-sigma {ref["Teff_1sig_range"][0]:.0f}-{ref["Teff_1sig_range"][1]:.0f}), '
             f'log g {ref["logg"]}, [M/H] {ref["mh"]:+.1f}, R {ref["R"]:.3f} Rsun, L {ref["L"]:.4f} Lsun, chi2 {ref["chi2"]:.1f} for {ref["ndof"]} dof\n'
             f'# all bands: Teff {a["Teff"]:.0f} K, R {a["R"]:.3f} Rsun, L {a["L"]:.4f} Lsun, chi2 {a["chi2"]:.1f} for {a["ndof"]} dof; '
             f'Monte Carlo (no u, v): Teff {mcv["Teff"][0]:.0f} +- {mcv["Teff"][1]:.0f} K, R {mcv["R"][0]:.3f} +- {mcv["R"][1]:.3f}, L {mcv["L"][0]:.4f} +- {mcv["L"][1]:.4f} (statistical only)\n'
             '# A_V checks (no u, v): ' + '; '.join(f'A_V {k[-4:]}: Teff {v["Teff"]:.0f} K, R {v["R"]:.3f}, L {v["L"]:.4f}, chi2 {v["chi2"]:.1f}'
                                                  for k, v in out.items() if k.startswith('no_uv_Av')) + '\n'
             '# The adopted parameters (README 3.1: Teff 4430 +- 100 K, L 0.190 +- 0.008 Lsun, R 0.74 +- 0.04 Rsun) combine this fit with other Teff estimates.\n'
             '# columns: band (SVO filter id), catalogue, mag, mag_err, system, fitted (reference fit), lambda_mean_A, model_mag (reference fit, scaled), '
             'residual (observed - model)\n')
    w = csv.writer(fo); w.writerow(['band', 'catalogue', 'mag', 'mag_err', 'system', 'fitted', 'lambda_mean_A', 'model_mag', 'residual'])
    for f, m, e, sy, cat in PHOT:
        wl, T, zp = FIL[f]; mm = synmag(lam, fl, f) + ref['off']
        w.writerow([f, cat, f'{m:.4f}', '' if e is None else f'{e:.4f}', sy, 'yes' if f in NOUV else 'no', f'{np.trapezoid(T * wl, wl) / np.trapezoid(T, wl):.0f}',
                    f'{mm:.4f}', f'{m - mm:+.4f}'])
print('written', fn)
