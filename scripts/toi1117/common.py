"""Shared definitions for the TOI-1117 c (TIC 295541511) scripts (README section 6.10).

Light-curve products (downloaded by download_toi1117.py into <work>/295541511/dl/lc/), labels PROV_SECTOR_KIND:
  SPOC      SPOC 2-min light curves, sectors 39, 100-103 (kind pdc = PDCSAP_FLUX, sap = SAP_FLUX)
  FAST      SPOC 20-s light curve, sector 39
  TSPOC     TESS-SPOC full-frame-image light curves, sectors 13 (1800 s) and 39 (600 s)
  QLP       QLP full-frame-image light curves, sectors 13, 39, 100-103 (ksp = KSPSAP_FLUX or DET_FLUX, sap = SAP_FLUX)
  ELE       eleanor-lite light curve, sector 13 (corr = CORR_FLUX, pca = PCA_FLUX)
  TARS      TARS light curve, sectors 13 and 39 (FLUX; no quality column)
Cadences with QUALITY = 0 (SPOC and TESS-SPOC: also finite PDCSAP_FLUX, so that SAP and PDCSAP use the same cadences),
finite positive flux, divided by the sector median. Time: BTJD = BJD_TDB - 2457000.

Planet b is masked in every light curve (+-(T14_b/2 + 1 h) around its transits). Gaps longer than 0.15 d are trimmed by
0.25 d on both sides. Detrending: wotan biweight, window 1.0 d, with the windows of c masked during the fit unless stated.

Transit model: batman (Kreidberg 2015), circular orbit, quadratic limb darkening, averaged over each exposure by sampling
it at the centres of NSUB(cadence) equal parts (midpoint sampling; batman's supersample_factor is not used).
"""
import glob, os, sys, types
import numpy as np

TIC = 295541511
GAIA = 6436995652638923392
RA, DEC = 273.60205521267, -66.41997959733            # Gaia DR3 (VizieR I/355/gaiadr3), epoch 2016.0
SECTORS = [13, 39, 100, 101, 102, 103]
SPOC2 = [39, 100, 101, 102, 103]                         # sectors with SPOC 2-min data
# planet b: linear ephemeris refitted on these TESS data in the search that selected the signal (Lockley et al. 2025:
# P = 2.22816 +- 0.00001 d, T0 = BTJD 2386.964 +- 0.004); T14_b = 2.015 h
PB, TB, DB = 2.2281822816, 2917.2777438537, 2.015 / 24
# planet c: starting ephemeris of the search (box fit); the adopted ephemeris is the MCMC fit (data/fit_295541511.csv)
PC, TC = 4.578660, 1670.168
TBOX, DBOX = TC - 0.0025, 1.625 / 24                     # box mid-time and duration used by the box measurements
PD = 8.665                                               # planet d (RV only; Lockley et al. 2025)
# host (Lockley et al. 2025, joint-fit posterior): M = 0.97 +- 0.02 Msun, R = 1.06 +- 0.03 Rsun
MS, MS_E, RS, RS_E = 0.97, 0.02, 1.06, 0.03
RHO_SUN = 1.40999                                        # g/cm^3
DATE = '2026-10-23'


def import_batman():
    """batman-package 2.5.x imports distutils, which Python 3.12 no longer ships; provide a stub."""
    if 'distutils' not in sys.modules:
        d = types.ModuleType('distutils'); c = types.ModuleType('distutils.ccompiler')
        c.new_compiler = lambda *a, **k: None; d.ccompiler = c
        sys.modules['distutils'] = d; sys.modules['distutils.ccompiler'] = c
    import batman
    return batman


def nsub(cad_d):
    """number of sub-samples per exposure (midpoint sampling)"""
    return max(3, int(np.ceil(cad_d * 86400 / 40.0)))


def subgrid(t, cad_d):
    """times at the centres of n equal parts of each exposure; returns (expanded times, n)"""
    n = nsub(cad_d); off = ((np.arange(n) + 0.5) / n - 0.5) * cad_d
    return (np.asarray(t)[:, None] + off[None, :]).ravel(), n


def workdir(work):
    d = os.path.join(work, str(TIC)); os.makedirs(d, exist_ok=True); return d


def dldir(work, sub=''):
    d = os.path.join(work, str(TIC), 'dl', sub); os.makedirs(d, exist_ok=True); return d


def mask_b(t, extra=1.0 / 24):
    n = np.round((t - TB) / PB); return np.abs(t - (TB + n * PB)) > DB / 2 + extra


def robust_sd(v):
    v = np.asarray(v); v = v[np.isfinite(v)]; return 1.4826 * np.median(np.abs(v - np.median(v)))


def phase(t, P, T0):
    return (t - T0 + 0.5 * P) % P - 0.5 * P


def load(work, label):
    """label PROV_SECTOR_KIND (see the module docstring); returns dict(t, f, bkg, c1, c2, cad, label, sec, file) or None"""
    from astropy.io import fits
    prov, sec, *kind = label.split('_'); sec = int(sec); kind = kind[0] if kind else ''
    D = dldir(work, 'lc')
    pat = {'SPOC': f'tess*-s{sec:04d}-0000000295541511-*-s_lc.fits', 'FAST': f'tess*-s{sec:04d}-0000000295541511-*-a_fast-lc.fits',
           'TSPOC': f'hlsp_tess-spoc*s{sec:04d}_tess_v1_lc.fits', 'QLP': f'hlsp_qlp*s{sec:04d}-*llc.fits',
           'ELE': f'hlsp_gsfc-eleanor*s{sec:04d}*.fits', 'TARS': f'hlsp_tars*s{sec:04d}*lc.fits'}[prov]
    fp = sorted(glob.glob(os.path.join(D, pat)))
    if not fp:
        return None
    with fits.open(fp[0]) as h:
        d = h[1].data; t = np.array(d['TIME'], float); n = len(t)
        bkg = np.full(n, np.nan); c1 = np.full(n, np.nan); c2 = np.full(n, np.nan); q = np.zeros(n, int)
        if prov in ('SPOC', 'FAST', 'TSPOC'):
            col = 'PDCSAP_FLUX' if kind == 'pdc' else 'SAP_FLUX'; f = np.array(d[col], float); q = np.array(d['QUALITY'])
            bkg = np.array(d['SAP_BKG'], float); c1 = np.array(d['MOM_CENTR1'], float); c2 = np.array(d['MOM_CENTR2'], float)
            good = (q == 0) & np.isfinite(np.array(d['PDCSAP_FLUX'], float))
        elif prov == 'QLP':
            col = 'SAP_FLUX' if kind == 'sap' else ('KSPSAP_FLUX' if 'KSPSAP_FLUX' in d.columns.names else 'DET_FLUX')
            f = np.array(d[col], float); q = np.array(d['QUALITY']); bkg = np.array(d['SAP_BKG'], float)
            c1 = np.array(d['SAP_X'], float); c2 = np.array(d['SAP_Y'], float); good = q == 0
        elif prov == 'ELE':
            col = 'PCA_FLUX' if kind == 'pca' else 'CORR_FLUX'; f = np.array(d[col], float); q = np.array(d['QUALITY']); good = q == 0
            bkg = np.array(d['FLUX_BKG'], float)
        else:
            col = 'FLUX'; f = np.array(d[col], float); good = np.ones(n, bool)
    good &= np.isfinite(t) & np.isfinite(f) & (f > 0)
    o = np.argsort(t[good]); t, f = t[good][o], f[good][o]
    return dict(t=t, f=f / np.median(f), bkg=bkg[good][o], c1=c1[good][o], c2=c2[good][o], q=q[good][o],
                cad=float(np.median(np.diff(t))), label=label, sec=sec, file=os.path.basename(fp[0]), col=col)


def detrend(t, f, window=1.0, mask=None):
    """f divided by a wotan biweight trend fitted to the points where mask is True (all points if None)"""
    from wotan import flatten
    if mask is None:
        mask = np.ones(len(t), bool)
    tr = flatten(t[mask], f[mask], method='biweight', window_length=window, edge_cutoff=0, break_tolerance=0.4, return_trend=True)[1]
    ok = np.isfinite(tr)
    return f / np.interp(t, t[mask][ok], tr[ok])


def trim_gaps(t, gap=0.15, cut=0.25):
    """True = keep; drops `cut` days after and before every gap longer than `gap` and at the start and end"""
    keep = np.ones(len(t), bool); g = np.where(np.diff(t) > gap)[0]
    for a in np.r_[t[0], t[g + 1]]:
        keep &= ~((t >= a) & (t < a + cut))
    for b in np.r_[t[g], t[-1]]:
        keep &= ~((t <= b) & (t > b - cut))
    return keep


def prep(work, lab, P=PC, T0=TBOX, D=DBOX, maskc=True, win=1.0, invert=False):
    """b masked, gaps trimmed, biweight-detrended residuals (flux - 1) with a two-sided 5-sigma clip"""
    L = load(work, lab)
    if L is None:
        return None
    t, f = L['t'], L['f']
    if invert:
        f = 2 - f
    k = mask_b(t) & trim_gaps(t); t, f = t[k], f[k]
    inc = np.abs(phase(t, P, T0)) < 0.75 * D + 0.5 / 24
    r = detrend(t, f, win, mask=(~inc) if maskc else None) - 1
    ok = np.isfinite(r); t, r = t[ok], r[ok]
    sd = robust_sd(r); cl = np.abs(r) < 5 * sd
    return dict(t=t[cl], r=r[cl], cad=L['cad'], lab=lab, sec=L['sec'], col=L['col'])


def box(t, r, P, T0, D):
    """box depth (ppm-free units), white error from the robust out-of-box scatter, number of in-box points"""
    m = np.abs(phase(t, P, T0)) < D / 2
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    sd = robust_sd(r[~m]); return float(-r[m].mean()), float(sd / np.sqrt(m.sum())), int(m.sum())


def emp_err(t, r, P, T0, D, n=300, seed=0):
    """standard deviation of box depths at n random phase shifts (|shift| > 2 D) of the same light curve"""
    rng = np.random.default_rng(seed); out = []
    while len(out) < n:
        s = rng.uniform(-0.5, 0.5) * P
        if abs(s) < 2 * D:
            continue
        d, e, k = box(t, r, P, T0 + s, D)
        if np.isfinite(d):
            out.append(d)
    return float(np.std(out)), float(np.mean(out))


def a_rs(rho, P):
    """a/R* from the stellar density (g/cm^3) and the period (d), circular orbit"""
    return (6.674e-8 * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def rho_prior():
    mu = MS / RS ** 3 * RHO_SUN
    return mu, mu * np.sqrt((MS_E / MS) ** 2 + (3 * RS_E / RS) ** 2)


def read_csv(path):
    """read a repository CSV (lines starting with '#' skipped) into a dict of arrays (float where possible)"""
    import csv
    with open(path) as f:
        rows = [r for r in csv.reader(l for l in f if not l.startswith('#'))]
    out = {}
    for i, h in enumerate(rows[0]):
        col = [r[i] for r in rows[1:]]
        try:
            out[h] = np.array([float(x) if x != '' else np.nan for x in col])
        except ValueError:
            out[h] = np.array(col)
    return out


def write_csv(path, comments, cols, rows):
    """comment lines (prefixed '# '), then the column line, then rows (sequences of formatted strings or numbers; non-finite
    floats are written empty; fields containing commas are quoted)"""
    import csv
    with open(path, 'w', newline='') as f:
        for l in comments:
            f.write('# ' + l + '\n')
        w = csv.writer(f, lineterminator='\n'); w.writerow(cols)
        for r in rows:
            w.writerow(['' if ((isinstance(x, float) and not np.isfinite(x)) or str(x) in ('nan', '+nan', '-nan')) else str(x) for x in r])
    print('wrote', path, len(rows), 'rows')
