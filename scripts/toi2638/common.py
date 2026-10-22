"""Shared definitions for the TOI-2638.01 (TIC 458342769) scripts (README section 6.9).

Light-curve products (downloaded by download_toi2638.py into <work>/458342769/dl/):
  SPOC     2-min light curves, sectors 36, 37, 63, 64, 90, 99 (PDCSAP_FLUX, SAP_FLUX)
  TESSSPOC TESS-SPOC full-frame-image light curves, sectors 36, 37, 63, 64 (PDCSAP_FLUX; 600 s and 200 s)
  QLP      QLP full-frame-image light curves, sectors 10, 36, 37, 63, 64, 90, 99 (KSPSAP_FLUX; 1800, 600 and 200 s)
  TGLC     TGLC light curve of sector 10 (cal_psf_flux)
Cadences with QUALITY = 0 (TGLC: TGLC_flags = TESS_flags = 0), finite positive flux, divided by the sector median.

Transit model: batman (Kreidberg 2015), circular orbit, quadratic limb darkening fixed at U, averaged over each exposure
by sampling the exposure at the centres of NSUB(cadence) equal parts (midpoint sampling; 3-4 sub-samples at approximately
120 s, 5-6 at 200 s, 15-16 at 600 s and 45 at 1800 s).
Windows: +-0.30 d around each predicted mid-transit; a quadratic polynomial in time fitted to the points more than
0.05 d from mid-transit is divided out; a window is used if it has at least 60 % of its cadences, at least half of the
cadences within +-0.03 d and at least 5 out-of-transit points on each side; the per-point error is the standard
deviation of the out-of-transit points of the window.
"""
import glob, os, sys, types
import numpy as np

TIC = 458342769
GAIA_TARGET = 5353676910584543872
GAIA_NB17 = 5353676910584550272          # TIC 458342741 = TYC 8605-330-1, 17" from the target
TIC_NB17 = 458342741
RA, DEC = 159.84255570616043, -55.491806335768544   # Gaia DR3, epoch 2016.0
SPOC_SECTORS = [36, 37, 63, 64, 90, 99]
QLP_SECTORS = [10, 36, 37, 63, 64, 90, 99]
TESSSPOC_SECTORS = [36, 37, 63, 64]
P_TOI, T0_TOI = 3.929934391848882, 2459282.149063 - 2457000   # ExoFOP TOI table (TESS project, 2025-07-09)
U = [0.17, 0.27]          # quadratic limb darkening, TESS band, Teff about 9000-9500 K, log g 4.2 (approximate, fixed)
HALF, MASK = 0.30, 0.05   # window half-width and out-of-transit mask (d)
REF_YEAR = 2023.0
DATE = '2026-10-22'


def import_batman():
    """batman-package 2.5.x imports distutils, which Python 3.12 no longer ships; provide a stub."""
    if 'distutils' not in sys.modules:
        d = types.ModuleType('distutils'); c = types.ModuleType('distutils.ccompiler')
        c.new_compiler = lambda *a, **k: None; d.ccompiler = c
        sys.modules['distutils'] = d; sys.modules['distutils.ccompiler'] = c
    import batman
    return batman


batman = import_batman()


def nsub(cad_d):
    """number of sub-samples per exposure (midpoint sampling)"""
    return max(3, int(np.ceil(cad_d * 86400 / 40.0)))


def subgrid(t, cad_d):
    """times at the centres of n equal parts of each exposure; returns (expanded times, n)"""
    n = nsub(cad_d); off = ((np.arange(n) + 0.5) / n - 0.5) * cad_d
    return (t[:, None] + off[None, :]).ravel(), n


def year(btjd):
    return 2000.0 + (np.asarray(btjd) + 2457000 - 2451545.0) / 365.25


def dldir(work):
    d = os.path.join(work, str(TIC), 'dl'); os.makedirs(d, exist_ok=True); return d


def load(work, prov, flux=None):
    """dict sector -> dict(t, f, e, cad, cam, ccd, hdr) for one product"""
    out = {}
    from astropy.io import fits
    pat = {'SPOC': 'tess*-s00??-0000000458342769-*-s_lc.fits', 'TESSSPOC': 'hlsp_tess-spoc_*458342769*_lc.fits',
           'QLP': 'hlsp_qlp_tess_ffi_s00??-0000000458342769_tess_v01_llc.fits', 'TGLC': 'hlsp_tglc_*5353676910584543872*.fits'}[prov]
    for fn in sorted(glob.glob(os.path.join(dldir(work), 'lc', pat))):
        with fits.open(fn) as h:
            d = h[1].data; h0 = h[0].header; h1 = h[1].header
            s = int(h0.get('SECTOR', 0)) or int(fn.split('-s00')[1][:2])
            if prov == 'TGLC':
                t = np.array(d['time'], float); fl = np.array(d[flux or 'cal_psf_flux'], float)
                ok = (np.array(d['TGLC_flags']) == 0) & (np.array(d['TESS_flags']) == 0)
                fe = np.full_like(fl, np.nan)
            else:
                t = np.array(d['TIME'], float); q = np.array(d['QUALITY'], int); ok = q == 0
                col = flux or ('KSPSAP_FLUX' if prov == 'QLP' else 'PDCSAP_FLUX')
                if prov == 'QLP' and col not in d.names:
                    col = 'DET_FLUX'
                fl = np.array(d[col], float)
                fe = np.array(d[col + '_ERR'], float) if col + '_ERR' in d.names else np.full_like(fl, np.nan)
            ok &= np.isfinite(t) & np.isfinite(fl) & (fl > 0)
            t, fl, fe = t[ok], fl[ok], fe[ok]; med = np.nanmedian(fl)
            out[s] = dict(t=t, f=fl / med, e=fe / med, cad=float(np.median(np.diff(t))), cam=h0.get('CAMERA'), ccd=h0.get('CCD'),
                          hdr={k: h1.get(k) for k in ('CROWDSAP', 'FLFRCSAP', 'PDCMETHD', 'TIMEDEL')}, file=os.path.basename(fn))
    return out


def crowd_correct(d):
    """SAP divided by its median, corrected for the CROWDSAP fraction of other stars' light: (f - (1 - c)) / c"""
    for s, x in d.items():
        c = x['hdr']['CROWDSAP']; x['f'] = (x['f'] - (1 - c)) / c; x['e'] = x['e'] / c
    return d


def windows(d, P, T0, half=HALF, mask=MASK):
    """detrended transit windows per sector: dict sector -> dict(t, f, e, ep, cad, ntr, year)"""
    W = {}
    for s, x in d.items():
        t, f = x['t'], x['f']; cad = x['cad']
        E = np.round((t - T0) / P); tt, ff, ee, ep = [], [], [], []
        for n in np.unique(E):
            tc = T0 + n * P; m = np.abs(t - tc) < half
            if m.sum() < 0.6 * 2 * half / cad:
                continue
            dt = t[m] - tc; oot = np.abs(dt) > mask
            if (np.abs(dt) < 0.03).sum() < 0.5 * 0.06 / cad or (dt[oot] < 0).sum() < 5 or (dt[oot] > 0).sum() < 5:
                continue
            c = np.polyfit(dt[oot], f[m][oot], 2); r = f[m] / np.polyval(c, dt)
            tt.append(t[m]); ff.append(r); ee.append(np.full(m.sum(), np.std(r[oot]))); ep.append(np.full(m.sum(), int(n)))
        if tt:
            t_ = np.concatenate(tt)
            W[s] = dict(t=t_, f=np.concatenate(ff), e=np.concatenate(ee), ep=np.concatenate(ep), cad=cad, ntr=len(tt), year=float(year(np.median(t_))))
    return W


class Model:
    """exposure-averaged batman model on the windows of one sector, phase relative to (P, T0)"""
    def __init__(self, t, cad, P, T0):
        self.ph = (t - T0 + 0.5 * P) % P - 0.5 * P
        tg, self.n = subgrid(self.ph, cad)
        pm = batman.TransitParams(); pm.t0 = 0.0; pm.per = P; pm.rp = 0.05; pm.a = 10.0; pm.inc = 85.0; pm.ecc = 0.0; pm.w = 90.0
        pm.limb_dark = 'quadratic'; pm.u = list(U); self.pm = pm; self.P = P
        self.m = batman.TransitModel(pm, tg)

    def depth(self, k, b, aR, dt0=0.0):
        """1 - flux (positive in transit), exposure-averaged"""
        pm = self.pm; pm.rp = k; pm.a = aR; pm.inc = np.degrees(np.arccos(min(b / aR, 0.9999))); pm.t0 = dt0
        return 1.0 - self.m.light_curve(pm).reshape(-1, self.n).mean(1)


def peak_depth(k, b, aR, P):
    """instantaneous model depth at mid-transit (no exposure averaging)"""
    pm = batman.TransitParams(); pm.t0 = 0.0; pm.per = P; pm.rp = k; pm.a = aR; pm.inc = np.degrees(np.arccos(min(b / aR, 0.9999)))
    pm.ecc = 0.0; pm.w = 90.0; pm.limb_dark = 'quadratic'; pm.u = list(U)
    return float(1.0 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0])


def model_t14(k, b, aR, P):
    """first-to-fourth-contact duration (min) of the circular model; 0 if no transit"""
    if b >= 1 + k:
        return 0.0
    inc = np.arccos(b / aR)
    return float(P / np.pi * np.arcsin(min(1.0, np.sqrt((1 + k) ** 2 - b ** 2) / aR / np.sin(inc))) * 1440)


def trapezoid(p, ph, cad):
    """trapezoid (depth, T14, ingress fraction f12 = 2 T12 / T14, centre offset), exposure-averaged (days)"""
    dep, T14, f12, dt0 = p; T12 = f12 * T14 / 2
    tg, n = subgrid(ph, cad); x = np.abs(tg - dt0); m = np.zeros_like(tg)
    m[x <= T14 / 2 - T12] = 1.0; ing = (x > T14 / 2 - T12) & (x < T14 / 2); m[ing] = (T14 / 2 - x[ing]) / T12
    return dep * m.reshape(-1, n).mean(1)


def wlinfit(x, y, e):
    """weighted linear fit y = a + b x; returns (a, b), (ea, eb), chi2"""
    A = np.c_[np.ones(len(x)), x]; W = 1 / np.asarray(e) ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None])); c = C @ (A.T @ (W * y)); chi = float(np.sum((y - A @ c) ** 2 * W))
    return c, np.sqrt(np.diag(C)), chi


def hdr(f, lines):
    """comment lines, then the column line (last element) without '#'"""
    for l in lines[:-1]:
        f.write('# ' + l + '\n')
    f.write(lines[-1] + '\n')
