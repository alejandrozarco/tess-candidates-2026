"""T11: own aperture photometry from SPOC 2-min target pixel files (no PDC). Own code.
compact <tpf_dir>         : store per-sector pixel cubes (FLUX+FLUX_BKG = raw calibrated, FLUX_ERR, QUALITY, aperture) as npz
phot <ap> <outdir> [inj.json]: own background (median of faint background-flagged pixels), aperture sum, flux-weighted
                            centroids; writes LC-like npz files (keys as t1_extract) so the window code runs unchanged.
                            inj.json: constant-geometry transit multiplied into the background-subtracted pixel fluxes."""
import sys, os, glob, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
PIXDIR = ROOT + '/out/tpfpix'
def compact(src):
    from astropy.io import fits
    os.makedirs(PIXDIR, exist_ok=True)
    for f in sorted(glob.glob(os.path.join(src, '*_tp.fits'))):
        with fits.open(f, memmap=True) as h:
            sec = h[0].header['SECTOR']; out = os.path.join(PIXDIR, 's%03d.npz' % sec)
            if os.path.exists(out): continue
            d = h[1].data; t = np.array(d['TIME'], float); ok = np.isfinite(t)
            np.savez(out, time=t[ok], raw=(np.array(d['FLUX'], np.float32) + np.array(d['FLUX_BKG'], np.float32))[ok],
                     ferr=np.array(d['FLUX_ERR'], np.float32)[ok], quality=np.array(d['QUALITY'], int)[ok], ap=np.array(h[2].data),
                     exptime=h[1].header['TIMEDEL'] * 86400, cam=h[0].header['CAMERA'], ccd=h[0].header['CCD'])
            print('compact', sec, ok.sum(), flush=True)
def model_inj(t, inj, exptime):
    import batman
    P, T0, _ = PL['d']
    p = batman.TransitParams(); p.t0 = T0 + inj['phase'] * P; p.per = P; p.rp = inj['k']; p.a = inj['aR']; p.inc = np.degrees(np.arccos(inj['b'] / inj['aR']))
    p.ecc = 0.; p.w = 90.; p.limb_dark = 'quadratic'; p.u = [inj['u1'], inj['u2']]
    return batman.TransitModel(p, t, supersample_factor=7, exp_time=exptime / 86400.).light_curve(p)
def phot(apmode, outdir, inj=None):
    os.makedirs(outdir, exist_ok=True)
    lcmeta = {int(z['meta'][0]['sector']): z['meta'][0] for z in (np.load(f, allow_pickle=True) for f in glob.glob(ROOT + '/out/lc/%d_s*_2m.npz' % TIC))}
    for f in sorted(glob.glob(PIXDIR + '/s*.npz')):
        z = np.load(f); sec = int(os.path.basename(f)[1:4])
        t, raw, q, ap, ferr = z['time'], z['raw'].astype(float), z['quality'], z['ap'], z['ferr'].astype(float)
        bgpix = ((ap & 4) > 0) & ((ap & 2) == 0)
        med = np.nanmedian(raw[::20], 0)
        faint = bgpix & (med <= np.nanpercentile(med[bgpix], 40))
        bkg = np.nanmedian(raw[:, faint], axis=1)
        img = raw - bkg[:, None, None]
        if inj is not None: img = img * model_inj(t, inj, float(z['exptime']))[:, None, None]
        medimg = np.nanmedian(img[::20], 0)
        ny, nx = medimg.shape; yy, xx = np.mgrid[:ny, :nx]
        spoc = (ap & 2) > 0
        # target centre: flux-weighted centroid of the SPOC aperture on the median image
        wsum = np.nansum(medimg * spoc); cy = np.nansum(medimg * spoc * yy) / wsum; cx = np.nansum(medimg * spoc * xx) / wsum
        rr = np.hypot(yy - cy, xx - cx)
        mask = spoc if apmode == 'spoc' else rr <= {'r15': 1.5, 'r25': 2.5}[apmode]
        flux = np.nansum(img[:, mask], axis=1); err = np.sqrt(np.nansum(ferr[:, mask] ** 2, axis=1))
        cm = rr <= 2.5; fw = np.nansum(img[:, cm], axis=1)
        c1 = np.nansum(img[:, cm] * xx[cm], axis=1) / fw; c2 = np.nansum(img[:, cm] * yy[cm], axis=1) / fw
        good = np.isfinite(flux) & (flux > 0)
        meta = dict(lcmeta[sec]); meta.update(dict(own_ap=apmode, ap_npix=int(mask.sum()), crowd=1.0, flfr=1.0))
        np.savez(os.path.join(outdir, '%d_s%03d_2m.npz' % (TIC, sec)), meta=np.array([meta], dtype=object), time=t[good],
                 pdcsap_flux=flux[good], pdcsap_flux_err=err[good], sap_flux=flux[good], sap_flux_err=err[good], sap_bkg=bkg[good] * mask.sum(),
                 pos_corr1=c1[good] - np.nanmedian(c1), pos_corr2=c2[good] - np.nanmedian(c2), mom_centr1=c1[good], mom_centr2=c2[good], quality=q[good])
if __name__ == '__main__':
    if sys.argv[1] == 'compact': compact(sys.argv[2])
    else: phot(sys.argv[2], sys.argv[3], json.load(open(sys.argv[4])) if len(sys.argv) > 4 else None)
