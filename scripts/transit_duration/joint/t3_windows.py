"""T3: per-transit windows of L 98-59 d with masks (b, c, contaminating EB, flares) and per-transit telemetry. Own code.
Usage: t3_windows.py <cad 2m|20s> <flux pdc|sap> <eb mask|model|none> <out.npy> [lcdir] [inject.json]
inject.json (optional): constant-geometry injection {'phase':..,'k':..,'aR':..,'b':..,'u1':..,'u2':..} multiplied into the
 quality-selected, normalised flux at T0_d + phase*P, before windows, masks and clipping (null tests)."""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import *
from scipy.ndimage import median_filter
EB = json.load(open(ROOT + '/out/eb_ephem.json'))
Pd, T0d, Dd = PL['d']; Dd /= 24.
HALF = 0.16
import batman
def inject(t, f, exptime, inj, T0x):
    p = batman.TransitParams(); p.t0 = T0x; p.per = Pd; p.rp = inj['k']; p.a = inj['aR']; p.inc = np.degrees(np.arccos(inj['b'] / inj['aR']))
    p.ecc = 0.; p.w = 90.; p.limb_dark = 'quadratic'; p.u = [inj['u1'], inj['u2']]
    m = batman.TransitModel(p, t, supersample_factor=7, exp_time=exptime / 86400.)
    return f * m.light_curve(p)
def make_windows(lcs, fluxk, ebmode, inj=None, regs=False, apply_inj=True):
    segs, flares = [], []
    if inj is not None:
        import batman
    T0x = T0d + inj['phase'] * Pd if inj is not None else None
    for lc in lcs:
        meta = lc['meta']; q = lc['quality']
        fl = lc['pdcsap_flux'] if fluxk == 'pdc' else lc['sap_flux']
        fe = lc['pdcsap_flux_err'] if fluxk == 'pdc' else lc['sap_flux_err']
        g = (q == 0) & np.isfinite(fl) & np.isfinite(fe)
        t = lc['time'][g]; med = np.nanmedian(fl[g]); f = fl[g] / med; e = fe[g] / med
        bkg = lc['sap_bkg'][g] / med; pc1 = lc['pos_corr1'][g]; pc2 = lc['pos_corr2'][g]; mc1 = lc['mom_centr1'][g]; mc2 = lc['mom_centr2'][g]
        if inj is not None and apply_inj: f = inject(t, f, meta['exptime'], inj, T0x)
        tall = lc['time']; desat = tall[(q & 32) > 0]; strayl = tall[(q & (2048 | 4096)) > 0]
        # masks of other planets
        other = np.zeros(len(t), bool)
        for pn in ('b', 'c'):
            P, T0, D = PL[pn]; ph = (t - T0 + 0.5 * P) % P - 0.5 * P; other |= np.abs(ph) < (0.5 * D / 24 + 0.025)
        if inj is not None:   # mask the real d when injecting elsewhere
            ph = (t - T0d + 0.5 * Pd) % Pd - 0.5 * Pd; other |= np.abs(ph) < (0.5 * Dd + 0.025)
        phe = (t - EB['T0'] + 0.5 * EB['P']) % EB['P'] - 0.5 * EB['P']
        ebin = np.abs(phe) < (0.5 * EB['D_h'] / 24 + 0.015)     # box half-width 0.59 h + 22 min (folded profile reaches 0 at ~0.9 h)
        if ebmode == 'model':   # divide out a trapezoid (total 1.7 h, flat 0.8 h) scaled to the per-sector depth
            dep = EB['per_sector'].get(str(meta['sector']), [EB['depth_ppm']])[0] * 1e-6
            x = np.abs(phe) * 24; shp = np.clip((0.85 - x) / 0.45, 0, 1)
            f = f / (1 - dep * shp)
        # flare detection on the full sector: running median over 2 h, planets masked
        Tc_self = T0x if inj is not None else T0d
        phs = (t - Tc_self + 0.5 * Pd) % Pd - 0.5 * Pd
        pm = other | (np.abs(phs) < 0.5 * Dd + 0.01)
        nwin = int(round((2 / 24.) / (meta['exptime'] / 86400.))) | 1
        trend = np.interp(t, t[~pm], median_filter(f[~pm], size=nwin, mode='nearest'))
        r = f - trend; s = 1.4826 * np.median(np.abs(r[~pm] - np.median(r[~pm])))
        hi = r > 3.0 * s
        cand = (hi & (np.roll(hi, 1) | np.roll(hi, -1))) | (r > 5 * s)
        flm = np.zeros(len(t), bool)
        ncad = meta['exptime'] / 86400.
        for j in np.where(cand)[0]:
            flm |= (t >= t[j] - 2 * ncad - 1e-6) & (t <= t[j] + max(10 * ncad, 0.0139))   # -2 cadences, +max(10 cad, 20 min)
        # flare list (groups)
        if cand.any():
            idx = np.where(cand)[0]; grp = np.split(idx, np.where(np.diff(idx) > 3)[0] + 1)
            for gg in grp: flares.append((int(meta['sector']), float(t[gg[0]]), float(r[gg].max() / s), float(r[gg].max())))
        n0, n1 = np.ceil((t.min() - Tc_self) / Pd), np.floor((t.max() - Tc_self) / Pd)
        for n in np.arange(n0, n1 + 1):
            tc = Tc_self + n * Pd
            sel = np.abs(t - tc) < HALF
            if sel.sum() < 10: continue
            drop = other[sel] | flm[sel]
            if ebmode == 'mask': drop |= ebin[sel]
            tt = t[sel]; intr = np.abs(tt - tc) < 0.5 * Dd
            cad_d = meta['exptime'] / 86400.
            cov_in = (intr & ~drop).sum() / (Dd / cad_d)
            pre = ((tt < tc - 0.5 * Dd - 0.01) & ~drop).sum() * cad_d; post = ((tt > tc + 0.5 * Dd + 0.01) & ~drop).sum() * cad_d
            ok = cov_in > 0.8 and pre > 0.04 and post > 0.04
            keep = ~drop
            if keep.sum() < 10: continue
            # telemetry for the window (all kept cadences) and near the transit
            near = keep & (np.abs(tt - tc) < Dd)
            dd = np.min(np.abs(desat - tc)) if len(desat) else 99.
            ds = np.min(np.abs(strayl - tc)) if len(strayl) else 99.
            gaps = np.where(np.diff(t) > 0.1)[0]
            edges = np.concatenate([[t[0]], t[gaps], t[gaps + 1], [t[-1]]])
            tel = dict(bkg_mean=float(np.mean(bkg[sel][keep])), bkg_in_minus_out=float(np.mean(bkg[sel][near]) - np.mean(bkg[sel][keep & ~near])) if near.any() else np.nan,
                       bkg_slope=float(np.polyfit(tt[keep] - tc, bkg[sel][keep], 1)[0]),
                       pc_rms=float(np.hypot(np.std(pc1[sel][keep]), np.std(pc2[sel][keep]))) if np.isfinite(pc1[sel]).any() else np.nan,
                       mc_rms=float(np.hypot(np.nanstd(mc1[sel][keep]), np.nanstd(mc2[sel][keep]))),
                       mc_in_minus_out=float(np.hypot(np.nanmean(mc1[sel][near]) - np.nanmean(mc1[sel][keep & ~near]), np.nanmean(mc2[sel][near]) - np.nanmean(mc2[sel][keep & ~near]))) if near.any() else np.nan,
                       n_desat_win=int(np.sum(np.abs(desat - tc) < HALF)), dt_desat=float(dd), dt_stray=float(ds),
                       dt_edge=float(np.min(np.abs(edges - tc))), n_flare_masked=int((flm[sel] & ~other[sel]).sum()),
                       eb_overlap=bool((ebin[sel] & intr).any()), eb_in_window=bool(ebin[sel].any()), cov_in=float(cov_in),
                       other_in_window=bool(other[sel].any()), rms_raw=float(np.std(np.diff(f[sel][keep])) / np.sqrt(2)))
            segs.append(dict(sector=int(meta['sector']), epoch=int(n), tc_lin=float(tc), t=tt[keep], f=f[sel][keep], e=e[sel][keep],
                             reg=np.vstack([pc1[sel][keep], pc2[sel][keep], bkg[sel][keep]]).T,
                             ebin=ebin[sel][keep], ok=bool(ok), cam=int(meta['cam']), ccd=int(meta['ccd']), crowd=float(meta['crowd']),
                             flfr=float(meta['flfr']), procver=meta['procver'], data_rel=int(meta['data_rel']), exptime=float(meta['exptime']),
                             ap_npix=meta['ap_npix'], tel=tel))

    good = [s for s in segs if s['ok']]
    return good, flares

if __name__ == '__main__':
    cad, fluxk, ebmode, outp = sys.argv[1:5]
    lcdir = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != '-' else None
    inj = json.load(open(sys.argv[6])) if len(sys.argv) > 6 else None
    good, flares = make_windows(load_lc(cad, lcdir), fluxk, ebmode, inj)
    print(cad, fluxk, ebmode, 'usable', len(good), 'with EB in transit', sum(s['tel']['eb_overlap'] for s in good),
          'with b/c points masked', sum(s['tel']['other_in_window'] for s in good), 'flare groups', len(flares))
    from collections import Counter
    print(sorted(Counter(s['sector'] for s in good).items()))
    np.save(outp, np.array(good, dtype=object), allow_pickle=True)
    np.save(outp.replace('.npy', '_flares.npy'), np.array(flares))
