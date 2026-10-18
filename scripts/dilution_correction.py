"""Target flux fraction in the circular apertures from Gaia DR3 fluxes, and the dilution-corrected depth and radius.

usage: python scripts/dilution_correction.py <tic> --fits fit_r1.5_<epochs>,fit_r2.0_<epochs>,fit_r3.0_<epochs> [--work ./work]
input:  <work>/<tic>/gaia_sources.json (all Gaia DR3 sources within 60" of the TIC position, written by
        difference_images.py), <work>/<tic>/diffimg.json (per-epoch Gaussian PSF sigma of the out-of-transit image),
        <work>/<tic>/cut/s<sector>.npz (cutout WCS, target pixel), <work>/<tic>/<fit>_shape.json (fit_transits.py,
        one per aperture; the aperture radius is read from the file)
output: printed fractions per epoch and aperture, and <work>/<tic>/dilution.json
Model: every Gaia DR3 source with a G magnitude is placed at its proper-motion-corrected position in the cutout with a
circular Gaussian PSF of the per-epoch sigma, and again with 1.5 x that sigma (a bracket for the PSF wings that the
Gaussian fit underestimates); the flux of each source is proportional to 10^(-0.4 G). The target flux fraction of a
circular aperture (same centre and radius as aperture_photometry.py) is averaged over the epochs and the two sigmas.
Dilution-corrected depth = k^2 (raw-aperture common-shape fit) / fraction; Rp = sqrt(corrected depth) x R* with the
TIC v8.2 radius (`star` in targets.py). Sources beyond 60" are not included.
"""
import argparse, json, os
import numpy as np
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.coordinates import SkyCoord
from targets import TARGETS, workdir, add_common_args, lab, sector_of

p = add_common_args(argparse.ArgumentParser(description=__doc__)); p.add_argument('--fits', required=True)
args = p.parse_args(); c = TARGETS[args.tic]; wd = workdir(args, args.tic); st = c['star']
G = [r for r in json.load(open(os.path.join(wd, 'gaia_sources.json'))) if r['phot_g_mean_mag'] is not None]
DJ = json.load(open(os.path.join(wd, 'diffimg.json'))); N = 21; yy, xx = np.indices((N, N)); RADII = [1.5, 2.0, 2.5, 3.0]
print(f'{len(G)} Gaia DR3 sources with G within the query radius')
frac = {}
for key, v in DJ.items():
    if key == 'gaia_psf_test_stacked':
        continue
    sec = sector_of(lab(key)); z = np.load(os.path.join(wd, 'cut', f's{sec}.npz')); w = WCS(Header.fromstring(str(z['wcs'])))
    yr = float(z['epoch']) if 'epoch' in z else 2016 + (v['tc'] + 2457000 - 2457389) / 365.25; tx, ty = float(z['tx']), float(z['ty']); rr = np.hypot(xx - tx, yy - ty)
    frac[key] = {}
    for sg in (v['psf_sigma_px'], 1.5 * v['psf_sigma_px']):
        tot = np.zeros((N, N)); tgt = None
        for r in G:
            ra = r['ra'] + (r['pmra'] or 0) / 3.6e6 * (yr - 2016) / np.cos(np.radians(r['dec'])); de = r['dec'] + (r['pmdec'] or 0) / 3.6e6 * (yr - 2016)
            x0, y0 = [float(q) for q in w.world_to_pixel(SkyCoord(ra, de, unit='deg'))]
            img = 10 ** (-0.4 * r['phot_g_mean_mag']) * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sg * sg)) / (2 * np.pi * sg * sg); tot += img
            if str(r['source_id']) == str(c['gaia_dr3']):
                tgt = img
        frac[key][f'sigma{sg:.2f}'] = {f'r{a}': float(tgt[rr < a].sum() / tot[rr < a].sum()) for a in RADII}
    print(f'epoch {key} (S{sec}): PSF sigma {v["psf_sigma_px"]:.2f} px; target fraction', {k: {a: round(x, 3) for a, x in d.items()} for k, d in frac[key].items()})
out = dict(per_epoch=frac, corrected={})
print(f'R* = {st["R"]} +- {st["eR"]} Rsun (TIC v8.2)')
for fn in args.fits.split(','):
    S = json.load(open(os.path.join(wd, fn.replace('.json', '') + '_shape.json'))); a = S['aperture']
    fr = np.mean([np.mean([d[f'r{a}'] for d in e.values()]) for e in frac.values()]); raw = S['k'] ** 2 * 1e6; dc = raw / fr; k = np.sqrt(dc * 1e-6); rp = k * st['R'] * 109.08
    out['corrected'][a] = dict(raw_depth_ppm=raw, target_fraction=fr, corrected_depth_ppm=dc, k=k, Rp_Rearth=rp, Rp_RJup=rp / 11.21)
    print(f'r = {a} px: raw depth {raw:.0f} ppm, target fraction {fr:.3f}, corrected depth {dc:.0f} ppm, k = {k:.3f}, Rp = {rp:.1f} Rearth = {rp / 11.21:.2f} RJup')
json.dump(out, open(os.path.join(wd, 'dilution.json'), 'w'), indent=1)
