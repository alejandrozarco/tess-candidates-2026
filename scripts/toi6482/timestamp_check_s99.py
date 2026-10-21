"""TOI-6482.01 (TIC 189300053): compare the QLP sector-99 time stamps with SPOC FFI time stamps referred to the target position
(README section 6.8; MAST reports a QLP barycentric-correction error from orbit 210).

usage: python scripts/toi6482/timestamp_check_s99.py [--work ./work]
input:  the QLP sector-99 light curve in <work>/189300053/dl/ (fit_toi6482.py downloads it); a 3x3 TESScut cutout of
        sector 99 is requested from MAST.
TESScut TIME carries the barycentric correction for the FFI reference position (TIMECORR); it is removed and the
target's own correction r_earth . n_target / c is added (Earth barycentric position from astropy; TESS is within about
0.003 au of the Earth, below 2 s). Prints QLP minus recomputed time per orbit (QUALITY = 0) and across the +-0.7-d
windows of the two sector-99 transits. Requires numpy, astropy, astroquery.
"""
import argparse, glob, os
import numpy as np
from astropy.io import fits
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord, get_body_barycentric
import astropy.units as u
from astropy.time import Time

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); args = p.parse_args()
c = SkyCoord(141.92963828650713, -24.96502818809507, unit='deg')          # Gaia DR3 5662277217084277888 (epoch 2016.0)
fq = glob.glob(os.path.join(args.work, '189300053', 'dl', '**', 'hlsp_qlp_tess_ffi_s0099-0000000189300053_tess_v01_llc.fits'), recursive=True)[0]
q = fits.getdata(fq, 1)
h = Tesscut.get_cutouts(coordinates=c, size=3, sector=99)[0][1].data
t, tc = h['TIME'], h['TIMECORR']; g = np.isfinite(t)
tt = Time(t[g] - tc[g] + 2457000, format='jd', scale='tdb'); r = get_body_barycentric('earth', tt).xyz.to(u.m).value
tb = t[g] - tc[g] + (c.cartesian.xyz.value @ r) / 2.99792458e8 / 86400
tq = q['TIME']; ok = np.isfinite(tq) & (q['QUALITY'] == 0)
i = np.clip(np.searchsorted(tb, tq), 1, len(tb) - 1); j = np.where(np.abs(tb[i] - tq) < np.abs(tb[i - 1] - tq), i, i - 1); dt = (tq - tb[j]) * 86400
for o in (209, 210):
    m = ok & (q['ORBITID'] == o); print(f'orbit {o}: QLP - TESScut(target) = {dt[m].min():+.2f} .. {dt[m].max():+.2f} s (N {m.sum()})')
for tcen in (4048.21961, 4068.27927):
    m = ok & (np.abs(tq - tcen) < 0.7); print(f'window around BTJD {tcen}: {dt[m].min():+.3f} .. {dt[m].max():+.3f} s')
