"""Predicted mid-transit times from the linear ephemeris and the TESS sectors that contain them.

usage: python scripts/predict_transits.py <tic> [--start 2026-10-18] [--months 12] [--out data]
output: data/predicted_transits_<tic>.csv and data/tess_sectors_<tic>.csv (formats in the file headers).
Uses tess-point (tess_stars2px) for the sectors in which the target falls on silicon and for the nominal sector
mid-times (window = mid-time +- 13.7 d). Errors are the propagated linear-ephemeris errors only; for TIC 402710280
an additional transit-timing-variation envelope of +-30 min applies (observed O-C up to +-29 min), for TIC 455017950
one of +-15 min (observed O-C up to +12 min); `ttv_min`, `ttv_note` in targets.py.
"""
import argparse, csv, os
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import tess_stars2px as TP
from targets import TARGETS

p = argparse.ArgumentParser(description=__doc__); p.add_argument('tic', type=int, choices=sorted(TARGETS))
p.add_argument('--start', default=Time.now().iso[:10]); p.add_argument('--months', type=float, default=12); p.add_argument('--out', default='data')
args = p.parse_args(); c = TARGETS[args.tic]; tic = args.tic; os.makedirs(args.out, exist_ok=True)
COORD = SkyCoord(c['ra'] * u.deg, c['dec'] * u.deg); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
def utc(bjd_tdb):
    """Geocentric UTC of a BJD_TDB instant: remove the barycentric light-travel time toward the target, then TDB -> UTC."""
    t = Time(bjd_tdb, format='jd', scale='tdb', location=GEO)
    for _ in range(2): t = Time(bjd_tdb, format='jd', scale='tdb', location=GEO) - (Time(t, location=GEO).light_travel_time(COORD, kind='barycentric'))
    return t.utc.iso[:16]
sc = TP.TESS_Spacecraft_Pointing_Data(); mid = {int(s): m for s, m in zip(sc.sectors, sc.midtimes)}
obs = sorted(set(int(x) for x in TP.tess_stars2px_function_entry(tic, c['ra'], c['dec'])[3]))
PS = f'{c["P"]:.10g}'; EPS = f'{c["eP"]:.{len(PS.split(".")[1])}f}'   # period and error with the decimals given in targets.py
t_start = Time(args.start).jd; t_end = t_start + args.months * 30.5
hdr = lambda lines: ''.join('# ' + l + '\n' for l in lines)
with open(os.path.join(args.out, f'predicted_transits_{tic}.csv'), 'w') as fo:
    env = f'; a TTV envelope of +-{c["ttv_min"]:.0f} min should be added ({c["ttv_note"]})' if c['ttv_min'] else ''
    fo.write(hdr([f'TIC {tic} = Gaia DR3 {c["gaia_dr3"]}: predicted mid-transit times for the {args.months:g} months from {args.start}, linear ephemeris T0 = BTJD {c["T0"]:.5f} +- {c["eT0"]:.5f}, P = {PS} +- {EPS} d',
                  f'err_min = 1-sigma propagated ephemeris error sqrt(eT0^2 + (n eP)^2){env}; *_utc = geocentric UTC (barycentric light-travel time toward the target removed; a topocentric observer differs by < 0.1 s)',
                  f'T14 = {c["T14h"]:.2f} h; tess_sector = sector whose nominal 27.4-d window contains the transit and in which tess-point places the target on silicon (empty = not observed by TESS)']))
    w = csv.writer(fo); w.writerow(['epoch', 'tc_btjd', 'tc_bjd_tdb', 'tc_utc', 'err_min', 'ingress_utc', 'egress_utc', 'tess_sector'])
    n = int(np.ceil((t_start - 2457000 - c['T0']) / c['P']))
    while c['T0'] + n * c['P'] + 2457000 <= t_end:
        tp = c['T0'] + n * c['P']; err = np.hypot(c['eT0'], n * c['eP']) * 1440; sec = [s for s, m in mid.items() if abs(tp - (m - 2457000)) < 13.7 and s in obs]
        w.writerow([n, f'{tp:.4f}', f'{tp + 2457000:.4f}', utc(tp + 2457000), f'{err:.1f}', utc(tp + 2457000 - c['T14h'] / 48), utc(tp + 2457000 + c['T14h'] / 48), ';'.join(map(str, sec))])
        print(f'n={n:3d} BTJD {tp:.4f} {utc(tp + 2457000)} UTC +-{err:.1f} min TESS {sec}'); n += 1
with open(os.path.join(args.out, f'tess_sectors_{tic}.csv'), 'w') as fo:
    fo.write(hdr([f'TIC {tic}: TESS sectors in which tess-point places the target on silicon, with nominal sector windows (midtime +- 13.7 d) and the linear-ephemeris transits inside them, {args.start}']))
    w = csv.writer(fo); w.writerow(['sector', 'window_start_utc', 'window_end_utc', 'predicted_transits_btjd', 'status'])
    for s in obs:
        m = mid[s] - 2457000; lo, hi = m - 13.7, m + 13.7; ns = np.arange(int((lo - c['T0']) / c['P']) - 1, int((hi - c['T0']) / c['P']) + 2); tp = c['T0'] + ns * c['P']; tp = tp[(tp > lo) & (tp < hi)]
        w.writerow([s, Time(lo + 2457000, format='jd').iso[:10], Time(hi + 2457000, format='jd').iso[:10], ';'.join(f'{x:.3f}' for x in tp), 'observed' if hi + 2457000 < t_start else 'future'])
print('sectors on silicon:', obs)
