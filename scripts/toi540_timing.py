"""TOI-540 b (TIC 200322593): audit of the sector-98 time stamps.

usage: python scripts/toi540_timing.py [--work ./work] [--out data]
input:  SPOC sector-98 products of TOI-540 (other_download.py 200322593) and, if present, of TOI-270 (--controls)
output: data/timing_audit_200322593.csv
Compares the cadence spacing and the spacecraft-frame times (TIME - TIMECORR) of the 2-min target pixel file, the
2-min and 20-s light curves of TOI-540 and the 2-min light curve of TOI-270 (same camera and CCD), and compares
TIMECORR with the Roemer delay computed from the JPL Horizons barycentric vectors of TESS (-95) (network access).
"""
import argparse, csv, glob, os
import numpy as np
from astropy.io import fits
import astropy.units as u
from targets import OTHER

p = argparse.ArgumentParser(description=__doc__); p.add_argument('--work', default='./work'); p.add_argument('--out', default='data')
args = p.parse_args(); TIC = 200322593; c = OTHER[TIC]; DL = os.path.join(args.work, str(TIC), 'dl', 'mastDownload', 'TESS')
DATE = __import__('datetime').datetime.now(__import__('datetime').timezone.utc).strftime('%Y-%m-%d')
hdr = lambda lines: ''.join('# ' + l + '\n' for l in lines)
rows_t = []
def tinfo(fn):
    with fits.open(fn, memmap=True) as h: d = h[1].data; return np.array(d['TIME'], float), np.array(d['TIMECORR'], float), np.array(d['CADENCENO']), h[0].header, h[1].header
fl = {'TPF 2-min': glob.glob(os.path.join(DL, f'*s0098-{TIC:016d}-*-s', '*-s_tp.fits')), 'LC 2-min': glob.glob(os.path.join(DL, f'*s0098-{TIC:016d}-*-s', '*-s_lc.fits')),
      'LC 20-s': glob.glob(os.path.join(DL, f'*s0098-{TIC:016d}-*', '*fast-lc.fits')), 'TOI-270 LC 2-min': glob.glob(os.path.join(args.work, 'controls', 'dl', 'mastDownload', 'TESS', '*s0098-0000000259377017*', '*-s_lc.fits'))}
ref = None
for k, f_ in fl.items():
    if not f_: continue
    t, tcor, cad, h0, h1 = tinfo(f_[0]); ok = np.isfinite(t); dtc = np.diff(t[ok]) * 86400 / np.diff(cad[ok])
    row = [k, h0['CAMERA'], h0['CCD'], h0.get('PROCVER'), h0.get('DATA_REL'), h1.get('TIMESYS'), h1.get('TIMEPIXR'), f'{np.median(np.diff(t[ok])) * 86400:.4f}', f'{dtc.min():.4f}', f'{dtc.max():.4f}', '', '']
    if ref is None: ref = (t, tcor, cad)
    else:
        cm, i1, i2 = np.intersect1d(ref[2], cad, return_indices=True); ff = np.isfinite(ref[0][i1]) & np.isfinite(t[i2])
        if ff.sum(): row[10] = f'{np.nanmax(np.abs((ref[0][i1] - ref[1][i1]) - (t[i2] - tcor[i2]))[ff]) * 86400:.3f}'; row[11] = f'{np.nanmedian((ref[0][i1] - t[i2])[ff]) * 86400:.1f}'
    rows_t.append(row)
try:
    import requests
    t, tcor, cad, h0, h1 = tinfo(fl['LC 2-min'][0]); ok = np.isfinite(t) & (tcor != 0); idx = np.where(ok)[0][np.linspace(0, ok.sum() - 1, 12).astype(int)]
    prm = {'format': 'json', 'COMMAND': "'-95'", 'EPHEM_TYPE': 'VECTORS', 'CENTER': "'@0'", 'REF_PLANE': 'FRAME', 'REF_SYSTEM': 'ICRF', 'OUT_UNITS': 'AU-D', 'VEC_TABLE': '1', 'CSV_FORMAT': 'YES',
           'TIME_TYPE': 'TDB', 'TLIST': ' '.join(f"'{x:.8f}'" for x in t[idx] - tcor[idx] + 2457000)}
    txt = requests.get('https://ssd.jpl.nasa.gov/api/horizons.api', params=prm, timeout=90).json()['result']
    pos = np.array([[float(x) for x in ln.split(',')[2:5]] for ln in txt[txt.index('$$SOE') + 5: txt.index('$$EOE')].strip().splitlines()]) * u.au.to(u.m)
    yr = 2016.0 + (np.nanmedian(t) + 2457000 - 2457388.5) / 365.25; ra = np.radians(c['ra'] + c['pmra'] * (yr - 2016) / 3.6e6 / np.cos(np.radians(c['dec']))); de = np.radians(c['dec'] + c['pmdec'] * (yr - 2016) / 3.6e6)
    roem = pos @ np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)]) / 2.99792458e8 / 86400
    rows_t.append(['TIMECORR minus Roemer delay from JPL Horizons TESS (-95) vectors, 12 cadences', '', '', '', '', '', '', '', '', '', f'{np.max(np.abs(tcor[idx] - roem)) * 86400:.3f}', ''])
except Exception as ex:
    rows_t.append([f'Horizons comparison not made: {repr(ex)[:80]}'] + [''] * 11)
with open(os.path.join(args.out, 'timing_audit_200322593.csv'), 'w') as fo:
    fo.write(hdr([f'TOI-540 b = TIC 200322593: sector-98 time stamps (scripts/toi540_analysis.py), {DATE}',
                  'cadence_s = median cadence spacing; step_min/max_s = spacing per cadence-number step; max_dt_vs_tpf_s = largest |(TIME - TIMECORR)| difference against the 2-min TPF on common cadence numbers (Horizons row: largest |TIMECORR - Roemer delay|)',
                  'median_time_diff_vs_tpf_s = median TIME difference on common cadences (non-zero for another star: different barycentric correction)']))
    w = csv.writer(fo); w.writerow(['product', 'camera', 'ccd', 'procver', 'data_rel', 'timesys', 'timepixr', 'cadence_s', 'step_min_s', 'step_max_s', 'max_dt_vs_tpf_s', 'median_time_diff_vs_tpf_s']); w.writerows(rows_t)

print('written', os.path.join(args.out, 'timing_audit_200322593.csv'))
