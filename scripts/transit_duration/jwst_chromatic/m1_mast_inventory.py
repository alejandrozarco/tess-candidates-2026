"""M1: complete MAST inventory of JWST observations of L 98-59 (TIC 307210830). Own code.
Queries CAOM (astroquery.mast.Observations) by position (L 98-59, 3 arcmin cone, all JWST), by target name
variants, and by proposal id 10498; prints one row per observation (obs_id, proposal, PI, instrument, filters,
t_min/t_max as UTC, dataRights, t_obs_release, calib_level). Positive controls (silent-failure rule): the same calls
for a known public L 98-59 product (jw01224-o001) and an unrelated known public program (WASP-39 b, ERS 1366),
both must return rows; an empty result for a control aborts with a non-zero exit.
Usage: m1_mast_inventory.py <out.csv>"""
import sys, json
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Observations
outp = sys.argv[1]
RA, DEC = 124.5317, -68.3130          # L 98-59, ICRS (J2000 position; proper motion ~0.3"/yr is small vs the cone)
def utc(mjd):
    try: return Time(float(mjd), format='mjd').utc.iso[:16]
    except Exception: return ''
def ctrl(name, **kw):
    t = Observations.query_criteria(**kw)
    print('CONTROL', name, kw, 'rows', len(t), flush=True)
    if len(t) == 0: raise SystemExit('positive control FAILED: ' + name)
    return t
# positive controls
c1 = Observations.query_criteria(obs_collection='JWST', proposal_id='1224')
print('CONTROL proposal 1224 rows', len(c1), flush=True)
if len(c1) == 0: raise SystemExit('positive control FAILED: proposal 1224')
c2 = Observations.query_criteria(obs_collection='JWST', proposal_id='1366')
print('CONTROL ERS 1366 (WASP-39) rows', len(c2), flush=True)
if len(c2) == 0: raise SystemExit('positive control FAILED: proposal 1366')
c3 = Observations.query_region(SkyCoord(217.3266, -3.4445, unit='deg'), radius=2 * u.arcmin)   # WASP-39 cone, all collections
c3 = c3[c3['obs_collection'] == 'JWST']; print('CONTROL WASP-39 cone JWST rows', len(c3), flush=True)
if len(c3) == 0: raise SystemExit('positive control FAILED: WASP-39 cone')
# science queries
tabs = []
cone = Observations.query_region(SkyCoord(RA, DEC, unit='deg'), radius=3 * u.arcmin)
print('L 98-59 cone, all collections:', len(cone), 'rows; collections', sorted(set(cone['obs_collection'])), flush=True)
cj = cone[cone['obs_collection'] == 'JWST']; print('L 98-59 cone JWST rows', len(cj), flush=True); tabs.append(cj)
for tn in ('L-98-59', 'L98-59', 'L 98-59', 'L-98-59-D', 'L-98-59-B', 'L-98-59-C', 'TOI-175', 'TOI-175.01', 'LHS-1966', 'L98-59D', 'L-98-59D', 'TIC-307210830'):
    try:
        t = Observations.query_criteria(obs_collection='JWST', target_name=tn)
    except Exception as e:
        print('target_name', tn, 'error', e); continue
    print('target_name', repr(tn), 'rows', len(t), flush=True)
    if len(t): tabs.append(t)
for pid in ('10498', '1201', '1224', '2512', '4098', '3942'):
    t = Observations.query_criteria(obs_collection='JWST', proposal_id=pid)
    t = t[np.array([abs(float(r['s_ra']) - RA) < 0.2 and abs(float(r['s_dec']) - DEC) < 0.1 if np.isfinite(float(r['s_ra'])) else False for r in t], bool)] if len(t) else t
    print('proposal', pid, 'rows near L 98-59', len(t), flush=True)
    if len(t): tabs.append(t)
from astropy.table import vstack, unique
T = unique(vstack(tabs, metadata_conflicts='silent'), keys='obsid')
T.sort('t_min')
cols = ['obs_id', 'proposal_id', 'proposal_pi', 'instrument_name', 'filters', 'target_name', 't_min', 't_max', 'dataRights', 't_obs_release', 'calib_level', 'obsid', 'dataproduct_type', 't_exptime']
rows = []
for r in T:
    d = {c: (str(r[c]) if c in T.colnames else '') for c in cols}
    d['start_utc'] = utc(r['t_min']); d['end_utc'] = utc(r['t_max']); d['release_utc'] = utc(r['t_obs_release'])
    rows.append(d)
import csv
with open(outp, 'w', newline='') as fo:
    w = csv.DictWriter(fo, fieldnames=cols + ['start_utc', 'end_utc', 'release_utc']); w.writeheader(); w.writerows(rows)
now = Time.now().mjd
for d in rows:
    pub = 'PUBLIC' if d['release_utc'] and float(T[T['obs_id'] == d['obs_id']]['t_obs_release'][0]) <= now else 'EXCLUSIVE'
    print('%-48s %6s %-22s %-10s %-22s %-14s %s..%s rights=%s release=%s %s lvl=%s' % (d['obs_id'], d['proposal_id'], d['proposal_pi'][:22], d['instrument_name'], d['filters'][:22], d['target_name'][:14], d['start_utc'], d['end_utc'][11:], d['dataRights'], d['release_utc'], pub, d['calib_level']))
print('n rows', len(rows))
