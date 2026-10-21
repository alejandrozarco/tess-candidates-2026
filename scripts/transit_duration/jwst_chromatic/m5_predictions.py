"""M5: predicted T14 of L 98-59 d for the GO 10498 transits (observation 4 executed 2026-08-01; eight visits planned
for 2027), on the scales of the JWST reductions of sections 6.6.10 and 6.6.11, for the TESS linear trend
(-1.185 +- 0.112 min/yr, section 6.6.10) and for no change. Own code.
Scales (T14 of the GO 4098 visit, 2024.195): joint NRS1 + NRS2 fit ns_main_v2, NRS1 only ns_nrs1, NRS2 only ns_nrs2,
pipeline NRS2 pipe_nrs2 ($TDJ_WORK/out/jwst/fits/, ../joint/run_joint.sh), and the detector average with the
systematic floor of c11_summary.py ($TDC_WORK/out/jwst_summary.json, Kurucz-prior LD). Mid-transit from the Cadieux et
al. (2025) linear ephemeris; the TTV of about 1.7 min shifts the mid-time, not T14. The visit dates are the STScI plan
windows of 2026-10-21 (GO 10498 visit-status page). Usage: m5_predictions.py <out.json>"""
import json, os, sys, numpy as np
from astropy.time import Time
F = os.path.join(os.environ.get('TDJ_WORK', 'work/transit_duration_joint'), 'out', 'jwst', 'fits')
H = os.environ.get('TDC_WORK', 'work/transit_duration_chromatic')
P, T0 = 7.450729, 1362.7400; S, SE = -1.185, 0.112
def yr(btjd): return 2000.0 + (btjd + 2457000.0 - 2451544.5) / 365.25
t_ref = yr(3381.8876)                         # GO 4098 obs. 2 mid-transit
sc = {}
for n, key in (('joint NRS1+NRS2 (ns_main_v2)', 'ns_main_v2'), ('NRS1 only (ns_nrs1)', 'ns_nrs1'), ('NRS2 only (ns_nrs2)', 'ns_nrs2'), ('NRS2 pipeline (pipe_nrs2)', 'pipe_nrs2')):
    d = json.load(open(F + '/%s.json' % key)); sc[n] = [float(x) for x in d['T14']['4098'][:2]]
va = json.load(open(H + '/out/jwst_summary.json'))['prior']['visit_avg']['4098']
sc['detector average with floor (c11_summary)'] = [float(va[0]), float(va[2])]
dates = [('obs 4 (executed)', '2026-08-01T02:56'), ('obs 2', '2027-03-20T02:17'), ('obs 5', '2027-04-03T23:55'), ('obs 7', '2027-05-26T03:38'), ('obs 3', '2027-06-10T01:16'),
         ('obs 9', '2027-06-17T12:06'), ('obs 6', '2027-07-02T09:44'), ('obs 8', '2027-07-09T20:33'), ('obs 1', '2027-08-01T05:00')]
print('reference epoch (GO 4098 obs 2) %.3f' % t_ref)
out = dict(reference_year=t_ref, slope=[S, SE], scales=sc, transits=[])
for lab, iso in dates:
    b = Time(iso).jd - 2457000; n = np.round((b - T0) / P); tc = T0 + n * P; ty = yr(tc); dt = ty - t_ref
    row = dict(label=lab, tc_btjd=float(tc), tc_utc=Time(tc + 2457000, format='jd', scale='tdb').utc.iso[:16], year=float(ty), dt_yr=float(dt), change=[S * dt, abs(SE * dt)], scales={})
    line = '%-18s %s  %.3f  dt %.3f yr  change %.2f +- %.2f |' % (lab, row['tc_utc'], ty, dt, S * dt, abs(SE * dt))
    for k, (T, e) in sc.items():
        row['scales'][k] = dict(trend=[float(T + S * dt), float(np.hypot(e, SE * dt))], none=[float(T), float(e)])
        line += ' %s: trend %.2f +- %.2f, none %.2f +- %.2f |' % (k.split(' (')[0], T + S * dt, np.hypot(e, SE * dt), T, e)
    out['transits'].append(row); print(line)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
