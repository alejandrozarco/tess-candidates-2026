"""Inject a constant-geometry copy of planet d at a phase offset into the real PDCSAP/SAP fluxes. Usage: v8_inject.py <phase_off> <b_inj> <outdir>
Input: the extracted light curves in $LCDIR (default $TD_WORK/lc); output: the same files with the injected transit in <outdir>."""
import sys, os, glob, numpy as np, batman
off, binj, outd = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
os.makedirs(outd, exist_ok=True)
P, T0 = 7.450729, 1362.74002 + off * 7.450729
k, aR, u = 0.0453, 34.80, [0.1895, 0.4604]
LCDIR = os.environ.get('LCDIR', os.path.join(os.environ.get('TD_WORK', os.path.join('work', 'transit_duration')), 'lc'))
for f in glob.glob(os.path.join(LCDIR, '307210830_s*_*.npz')):
    z = dict(np.load(f, allow_pickle=True))
    exp = 120. if f.endswith('_2m.npz') else 20.
    p = batman.TransitParams(); p.t0 = T0; p.per = P; p.rp = k; p.a = aR; p.inc = np.degrees(np.arccos(binj / aR)); p.ecc = 0; p.w = 90; p.limb_dark = 'quadratic'; p.u = u
    m = batman.TransitModel(p, z['t'], supersample_factor=7 if exp > 30 else 1, exp_time=exp / 86400. * 6 / 7).light_curve(p)   # midpoint sub-exposures
    z['pdc'] = z['pdc'] * m; z['sap'] = z['sap'] * m
    np.savez(os.path.join(outd, os.path.basename(f)), **z)
print('injected', T0)
