"""T12: fit the real d transits in an alternative light-curve set (own TPF apertures, SAP, 20-s, EB modelled) with the
same chain (windows, prepare, const/lin_b least squares); optional pointing/background regressors in the baseline.
Usage: t12_ownphot.py <label> <cad> <flux pdc|sap> <eb mask|model> <lcdir or -> <regs 0|1> <out.json>"""
import sys, json, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from t5lib import *; from t3_windows import make_windows
lab, cad, fl, eb, lcdir, regs, outp = sys.argv[1:8]
segs, flares = make_windows(load_lc(cad, None if lcdir == '-' else lcdir), fl, eb)
D = prepare(segs, regs=bool(int(regs))); r = fit_pair(D)
r.update(label=lab, ntr=len(segs), npts=int(D.npts))
json.dump(r, open(outp, 'w'), indent=1)
print(lab, 'ntr %d  db/dt %.5f +- %.5f  dT14/dt %.3f +- %.3f min/yr  dchi2 %.1f  k %.4f b0 %.4f lnrho %.3f' % (len(segs), r['bdot'], r['bdot_err'], r['dT14'], r['dT14_err'], r['dchi2'], r['k'], r['b0'], r['lnrho']))
