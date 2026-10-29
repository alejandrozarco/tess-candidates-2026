"""TOI-421 b (TIC 94986319), JWST visit 1 (2023-09-10, NIRSpec G395M, NRS1, program 1935).
Convert the Davenport+2025 Eureka! Stage-4 white-light light curve (Zenodo 10.5281/zenodo.15176438,
TOI421b/NIRSpec/Stage4/White_Light/S4_toi421b_ap3_bg8_LCData.h5; 2.87-5.168 um) into a plain npz
for fit_wlc.py. Reads the h5 as data only. Time: BMJD_TDB (Eureka! convention, mid-integration).
Usage: python -I prep_eureka_nrs1.py <LCData.h5> <out.npz>
"""
import sys
import numpy as np
import h5py

f = h5py.File(sys.argv[1], "r")
t = f["time"][:]
fl = f["flux_white"][:]
er = f["err_white"][:]
m = f["mask_white"][:].astype(bool)
ok = (~m) & np.isfinite(fl) & np.isfinite(er)
t, fl, er = t[ok], fl[ok], er[ok]
med = np.median(fl)
dt = float(np.median(np.diff(t)))
meta = ("planet=TOI-421 b (TIC 94986319); visit=2023-09-10 JWST NIRSpec G395M (program 1935, obs 2); "
        "detector=NRS1; band=2.87-5.168 um white light; source=Davenport+2025 Eureka! S4 white-light LC "
        "(Zenodo 15176438, S4_toi421b_ap3_bg8_LCData.h5, ap3_bg8); time=BMJD_TDB mid-integration; "
        "flux normalised by median; trimming=none")
np.savez(sys.argv[2], t=t, f=fl / med, e=er / med, exptime=dt, meta=meta,
         label="V1_NIRSpecG395M_NRS1_eureka", ld_band="NIRSpec_G395M_NRS1_2.87-5.17um")
print(len(t), "points; cadence", dt * 86400, "s; span", (t[-1] - t[0]) * 24, "h")
