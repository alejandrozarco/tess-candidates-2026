# TOI-421 b (TIC 94986319) JWST white-light theoretical quadratic LD.
# Star: Teff 5291 K, logg 4.48, [Fe/H] -0.02. Source: ExoTiC-LD 3.2.0 (Grant & Wakeford 2024), on-the-fly grids.
# usage: python -I ld_theory.py <exotic_ld data directory> <out.json>
import sys, json, numpy as np
from exotic_ld import StellarLimbDarkening
out = {"planet": "TOI-421 b (TIC 94986319)", "star": {"Teff": 5291, "logg": 4.48, "M_H": -0.02},
       "source": "ExoTiC-LD 3.2.0, quadratic law, trilinear interpolation", "bands": {}}
bands = {"NIRSpec_G395M_NRS1_2.87-5.17um": ("JWST_NIRSpec_G395M", [28700., 51680.]),
         "NIRISS_SOSS_order1_0.85-2.83um": ("JWST_NIRISS_SOSSo1", [8500., 28300.]),
         "TESS": ("TESS", [5800., 11000.])}
for model in ["mps1", "stagger"]:
    sld = StellarLimbDarkening(M_H=-0.02, Teff=5291, logg=4.48, ld_model=model,
                               ld_data_path=sys.argv[1], interpolate_type="trilinear", verbose=0)
    for name, (mode, wr) in bands.items():
        u = sld.compute_quadratic_ld_coeffs(wavelength_range=wr, mode=mode)
        out["bands"].setdefault(name, {})[model] = [float(u[0]), float(u[1])]
        print(model, name, u)
json.dump(out, open(sys.argv[2], "w"), indent=1)
