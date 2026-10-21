"""Wavelength bins (um) used for the chromatic test; identical edges for every visit of a detector. Own code."""
import numpy as np
def _eq(lo, hi, n): e = np.linspace(lo, hi, n + 1); return [(float(round(a, 4)), float(round(b, 4))) for a, b in zip(e[:-1], e[1:])]
BINS_UM = {'NRS1': _eq(2.87, 3.71, 8), 'NRS2': _eq(3.83, 5.15, 8),
           'SOSS1': [(0.90, 1.10), (1.10, 1.30), (1.30, 1.50), (1.50, 1.70), (1.70, 1.90), (1.90, 2.10), (2.10, 2.40), (2.40, 2.80)],
           'SOSS2': [(0.63, 0.75), (0.75, 0.85)]}
