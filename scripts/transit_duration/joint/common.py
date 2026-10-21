"""Shared definitions for the joint analysis of L 98-59 d (section 6.6.10). Own code.
All inputs and outputs live under ROOT = $TDJ_WORK (default work/transit_duration_joint, relative to the repository
root); outputs go to ROOT/out."""
import numpy as np, glob, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('TDJ_WORK', os.path.join('work', 'transit_duration_joint'))
TIC = 307210830
# Cadieux et al. 2025 (AJ 170, 154; arXiv:2507.09343) linear ephemerides, BTJD; T14 in hours
PL = {'b': (2.2531140, 1366.1706, 1.01), 'c': (3.6906764, 1367.2730, 1.28), 'd': (7.450729, 1362.7400, 0.84)}
PTTV = 396.32   # c-d super-period, Cadieux et al. 2025
# stellar density prior: ln rho[SI] ~ N(9.50, 0.09) (Cadieux et al. 2025, from Mann et al. 2015/2019 M, R)
LNRHO = (9.50, 0.09)
G = 6.674e-11
def aR_from_lnrho(lnrho, P):
    return (G * np.exp(lnrho) * (P * 86400.) ** 2 / (3 * np.pi)) ** (1 / 3.)
# Claret (2017) TESS quadratic LD: out/claret_ld.json (copied from scripts/transit_duration/joint/claret_ld.json by run_joint.sh)
def year(btjd):
    return 2000.0 + (btjd + 2457000.0 - 2451544.5) / 365.25
def campaign(sector):
    return 0 if sector < 20 else 1 if sector < 50 else 2 if sector < 80 else 3
def load_lc(cad='2m', lcdir=None):
    lcdir = lcdir or os.path.join(ROOT, 'out', 'lc')
    out = []
    for f in sorted(glob.glob(os.path.join(lcdir, '%d_s*_%s.npz' % (TIC, cad)))):
        z = np.load(f, allow_pickle=True)
        d = {k: z[k] for k in z.files if k != 'meta'}; d['meta'] = z['meta'][0]
        out.append(d)
    return out
