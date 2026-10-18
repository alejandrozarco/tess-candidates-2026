"""Target definitions shared by the scripts in this folder.

Candidates: TIC 402710280 (Gaia DR3 5828441289812582400), TIC 324429952 (Gaia DR3 6416596855786566528) and
TIC 89247735 (Gaia DR3 6708071581514074368).
Positions are Gaia DR3 (epoch 2016.0) unless noted; `tic_ra/tic_dec` are the TIC v8.2 coordinates used for the
TESScut requests (`cut_center='gaia'`: the cutout is requested at the Gaia DR3 position and the target pixel is the
proper-motion-corrected Gaia position at the median time of each sector). Ephemerides are the final linear
ephemerides of README.md; `tc_guess` are the starting values of the per-epoch fits (from the detection stage); `ncyc`
are the cycle numbers relative to T0. Transit epochs are labelled by their sector number; when a sector contains two
transits they are labelled '<sector>a', '<sector>b' (TIC 89247735, sector 93); `sector_of` returns the sector.
`diff` holds per-target settings of difference_images.py (defaults there when absent); `star` the TIC v8.2 stellar
parameters used by dilution_correction.py and gaia_companion_limits.py.
Work directory layout (set with --work, default ./work): <work>/<tic>/cut/s<sector>.npz (cutouts),
<work>/<tic>/cut/lc_s<sector>.npz (photometry), <work>/<tic>/dl/ (MAST products), <work>/<tic>/fit_*.json (fits).
"""
import os, re

TARGETS = {
    402710280: dict(
        gaia_dr3='5828441289812582400',
        ra=248.5242026111537, dec=-63.54961908541294, pmra=59.909, pmdec=13.685,
        tic_ra=248.52360473, tic_dec=-63.54967989,
        sectors=[12, 39, 66, 93, 100, 101, 102, 103],
        P=36.34005, eP=0.00006, T0=4128.92695, eT0=0.00140,
        ncyc={39: -48, 66: -28, 93: -8, 100: -1, 102: 0},
        tc_guess={39: 2384.6092, 66: 3111.3942, 93: 3838.1853, 100: 4092.5841, 102: 4128.9448},
        epochs_fit=[39, 93, 100, 102],      # S66 is in scattered light and is fitted separately
        u=[0.46, 0.19],                     # quadratic limb darkening, TESS band, Teff ~5100 K, log g 4.5
        k0=0.043, b0=0.3, aR0=50.0, fix=None,
        T14h=5.27, spoc_2min=None, ttv_min=30.0,
        sources={'target': (248.5242026111537, -63.54961908541294, 59.909, 13.685),
                 'Gaia DR3 5828441289812582656 (G 16.15)': (248.52270239, -63.55016350, 0, 0),
                 'G 14.46 (39")': (248.53744129, -63.55879284, 0, 0),
                 'G 12.12 (54")': (248.55692634, -63.55335536, 0, 0),
                 'G 12.38 (60")': (248.55824678, -63.54260739, 0, 0)},
    ),
    324429952: dict(
        gaia_dr3='6416596855786566528',
        ra=285.95905536343747, dec=-73.11401715610884, pmra=25.209, pmdec=-66.783,
        tic_ra=285.95867, tic_dec=-73.11372,
        sectors=[13, 27, 39, 66, 67, 93, 94, 100, 101, 102, 103],
        P=42.37574, eP=0.00008, T0=4174.88770, eT0=0.00290,
        ncyc={13: -59, 27: -50, 66: -25, 93: -8, 94: -7, 100: -2, 102: -1, 103: 0},
        tc_guess={s: 4174.89100 + n * 42.375818 for s, n in
                  {13: -59, 27: -50, 66: -25, 93: -8, 94: -7, 100: -2, 102: -1, 103: 0}.items()},
        epochs_fit=[13, 27, 66, 93, 94, 100, 102, 103],
        u=[0.40, 0.21],                     # Teff ~5750 K, log g 4.5
        k0=0.038, b0=0.6, aR0=50.0, fix=(0.0386, 0.80, 55.9),   # shape fixed at the stellar-density solution (README)
        T14h=3.84, spoc_2min=13, ttv_min=None,
        sources={'target': (285.95905536, -73.11401716, 25.209, -66.783),
                 'Gaia DR3 6416597615998643200 (M dwarf, G 15.17, 43.6")': (285.96930769, -73.10227376, 24.36, -65.56),
                 'G 15.15 (52")': (286.00427087, -73.10800433, 0, 0),
                 'G 15.07 (77")': (285.93786990, -73.13440608, 0, 0),
                 'G 14.33 (86")': (285.95625500, -73.09025762, 0, 0)},
    ),
    89247735: dict(
        gaia_dr3='6708071581514074368',
        ra=275.2234460704073, dec=-45.56163775820568, pmra=-39.073, pmdec=-105.840,
        tic_ra=275.22369241, tic_dec=-45.56116720, cut_center='gaia',
        sectors=[13, 66, 93, 104],
        P=18.323173, eP=0.000031, T0=4196.80899, eT0=0.00100,
        ncyc={66: -59, '93a': -20, '93b': -19, 104: 0},
        tc_guess={k: 4196.8100 + n * 18.323343 for k, n in {66: -59, '93a': -20, '93b': -19, 104: 0}.items()},
        epochs_fit=[66, '93a', '93b', 104],
        u=[0.52, 0.18],                     # Teff ~4400 K, log g 4.5
        k0=0.1, b0=0.3, aR0=30.0, fix=None,
        T14h=4.50, spoc_2min=None, ttv_min=None,
        star=dict(R=0.785, eR=0.066, M=0.69, eM=0.08, rho=1.43, erho=0.43, d_pc=105.66),   # TIC v8.2
        diff=dict(hw=0.07, oot_min=0.13, oot_order=1, cen_radius=4.0, sigma0=0.8, sigma_min=0.3, seed=1, fake_excl=1.0, fake_min_frac=0.7,
                  gaia_radius=60.0, gaia_gmax=19.0),
        sources={'target': (275.2234460704073, -45.56163775820568, -39.073, -105.840),
                 'Gaia DR3 6708071581515491456 (G 15.49, 9.5")': (275.22077987, -45.56282528, 1.694, -6.031),
                 'Gaia DR3 6708071581515486720 (G 17.66, 12.4")': (275.22804007, -45.56279934, -0.457, 0.556),
                 'Gaia DR3 6708071207865931392 (G 14.84, 27.9")': (275.21905155, -45.56818791, -1.346, 0.273),
                 'Gaia DR3 6708071203558361728 (G 15.39, 27.8")': (275.21553903, -45.56637706, -15.105, -21.808)},
    ),
}


def lab(x):
    """epoch label: int for a plain sector number, str for '<sector>a/b'"""
    return int(x) if str(x).isdigit() else str(x)


def sector_of(label):
    return int(re.match(r'\d+', str(label)).group())


def workdir(args, tic):
    d = os.path.join(args.work, str(tic))
    os.makedirs(os.path.join(d, 'cut'), exist_ok=True)
    return d


def add_common_args(p):
    p.add_argument('tic', type=int, choices=sorted(TARGETS))
    p.add_argument('--work', default='./work', help='work directory (default ./work)')
    return p


def import_batman():
    """batman-package 2.5.x imports distutils, which Python 3.12 no longer ships; provide a stub."""
    import sys, types
    if 'distutils' not in sys.modules:
        d = types.ModuleType('distutils'); c = types.ModuleType('distutils.ccompiler')
        c.new_compiler = lambda *a, **k: None; d.ccompiler = c
        sys.modules['distutils'] = d; sys.modules['distutils.ccompiler'] = c
    import batman
    return batman
