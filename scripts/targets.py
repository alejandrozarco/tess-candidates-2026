"""Target definitions shared by the scripts in this folder.

Candidates: TIC 402710280 (Gaia DR3 5828441289812582400), TIC 324429952 (Gaia DR3 6416596855786566528),
TIC 89247735 (Gaia DR3 6708071581514074368) and TIC 455017950 (HD 134147, Gaia DR3 5825217968423162752).
Positions are Gaia DR3 (epoch 2016.0) unless noted; `tic_ra/tic_dec` are the TIC v8.2 coordinates used for the
TESScut requests (`cut_center='gaia'`: the cutout is requested at the Gaia DR3 position and the target pixel is the
proper-motion-corrected Gaia position at the median time of each sector). Ephemerides are the final linear
ephemerides of README.md; `tc_guess` are the starting values of the per-epoch fits (from the detection stage); `ncyc`
are the cycle numbers relative to T0. Transit epochs are labelled by their sector number; when a sector contains two
transits they are labelled '<sector>a', '<sector>b' (TIC 89247735, sector 93); `sector_of` returns the sector.
`diff` holds per-target settings of difference_images.py (defaults there when absent); `star` the TIC v8.2 stellar
parameters used by dilution_correction.py and gaia_companion_limits.py; `star_sed` (TIC 89247735, TIC 455017950) the
adopted stellar parameters of README sections 3.1 and 4.1 (SED fit, sed_fit.py), used by fit_free_dilution.py;
`free_dilution` the priors and starting point of fit_free_dilution.py (f_t prior, optional stellar-density prior).
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
        T14h=5.27, spoc_2min=None, ttv_min=30.0, ttv_note='observed O-C up to +-29 min',
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
        star_sed=dict(Teff=4430, eTeff=100, L=0.190, eL=0.008, R=0.74, eR=0.04, M=0.69, eM=0.03, rho=2.40, erho=0.40),   # README 3.1
        free_dilution=dict(ft=(0.735, 0.03), p0=(0.106, 0.2, 2.3, 0.4, 0.35), rho_prior=None),   # README 3.2
        diff=dict(hw=0.07, oot_min=0.13, oot_order=1, cen_radius=4.0, sigma0=0.8, sigma_min=0.3, seed=1, fake_excl=1.0, fake_min_frac=0.7,
                  gaia_radius=60.0, gaia_gmax=19.0),
        sources={'target': (275.2234460704073, -45.56163775820568, -39.073, -105.840),
                 'Gaia DR3 6708071581515491456 (G 15.49, 9.5")': (275.22077987, -45.56282528, 1.694, -6.031),
                 'Gaia DR3 6708071581515486720 (G 17.66, 12.4")': (275.22804007, -45.56279934, -0.457, 0.556),
                 'Gaia DR3 6708071207865931392 (G 14.84, 27.9")': (275.21905155, -45.56818791, -1.346, 0.273),
                 'Gaia DR3 6708071203558361728 (G 15.39, 27.8")': (275.21553903, -45.56637706, -15.105, -21.808)},
    ),
    455017950: dict(
        gaia_dr3='5825217968423162752',
        ra=227.9087992432127, dec=-64.9369219229812, pmra=-38.134, pmdec=-46.527,
        tic_ra=227.90881181, tic_dec=-64.93691547, cut_center='gaia',
        sectors=[12, 38, 39, 65, 99, 100, 101, 102],
        P=40.605660, eP=0.000074, T0=4055.48930, eT0=0.00179,
        ncyc={38: -42, 65: -24, 99: 0, 100: 1, 102: 2},
        tc_guess={s: 4055.48889 + n * 40.6056453 for s, n in {38: -42, 65: -24, 99: 0, 100: 1, 102: 2}.items()},
        epochs_fit=[38, 65, 99, 100, 102],
        u=[0.42, 0.20],                     # Teff ~5300 K, log g 4.6
        k0=0.035, b0=0.3, aR0=65.0, fix=None,
        T14h=4.74, spoc_2min=None, ttv_min=15.0, ttv_note='observed O-C up to +12 min',
        star=dict(R=0.794, eR=0.042, M=0.91, eM=0.11, rho=2.56, erho=0.58, d_pc=65.09),   # TIC v8.2
        star_sed=dict(Teff=5300, eTeff=60, L=0.455, eL=0.02, R=0.800, eR=0.024, M=0.88, eM=0.05, rho=2.42, erho=0.26),   # README 4.1
        free_dilution=dict(ft=(0.966, 0.010), p0=(0.0350, 0.3, 2.5, 0.38, 0.34), rho_prior=(2.40, 0.30)),   # README 4.2
        diff=dict(hw=0.07, oot_min=0.13, oot_order=2, cen_radius=4.0, sigma0=0.8, sigma_min=0.3, seed=1, fake_excl=1.0, fake_min_frac=0.7,
                  gaia_radius=60.0, gaia_gmax=19.0),
        sources={'target': (227.9087992432127, -64.9369219229812, -38.134, -46.527),
                 'Gaia DR3 5825217964086439552 (G 14.91, 10.0")': (227.91094627, -64.93428973, -6.808, -6.405),
                 'Gaia DR3 5825217968399119488 (G 17.77, 6.1")': (227.90566857, -64.93588642, -5.635, -9.088),
                 'Gaia DR3 5825217762264810752 (G 16.27, 27.5")': (227.91292922, -64.94436290, -5.035, -6.333),
                 'Gaia DR3 5825217762264718976 (G 14.92, 45.3")': (227.90302521, -64.94926701, -8.561, -5.465),
                 'Gaia DR3 5825218002782625024 (G 14.49, 51.4")': (227.90523379, -64.92272793, -6.597, -8.488),
                 'Gaia DR3 5825217796624167808 (G 13.91, 53.1")': (227.94087808, -64.94269542, -10.296, -10.438)},
    ),
    311179599: dict(
        gaia_dr3='5801882036668103680',
        ra=257.6466511014627, dec=-75.67187134421917, pmra=30.842, pmdec=-69.334,
        tic_ra=257.646094706738, tic_dec=-75.6715631372543, cut_center='gaia',
        sectors=[12, 13, 39, 66, 93, 100, 101, 102, 103],
        P=24.826859, eP=0.000044, T0=4082.33414, eT0=0.00154,
        ncyc={12: -98, 13: -97, 39: -69, 66: -39, 93: -10, 100: 0, 102: 2, 103: 3},
        tc_guess={12: 1649.3224, 13: 1674.1354, 39: 2369.1873, 66: 3114.0528, 93: 3834.0824, 100: 4082.3378, 102: 4131.9854, 103: 4156.8144},
        epochs_fit=[12, 13, 39, 66, 93, 100, 102, 103],
        ephem_exclude={39: 'sector 39 (-2.2 h); not used in the adopted ephemeris'},   # README 5.3; the all-epoch ephemeris is also written
        u=[0.33, 0.27],                     # Teff ~3950 K, log g 4.65
        k0=0.036, b0=0.3, aR0=60.0, fix=None,
        T14h=3.14, spoc_2min=None, ttv_min=60.0, ttv_note='O-C -50 to +28 min in the epochs of the ephemeris, -2.2 h in sector 39',
        star=dict(R=0.599, eR=0.018, M=0.587, eM=0.020, rho=3.86, erho=0.21, d_pc=47.49),   # TIC v8.2
        star_sed=dict(R=0.599, eR=0.018, M=0.587, eM=0.012, rho=3.85, erho=0.35),   # README 5.1: Mann et al. (2015, 2019) M_Ks relations
        free_dilution=dict(ft=(0.955, 0.010), p0=(0.038, 0.6, 3.9, 0.36, 0.30), rho_prior=(3.85, 0.35)),   # README 5.2
        diff=dict(hw=0.05, oot_min=0.10, oot_order=2, cen_radius=4.0, sigma0=0.8, sigma_min=0.3, seed=1, fake_excl=1.0, fake_min_frac=0.7,
                  gaia_radius=60.0, gaia_gmax=19.0),
        sources={'target': (257.6466511014627, -75.67187134421917, 30.842, -69.334),
                 'Gaia DR3 5801882032368601856 (G 17.77, 16.8")': (257.66259555, -75.66938311, -10.084, -3.626),
                 'Gaia DR3 5801883514136854528 (G 17.43, 26.5")': (257.65336001, -75.66470274, -2.160, 1.523),
                 'Gaia DR3 5801882002308365184 (G 15.99, 27.9")': (257.63075376, -75.67855605, 1.488, -10.157),
                 'Gaia DR3 5801883509837351936 (G 16.90, 32.1")': (257.67276323, -75.66573075, -1.049, -7.239)},
        # candidate second transit signal (README section 5.10): three-event linear ephemerides for the two allowed
        # periods (n = 1: P = (t_B - t_A); n = 2: half of it) and the joint-fit duration (fit_second_signal.py, r = 2.0 px,
        # data/transit_times_311179599c.csv); used by predict_transits.py --signal c
        # `events` = (sector, starting mid-time) of the three events; `mask_b` = mid-times of the 24.8-d transits (and the
        # sector-101 prediction) masked in the fit; fit_second_signal.py
        second=dict(name='c', T14h=5.38, T0=4088.51641, eT0=0.00664, P_ref=73.87,
                    events={'A': (100, 4088.5007), 'B': (103, 4162.4007), 'C': (12, 1650.8937)},
                    camera_ccd={'A': '4/3', 'B': '2/1', 'C': '2/2'},
                    mask_b=[1649.32236, 1674.13540, 2369.18732, 3114.05281, 3834.08235, 4082.33775, 4107.16, 4131.98535, 4156.81437],
                    aliases=[(1, 73.867214, 0.000455), (2, 36.933607, 0.000228)],
                    ttv_note='O-C -25 and +26 min of the sector-100 and sector-103 events'),
    ),
}

# Other TESS results (README section 6); used by eb_*.py and toi540_*.py only.
# 11613676: eclipsing binary CD-30 448. Cutouts requested at the Gaia DR3 position; target pixel = proper-motion-
#   corrected Gaia position at the sector median time. `eph0` = starting ephemeris of the per-eclipse fits (detection
#   stage), `ph2` = starting phase of the secondary eclipse; `neighbours` = the six Gaia DR3 sources nearest to the
#   target (30-97"; ra, dec, G, BP-RP) for the contamination estimate.
# 200322593: TOI-540 (Ment et al. 2021). `star` = mass and radius of Ment et al. (2021); `toi` = ExoFOP TOI-table
#   ephemeris (2026-10-18; starting values only); `controls` = ExoFOP TOI-table values (2026-10-18) of the control
#   planets observed in sector 98 at 2-min cadence: (TOI, TIC, P, eP, T0 BJD_TDB, eT0, T14 h, depth ppm).
# 124751008: eclipsing binary 2MASS J06162084-3509125 (README section 6.4). Cutouts requested at the Gaia DR3
#   position; `eph0` = starting ephemeris of the per-eclipse fits (detection stage), `ph2` = phase of the secondary
#   eclipse (2019); `bleed` = image level (e-/s) above which pixels belong to the saturation bleed column of the G = 4.1
#   star kappa Col (158"); `neighbours` = Gaia DR3 sources within 60" (ra, dec, G, BP-RP) for the contamination estimate.
OTHER = {
    11613676: dict(
        gaia_dr3='5029330944688623360', name='CD-30 448',
        ra=20.632566548307413, dec=-29.43755526280964, pmra=6.2715, pmdec=-32.8469, cut_center='gaia',
        sectors=[3, 30, 97, 104, 106], spoc_2min=[3, 30],
        eph0=(3979.60469, 12.115919), ph2=0.557,
        neighbours=[(20.642006860, -29.435474122, 21.077, 2.393), (20.622196859, -29.437330474, 17.128, 1.069),
                    (20.652595477, -29.442875136, 20.288, 0.721), (20.618010897, -29.454294068, 20.570, 2.161),
                    (20.646714805, -29.455664264, 15.626, 0.733), (20.632790522, -29.464228349, 18.994, 0.921)],
    ),
    124751008: dict(
        gaia_dr3='2885572716963709952', name='2MASS J06162084-3509125',
        ra=94.08686139415221, dec=-35.15313837719891, pmra=-3.1427, pmdec=83.5702, cut_center='gaia',
        sectors=[6, 7, 33, 87, 98, 108, 109], spoc_2min=[],
        eph0=(3997.3156, 12.61176), ph2=0.4986, bleed=20000.0, G=11.191, bprp=2.111,
        neighbours=[(94.08236243, -35.15061433, 17.964, 0.805), (94.09118199, -35.15706361, 16.210, 1.335),
                    (94.08211750, -35.14919168, 19.720, 2.328), (94.09254373, -35.14365280, 21.080, 0.613),
                    (94.08443369, -35.16439934, 20.750, 1.558), (94.09932185, -35.15799980, 15.300, 0.864),
                    (94.07333094, -35.15482260, 19.150, 1.367), (94.07517045, -35.16153183, 20.060, 1.876),
                    (94.09030819, -35.13828098, 16.700, 2.346)],
    ),
    200322593: dict(
        gaia_dr3='4785886941312921344', name='TOI-540',
        ra=76.30968124678876, dec=-47.937535542918475, pmra=-65.901, pmdec=24.895,
        sectors=[4, 5, 6, 31, 32, 98], pre=[4, 5, 6, 31, 32],
        star=dict(M=0.159, eM=0.014, R=0.1895, eR=0.0079), u=[0.25, 0.35],
        toi=(2459175.141647, 1.23915336),
        controls=[('269.01', 220479565, 3.69771670923474, 2.2476602e-06, 2458381.845549, 0.00024341463, 0.930826340882136, 4755.67),
                  ('175.02', 307210830, 7.4507305, 3.2e-06, 2460039.154503, 0.0004444, 0.475, 1542.0),
                  ('175.01', 307210830, 3.6906758578021837, 1.008001358872935e-06, 2460688.70837, 0.000380081725349, 1.107, 1657.0),
                  ('696.01', 77156829, 0.860240578952471, 1.0231897e-05, 2458411.473627, 0.00040476615, 0.512434461523044, 719.54),
                  ('448.01', 1528696, 0.882062425827087, 5.160957e-06, 2458438.402842, 0.00021536286, 0.575967222552898, 11755.45),
                  ('206.01', 55650590, 0.736313219931196, 4.7654476e-06, 2458325.542536, 0.0002835781, 0.673685332264972, 1322.55),
                  ('4508.01', 170789802, 3.46727789358939, 9.109895e-06, 2458438.012099, 0.00049501995, 0.405592768803134, 3657.06),
                  ('2224.01', 388198242, 0.505510515897485, 2.8973066e-06, 2459988.772614, 0.0012692922, 0.648232028256928, 2269.20),
                  ('4981.01', 14029793, 0.2097217, 4e-07, 2459199.878456, 0.0008661, 0.25, 1250.0),
                  ('270.03', 259377017, 3.3601670547055087, 3.522546590899222e-06, 2460987.858405, 0.002242517381115937, 1.3266673034451342, 976.69)],
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
