"""TOI-421 b/c (TIC 94986319): nodal precession of the b-c pair in short N-body (REBOUND WHFast) runs.
Line of sight = +z; sky plane = x-y. Impact parameter b = (r_transit/R*) |cos i| with i the inclination to the z axis.
Inputs: M* 0.852 Msun, R* 0.871 Rsun (Carleo+2020); m_b 7.17, m_c 16.42 M_earth (Carleo+2020; mass range scanned);
P_b 5.19760 d, P_c 16.06750 d; e_b 0.163 (w 129 deg; Davenport+2025 SOSS config), e_c 0.15 (w free, 90 deg here; e has
little effect on the node rate); sky inclinations from the TESS 2018-20 b values (b_b 0.944 at a/R* 13.76; b_c 0.766 at
a/R* 31.91); the sky-plane node difference dOmega is set so that the mutual inclination equals i_bc (scan).
Output: d(cos i)/dt for both planets (secular, linear fit over 40 yr), the precession period, and the L_c/L_b check.
usage: python scripts/toi421/nbody_toi421.py  -> work/toi421/out/nodal_nbody.json"""
import json, os, numpy as np, rebound
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "work", "toi421", "out")

MS, RS = 0.852, 0.871
AU_RS = 215.032
ME = 3.003e-6
PB, PC = 5.19760 / 365.25, 16.06750 / 365.25


def run(i_bc_deg, mb=7.17, mc=16.42, same_hemi=True, sgn=+1, years=40., eb=0.163, ec=0.15, wb=129., wc=90.):
    ab = (MS * PB ** 2) ** (1 / 3.); ac = (MS * PC ** 2) ** (1 / 3.)
    cos_ic = 0.766 / 31.91; cos_ib = 0.944 / 13.76 * (1 if same_hemi else -1)
    ic, ib = np.arccos(cos_ic), np.arccos(cos_ib)
    cd = (np.cos(np.radians(i_bc_deg)) - np.cos(ib) * np.cos(ic)) / (np.sin(ib) * np.sin(ic))
    if abs(cd) > 1: return None
    dO = sgn * np.arccos(cd)
    sim = rebound.Simulation(); sim.units = ('yr', 'AU', 'Msun'); sim.G = 4 * np.pi ** 2
    sim.add(m=MS)
    sim.add(m=mb * ME, P=PB, e=eb, inc=ib, Omega=0.0, omega=np.radians(wb) - 0.0, M=0.3)
    sim.add(m=mc * ME, P=PC, e=ec, inc=ic, Omega=dO, omega=np.radians(wc), M=2.1)
    sim.move_to_com(); sim.integrator = 'whfast'; sim.dt = PB / 40.
    ts = np.linspace(0, years, 801); cb, cc, nb, nc = [], [], [], []
    for t in ts:
        sim.integrate(t, exact_finish_time=0)
        o = sim.particles[1].orbit(primary=sim.particles[0]); p = sim.particles[2].orbit(primary=sim.particles[0])
        cb.append(np.cos(o.inc)); cc.append(np.cos(p.inc))
    cb, cc = np.array(cb), np.array(cc)
    sb = np.polyfit(ts, cb, 1)[0]; sc = np.polyfit(ts, cc, 1)[0]
    # mutual inclination check
    Lb = mb * np.sqrt(ab); Lc = mc * np.sqrt(ac)
    # precession frequency from the relative node rate over a longer run is not needed: |dn/dt| = w_p sin(theta)
    return dict(i_bc=i_bc_deg, dOmega_deg=float(np.degrees(dO)), same_hemi=same_hemi, sgn=sgn, mb=mb, mc=mc,
                dcosi_b_dt=float(sb), dcosi_c_dt=float(sc), ratio=float(-sb / sc) if sc != 0 else None, Lc_over_Lb=float(Lc / Lb),
                # b-rate in the measurement conventions (circular; a/R* 13.76 for b, 31.91 for c) and physical a/R*
                db_b_dt_conv=float(13.76 * sb * np.sign(cos_ib)), db_c_dt_conv=float(31.91 * sc),
                db_c_dt_phys=float(ac * AU_RS / RS * sc), aR_c_phys=float(ac * AU_RS / RS), aR_b_phys=float(ab * AU_RS / RS))


def precession_period(i_bc_deg, mb=7.17, mc=16.42):
    """period of the relative node (Omega_c - Omega_b) in the invariable frame: long run, sample, fit."""
    sim = rebound.Simulation(); sim.units = ('yr', 'AU', 'Msun'); sim.G = 4 * np.pi ** 2
    sim.add(m=MS); sim.add(m=mb * ME, P=PB, e=0.163, inc=0.0, M=0.3)
    sim.add(m=mc * ME, P=PC, e=0.15, inc=np.radians(i_bc_deg), Omega=0.0, M=2.1)
    sim.move_to_com(); sim.integrator = 'whfast'; sim.dt = PB / 30.
    # rotate to invariable plane
    L = np.zeros(3)
    for p in sim.particles[1:]:
        r = np.array(p.xyz); v = np.array(p.vxyz); L += p.m * np.cross(r, v)
    Lh = L / np.linalg.norm(L)
    ts = np.linspace(0, 400, 401); ang = []
    for t in ts:
        sim.integrate(t, exact_finish_time=0)
        p = sim.particles[2]; r = np.array(p.xyz) - np.array(sim.particles[0].xyz); v = np.array(p.vxyz) - np.array(sim.particles[0].vxyz)
        n = np.cross(r, v); n /= np.linalg.norm(n)
        # azimuth of n around Lh
        e1 = np.cross(Lh, [0, 0, 1.0]); e1 = e1 / np.linalg.norm(e1) if np.linalg.norm(e1) > 1e-9 else np.array([1., 0, 0]); e2 = np.cross(Lh, e1)
        ang.append(np.arctan2(n @ e2, n @ e1))
    ang = np.unwrap(ang); rate = np.polyfit(ts, ang, 1)[0]
    return float(2 * np.pi / abs(rate)), float(rate)


if __name__ == '__main__':
    out = dict(meta=__doc__, runs=[], periods={})
    for ibc in (5, 8, 10, 12, 15, 20, 25, 30, 35, 42, 49):
        out['periods'][str(ibc)] = precession_period(ibc)
        for hemi in (True, False):
            for sgn in (+1, -1):
                r = run(ibc, same_hemi=hemi, sgn=sgn)
                if r: out['runs'].append(r)
        print(ibc, out['periods'][str(ibc)], [(r['same_hemi'], r['sgn'], round(r['db_c_dt_conv'], 4), round(r['db_b_dt_conv'], 4)) for r in out['runs'] if r['i_bc'] == ibc], flush=True)
    # mass range
    out['mass_scan'] = []
    for mb, mc in ((6.4, 14.0), (7.17, 16.42), (7.8, 17.5)):
        for ibc in (10, 12, 15, 35):
            r = run(ibc, mb=mb, mc=mc, same_hemi=True, sgn=+1); out['mass_scan'].append(r)
    json.dump(out, open(os.path.join(OUT, 'nodal_nbody.json'), 'w'), indent=1)
