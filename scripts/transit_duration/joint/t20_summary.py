"""T20: collect all outputs of the joint analysis into out/summary.json and print the tables of section 6.6.10."""
import json, glob, os, numpy as np, sys
from scipy import stats
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from common import ROOT
O = ROOT + '/out'; S = {}
def J(p):
    try: return json.load(open(p))
    except Exception: return None
# --- TESS MCMC variants
mc = {}
for f in sorted(glob.glob(O + '/mcmc/*.json')):
    d = J(f); n = os.path.basename(f)[:-5]
    if d is None: continue
    e = dict(chi2=d['chi2_map'], bic=d['bic'], lnZ=d['lnZ_is'], ess=d['is_ess'], lnZ_lap=d['lnZ_laplace'], npar=len(d['names']),
             rhat_max=max(d['rhat'].values()) if d['rhat'] else None, tau_max=max(max(t) for t in d['tau']), nsteps=d['nsteps'], nens=d['nens'],
             n_over_tau=(d['nsteps'] - d['nsteps'] // 3) / max(max(t) for t in d['tau']), post=d['post'], derived=d['derived'])
    mc[n] = e
S['mcmc'] = mc
print('%-18s %5s %10s %9s %9s %6s %6s %7s   %-24s %-22s' % ('variant', 'npar', 'chi2', 'BIC', 'lnZ', 'Rhat', 'N/tau', 'ESS_IS', 'dT14/dt (min/yr)', 'db/dt (/yr)'))
for n, e in mc.items():
    dt = e['derived'].get('dT14dt_min_per_yr'); db = e['derived'].get('dbdt_per_yr')
    fm = lambda q, f: (f % q[0]) + ' +%s/-%s' % ((f % (q[2] - q[0])).strip(), (f % (q[0] - q[1])).strip()) if q else '-'
    print('%-18s %5d %10.2f %9.2f %9.2f %6.3f %6.0f %7.0f   %-24s %-22s' % (n, e['npar'], e['chi2'], e['bic'], e['lnZ'], e['rhat_max'] or 0, e['n_over_tau'], e['ess'], fm(dt, '%.3f'), fm(db, '%.4f')))
# model comparisons (evolving vs constant), prior on bdot uniform width 0.1/yr -> +ln(10); on idot width 1 deg/yr -> +ln(1); the i0 prior
# width of lin_i is not normalised; key lnBF = profile-evidence log contrast (baselines profiled, not marginalised)
pairs = [('main_lin_b', 'main_const', np.log(10)), ('v_lin_b_kcamp', 'v_const_kcamp', np.log(10)), ('v_lin_b_ldfree', 'v_const_ldfree', np.log(10)),
         ('v_lin_b_ecc', 'v_const_ecc', np.log(10)), ('v_lin_b_fixtc', 'v_const_fixtc', np.log(10)), ('v_lin_i', 'main_const', 0.0),
         ('v_camp_slope', 'v_camp', np.log(10)), ('v_camp', 'main_lin_b', None), ('v_lin_b_rhofree', 'v_const_rhofree', np.log(10))]
cmp_ = {}
print('\nmodel comparison: dchi2, dBIC (positive favours first), profile-evidence log contrast (importance sampling of the profiled likelihood; bdot prior width added; not a normalised Bayes factor)')
for a, b, lnw in pairs:
    if a in mc and b in mc:
        dc = mc[b]['chi2'] - mc[a]['chi2']; dB = mc[b]['bic'] - mc[a]['bic']; dk = mc[a]['npar'] - mc[b]['npar']
        lnBF = (mc[a]['lnZ'] + lnw - mc[b]['lnZ']) if lnw is not None else None
        nested = (a, b) != ('v_camp', 'main_lin_b')   # one b per campaign vs b linear are not nested: no LR calibration
        p = stats.chi2.sf(dc, dk) if nested and dk > 0 and dc > 0 else None
        z = stats.norm.isf(p / 2) if p else None
        cmp_[a + '_vs_' + b] = dict(dchi2=dc, ddof=dk, dBIC=dB, lnBF=lnBF, p_LR=p, z_LR=z)
        print('  %-16s vs %-16s dchi2 %7.2f (%d dof, LR z %s)  dBIC %7.2f  profile-evidence log contrast %s' % (a, b, dc, dk, '%.1f' % z if z else '-', dB, '%.1f' % lnBF if lnBF is not None else '-'))
S['comparisons'] = cmp_
# --- bootstraps
bt = {}
for m in ('transit', 'sector', 'jack'):
    f = O + '/boot_%s.npz' % m
    if os.path.exists(f):
        z = np.load(f, allow_pickle=True); bd, dt = z['bdot'], z['dT14']
        bt[m] = dict(n=len(bd), bdot_med=float(np.median(bd)), bdot_sd=float(np.std(bd)), dT14_med=float(np.median(dt)), dT14_sd=float(np.std(dt)),
                     frac_dT14_ge0=float(np.mean(dt >= 0)), labels=[str(x) for x in z['label']] if m == 'jack' else None, dT14_all=dt.tolist() if m == 'jack' else None)
        print('bootstrap %-8s n=%d  dT14/dt %.3f +- %.3f  db/dt %.5f +- %.5f  frac(dT14>=0) %.4f' % (m, len(bd), np.median(dt), np.std(dt), np.median(bd), np.std(bd), np.mean(dt >= 0)))
S['bootstrap'] = bt
# --- injections
for m in ('lc', 'pix2'):
    f = O + '/inj_%s.jsonl' % m
    if os.path.exists(f):
        R = [json.loads(l) for l in open(f)]
        dt = np.array([r['dT14'] for r in R]); de = np.array([r['dT14_err'] for r in R]); zz = dt / de; dc = np.array([r['dchi2'] for r in R])
        S['inj_' + m] = dict(n=len(R), dT14_mean=float(dt.mean()), dT14_sd=float(dt.std()), z_mean=float(zz.mean()), z_sd=float(zz.std()), dchi2_med=float(np.median(dc)),
                             dchi2_max=float(dc.max()), frac_dchi2_gt9=float(np.mean(dc > 9)), frac_z_lt_m3=float(np.mean(zz < -3)), z_min=float(zz.min()), z_max=float(zz.max()), ks_dchi2_vs_chi2_1=float(stats.kstest(dc, 'chi2', args=(1,)).pvalue))
        print('injections %-3s n=%d  dT14/dt mean %.3f sd %.3f min/yr; z mean %.2f sd %.2f (min %.2f max %.2f); dchi2 median %.2f max %.1f; frac dchi2>9 %.3f; KS vs chi2(1) p=%.3f' % (
            m, len(R), dt.mean(), dt.std(), zz.mean(), zz.std(), zz.min(), zz.max(), np.median(dc), dc.max(), np.mean(dc > 9), S['inj_' + m]['ks_dchi2_vs_chi2_1']))
# --- alternative data sets
alt = {}
for f in sorted(glob.glob(O + '/alt_*.json')):
    d = J(f); alt[d['label']] = {k: d[k] for k in ('ntr', 'bdot', 'bdot_err', 'dT14', 'dT14_err', 'dchi2', 'k', 'b0', 'lnrho')}
    print('alt %-18s ntr %2d dT14/dt %.3f +- %.3f  db/dt %.5f +- %.5f  dchi2 %.1f  k %.4f' % (d['label'], d['ntr'], d['dT14'], d['dT14_err'], d['bdot'], d['bdot_err'], d['dchi2'], d['k']))
S['alt'] = alt
S['audit'] = J(O + '/audit_2m_pdc.json')
# --- JWST
jw = {}
for f in sorted(glob.glob(O + '/jwst/fits/*.json')):
    d = J(f); n = os.path.basename(f)[:-5]
    if d is None or 'T14' not in d: continue
    jw[n] = dict(T14=d['T14'], b=d['b'], aR=d['aR'], dT14=d.get('dT14'), db=d.get('db'), chi2=d['chi2_map'], npts=d['npts'], npar=len(d['names']), noise=d['noise'], n_over_tau=d['n_over_tau'], rhat_max=max(d['rhat'].values()))
    s = '  '.join('%s %.2f+-%.2f' % (v, q[0], q[1]) for v, q in d['T14'].items())
    print('jwst %-16s %s  aR %.2f+-%.2f %s chi2 %.1f N/tau %.0f Rhat %.3f' % (n, s, d['aR'][0], d['aR'][1], ('dT14 %.3f+-%.3f db %.4f+-%.4f' % (d['dT14'][0], d['dT14'][1], d['db'][0], d['db'][1])) if d.get('dT14') else '', d['chi2_map'], d['n_over_tau'], max(d['rhat'].values())))
S['jwst'] = jw
json.dump(S, open(O + '/summary.json', 'w'), indent=1, default=float)
