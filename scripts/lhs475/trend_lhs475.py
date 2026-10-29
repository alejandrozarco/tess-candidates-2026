"""LHS 475 b (TIC 369327947): trends of the per-sector impact parameter and transit duration, comparison with the JWST
epoch, and the forecast table (README section 6.11).

usage: python scripts/lhs475/trend_lhs475.py [--data data]
input:  <data>/sectors_369327947.csv (per-sector fits per light-curve product), <data>/jwst_369327947.csv (JWST fits)
output: <data>/trend_369327947.csv, <data>/forecast_369327947.csv (formats in the file headers)
Steps:
 1. weighted linear fit of b against year (reference 2022.0) per product, error scaled by sqrt(chi2/dof) when it
    exceeds 1; for SPOC 2-min PDCSAP also with the empirical (injection-calibrated) errors, and leaving out one sector.
 2. weighted means of b in the epoch groups 2019-2021, 2023 and 2025-2026, converted to T14 at the fixed template, and
    that T14 converted back to b at the a/R* of the stellar density (15.02 +- 0.46).
 3. b and T14 of the linear trend at the JWST epoch (2022.67); offsets of the JWST b (a/R* fixed at the template)
    from the trend and from constant-b hypotheses, with an added inter-instrument allowance of 0, 0.005 and 0.010.
 4. forecast: b(t) = 0.780 + 0.0185 (t - 2022.67) with errors 0.005 (anchor) and 0.0045 /yr (slope) combined in
    quadrature, i.e. the uncalibrated PDCSAP trend rounded as recorded on 2026-10-28; T14 and T23 at the template;
    and the same with the slope and anchor of the injection-calibrated fit (full covariance).
Requires numpy only.
"""
import argparse, csv, os
import numpy as np

P, K, AR = 2.029107654, 0.03124, 15.77          # template of sectors_369327947.csv (Ment et al. 2024, arXiv version)
ARS, EARS = 15.02, 0.46                          # a/R* from the stellar density (Gaia DR3 parallax, 2MASS Ks, Mann et al. 2015, 2019)
X0, YJ = 2022.0, 2022.67
PDC = 'SPOC 2-min PDCSAP'
GROUPS = (('2019-2021', 2018, 2022), ('2023', 2022, 2024.5), ('2025-2026', 2024.5, 2027))
YEARS = (2026.9, 2027.0, 2027.25, 2027.5, 2027.75, 2028.0, 2029.0, 2030.0)


def t14(b, k=K, a=AR, inner=False):
    v = ((1 - k) if inner else (1 + k)) ** 2 - b * b
    return P / np.pi * np.arcsin(np.sqrt(max(v, 0.0)) / np.sqrt(a * a - b * b)) * 1440


def t14err(b, s):
    return abs(t14(b + s) - t14(b - s)) / 2


def b_at(T, ar):
    """impact parameter that gives T14 = T (min) at a/R* = ar and k = K (bisection; T14 decreases with b)"""
    lo, hi = 0.0, 1 + K - 1e-9
    for _ in range(100):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if t14(m, a=ar) > T else (lo, m)
    return 0.5 * (lo + hi)


def read(path):
    with open(path) as f:
        return list(csv.DictReader(l for l in f if not l.startswith('#')))


def linfit(x, y, e):
    w = 1 / e ** 2; A = np.vstack([np.ones_like(x), x - X0]).T
    C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ A.T @ (w * y)
    chi2 = float(np.sum(w * (y - A @ p) ** 2)); dof = len(y) - 2; s = max(1.0, chi2 / dof)
    return p, C * s, chi2, dof, np.sqrt(s)


def wmean(y, e):
    w = 1 / e ** 2; return float(np.sum(w * y) / np.sum(w)), float(1 / np.sqrt(np.sum(w)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--data', default='data'); a = ap.parse_args()
    S = read(os.path.join(a.data, 'sectors_369327947.csv')); J = read(os.path.join(a.data, 'jwst_369327947.csv'))
    out = []

    def row(q, var, val, err, unit, note, fmt='%.4f'):
        out.append([q, var, fmt % val, '' if err is None else fmt % err, '' if err in (None, 0) else '%.2f' % (val / err), unit, note])

    variants = list(dict.fromkeys(r['variant'] for r in S))
    for v in variants:
        R = [r for r in S if r['variant'] == v]
        x = np.array([float(r['year']) for r in R]); y = np.array([float(r['b']) for r in R]); e = np.array([float(r['b_err']) for r in R])
        p, C, chi2, dof, sc = linfit(x, y, e)
        row('b slope', v, p[1], np.sqrt(C[1, 1]), '1/yr', f'{len(R)} sectors; chi2 = {chi2:.1f} for {dof} dof; error scaled by {sc:.2f}; b(2022.0) = {p[0]:.3f}')
        if v == PDC:
            ee = np.array([float(r['b_err_empirical']) for r in R])
            pc, Cc, chi2c, dofc, scc = linfit(x, y, ee)
            row('b slope', v + ', empirical errors', pc[1], np.sqrt(Cc[1, 1]), '1/yr', f'{len(R)} sectors; chi2 = {chi2c:.1f} for {dofc} dof; error scaled by {scc:.2f}; b(2022.0) = {pc[0]:.3f}')
            for i, r in enumerate(R):
                m = np.arange(len(R)) != i
                q = linfit(x[m], y[m], e[m]); qc = linfit(x[m], y[m], ee[m])
                row('b slope, leaving out one sector', f'{v}, without sector {r["sector"]}', q[0][1], np.sqrt(q[1][1, 1]), '1/yr',
                    f'empirical errors: {qc[0][1]:+.4f} +- {np.sqrt(qc[1][1, 1]):.4f} ({qc[0][1] / np.sqrt(qc[1][1, 1]):.1f} sigma)')
            G = {}
            for lab, lo, hi in GROUPS:
                m = (x > lo) & (x < hi); mb, sb = wmean(y[m], e[m]); mbc, sbc = wmean(y[m], ee[m]); G[lab] = (mb, sb, sbc)
                secs = ', '.join(R[i]['sector'] for i in np.where(m)[0])
                row('b, epoch group', f'{v}, {lab}', mb, sb, '', f'sectors {secs}; weighted mean; empirical errors: {mbc:.3f} +- {sbc:.3f}')
                row('T14, epoch group', f'{v}, {lab}', t14(mb), t14err(mb, sb), 'min', f'T14 of the group b at the template; with the empirical error +- {t14err(mb, sbc):.2f} min', fmt='%.2f')
            for lab in G:
                T = t14(G[lab][0]); bs = [b_at(T, ar) for ar in (ARS - EARS, ARS, ARS + EARS)]
                row('b at the stellar-density a/R*, epoch group', f'{v}, {lab}', bs[1], None, '',
                    f'b with T14 = {T:.2f} min at a/R* = {ARS} (k fixed); {bs[0]:.3f} at a/R* = {ARS - EARS:.2f}, {bs[2]:.3f} at {ARS + EARS:.2f}; a conversion, not a fit', fmt='%.3f')
            d = G['2025-2026'][0] - G['2019-2021'][0]
            row('b difference, 2025-2026 minus 2019-2021', v, d, np.hypot(G['2025-2026'][1], G['2019-2021'][1]), '',
                f'empirical errors: +- {np.hypot(G["2025-2026"][2], G["2019-2021"][2]):.3f}')
            # trend at the JWST epoch, and the JWST comparison
            g = np.array([1.0, YJ - X0]); bt = float(g @ p); st = float(np.sqrt(g @ C @ g))
            btc = float(g @ pc); stc = float(np.sqrt(g @ Cc @ g))
            row('b of the linear trend at 2022.67', v, bt, st, '', f'empirical errors: {btc:.4f} +- {stc:.4f}')
            row('T14 of the linear trend at 2022.67', v, t14(bt), t14err(bt, st), 'min', '', fmt='%.2f')
            j = [r for r in J if r['config'] == 'aRfix'][0]; bJ = float(j['b']); sJ = 0.5 * (float(j['b_minus']) + float(j['b_plus']))
            for sysv in (0.0, 0.005, 0.010):
                sj = np.hypot(sJ, sysv)
                for lab, (bh, sh) in (('linear trend', (bt, st)), ('constant at 2019-2021', G['2019-2021'][:2]), ('constant at 2025-2026', G['2025-2026'][:2])):
                    row('JWST b minus hypothesis', f'{lab}; allowance {sysv:.3f}', bJ - bh, np.hypot(sh, sj), '',
                        f'JWST b = {bJ:.4f} +- {sJ:.4f} (a/R* = 15.77); hypothesis b = {bh:.4f} +- {sh:.4f}')
            sh = G['2019-2021'][2]
            row('JWST b minus hypothesis', 'constant at 2019-2021, empirical group error; allowance 0.000', bJ - G['2019-2021'][0], np.hypot(sh, sJ), '',
                f'group error {sh:.4f} from the empirical per-sector errors')
            fit_c = (pc, Cc)
    hdr = ('# LHS 475 b = TIC 369327947: trends of the per-sector impact parameter b and duration T14, epoch groups, and comparison with the JWST epoch; '
           'scripts/lhs475/trend_lhs475.py from sectors_369327947.csv and jwst_369327947.csv\n'
           '# slopes: weighted linear fits against year (reference 2022.0), errors scaled by sqrt(chi2/dof) when it exceeds 1; z = value / err; '
           'empirical errors = b_err_empirical of sectors_369327947.csv; T14 at the template k = 0.03124, a/R* = 15.77\n'
           '# all values depend on the fixed a/R* (section 6.11.4); the significances are local (the planet was selected from a screen of 271 planets)\n')
    with open(os.path.join(a.data, 'trend_369327947.csv'), 'w') as f:
        f.write(hdr); w = csv.writer(f, lineterminator='\n'); w.writerow(['quantity', 'variant', 'value', 'err', 'z', 'unit', 'note']); w.writerows(out)

    # forecast
    pc, Cc = fit_c
    rows = []
    for kind, years in (('forecast', YEARS), ('curve', np.round(np.arange(2019.0, 2030.01, 0.25), 2))):
        for y in years:
            b = 0.780 + 0.0185 * (y - YJ); s = float(np.hypot(0.005, 0.0045 * (y - YJ)))
            g = np.array([1.0, y - X0]); bc = float(g @ pc); sc = float(np.sqrt(g @ Cc @ g))
            rows.append([kind, '%.2f' % y, '%.4f' % b, '%.4f' % s, '%.1f' % t14(b), '%.1f' % t14(b + s), '%.1f' % t14(b - s), '%.1f' % t14(b, inner=True),
                         '%.4f' % bc, '%.4f' % sc, '%.1f' % t14(bc), '%.1f' % t14(bc + sc), '%.1f' % t14(bc - sc),
                         0.736, '%.1f' % t14(0.736), 0.79, '%.1f' % t14(0.79)])
    with open(os.path.join(a.data, 'forecast_369327947.csv'), 'w') as f:
        f.write('# LHS 475 b = TIC 369327947: forecast of the impact parameter b and transit duration T14 (first to fourth contact; T23 second to third) at the template k = 0.03124, a/R* = 15.77; '
                'scripts/lhs475/trend_lhs475.py\n'
                '# trend_*: b(t) = 0.780 + 0.0185 (t - 2022.67), error sqrt(0.005^2 + (0.0045 (t - 2022.67))^2): the uncalibrated SPOC 2-min PDCSAP trend as recorded on 2026-10-28, '
                'computed from the TESS fits alone; T14_lo/T14_hi = T14 at b + err / b - err\n'
                '# cal_*: the same from the fit with the empirical (injection-calibrated) errors, full covariance; const_early: b constant at the 2019-2021 value; '
                'const_mean: b constant at about the TESS mean (0.79)\n'
                '# kind = forecast: the epochs of the recorded forecast; kind = curve: grid for the figure. Predictions, not measurements; they scale with the adopted a/R* (section 6.11.4)\n')
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['kind', 'year', 'trend_b', 'trend_b_err', 'trend_T14_min', 'trend_T14_lo_min', 'trend_T14_hi_min', 'trend_T23_min',
                    'cal_b', 'cal_b_err', 'cal_T14_min', 'cal_T14_lo_min', 'cal_T14_hi_min', 'const_early_b', 'const_early_T14_min', 'const_mean_b', 'const_mean_T14_min'])
        w.writerows(rows)
    print('wrote trend_369327947.csv and forecast_369327947.csv')


if __name__ == '__main__':
    main()
