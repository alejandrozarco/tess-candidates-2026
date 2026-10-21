"""T8: per-transit impact parameter, mid-time and T14 with the transit shape (k, a/R*, LD) fixed at the main posterior
median; quadratic baseline per transit; error from the least-squares covariance scaled by sqrt(chi2_red) if > 1.
Also a stretch-factor width w (template time axis scaled, free depth) as a shape-agnostic measure. Own code.
Usage: t8_pertransit.py <prep.npz> <main_posterior.json> <out.json>"""
import sys, json, numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, __file__.rsplit('/', 1)[0]); from t5lib import *
prep, postf, outp = sys.argv[1:4]
pr = np.load(prep, allow_pickle=True); segs = list(np.load(str(pr['win']), allow_pickle=True))
post = json.load(open(postf)); P = post['post']
k = P['k'][0]; aR = float(aR_from_lnrho(P['lnrho'][0], P_D)); u1, u2 = ld_q2u(P['q1'][0], P['q2'][0])
D = TData(segs); D.build(list(pr['keep'])); D.set_sig(pr['sig'])
S = Spec('lin_b', D=D, ld0=LD0); th = np.array([post['map'][n] for n in S.names]); pm = S.unpack(th)
rows = []
for i in range(D.n):
    Di = TData([segs[i]], sig=[D.sig[i]]); Di.build([D.keep[i]]); Di.set_sig([D.sig[i]])
    tc0 = pm['tc'][i]
    def res_b(x):
        b, dtc = x
        if not (0 <= b < 1 + k - 0.002): return np.full(Di.npts, 1e3)
        m = Di.model(np.array([tc0 + dtc / 1440.]), np.array([b]), k, aR, u1, u2); return Di.chi2(m, ret='resid')
    r = least_squares(res_b, [pm['b'][i], 0.], x_scale=[0.01, 0.5], diff_step=1e-4)
    chi2 = np.sum(r.fun ** 2); dof = Di.npts - 2 - 3; s = max(1., np.sqrt(chi2 / dof))
    cov = np.linalg.inv(r.jac.T @ r.jac) * s * s
    b, dtc = r.x; eb = np.sqrt(cov[0, 0])
    T = t14_from(P_D, aR, b, k) * 1440; dTdb = (t14_from(P_D, aR, b + 1e-5, k) - t14_from(P_D, aR, b - 1e-5, k)) / 2e-5 * 1440
    # stretch width: model evaluated at tc + (t - tc)/w with free depth scale A (k, aR, LD, b fixed at the main fit's b(t))
    bfix = pm['b'][i]
    def res_w(x):
        w, A, dtc2 = x
        if not (0.3 < w < 3 and 0 < A < 3): return np.full(Di.npts, 1e3)
        Di_t = Di.tss.copy()
        tcw = tc0 + dtc2 / 1440.
        Di.tss = tcw + (Di_t - tcw) / w
        m = Di.model(np.array([tcw]), np.array([bfix]), k, aR, u1, u2); Di.tss = Di_t
        return Di.chi2(1 - A * (1 - m), ret='resid')
    rw = least_squares(res_w, [1., 1., dtc], x_scale=[0.02, 0.05, 0.5], diff_step=1e-4)
    covw = np.linalg.inv(rw.jac.T @ rw.jac) * max(1., np.sum(rw.fun ** 2) / (Di.npts - 6))
    sg = segs[i]
    rows.append(dict(i=i, sector=int(sg['sector']), epoch=int(sg['epoch']), year=float(D.yr[i]), camp=int(D.camp[i]), b=float(b), b_err=float(eb),
                     T14=float(T), T14_err=float(abs(dTdb) * eb), dtc_min=float(dtc), dtc_err=float(np.sqrt(cov[1, 1])), chi2=float(chi2), dof=int(dof),
                     w=float(rw.x[0]), w_err=float(np.sqrt(covw[0, 0])), A=float(rw.x[1]), A_err=float(np.sqrt(covw[1, 1])),
                     sig=float(D.sig[i]), cam=sg['cam'], ccd=sg['ccd'], crowd=sg['crowd'], flfr=sg['flfr'], data_rel=sg['data_rel'],
                     procver=sg['procver'], ap_npix=sg['ap_npix'], exptime=sg['exptime'], **sg['tel']))
    print('%2d S%03d %.3f b %.4f+-%.4f T14 %.2f+-%.2f w %.3f+-%.3f A %.2f' % (i, sg['sector'], D.yr[i], b, eb, T, rows[-1]['T14_err'], rw.x[0], rows[-1]['w_err'], rw.x[1]), flush=True)
json.dump(dict(k=k, aR=aR, u1=float(u1), u2=float(u2), rows=rows), open(outp, 'w'), indent=1)
y = np.array([r['year'] for r in rows]); T = np.array([r['T14'] for r in rows]); e = np.array([r['T14_err'] for r in rows])
A = np.vstack([np.ones_like(y), y - TREF]).T; W = 1 / e ** 2
c = np.linalg.solve(A.T @ (A * W[:, None]), A.T @ (W * T)); cv = np.linalg.inv(A.T @ (A * W[:, None])); chi2 = np.sum(W * (T - A @ c) ** 2)
print('per-transit T14 slope %.3f +- %.3f min/yr (chi2 %.1f / %d; scaled err %.3f)' % (c[1], np.sqrt(cv[1, 1]), chi2, len(T) - 2, np.sqrt(cv[1, 1] * max(1, chi2 / (len(T) - 2)))))
