"""Two-star Gaussian PSF fit to median image; decompose difference image into target/neighbour amplitudes.
Usage: p2_psf.py s0008 s0035 ... (reads $TD_WORK/tpf620/<s>_med.npy, _diff.npy, .json written by p1_tpf.py)"""
import sys, os, json, numpy as np
OUT = {}
W = os.path.join(os.environ.get('TD_WORK', os.path.join('work', 'transit_duration')), 'tpf620')
from scipy.optimize import least_squares
for s in sys.argv[1:]:
    med = np.load(os.path.join(W, '%s_med.npy' % s)); dimg = np.load(os.path.join(W, '%s_diff.npy' % s)); J = json.load(open(os.path.join(W, '%s.json' % s)))
    ny, nx = med.shape; yy, xx = np.mgrid[:ny, :nx]
    def g(x0, y0, sx, sy): return np.exp(-0.5 * (((xx - x0) / sx) ** 2 + ((yy - y0) / sy) ** 2))
    bp = J['brightest_px']; tx, ty = J['target_xy']
    # initial: target = brightest px within 1.5 px of header pos, neighbour = brightest px > 2 px away
    tgt = min(bp, key=lambda p: (p[0] - tx) ** 2 + (p[1] - ty) ** 2)
    nb = max([p for p in bp if (p[0] - tgt[0]) ** 2 + (p[1] - tgt[1]) ** 2 > 4], key=lambda p: p[2])
    def model(p):
        a1, x1, y1, a2, x2, y2, sx, sy, c = p
        return a1 * g(x1, y1, sx, sy) + a2 * g(x2, y2, sx, sy) + c
    p0 = [tgt[2], tgt[0], tgt[1], nb[2], nb[0], nb[1], 0.8, 0.8, np.percentile(med, 10)]
    ok = np.isfinite(med)
    r = least_squares(lambda p: (model(p) - med)[ok], p0)
    a1, x1, y1, a2, x2, y2, sx, sy, c = r.x
    G1, G2 = g(x1, y1, sx, sy), g(x2, y2, sx, sy)
    okd = np.isfinite(dimg)
    A = np.vstack([G1[okd], G2[okd], np.ones(okd.sum())]).T
    sol, res, *_ = np.linalg.lstsq(A, dimg[okd], rcond=None)
    resid = dimg[okd] - A @ sol; s2 = np.sum(resid ** 2) / (okd.sum() - 3)
    cov = np.linalg.inv(A.T @ A) * s2
    # also: free single-Gaussian position fit to difference image
    rd = least_squares(lambda p: (p[0] * g(p[1], p[2], sx, sy) + p[3] - dimg)[okd], [dimg[okd].max(), x1, y1, 0])
    sep = np.hypot(x2 - x1, y2 - y1)
    frac_t = sol[0] * G1.sum(); frac_n = sol[1] * G2.sum()
    OUT[s] = dict(target_xy=[float(x1), float(y1)], neighbour_xy=[float(x2), float(y2)], sep_px=float(sep), flux_ratio_n_t=float(a2 / a1),
                  amp_target=float(sol[0]), amp_target_err=float(np.sqrt(cov[0, 0])), amp_neighbour=float(sol[1]), amp_neighbour_err=float(np.sqrt(cov[1, 1])),
                  share_neighbour=float(frac_n / (frac_t + frac_n)), free_xy=[float(rd.x[1]), float(rd.x[2])], free_offset_px=float(np.hypot(rd.x[1] - x1, rd.x[2] - y1)))
    print('S%s target (%.2f,%.2f) F=%.0f  neighbour (%.2f,%.2f) F=%.0f sep %.2f px flux ratio n/t %.2f | diff: amp_t %.1f+-%.1f amp_n %.1f+-%.1f -> share on neighbour %.2f | free diff pos (%.2f,%.2f) offset from target %.2f px'
          % (s, x1, y1, a1 * G1.sum(), x2, y2, a2 * G2.sum(), sep, a2 / a1, sol[0], np.sqrt(cov[0, 0]), sol[1], np.sqrt(cov[1, 1]), frac_n / (frac_t + frac_n), rd.x[1], rd.x[2], np.hypot(rd.x[1] - x1, rd.x[2] - y1)))
json.dump(OUT, open(os.path.join(W, 'psf_decomposition.json'), 'w'), indent=1)
