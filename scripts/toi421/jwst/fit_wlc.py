"""TOI-421 b (TIC 94986319) JWST white-light transit fit with fixed TESS-comparable conventions.

Conventions (fixed; README section 6.12):
  P = 5.19759807 d (fixed), reference T0 = BTJD 1441.283952 (BJD_TDB = BTJD + 2457000); tc fitted free.
  Circular orbit. batman quadratic LD, model integrated over the integration time (supersample 5).
  Free: tc, k=Rp/R*, b, systematics = quadratic in time + exponential ramp A*exp(-(t-t_start)/tau),
  jitter term added in quadrature to the per-point errors.
  Variants:
    a  : a/R* fixed 13.76; u1,u2 free via Kipping q1,q2 ~ U(0,1)
    b  : a/R* fixed 13.76; u1,u2 fixed to ExoTiC-LD 3.2.0 MPS1 values (ld_theory.json)
    bs : as b but ExoTiC-LD Stagger-grid values (sensitivity)
    c  : a/R* free U(3,40), no density prior; q1,q2 free
    bl : as b but systematics linear in time only, no ramp (systematics sensitivity)
Priors: b ~ U(0, 1+k); k ~ U(0.005, 0.1); tc ~ U(pred-0.04, pred+0.04) d.
Usage: python -I fit_wlc.py <lc.npz> <variant> <ld_theory.json> <outdir> [nsteps] [trim_min]
Input npz: t (BMJD_TDB), f, e (normalised), exptime (d), meta, label, ld_band.
"""
import sys, json, os
import numpy as np
import batman, emcee
from scipy.optimize import minimize

P = 5.19759807
T0_BMJD = 1441.283952 + 2457000.0 - 2400000.5
AR_FIX = 13.76

lcfile, variant, ldfile, outdir = sys.argv[1:5]
nsteps = int(sys.argv[5]) if len(sys.argv) > 5 else 12000
trim_min = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
d = np.load(lcfile, allow_pickle=False)
t, f, e = d["t"].astype(float), d["f"].astype(float), d["e"].astype(float)
exptime = float(d["exptime"])
meta, label, ld_band = str(d["meta"]), str(d["label"]), str(d["ld_band"])
if trim_min > 0:
    k_ = t >= t[0] + trim_min / 1440.0
    t, f, e = t[k_], f[k_], e[k_]
ldth = json.load(open(ldfile))["bands"][ld_band]
epoch = np.round((np.mean(t) - T0_BMJD) / P)
tc_pred = T0_BMJD + epoch * P
assert t[0] < tc_pred < t[-1], "predicted mid-time outside visit"
tmid = 0.5 * (t[0] + t[-1])
x = t - tmid
tstart = t[0]

free_ld = variant in ("a", "c")
free_ar = variant == "c"
ramp = variant != "bl"
quad = variant != "bl"
if variant == "b" or variant == "bl":
    u_fix = ldth["mps1"]
elif variant == "bs":
    u_fix = ldth["stagger"]
else:
    u_fix = None

names = ["tc", "k", "b"]
if free_ld: names += ["q1", "q2"]
if free_ar: names += ["aR"]
names += ["c0", "c1"]
if quad: names += ["c2"]
if ramp: names += ["A", "lntau"]
names += ["lnjit"]
ix = {n: i for i, n in enumerate(names)}

pm = batman.TransitParams()
pm.t0 = tc_pred; pm.per = P; pm.rp = 0.028; pm.a = AR_FIX; pm.inc = 86.0
pm.ecc = 0.0; pm.w = 90.0; pm.limb_dark = "quadratic"; pm.u = [0.1, 0.1]
bm = batman.TransitModel(pm, t, supersample_factor=5, exp_time=exptime)


def unpack(p):
    k = p[ix["k"]]; b = p[ix["b"]]
    if free_ld:
        q1, q2 = p[ix["q1"]], p[ix["q2"]]
        u = [2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)]
    else:
        u = list(u_fix)
    aR = p[ix["aR"]] if free_ar else AR_FIX
    return k, b, u, aR


def model(p):
    k, b, u, aR = unpack(p)
    pm.t0 = p[ix["tc"]]; pm.rp = k; pm.a = aR
    pm.inc = np.degrees(np.arccos(b / aR)); pm.u = u
    lc = bm.light_curve(pm)
    s = p[ix["c0"]] + p[ix["c1"]] * x
    if quad: s = s + p[ix["c2"]] * x ** 2
    if ramp: s = s + p[ix["A"]] * np.exp(-(t - tstart) / np.exp(p[ix["lntau"]]))
    return lc * s


def lnprior(p):
    k = p[ix["k"]]; b = p[ix["b"]]
    if not (0.005 < k < 0.1 and 0 <= b < 1 + k): return -np.inf
    if abs(p[ix["tc"]] - tc_pred) > 0.04: return -np.inf
    if free_ld and not (0 < p[ix["q1"]] < 1 and 0 < p[ix["q2"]] < 1): return -np.inf
    if free_ar and not (3 < p[ix["aR"]] < 40): return -np.inf
    if ramp and not (-0.005 < p[ix["A"]] < 0.005 and np.log(0.5 / 1440) < p[ix["lntau"]] < np.log(0.08)):
        return -np.inf
    if not (-15 < p[ix["lnjit"]] < -4): return -np.inf
    return 0.0


def lnprob(p):
    lp = lnprior(p)
    if not np.isfinite(lp): return -np.inf
    m = model(p)
    s2 = e ** 2 + np.exp(2 * p[ix["lnjit"]])
    return lp - 0.5 * np.sum((f - m) ** 2 / s2 + np.log(s2))


p0 = np.zeros(len(names))
p0[ix["tc"]] = tc_pred; p0[ix["k"]] = 0.029; p0[ix["b"]] = 0.95
if free_ld: p0[ix["q1"]] = 0.1; p0[ix["q2"]] = 0.3
if free_ar: p0[ix["aR"]] = AR_FIX
p0[ix["c0"]] = 1.0
if ramp: p0[ix["A"]] = 1e-4; p0[ix["lntau"]] = np.log(15 / 1440)
p0[ix["lnjit"]] = np.log(1e-5)
best = None
for b0 in (0.90, 0.95, 0.97):
    for tc0 in (tc_pred - 0.005, tc_pred, tc_pred + 0.005):
        q = p0.copy(); q[ix["b"]] = b0; q[ix["tc"]] = tc0
        r = minimize(lambda z: -lnprob(z), q, method="Nelder-Mead",
                     options={"maxiter": 40000, "maxfev": 40000, "xatol": 1e-9, "fatol": 1e-4})
        if best is None or r.fun < best.fun: best = r
pbest = best.x
nw = max(48, 4 * len(names))
rng = np.random.default_rng(42)
scale = np.full(len(names), 1e-5)
scale[ix["b"]] = 3e-3; scale[ix["k"]] = 3e-4; scale[ix["tc"]] = 2e-4
if free_ld: scale[ix["q1"]] = 0.02; scale[ix["q2"]] = 0.02
if free_ar: scale[ix["aR"]] = 0.1
if ramp: scale[ix["lntau"]] = 0.1
scale[ix["lnjit"]] = 0.1
start = []
while len(start) < nw:
    q = pbest + scale * rng.standard_normal(len(names))
    if np.isfinite(lnprob(q)): start.append(q)
sampler = emcee.EnsembleSampler(nw, len(names), lnprob)
sampler.run_mcmc(np.array(start), nsteps, progress=False)
burn = nsteps // 3
ch = sampler.get_chain(discard=burn, flat=True, thin=5)
try:
    tau = sampler.get_autocorr_time(discard=burn, quiet=True)
except Exception:
    tau = np.full(len(names), np.nan)
lnp = sampler.get_log_prob(discard=burn, flat=True, thin=5)
pmap = ch[np.argmax(lnp)]

# derived
k_s, b_s = ch[:, ix["k"]], ch[:, ix["b"]]
aR_s = ch[:, ix["aR"]] if free_ar else np.full(len(ch), AR_FIX)
sini = np.sqrt(1 - (b_s / aR_s) ** 2)
a14 = np.clip(np.sqrt(np.clip((1 + k_s) ** 2 - b_s ** 2, 0, None)) / (aR_s * sini), 0, 1)
a23 = np.clip(np.sqrt(np.clip((1 - k_s) ** 2 - b_s ** 2, 0, None)) / (aR_s * sini), 0, 1)
T14 = P / np.pi * np.arcsin(a14) * 24
T23 = P / np.pi * np.arcsin(a23) * 24
inc = np.degrees(np.arccos(b_s / aR_s))
if free_ld:
    q1, q2 = ch[:, ix["q1"]], ch[:, ix["q2"]]
    u1s, u2s = 2 * np.sqrt(q1) * q2, np.sqrt(q1) * (1 - 2 * q2)
else:
    u1s, u2s = np.full(len(ch), u_fix[0]), np.full(len(ch), u_fix[1])


def st(a):
    q = np.percentile(a, [16, 50, 84])
    return {"median": float(q[1]), "minus": float(q[1] - q[0]), "plus": float(q[2] - q[1])}


res = {"planet": "TOI-421 b (TIC 94986319)", "label": label, "data_meta": meta, "variant": variant,
       "conventions": ("P=5.19759807 d fixed; ref T0=BTJD 1441.283952; circular; batman quadratic LD, "
                       "supersample 5 over exptime; " +
                       {"a": "a/R*=13.76 fixed, q1,q2 free U(0,1)",
                        "b": "a/R*=13.76 fixed, u1,u2 fixed ExoTiC-LD 3.2.0 MPS1",
                        "bs": "a/R*=13.76 fixed, u1,u2 fixed ExoTiC-LD 3.2.0 Stagger",
                        "c": "a/R* free U(3,40), q1,q2 free U(0,1)",
                        "bl": "a/R*=13.76 fixed, u1,u2 fixed ExoTiC-LD MPS1, linear baseline, no ramp"}[variant]
                       + ("; systematics quadratic+exp ramp" if ramp else "; systematics linear") +
                       "; jitter in quadrature; emcee"),
       "trim_first_min": trim_min, "n_points": int(len(t)), "exptime_s": exptime * 86400,
       "epoch": int(epoch), "tc_pred_BMJD_TDB": tc_pred, "u_fixed": u_fix,
       "nwalkers": nw, "nsteps": nsteps, "burn": burn, "autocorr": dict(zip(names, map(float, tau))),
       "params": {n: st(ch[:, i]) for i, n in enumerate(names)},
       "derived": {"T14_h": st(T14), "T23_h": st(T23), "inc_deg": st(inc), "aR": st(aR_s),
                   "u1": st(u1s), "u2": st(u2s), "depth_k2_ppm": st(k_s ** 2 * 1e6),
                   "frac_grazing_T23_zero": float(np.mean(T23 == 0)),
                   "tc_minus_pred_min": st((ch[:, ix["tc"]] - tc_pred) * 1440)},
       "corr_b_k": float(np.corrcoef(b_s, k_s)[0, 1]),
       "corr_b_aR": float(np.corrcoef(b_s, aR_s)[0, 1]) if free_ar else None,
       "map": dict(zip(names, map(float, pmap)))}
# residual diagnostics at MAP
m = model(pmap)
r = f - m
sig = np.sqrt(e ** 2 + np.exp(2 * pmap[ix["lnjit"]]))
res["diag"] = {"rms_ppm": float(np.std(r) * 1e6), "median_err_ppm": float(np.median(e) * 1e6),
               "jitter_ppm": float(np.exp(pmap[ix["lnjit"]]) * 1e6),
               "redchi2_raw_err": float(np.sum((r / e) ** 2) / (len(t) - len(names))),
               "redchi2_with_jit": float(np.sum((r / sig) ** 2) / (len(t) - len(names)))}
# time-averaging beta and in-transit binned residual scan (spot-crossing check)
betas = {}
for nb in (10, 30, 60, 120):
    nbin = len(r) // nb
    rb = r[:nbin * nb].reshape(nbin, nb).mean(1)
    betas[str(nb)] = float(np.std(rb) / (np.std(r) / np.sqrt(nb) * np.sqrt(nbin / (nbin - 1))))
res["diag"]["beta_time_avg"] = betas
pm.t0 = pmap[ix["tc"]]
kk, bb, uu, aa = unpack(pmap)
half = res["derived"]["T14_h"]["median"] / 48.0
intr = np.abs(t - pmap[ix["tc"]]) < half
binw = 5.0 / 1440
edges = np.arange(t[0], t[-1] + binw, binw)
idx = np.digitize(t, edges)
rows = []
for j in np.unique(idx):
    s = idx == j
    if s.sum() < 5: continue
    rows.append([float(t[s].mean()), float(r[s].mean() * 1e6), float(np.std(r) / np.sqrt(s.sum()) * 1e6),
                 bool(intr[s].mean() > 0.5)])
rows = np.array(rows, dtype=object)
zin = np.array([row[1] / row[2] for row in rows if row[3]], float)
zout = np.array([row[1] / row[2] for row in rows if not row[3]], float)
res["diag"]["inTransit_5min_bins"] = {"n": int(len(zin)), "max_abs_z": float(np.max(np.abs(zin))) if len(zin) else None,
                                     "chi2_per_bin": float(np.mean(zin ** 2)) if len(zin) else None,
                                     "max_pos_resid_ppm": float(max(row[1] for row in rows if row[3])) if len(zin) else None}
res["diag"]["outOfTransit_5min_bins"] = {"n": int(len(zout)), "max_abs_z": float(np.max(np.abs(zout))),
                                        "chi2_per_bin": float(np.mean(zout ** 2))}
os.makedirs(outdir, exist_ok=True)
tag = f"{label}_var{variant}" + (f"_trim{trim_min:g}" if trim_min > 0 else "")
json.dump(res, open(os.path.join(outdir, tag + ".json"), "w"), indent=1)
np.savez(os.path.join(outdir, tag + "_resid.npz"), t=t, f=f, model=m, resid=r, meta=meta + " | fit " + res["conventions"])
np.savez(os.path.join(outdir, tag + "_chain_b_k.npz"), b=b_s, k=k_s, aR=aR_s,
         meta=res["planet"] + " | " + meta + " | variant " + variant + ": " + res["conventions"])
print(tag, "b=", res["params"]["b"], "k=", res["params"]["k"], "corr", res["corr_b_k"], "diag", res["diag"]["redchi2_raw_err"], betas)
