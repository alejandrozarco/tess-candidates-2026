#!/usr/bin/env python
"""TOI-421 CHEOPS impact-parameter fits: build the run list (planet b and c, per visit and joint, PIPE and DRP,
k/LD variants) and execute it with a process pool. Conventions: see cheops_fit.py docstring.
Usage: python scripts/toi421/cheops/run_all.py OUTDIR NPROC [filter-substring] [nsteps] [burn] [variants]
Section 6.12 used OUTDIR = work/toi421/cheops/runs, nsteps 10000, burn 3000 (all 126 runs)."""
import sys, os, json, zlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import cheops_fit as cf

VAR = {  # name: (kfree, ld)
    'kfix_LDatlas': (False, 'atlas'), 'kfix_LDfree': (False, 'free'),
    'kfree_LDatlas': (True, 'atlas'), 'kfree_LDfree': (True, 'free'),
    'kfix_LDphoenix': (False, 'phoenix'), 'kfix_LDldcu': (False, 'ldcu'), 'kfix_LDtess': (False, 'tess'),
    'kfree_LDphoenix': (True, 'phoenix'), 'kfree_LDldcu': (True, 'ldcu'), 'kfree_LDtess': (True, 'tess'),
}
VARS_B = ['kfix_LDatlas', 'kfix_LDfree', 'kfree_LDatlas', 'kfree_LDfree', 'kfix_LDphoenix', 'kfix_LDldcu', 'kfix_LDtess']
VARS_C = ['kfree_LDatlas', 'kfree_LDfree', 'kfix_LDatlas', 'kfix_LDfree', 'kfree_LDphoenix', 'kfree_LDldcu', 'kfree_LDtess']
GROUPS = {'v1': ('b', [1]), 'v2': ('b', [2]), 'v3': ('b', [3]), 'v4': ('b', [4]),
          'jb2021': ('b', [1, 2, 3]), 'jball': ('b', [1, 2, 3, 4]),
          'v5': ('c', [5]), 'v6': ('c', [6]), 'jcall': ('c', [5, 6])}


def specs(nsteps, burn):
    out = []
    for prod in ('PIPE', 'DRP'):
        for g, (pl, vis) in GROUPS.items():
            for vn in (VARS_B if pl == 'b' else VARS_C):
                kf, ld = VAR[vn]
                nm = f"{pl}_{prod}_{g}_{vn}"
                out.append(dict(name=nm, planet=pl, visits=vis, product=prod, kfree=kf, ld=ld,
                                nsteps=nsteps, burn=burn,
                                seed=zlib.crc32(nm.encode())))
    return out


def work(args):
    spec, outdir = args
    if os.path.exists(os.path.join(outdir, f"fit_{spec['name']}.json")):
        return spec['name'] + ' exists'
    try:
        o = cf.run(spec, outdir)
        r = o['results']
        return f"{spec['name']} b={r['b'][0]:.4f} -{r['b'][1]:.4f} +{r['b'][2]:.4f} {o['runtime_s']:.0f}s"
    except Exception as e:
        import traceback
        return f"{spec['name']} FAILED {e!r} {traceback.format_exc()[-800:]}"


if __name__ == '__main__':
    outdir = sys.argv[1]; npr = int(sys.argv[2])
    flt = sys.argv[3] if len(sys.argv) > 3 else ''
    nsteps = int(sys.argv[4]) if len(sys.argv) > 4 else 6000
    burn = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    allow = sys.argv[6].split(',') if len(sys.argv) > 6 else None   # optional variant whitelist
    S = [s for s in specs(nsteps, burn) if all(f in s['name'] for f in flt.split(','))
         and (allow is None or any(s['name'].endswith('_' + a) for a in allow))]
    # heavy (joint) runs first
    S.sort(key=lambda s: -len(s['visits']))
    print(len(S), 'runs', flush=True)
    with Pool(npr) as p:
        for msg in p.imap_unordered(work, [(s, outdir) for s in S]):
            print(msg, flush=True)
