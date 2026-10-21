"""Keep only the part of a prepared JWST white-light curve after a given time (used to give the GO 4098 visit the same
pre-ingress coverage as the GTO 1224 visit). Usage: j3_cut.py <in.npz> <out.npz> <t_min_btjd>"""
import sys, numpy as np
z = np.load(sys.argv[1]); tmin = float(sys.argv[3]); m = z['t'] >= tmin
np.savez(sys.argv[2], t=z['t'][m], f=z['f'][m], e=z['e'][m])
print('kept', m.sum(), 'of', len(m), 'points from BTJD %.5f' % tmin)
