"""Filter windows+timing to drop transits with other-planet points in window. Usage: v7_filter.py win.npy tim.npz outwin.npy outtim.npz"""
import sys, numpy as np
w = np.load(sys.argv[1], allow_pickle=True); T = dict(np.load(sys.argv[2], allow_pickle=True))
keep = np.array([not s['other_in_window'] for s in w])
np.save(sys.argv[3], w[keep], allow_pickle=True)
T2 = dict(rows=T['rows'][keep], pars=T['pars'], tc=T['tc'][keep], sig=T['sig'][keep], masks=np.array(list(T['masks'][keep]) , dtype=object))
np.savez(sys.argv[4], **T2, allow_pickle=True)
print('kept', keep.sum(), 'of', len(w))
