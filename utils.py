import numpy as np
import h5py

vd = 1.62 # mm/us
drift_dist = 310. # mm
clock_interval = 0.1 # us

def pairWiseDist(posPairArray):
    return np.sqrt(np.sum(np.power(posPairArray[:,1,:] - posPairArray[:,0,:], 2), axis = -1))

def hit_to_3d(hits, event):
    t0 = event['ts_start']

    x = hits['px']
    y = hits['py']
    t = clock_interval*(hits['ts'] - t0)
    grp = hits['iogroup']
    #parity = np.power(-1, grp)
    parity = np.power(-1, grp + 1) 
    z = parity*(drift_dist - t*vd)

    pos3d = np.array([x, y, z])
    return pos3d
