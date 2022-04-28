import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

def draw_boundaries(ax):
    bounds = [[-300, 300], [-600, 600], [-310, 310]]
    cathodeZ = 0

    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][0]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][0]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][1], bounds[2][1]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][1], bounds[2][1]],
            color = 'black', ls = '--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]],
            color = 'black', ls = '--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]],
            color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]],
            color = 'black', ls = '--')

    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [cathodeZ, cathodeZ],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [cathodeZ, cathodeZ],
            color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [cathodeZ, cathodeZ],
            color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [cathodeZ, cathodeZ],
            color = 'black', ls = '--')


def plot_track(ax, track, f, plotHits = False):

    draw_boundaries(ax)

    event = f['events'][track['event_ref']]
    t0 = event['ts_start']

    hits = f['hits'][track['hit_ref']]

    x = hits['px']
    y = hits['py']
    t = 0.1*(hits['ts'] - t0)
    v = 1.62

    q = hits['q']
    grp = hits['iogroup']
    parity = np.power(-1, grp)
    z = parity*(310. - t*v) 

    if plotHits:
        ax.scatter(x, y, z, c = q)
    
    start = np.array([track['start'][0],
                      track['start'][1],
                      track['start'][2]])
    end = np.array([track['end'][0],
                    track['end'][1],
                    track['end'][2]])
    ax.scatter(*start,
               c = 'r')
    ax.scatter(*end,
               c = 'b')
    ax.plot(*zip(start, end),
            c = 'g')
