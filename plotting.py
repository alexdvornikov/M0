import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

from utils import *

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
    q = hits['q']

    if plotHits:
        pos3d = hit_to_3d(hits, event)
        ax.scatter(*pos3d, c = q)
    
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
    
    
    
#---------------------------------------------------------------------------------------------------------#
# Plot 10 tracks that are > 30 cm.
# Assuming you downloaded a file from
# https://portal.nersc.gov/project/dune/data/Module0/TPC1+2/dataRuns/tracksData/
# and put it in the same directory.
#---------------------------------------------------------------------------------------------------------#
f = h5py.File(file_name, 'r') # 

fig = plt.figure() 
ax = fig.add_subplot(111, projection = '3d') 
i = 0
for thisTrack in f['tracks']:
        if thisTrack['length'] > 300 and i < 10: 
                i += 1
                plot_track(ax, thisTrack, f, plotHits = True) 

plt.show() 
#---------------------------------------------------------------------------------------------------------#
