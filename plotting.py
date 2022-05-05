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

def main(args):
    f = h5py.File(args.infile, 'r')

    rawTracks = np.array(f['tracks'])

    trackMask = np.logical_and(rawTracks['length'] > 300,
                               rawTracks['nhit'] > 0) # more variables to cut on here

    tracks = rawTracks[trackMask]
    
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection = '3d') 

    for thisTrack in tracks[:args.n]:
        plot_track(ax, thisTrack, f, plotHits = True) 

    if args.o:
        plt.savefig(args.o, dpi = 300)
    else:
        plt.show() 

if __name__ == '__main__': 
    # Plot 10 tracks that are > 30 cm.
    # Assuming you downloaded a file from
    # https://portal.nersc.gov/project/dune/data/Module0/TPC1+2/dataRuns/tracksData/
    import argparse

    parser = argparse.ArgumentParser(description='Plot the first N tracks from a given file')
    parser.add_argument('infile',
                        help = 'intput larpix data with track reconstruction')
    parser.add_argument('-n',
                        default = 10,
                        type = int,
                        help = 'plot the first n tracks (default 10)')
    parser.add_argument('-o',
                        default = '',
                        type = str,
                        help = 'save the output plot to a file')

    args = parser.parse_args()

    main(args)
