import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

from utils import *

debug = False

def draw_boundaries(ax):
    for ix in range(detector.TPC_BORDERS.shape[0]):
        bounds = detector.TPC_BORDERS[ix]*10
                 
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

        
def plot_hits(ax, hits, track = None):
    if track:
        t0 = track['t0']
    else:
        event = f['events'][hits[0]['event_ref']]
        t0 = event['ts_start']
    
    pos3d = hit_to_3d(my_geometry, hits, t0)

    q = hits['q']
    if debug:
        color = ['blue' if hit['iogroup'] == 1 else 'red' if hit['iogroup'] == 2 else 'yellow'
                 for hit in hits]
    else:
        color = q
        
    ax.scatter(*pos3d, c = color)
    
def plot_track(ax, track, f):
    start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                      track['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                      track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                    track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                    track['end'][2] + my_geometry.tpc_offsets[0][2]*10])

    ax.scatter(*start,
               c = 'r')
    ax.scatter(*end,
               c = 'b')
    ax.plot(*zip(start, end),
            c = 'g')

def main(args):
    global my_geometry
    my_geometry = DetectorGeometry(args.d, args.g)

    global f
    f = h5py.File(args.infile, 'r')
    
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection = '3d') 

    ax.set_xlabel(r'x (horizontal) [mm]')
    ax.set_ylabel(r'y (vertical) [mm]')
    ax.set_zlabel(r'z (drift) [mm]')
    
    draw_boundaries(ax)

    if args.e > 0:
        # events = np.array(f['events'])
        # eventMask = (events['evid'] == args.e)
        # print (args.e)
        # thisEvent = events[eventMask]
        thisEvent = np.array(f['events'])[args.e]

        print ("this event has " + str(thisEvent['n_ext_trigs']) + " external triggers") 

        hits = f['hits'][thisEvent['hit_ref']]

        plot_hits(ax, hits)

        for ti in range(thisEvent['ntracks']):
            thisTrack = f['tracks'][thisEvent['track_ref']][ti]
            plot_track(ax, thisTrack, f)
    
    else:
        rawTracks = np.array(f['tracks'])
        trackMask = np.logical_and(rawTracks['length'] > 300,
                                   rawTracks['nhit'] > 0) # more variables to cut on here
        tracks = rawTracks[trackMask]


        for thisTrack in tracks[:args.n]:
            plot_track(ax, thisTrack, f)

            hits = f['hits'][thisTrack['hit_ref']]
            plot_hits(ax, hits, thisTrack)

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
    parser.add_argument('-e',
                        default = -1,
                        type = int,
                        help = 'plot the event with this evid (default, plot tracks)')
    parser.add_argument('-n',
                        default = 10,
                        type = int,
                        help = 'plot the first n tracks (default 10)')
    parser.add_argument('-g',
                        default = './pixel_layouts/multi_tile_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-o',
                        default = '',
                        type = str,
                        help = 'save the output plot to a file')

    args = parser.parse_args()

    main(args)
