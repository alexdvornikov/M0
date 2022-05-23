import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

from utils import *
from plotting import *

def track_selection(track, f):

    #TPC dims
    for ix in range(detector.TPC_BORDERS.shape[0]):
        bounds = detector.TPC_BORDERS[ix]*10 #in mm


    start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                      track['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                      track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                    track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                    track['end'][2] + my_geometry.tpc_offsets[0][2]*10])

    p1_x,p1_y,p1_z = start[0],start[1],start[2]
    p2_x,p2_y,p2_z = end[0],end[1],end[2]
    
    # Distances in mm
    epsilon = 10
    z_bound = bounds[0][1]
    y_top, y_bottom = bounds[1][1], bounds[1][0]
    
    if (
        ( (abs(abs(p1_z) - z_bound) < epsilon) ) and # Near either anode (in z)
        ( (abs(p2_y - y_top) < epsilon) ) and # Near top face (in y)
        ( abs(p1_z) - abs(p2_z) > 2*epsilon ) and # Give it some angle (avoid endpoints with nearby zs, parlle with the anode)
        ( abs(p1_z) < (z_bound + 2*epsilon) and abs(p2_z) < (z_bound +2*epsilon) ) # Make sure thhat both z endpoints are within the TPC. 
    ):
        return True 

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

    rawTracks = np.array(f['tracks'])
    trackMask = np.logical_and(rawTracks['length'] > 300,
                               rawTracks['nhit'] > 0) # more variables to cut on here
    tracks = rawTracks[trackMask]
    
    for thisTrack in tracks[:args.n]:
        if track_selection(thisTrack, f):
        #     plot_track(ax, thisTrack, f)
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

