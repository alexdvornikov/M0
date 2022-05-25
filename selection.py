import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

from utils import *
from plotting import *

# Distances in mm
epsilon = 10
    
def is_good_track(track, f):
    # TODO: implement global track cuts (remove noise tracks)
    return True

def is_anode_piercer(track, f):
    # TODO: compare start and end points to anode z values
    start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                      track['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                      track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                    track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                    track['end'][2] + my_geometry.tpc_offsets[0][2]*10])
    return True

def is_cathode_piercer(track, f):
    # TODO: compare start and end points to cathode z values
    start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                      track['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                      track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                    track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                    track['end'][2] + my_geometry.tpc_offsets[0][2]*10])
    return True

def is_side_piercer(track, f):
    start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                      track['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                      track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                    track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                    track['end'][2] + my_geometry.tpc_offsets[0][2]*10])

    print (detector.TPC_BORDERS*10)
    
    start_near_xLow_wall = abs(start[0] - bounds[0][0]) < epsilon
    start_near_xHigh_wall = abs(start[0] - bounds[0][1]) < epsilon

    start_near_yLow_wall = abs(start[1] - bounds[1][0]) < epsilon
    start_near_yHigh_wall = abs(start[1] - bounds[1][1]) < epsilon

    start_near_wall = (start_near_xLow_wall or
                       start_near_xHigh_wall or
                       start_near_yLow_wall or
                       start_near_yHigh_wall)
    
    end_near_xLow_wall = abs(end[0] - bounds[0][0]) < epsilon
    end_near_xHigh_wall = abs(end[0] - bounds[0][1]) < epsilon

    end_near_yLow_wall = abs(end[1] - bounds[1][0]) < epsilon
    end_near_yHigh_wall = abs(end[1] - bounds[1][1]) < epsilon

    end_near_wall = (end_near_xLow_wall or
                     end_near_xHigh_wall or
                     end_near_yLow_wall or
                     end_near_yHigh_wall)

    return 

def is_wall_piercer(track, f):
    return bool(np.sum([is_side_piercer(track, f),
                        is_anode_piercer(track, f),
                        is_cathode_piercer(track, f)], axis = -1))

def track_selection(track, f):
    # TODO: this function is the ultimate filter
    return is_wall_piercer(track, f)

def main(args):
    global my_geometry
    my_geometry = DetectorGeometry(args.d, args.g)

    global f
    f = h5py.File(args.infile, 'r')
    
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection = '3d') 

    draw_boundaries(ax)

    rawTracks = np.array(f['tracks'])
    trackMask = np.logical_and(rawTracks['length'] > 300,
                               rawTracks['nhit'] > 0) # more variables to cut on here
    tracks = rawTracks[trackMask]
    
    for thisTrack in tracks[:args.n]:
        # if track_selection(thisTrack, f):
        #     plot_track(ax, thisTrack, f)
        hits = f['hits'][thisTrack['hit_ref']]
        plot_hits(ax, hits, my_geometry, thisTrack)
        print ("this track passes the cuts: " + str(track_selection(thisTrack, f)))
        
    if args.o:
        plt.savefig(args.o, dpi = 300)
    else:
        plt.show() 

if __name__ == '__main__': 
    import argparse

    parser = argparse.ArgumentParser(description='Plot the first N tracks from a given file')
    parser.add_argument('infile',
                        help = 'intput larpix data with track reconstruction')
    parser.add_argument('-n',
                        default = 10,
                        type = int,
                        help = 'evaluate the selection criteria over the first N tracks')
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
                        help = 'save the data which passes the selection to a file')

    args = parser.parse_args()

    main(args)

