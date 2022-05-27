# python3 selection.py ~/Desktop/Desktop/datalog_2021_04_03_20_01_46_CESTevd.h5 -n 1000

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d

import h5py

from utils import *
from plotting import *


#------------------------------------------------------------------------------#
# Helper functions
#------------------------------------------------------------------------------#

def get_TPC_bounds():
    bounds = []
    for ix in range(detector.TPC_BORDERS.shape[0]):
        bounds.append(detector.TPC_BORDERS[ix]*10)
    return np.array( bounds )


def get_track_ends(track, my_geometry):
    track_start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                            track['start'][1] +my_geometry.tpc_offsets[0][1]*10,
                            track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    track_end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                          track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                          track['end'][2] + my_geometry.tpc_offsets[0][2]*10])
    return track_start, track_end


def shift_track_to_anode(track_start,track_end):

    p1_z = track_start[2]
    p2_z = track_end[2]

    if p1_z > 0: #If first TPC (anode at positive z)
        shift = anode_z - p1_z
        p1_z_shifted = anode_z
        p2_z_shifted = p2_z + shift

    if p1_z < 0: #If second TPC (anode at negative z)
        shift = anode_z - abs(p1_z) 
        p1_z_shifted = -anode_z
        p2_z_shifted = p2_z - shift

    return p1_z_shifted, p2_z_shifted

#------------------------------------------------------------------------------#
# Selection functions
#------------------------------------------------------------------------------#

def is_good_track(track_start,track_end):
    # Drift direction containment.
    # Parallel avoidance (ensure some separation in drift direction z).
    p1_z,p2_z = track_start[2], track_end[2]
    if ( ( abs(p1_z) < (anode_z + epsilon) and abs(p2_z) < (anode_z + epsilon) ) and
        ( (abs(p1_z) - abs(p2_z) > epsilon) ) ):
        return True

def is_cathode_piercer(track):
    # TODO: compare start and end points to cathode z values
    # track_start,track_end = get_track_ends(track,f)
    return 

def is_anode_piercer(p1_z,p2_z):
    # Using two abs() successively checks both anodes. 
    if ( abs(abs(p1_z) - anode_z) < epsilon ):
        return True

def is_top_piercer(p1_y,p2_y):
    if ( abs(p2_y - top) < epsilon ):
        return True

def is_bottom_piercer(p1_y,p2_y):
    if ( abs(p2_y - bottom) < epsilon ):
        return True

def is_upstream_piercer(p1_x,p2_x):
    if ( abs(p2_x - upstream) < epsilon ):
        return True

def is_downstream_piercer(p1_x,p2_x):
    if ( abs(p2_x - downstream) < epsilon ):
        return True

def is_side_piercer(track_start, track_end):
    p1_x, p1_y, p1_z = track_start
    p2_x, p2_y, p2_z = track_end
    if ( is_downstream_piercer(p1_x,p2_x) or
        is_upstream_piercer(p1_x,p2_x) or
        is_bottom_piercer(p1_y,p2_y) or
        is_top_piercer(p1_y,p2_y) ):
        return True

# def is_wall_piercer(track_start,track_end):
#     p1_x, p1_y, p1_z = track_start
#     p2_x, p2_y, p2_z = track_end
#     return bool(np.sum([is_side_piercer(track_start,track_end),
#                         is_anode_piercer(p1_z,p2_z),
#                         is_cathode_piercer(track)], axis = -1))

def track_selection(track_start,track_end):
    # This function is the ultimate filter.
    # Toggle desirable/undesirable conditions below. 
    p1_x, p1_y, p1_z = track_start
    p2_x, p2_y, p2_z = track_end
    if ( is_good_track(track_start,track_end) and 
        is_anode_piercer(p1_z,p2_z) and 
        # is_top_piercer(p1_y,p2_y) and 
        # is_bottom_piercer(p1_y,p2_y) and 
        # is_upstream_piercer(p1_x,p2_x) and 
        is_downstream_piercer(p1_x,p2_x)):
        return True

#------------------------------------------------------------------------------#
# Plotting
#------------------------------------------------------------------------------#
def plot_selected_track(ax, track_start,track_end):
    ax.plot(*zip(track_start, track_end),c='k',alpha = 0.5)




#------------------------------------------------------------------------------#
# Where the magic happens.
#------------------------------------------------------------------------------#
def main(args):
    global my_geometry
    global TPC_bounds, anode_z, top, bottom, upstream, downstream
    global passing_counter
    global length_cut
    global epsilon

    my_geometry = DetectorGeometry(args.d, args.g)
    TPC_bounds = get_TPC_bounds()
    anode_z = TPC_bounds[0][0][1]
    top = TPC_bounds[0][1][1]
    bottom = TPC_bounds[0][1][0]
    upstream = TPC_bounds[0][0][1]
    downstream = TPC_bounds[0][0][0]
    length_cut = 100 #mm
    epsilon = 50 #mm

    global f
    f = h5py.File(args.infile, 'r')
    
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection = '3d') 
    ax.set_xlabel(r'x (horizontal) [mm]')
    ax.set_ylabel(r'y (vertical) [mm]')
    ax.set_zlabel(r'z (drift) [mm]')
    # Set background to white
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # draw_boundaries(ax,TPC_bounds)
    draw_boundaries(ax)

    rawTracks = np.array(f['tracks'])
    # print(rawTracks.size)
    trackMask = np.logical_and(rawTracks['length'] > length_cut,
                               rawTracks['nhit'] > 0) # more variables to cut on here
    tracks = rawTracks[trackMask]
    # print(tracks.size)
    # print(100*(tracks.size/rawTracks.size))
    
    passing_counter = 0
    for thisTrack in tracks[:args.n]:
        track_start,track_end = get_track_ends(thisTrack,my_geometry)
        p1_z_shifted, p2_z_shifted = shift_track_to_anode(track_start,track_end)
        shifted_start = [ track_start[0], track_start[1], p1_z_shifted ]
        shifted_end = [ track_end[0], track_end[1], p2_z_shifted ]

        if track_selection(shifted_start,shifted_end):
            passing_counter += 1
            plot_selected_track(ax, shifted_start,shifted_end)
        # hits = f['hits'][thisTrack['hit_ref']]
        # plot_hits(ax, hits, my_geometry, thisTrack)
        # print ("this track passes the cuts: " + str(track_selection(thisTrack, f)))
        
    # print('Percentage of tracks over ' + str(length_cut) + '[mm] passing the cuts')
    # print(100*(passing_counter/tracks.size))

    passing_fraction = round( 100*(passing_counter/tracks.size), 2 )

    # Some helpful stats to print on the figure. 
    text = 'Number of tracks over ' + str(length_cut/10) + ' [cm] : ' + str(tracks.size)
    text2 = 'Percent of tracks over ' + str(length_cut/10) + ' [cm] passing the cuts: ' + str(passing_fraction) +'%'
    text3 = 'Anode + Downstream Crossers' #Change this to whichever selection is made
    ax.text2D(0.05, 0.95, text, transform=ax.transAxes)
    ax.text2D(0.05, 0.90, text2, transform=ax.transAxes)
    ax.text2D(0.05, 0.0, text3, transform=ax.transAxes)


    if args.o:
        plt.savefig(args.o, dpi = 300)
    else:
        plt.show() 



#------------------------------------------------------------------------------#
# Argument parsing
#------------------------------------------------------------------------------#
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