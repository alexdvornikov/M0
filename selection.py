# python3 selection.py ~/Desktop/Desktop/datalog_2021_04_03_20_01_46_CESTevd.h5 -n 1000

import numpy as np

import h5py

from utils import *

def is_good_track(track_start, track_end):
    # Drift direction containment.
    # Parallel avoidance (ensure some separation in drift direction z).
    p1_z,p2_z = track_start[2], track_end[2]
    return ( ( abs(p1_z) < (anode_z + epsilon) and abs(p2_z) < (anode_z + epsilon) ) and
             ( (abs(p1_z) - abs(p2_z) > epsilon) ) )

def is_cathode_piercer(track_start, track_end):
    p1_z = track_start[2]
    p2_z = track_end[2]

    p1Cross = approx_equals(abs(p1_z), cathode_z, epsilon)
    p2Cross = approx_equals(abs(p2_z), cathode_z, epsilon)
    return p1Cross or p2Cross
    
def is_anode_piercer(track_start, track_end, anode_z, epsilon):
    p1_z = track_start[2]
    p2_z = track_end[2]

    # Using two abs() successively checks both anodes. 
    p1Cross = approx_equals(abs(p1_z), anode_z, epsilon)
    p2Cross = approx_equals(abs(p2_z), anode_z, epsilon)
    return p1Cross or p2Cross


# If have two TPCs and two anodes can check if cross both. 
def is_both_anodes_piercer(track_start, track_end, anode_z, epsilon):
    p1_z = track_start[2]
    p2_z = track_end[2]

    p1Cross = approx_equals(abs(p1_z), anode_z, epsilon)
    p2Cross = approx_equals(abs(p2_z), anode_z, epsilon)
    return p1Cross and p2Cross


def is_top_piercer(track_start, track_end):
    p1_y = track_start[1]
    p2_y = track_end[1]

    p1Cross = approx_equals(p1_y, top, epsilon)
    p2Cross = approx_equals(p2_y, top, epsilon)
    return p1Cross or p2Cross

def is_bottom_piercer(track_start, track_end):
    p1_y = track_start[1]
    p2_y = track_end[1]

    p1Cross = approx_equals(p1_y, bottom, epsilon)
    p2Cross = approx_equals(p2_y, bottom, epsilon)
    return p1Cross or p2Cross
        
def is_upstream_piercer(track_start, track_end):
    p1_x = track_start[0]
    p2_x = track_end[0]

    p1Cross = approx_equals(p1_x, upstream, epsilon)
    p2Cross = approx_equals(p2_x, upstream, epsilon)
    return p1Cross or p2Cross
        
def is_downstream_piercer(track_start, track_end):
    p1_x = track_start[0]
    p2_x = track_end[0]

    p1Cross = approx_equals(p1_x, downstream, epsilon)
    p2Cross = approx_equals(p2_x, downstream, epsilon)
    return p1Cross or p2Cross
        
def is_side_piercer(track_start, track_end):
    p1_x, p1_y, p1_z = track_start
    p2_x, p2_y, p2_z = track_end

    return ( is_downstream_piercer(track_start, track_end) or
             is_upstream_piercer(track_start, track_end) or
             is_bottom_piercer(track_start, track_end) or
             is_top_piercer(track_start, track_end) )

def track_selection(track_start, track_end, cutString):
    # This function is the ultimate filter.
    # Toggle desirable/undesirable conditions below.
    conditional_funcs = [is_good_track]

    if 'upstream' in cutString:
        conditional_funcs.append(is_upstream_piercer)
    if 'downstream' in cutString:
        conditional_funcs.append(is_downstream_piercer)
    if 'top' in cutString:
        conditional_funcs.append(is_top_piercer)
    if 'bottom' in cutString:
        conditional_funcs.append(is_bottom_piercer)
        
    return all(conditional(track_start, track_end)
               for conditional in conditional_funcs)
    
def main(args):
    global my_geometry
    global TPC_bounds, anode_z, cathode_z, top, bottom, upstream, downstream
    global length_cut
    global epsilon

    my_geometry = DetectorGeometry(args.detector, args.geometry)
    TPC_bounds = get_TPC_bounds()
    anode_z = TPC_bounds[1][2][0]
    cathode_z = TPC_bounds[1][2][1]
    top = TPC_bounds[0][1][1]
    bottom = TPC_bounds[0][1][0]
    upstream = TPC_bounds[0][0][1]
    downstream = TPC_bounds[0][0][0]
    length_cut = 100 #mm
    epsilon = 50 #mm
    # epsilon = 10 #mm

    global f
    f = h5py.File(args.infile, 'r')
    
    if args.plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits import mplot3d as plt3d
        from plotting import draw_boundaries, plot_selected_track
        
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
    events = f['events']

    if args.n >= 0:
        N = args.n
    else:
        N = len(tracks)

    corrected_endpoints = []
    
    for thisTrack in tracks[:N]:
        thisEvent = events[thisTrack['event_ref']]

        if thisEvent['n_ext_trigs'] < 2:
            continue
            # track_start,track_end = get_track_ends(thisTrack, my_geometry)
            # p1_z_shifted, p2_z_shifted = shift_track_to_anode(track_start,
            #                                                   track_end,
            #                                                   anode_z)
            # corrected_start = [track_start[0],
            #                    track_start[1],
            #                    p1_z_shifted]
            # corrected_end = [track_end[0],
            #                  track_end[1],
            #                  p2_z_shifted]
        else:
            corrected_start, corrected_end = get_track_ends(thisTrack, my_geometry)
            
        if track_selection(corrected_start, corrected_end, args.cut):
            if args.plot:
                plot_selected_track(ax,
                                    corrected_start,
                                    corrected_end)

            theseHits = f['hits'][thisTrack['hit_ref']]
            startHitPos, endHitPos = get_extreme_hit_pos(thisTrack,
                                                         theseHits,
                                                         my_geometry)

            corrected_endpoints.append(startHitPos)
            corrected_endpoints.append(endHitPos)
        # hits = f['hits'][thisTrack['hit_ref']]
        # plot_hits(ax, hits, my_geometry, thisTrack)
        # print ("this track passes the cuts: " + str(track_selection(thisTrack, f)))
        
    if args.plot:
        # Some helpful stats to print on the figure. 
        text = 'Anode + Downstream Crossers' #Change this to whichever selection is made
        ax.text2D(0.05, 0.0, text, transform=ax.transAxes)

    if args.output:
        np.save(args.output, corrected_endpoints)
    if args.plot:
        plt.show()

    f.close()

if __name__ == '__main__': 
    import argparse

    parser = argparse.ArgumentParser(description='Plot the first N tracks from a given file')
    parser.add_argument('infile',
                        help = 'input larpix data with track reconstruction')
    parser.add_argument('-n',
                        default = -1,
                        type = int,
                        help = 'evaluate the selection criteria over the first N tracks')
    parser.add_argument('-g', '--geometry',
                        default = './pixel_layouts/multi_tile_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d', '--detector',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-o', '--output',
                        default = '',
                        type = str,
                        help = 'save the data which passes the selection to a file')
    parser.add_argument('-p', '--plot',
                        default = False,
                        type = bool,
                        help = 'show a plot of the selected tracks')
    parser.add_argument('-c', '--cut',
                        default = 'upstream',
                        type = str,
                        help = 'faces which tracks must intersect (upstream, downstream, top, bottom, cathode)')

    args = parser.parse_args()

    main(args)
