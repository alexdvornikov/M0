import numpy as np
import h5py
from utils_m1 import *

#Peter's library (get from GitHub)
from h5flow.data import dereference
    
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

    global f
    f = h5py.File(args.infile, 'r')
    data = Data(f)
    #----------------------------------------------------------------------------#
    
    if args.plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits import mplot3d as plt3d
        from plotting import draw_boundaries
        
        fig = plt.figure() 
        ax = fig.add_subplot(111, projection = '3d') 
        ax.set_xlabel(r'x (horizontal) [mm]')
        ax.set_ylabel(r'y (vertical) [mm]')
        ax.set_zlabel(r'z (drift) [mm]')
        # Set background to white
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        draw_boundaries(ax)

    trackMask = (data.rawTracks['length'] > length_cut) & (data.rawTracks['nhit'] > 0)
    track_idx = np.where(trackMask)[0]

    if args.n >= 0:
        N = args.n
    else:
        N = len(track_idx)

    corrected_endpoints = []
    
    for thisTrack_idx in track_idx[:N]:
        thisTrack = data.rawTracks[thisTrack_idx]
        thisEvent = dereference(thisTrack_idx, data.track_ref, data.rawEvents, region=data.track_reg, ref_direction=(1,0))
        evt_idx = thisEvent['id']
        t0_dtype = dereference(evt_idx, data.t0_ref, data.rawT0, region=data.t0_reg, ref_direction=(1,0))
        t0_unflattened = t0_dtype['ts'] #Light trigger
        t0 = t0_unflattened[0][0]

        if thisEvent['n_ext_trigs'] < 2:
            continue
        else:
            corrected_start, corrected_end = get_track_ends(thisTrack, my_geometry)
        #----------------------------------------------------------------------------#
        # Get all of the endpoints.
        #----------------------------------------------------------------------------#
        theseHits = dereference(thisTrack_idx, data.hit_ref, data.rawHits, region=data.hit_reg, ref_direction=(0,1))
        startHitPos, endHitPos = get_extreme_hit_pos(t0,thisTrack,theseHits[0],my_geometry)
        corrected_endpoints.append(startHitPos)
        corrected_endpoints.append(endHitPos)
        #----------------------------------------------------------------------------#
        
    corrected_endpoints = delete_nans(corrected_endpoints)

    if args.plot:
        ax.scatter(*corrected_endpoints.T)

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

    args = parser.parse_args()

    main(args)
