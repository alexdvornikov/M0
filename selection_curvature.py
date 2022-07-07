# Example usage
# python3 selection_curvature.py events_2022_02_08_07_36_25_CET.gz.h5


# If interested in timiming the code...
# from datetime import datetime
# startTime = datetime.now()
# Also need to uncomment a datetime line above __main__ at the bottom of the code. 

from utils_m1 import *
    
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
    epsilon = 10 #mm

    global f
    f = h5py.File(args.infile, 'r')
    data = Data(f)
    #----------------------------------------------------------------------------#            
    trackMask = (data.rawTracks['length'] > length_cut) & (data.rawTracks['nhit'] > 0)
    track_idx = np.where(trackMask)[0]

    if args.n >= 0:
        N = args.n
    else:
        N = len(track_idx)

    output = []

    for thisTrack_idx in track_idx[:N]:
        thisTrack = data.rawTracks[thisTrack_idx]
        thisEvent = dereference(thisTrack_idx, data.track_ref, data.rawEvents, region=data.track_reg, ref_direction=(1,0))
        evt_idx = thisEvent['id']
        t0_dtype = dereference(evt_idx, data.t0_ref, data.rawT0, region=data.t0_reg, ref_direction=(1,0))
        t0_unflattened = t0_dtype['ts'] #Light trigger
        t0 = t0_unflattened[0][0]

        if thisEvent['n_ext_trigs'] < 2:
            continue
        #----------------------------------------------------------------------------#
        # Get hits
        #----------------------------------------------------------------------------#
        theseHits = dereference(thisTrack_idx, data.hit_ref, data.rawHits, region=data.hit_reg, ref_direction=(0,1))
        pos3d = hit_to_3d(my_geometry, theseHits[0], t0)
        pos3d = np.array([np.array(entry) for entry in pos3d]) #Nested list to nested array
        #----------------------------------------------------------------------------#
        # Get endpoints
        #----------------------------------------------------------------------------#
        startHitPos, endHitPos = get_extreme_hit_pos(t0,thisTrack,theseHits[0],my_geometry)
        # Focus on a single TPC for now (avoid positive drift coordinates)
        if startHitPos[2] > 0 or endHitPos[2] > 0:
            continue
        #----------------------------------------------------------------------------#
        # Check anode crossing and plot
        #----------------------------------------------------------------------------#
        if is_anode_piercer(startHitPos, endHitPos, anode_z, epsilon):

            # For more on get_pca_endpts see utils_m1.py
            ds = distortions(t0, my_geometry, theseHits[0], pos3d, near_anode = True, nhit = 10)
            output.append(ds)

    # if args.output:
    #     np.save(args.output, output)

    #----------------------------------------------------------------------------#
    # Get reco and true coords
    #----------------------------------------------------------------------------#
    # print( len(output) )
    concat = lambda key : np.concatenate([out[key] for out in output])
    reco_coords = concat('reco')
    true_coords = concat('true')

    x, y, z = reco_coords.T
    dx, dy, dz = (reco_coords - true_coords).T

    #----------------------------------------------------------------------------#
    # Bin the offsets and save the counts and bin centers (hexbin)
    #----------------------------------------------------------------------------#
    global plt
    offset_range = 20
    hb_dx = plt.hexbin(z/10, dx/10, extent = (downstream/10, 0, -offset_range, offset_range), mincnt=1, gridsize = 50, bins='log')
    hb_dy = plt.hexbin(z/10, dy/10, extent = (downstream/10, 0, -offset_range, offset_range), mincnt=1, gridsize = 50, bins='log')

    hb_counts = []
    hb_pos = []
    # Centers of hexbins for dx offsets 
    xs_dx = hb_dx.get_offsets()[:, 0] 
    ys_dx = hb_dx.get_offsets()[:, 1]
    # Centers of hexbins for dy offsets
    xs_dy = hb_dy.get_offsets()[:, 0] 
    ys_dy = hb_dy.get_offsets()[:, 1]
    # Get dx vals
    arr_dx = hb_dx.get_array() # Masked array of counts
    counts_dx = arr_dx.data # Grab the data, assuming no mask
    # Get dy vals
    arr_dy = hb_dy.get_array() 
    counts_dy = arr_dy.data 

    hb_counts.append(np.asarray(counts_dx))
    hb_counts.append(np.asarray(counts_dy))

    hb_pos.append(xs_dx)
    hb_pos.append(ys_dx)

    hb_pos.append(xs_dy)
    hb_pos.append(ys_dy)

    hb_counts = np.array( hb_counts, dtype=object )

    if args.output: 
        np.save('hb_counts.npy', hb_counts)
        np.save('hb_pos.npy', hb_pos)

    f.close()

    # Show time elapsed for running the code
    # print(datetime.now() - startTime)

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
                        default = './pixel_layouts/module1_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d', '--detector',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-o', '--output',
                        default = True,
                        type = bool,
                        help = 'save the data which passes the selection to a file')

    args = parser.parse_args()

    main(args)
