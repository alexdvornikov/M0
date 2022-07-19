# Example use
# python3 selection_curvature.py events_2022_02_08_07_36_25_CET.gz.h5 -o1 hb_counts.npy -o2 hist2d_zx.npy -o3 hist2d_zy.npy
# python3 selection_curvature.py events_2022_02_08_07_36_25_CET.gz.h5 -o1 hb_counts_2anodes.npy -o2 hist2d_zx_2anodes.npy -o3 hist2d_zy_2anodes.npy

# python3 selection_curvature.py events_2022_02_08_07_36_25_CET.gz.h5 -o4 hist3d.npy
# python3 selection_curvature.py ~/Desktop/events_2022_02_09_17_23_09_CET.gz.h5 -o4 hist3d.npy


# If interested in timiming the code...
from datetime import datetime
startTime = datetime.now()
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
    cm = 0.1 #Conversion from mm to cm

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


    

    print('Number of tracks in file: ' + str( N ))
    i = 0
    n_tracks = 0
    n_selected_tracks = 0
    check_points = np.linspace(0,N,10) #for progress bar

    for thisTrack_idx in track_idx[:N]:

        n_tracks += 1

        #----------------------------------------------
        # Progress bar (10% increments)
        #----------------------------------------------
        if n_tracks % round( check_points[1] ) == 0:
            i += 1
            print('Completed ' + str(10*i) + '%')
            n_tracks = 0
        #----------------------------------------------


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

        #----------------------------------------------------------------------------#
        # July 14, 2022
        #----------------------------------------------------------------------------#
        if is_both_anodes_piercer(startHitPos, endHitPos, anode_z, epsilon):

            # For more on get_pca_endpts see utils_m1.py
            ds = distortions_2anodes(t0, my_geometry, pos3d, startHitPos, endHitPos)
            # ds = distortions(t0, my_geometry, theseHits[0], pos3d, near_anode = True, nhit = 10)
            output.append(ds)
            n_selected_tracks += 1

        # #----------------------------------------------------------------------------#
        # # Focus on a single TPC for now (left of the cathode)
        # if startHitPos[2] > cathode_z or endHitPos[2] > cathode_z:
        #     continue
        # #----------------------------------------------------------------------------#
        # # Check anode crossing and plot
        # #----------------------------------------------------------------------------#
        # if is_anode_piercer(startHitPos, endHitPos, anode_z, epsilon):

        #     # For more on get_pca_endpts see utils_m1.py
        #     ds = distortions(t0, my_geometry, theseHits[0], pos3d, near_anode = True, nhit = 10)
        #     output.append(ds)
        #     n_selected_tracks += 1
        # #----------------------------------------------------------------------------#
            

    print('Number of selected tracks: ' + str(n_selected_tracks))
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
    offset_range = 5 # [cm]
    # offset_range = 20 # [cm]
    # hex_bins = ( round(2*anode_z*cm), round( 2*offset_range*cm) )
    grdsize = 2*round(2*anode_z*cm) #If have ~2 pixels per cm can allow for such binning
    hb_dx = plt.hexbin(z*cm, dx*cm, extent = (-anode_z*cm, anode_z*cm, -offset_range, offset_range), mincnt=0, gridsize = grdsize, bins='log')
    hb_dy = plt.hexbin(z*cm, dy*cm, extent = (-anode_z*cm, anode_z*cm, -offset_range, offset_range), mincnt=0, gridsize = grdsize, bins='log')
    # hb_dx = plt.hexbin(z*cm, dx*cm, extent = (-anode_z*cm, cathode_z*cm, -offset_range, offset_range), mincnt=0, gridsize = grdsize, bins='log')
    # hb_dy = plt.hexbin(z*cm, dy*cm, extent = (-anode_z*cm, cathode_z*cm, -offset_range, offset_range), mincnt=0, gridsize = grdsize, bins='log')

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
    #----------------------------------------------------------------------------#

    #----------------------------------------------------------------------------#
    # 2d histogram for dx, dy offsets in the zx plane (also need yz plane)
    #----------------------------------------------------------------------------#
    zbins = round( (2*anode_z)*cm) #If crossing anode to anode
    # zbins = round( (anode_z - cathode_z)*cm) #If using single TPC
    xbins = round( (upstream + abs(downstream))*cm )
    ybins = round( (top + abs(bottom))*cm )
    # binned_statistic_2d by default gives the mean in each bin (different from hist2d)
     #----------------------------------------------------------------------------#
    def get_hist(drift_coord, transverse_coord, dx, dy, extent, bns):
        dx_hist = binned_statistic_2d(drift_coord, transverse_coord, dx, bins=bns, range=extent)
        dy_hist = binned_statistic_2d(drift_coord, transverse_coord, dy, bins=bns, range=extent)
        bin_means = np.array( [dx_hist.statistic, dy_hist.statistic] )

        # bin_edges_dx = np.array( [dx_hist.x_edge,dx_hist.y_edge], dtype=object )
        # bin_edges_dy = np.array( [dy_hist.x_edge,dy_hist.y_edge], dtype=object )

        # Bin edges should be the same for all of the files. 
        # Only need to get these once. 
        # Also the dx and dy edges should be the same just being extra looking at both. 
        # return bin_means, [bin_edges_dx, bin_edges_dy]
        return bin_means
    #----------------------------------------------------------------------------#
    bns = [zbins,xbins]
    extent = [(-anode_z*cm,anode_z*cm), (downstream*cm,upstream*cm)]
    # extent = [(-anode_z*cm,cathode_z*cm), (downstream*cm,upstream*cm)]
    zxmeans = get_hist(z*cm,x*cm,dx*cm,dy*cm, extent, bns)

    bns = [zbins,ybins]
    extent = [(-anode_z*cm,anode_z*cm), (bottom*cm,top*cm)]
    # extent = [(-anode_z*cm,cathode_z*cm), (bottom*cm,top*cm)]
    zymeans = get_hist(z*cm,y*cm,dx*cm,dy*cm, extent, bns)

    #----------------------------------------------------------------------------#
    # 3d histograms (for quiver plot)
    #----------------------------------------------------------------------------#
    def get_3dhist(x,y,z, dx, dy, dz, extent, bns):
        s = (x,y,z)
        dx_hist = binned_statistic_dd(sample = s, values = dx, bins=bns, range=extent)
        dy_hist = binned_statistic_dd(sample = s, values = dy, bins=bns, range=extent)
        dz_hist = binned_statistic_dd(sample = s, values = dz, bins=bns, range=extent)
        bin_means = np.array( [dx_hist.statistic, dy_hist.statistic, dz_hist.statistic] )
        return bin_means


    bns = [xbins, ybins, zbins]
    print( bns )
    extent = [ (downstream*cm,upstream*cm), (bottom*cm,top*cm), (-anode_z*cm,anode_z*cm)]
    xyzmeans = get_3dhist(x,y,z,dx,dy,dz, extent, bns)
    #----------------------------------------------------------------------------#

    if args.output4: 

        np.save(args.output4, xyzmeans)

        # if args.edges:
        #     np.save('hb_pos_2anodes.npy', hb_pos) #hexbin edges 

    if args.output1: 

        # # Save hexbin outputs
        # np.save('hb_counts.npy', hb_counts)
        # # np.save('hb_pos.npy', hb_pos)

        # # Save binned_statistic_2d outputs
        # np.save('hist2d_zx.npy', zxmeans)
        # np.save('hist2d_zy.npy', zymeans)
        # # np.save('hist2d_edges.npy', bin_edges_dx)

        np.save(args.output1, hb_counts)
        np.save(args.output2, zxmeans)
        np.save(args.output3, zymeans)

        np.save(args.output4, xyzmeans)

        if args.edges:
            np.save('hb_pos_2anodes.npy', hb_pos) #hexbin edges 
            # Don't need the binned_statistic_2d edges since using imshow extent. 
            # Same for all files. 

    f.close()


    # Show time elapsed for running the code
    print(datetime.now() - startTime)

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
    parser.add_argument('-o1', '--output1',
                        default = '',
                        type = str,
                        help = '1st output file')
    parser.add_argument('-o2', '--output2',
                        default = '',
                        type = str,
                        help = '2nd output file')
    parser.add_argument('-o3', '--output3',
                        default = '',
                        type = str,
                        help = '3rd output file')
    parser.add_argument('-o4', '--output4',
                        default = '',
                        type = str,
                        help = '4th output file')
    parser.add_argument('-e', '--edges',
                        default = False,
                        type = bool,
                        help = 'Save histogram bin edges')

    args = parser.parse_args()

    main(args)
