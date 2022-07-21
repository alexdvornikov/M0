# python3 get_selection.py events_2022_02_08_07_36_25_CET.gz.h5 -o selection_test.npz

# If interested in timiming the code...
from datetime import datetime
startTime = datetime.now()
# Also need to uncomment a datetime line above __main__ at the bottom of the code. 

# from tqdm import tqdm # If want progress bar use tqdm
from collections import defaultdict
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
    #-----------------------------------------------------------------------------------------#         
    trackMask = (data.rawTracks['length'] > length_cut) & (data.rawTracks['nhit'] > 0)
    track_idx = np.where(trackMask)[0]

    if args.n >= 0:
        N = args.n
    else:
        N = len(track_idx)


    output = defaultdict(list) # Initialize output dictionary 
    print('Number of tracks in file: ' + str( N ))
    i = 0
    n_selected_tracks = 0

    for thisTrack_idx in track_idx[:N]:
    # for thisTrack_idx in tqdm(track_idx[:N]): #If want progress bar use tqdm

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
        # Check if crossing both anodes
        #----------------------------------------------------------------------------#
        if is_both_anodes_piercer(startHitPos, endHitPos, anode_z, epsilon):

            output['pos'].append( pos3d )
            output['start'].append( startHitPos )
            output['end'].append( endHitPos )

            n_selected_tracks += 1


    print('Number of selected tracks: ' + str(n_selected_tracks))
    if args.output:
        # np.save(args.output, dict(output) ) 
        np.savez_compressed(args.output, data= dict(output) )


    # Close the data file
    f.close() 
    # Show time elapsed for running the code
    print(datetime.now() - startTime)
#-----------------------------------------------------------------------------------------------#
# Arguments through the command line 
#-----------------------------------------------------------------------------------------------#
if __name__ == '__main__': 
    import argparse

    parser = argparse.ArgumentParser(description='Get track selection.')
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
                        default = '',
                        type = str,
                        help = 'output file')
    args = parser.parse_args()

    main(args)
#-----------------------------------------------------------------------------------------------#
