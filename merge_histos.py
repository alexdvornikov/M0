# Output three histograms. 
# One hexbin and two binned_statistic_2d. 

#python merge_histos.py m1_histos -o1 hb_counts.npy -o2 hist2d_zx.npy -o3 hist2d_zy.npy

import os 
import argparse
import numpy as np

def main(indir,output1,output2,output3):

    #----------------------------------------------------------------------------#
    # Get array shape from the first file 
    # Need this to initialize an empty array below
    #----------------------------------------------------------------------------#
    i = 0
    array_shape = 0
    for file in os.listdir(indir):
        hb_counts = np.load(indir + '/' + file, allow_pickle=True)
        array_shape = hb_counts.shape
        i += 1
        if i >= 1:
            break
    #----------------------------------------------------------------------------#

    #----------------------------------------------------------------------------#
    # Initialize empty array and iteratively add to it.
    #----------------------------------------------------------------------------#
    n_files = 0
    counts = np.zeros(array_shape)
    for file in os.listdir(indir):
        thisCount = np.load(indir + '/' + file, allow_pickle=True)
        counts = counts + thisCount #np.add() if same shape
        n_files += 1
    #----------------------------------------------------------------------------#


    if output1: #hexbin merging 
        np.save(args.output1, counts)


    # binned_statistic_2d merging for one face
    if output2: 
        # Now counts are actually means. Need to do a final mean as well. 
        mean_counts = counts/n_files
        np.save(args.output2, mean_counts)

    # binned_statistic_2d merging for another face
    if output3: 
        # Now counts are actually means. Need to do a final mean as well. 
        mean_counts = counts/n_files
        np.save(args.output3, mean_counts)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge numpy arrays')

    parser.add_argument('--indir','-i',required=True,type=str)
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

    args = parser.parse_args()
    main(**vars(args))