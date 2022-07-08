# Output three histograms. 
# One hexbin and two binned_statistic_2d. 

# python3 merge_histos.py -i /Users/alex/Desktop/hists/hb -o1 hb_counts_merged.npy 
# python3 merge_histos.py -i /Users/alex/Desktop/hists/zx -o2 hist2d_zx_merged.npy 
# python3 merge_histos.py -i /Users/alex/Desktop/hists/zy -o3 hist2d_zy_merged.npy

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
        counts = np.load(indir + '/' + file, allow_pickle=True)
        array_shape = counts.shape
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
        print(thisCount.shape)
        print(counts.shape)
        counts = counts + thisCount #np.add() if same shape
        n_files += 1
    #----------------------------------------------------------------------------#


    if output1: #hexbin merging 
        np.save(args.output1, counts)


    # binned_statistic_2d merging for zx face
    if output2: 
        # Now counts are actually means. Need to do a final mean as well. 
        mean_counts = counts/n_files
        np.save(args.output2, mean_counts)

    # binned_statistic_2d merging for zy face
    if output3: 
        # Now counts are actually means. Need to do a final mean as well. 
        mean_counts = counts/n_files
        np.save(args.output3, mean_counts)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge numpy arrays')

    parser.add_argument('--indir','-i',
                        required=True,
                        type=str)
    parser.add_argument('-o1', '--output1',
                        default = '', # hb_counts.npy
                        type = str,
                        help = '1st output file')
    parser.add_argument('-o2', '--output2',
                        default = '', #hist2d_zx.npy
                        type = str,
                        help = '2nd output file')
    parser.add_argument('-o3', '--output3',
                        default = '', #hist2d_zy.npy
                        type = str,
                        help = '3rd output file')

    args = parser.parse_args()
    main(**vars(args))

#python merge_histos.py m1_histos -o1 hb_counts.npy -o2 hist2d_zx.npy -o3 hist2d_zy.npy