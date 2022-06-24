import os 
import argparse
import numpy as np

def main(indir,outfile):
    merged_array = np.concatenate([np.load(indir + '/' + i) for i in os.listdir(indir)])
    np.save(outfile, merged_array)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge numpy arrays')

    parser.add_argument('--indir','-i',
			required=True,
			type=str)
    parser.add_argument('--outfile','-o',
			required=True,
			type=str,
            help = 'Needs to be .npy extension')

    args = parser.parse_args()
    main(**vars(args))