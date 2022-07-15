# python curvature_selection_batch_manager_m1.py -i /global/project/projectdirs/dune/www/data/Module1/reco/charge_only/ -r /global/project/projectdirs/dune/www/data/Module1/runlist.txt -o /global/project/projectdirs/dune/users/olexiy/m1_2anode_histos -n 10

import os
import argparse

def main(indir, runlist, nfiles, outdir):
    print ("input directory: " + indir)
    print ("runlist file: " + runlist)
    print ("to be processed: " + str(nfiles))
    print ("output directory: " + outdir)

    print ("reading runlist...")
    with open(runlist) as runlistText:
        lines = runlistText.read().split('\n')
        #M1 runlist is all strings
        runMetaData = [{"e_field": line.split(' ')[0],
                        "charge_filename": line.split(' ')[1],
                        "light_filename": line.split(' ')[2],
                        "charge_thresholds": line.split(' ')[3],
                        "light_samples": line.split(' ')[4]}
                       for line in lines[1:] if line]

    runConditions = {"e_field": "500",
                     "charge_thresholds": "low"}

    if nfiles < 0:
        nfiles = len(runMetaData)
        
    n_launched = 0
    for thisRun in runMetaData:
        # print(thisRun)
        conditions_met = [thisRun[key] == value
                          for key, value in runConditions.items()]
        if all(conditions_met):
            if n_launched < nfiles:
                    
                rel_infilename = 'events_'+thisRun['charge_filename']+'.gz.h5' #module1
                # rel_infilename = 'datalog_'+thisRun['charge_filename']+'evd.h5' #module0 
                infileName = os.path.join(indir, rel_infilename)

                rel_outfileName1 = 'hb_counts_2anodes_'+thisRun['charge_filename']+'.npy'
                rel_outfileName2 = 'hist2d_zx_2_anodes_'+thisRun['charge_filename']+'.npy'
                rel_outfileName3 = 'hist2d_zy_2anodes_'+thisRun['charge_filename']+'.npy'
                outfileName1 = os.path.join(outdir, rel_outfileName1)
                outfileName2 = os.path.join(outdir, rel_outfileName2)
                outfileName3 = os.path.join(outdir, rel_outfileName3)

                # Run sbatch script for module1. 
                sbatch_cmd = " ".join(["sbatch ./curvature_selection_batch_m1.sh",
                                       infileName,
                                       outfileName1,
                                       outfileName2,
                                       outfileName3])

                os.system(sbatch_cmd)

                n_launched += 1

                print ("input: " + infileName)
                print ("output1: " + outfileName1)
                print ("output2: " + outfileName2)
                print ("output3: " + outfileName3)
                print ("sbatch command: " + sbatch_cmd)
                print ()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--indir','-i',
			required=True,
			type=str)
    parser.add_argument('--runlist','-r',
			required=True,
			type=str)
    parser.add_argument('--nfiles','-n',
			default = -1,
			type=int)
    parser.add_argument('--outdir','-o',
			required=True,
			type=str)

    args = parser.parse_args()

    main(**vars(args))
