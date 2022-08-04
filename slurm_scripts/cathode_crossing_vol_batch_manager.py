import os
import argparse

def main(indir, runlist, nfiles, outdir):
# def main(indir, runlist, cut, nfiles, outdir):
    print ("input directory: " + indir)
    print ("runlist file: " + runlist)
    # print ("cut string: " + cut)
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

        # Module 0
        # runMetaData = [{"e_field": float(line.split(' ')[0]),
        #                 "charge_filename": line.split(' ')[1],
        #                 "light_filename": line.split(' ')[2],
        #                 "charge_thresholds": line.split(' ')[3],
        #                 "light_samples": int(line.split(' ')[4])}
        #                for line in lines[1:] if line]

    # runConditions = {"e_field": 500}
    #M1 data only has low thresholds

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

                rel_outfileName = 'selected_'+thisRun['charge_filename']+'.npz' #Compressed
                # rel_outfileName = 'selected_'+thisRun['charge_filename']+'.npy'
                outfileName = os.path.join(outdir, rel_outfileName)

                # Run sbatch script for module1. 
                sbatch_cmd = " ".join(["sbatch ./cathode_crossing_vol_batch.sh",
                                       infileName,
                                       outfileName])

                # Run sbatch script for module0. 
                # sbatch_cmd = " ".join(["sbatch ./selection_batch.sh",
                #                        infileName,
                #                        outfileName,
                #                        cut])

                os.system(sbatch_cmd)

                n_launched += 1

                print ("input: " + infileName)
                print ("output: " + outfileName)
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
    # parser.add_argument('--cut','-c',
	# 		required=True,
	# 		type=str)
    parser.add_argument('--nfiles','-n',
			default = -1,
			type=int)
    parser.add_argument('--outdir','-o',
			required=True,
			type=str)

    args = parser.parse_args()

    main(**vars(args))
