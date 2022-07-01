# python3 selection_batch_manager.py -i /global/project/projectdirs/dune/www/data/Module0/TPC1+2/dataRuns/tracksData -r /global/project/projectdirs/dune/www/data/Module0/runlist.txt -o /global/project/projectdirs/dune/users/olexiy/M0/crossers/upstream -c upstream -n 5

import os
import argparse

def main(indir, runlist, nJobs, outdir):
    print ("input directory: " + indir)
    print ("runlist file: " + runlist)
    print ("to be processed: " + str(nJobs))
    print ("output directory: " + outdir)

    print ("reading runlist...")
    with open(runlist) as runlistText:
        lines = runlistText.read().split('\n')
        runMetaData = [{"e_field": float(line.split(' ')[0]),
                        "charge_filename": line.split(' ')[1],
                        "light_filename": line.split(' ')[2],
                        "charge_thresholds": line.split(' ')[3],
                        "light_samples": int(line.split(' ')[4])}
                       for line in lines[1:] if line]

    runConditions = {"e_field": 500,
                     "charge_thresholds": "high"}

    goodRuns = []
    
    if nJobs < 0:
        nJobs = 999999999999
        
    n_launched = 0
    for thisRun in runMetaData:
        conditions_met = [thisRun[key] == value
                          for key, value in runConditions.items()]
        if all(conditions_met) and os.path.exists(thisRun['charge_filename']):
            goodRuns.append(thisRun)

    for thisRunA in goodRuns:
        for thisRunB in goodRuns:
            if hash(thisRunA['charge_filename']) >= hash(thisRunB['charge_filename']):
                if n_launched < nJobs:
                    
                    rel_infilenameA = 'datalog_'+thisRunA['charge_filename']+'evd.h5'
                    infileNameA = os.path.join(indir, rel_infilenameA)

                    rel_infilenameB = 'datalog_'+thisRunB['charge_filename']+'evd.h5'
                    infileNameB = os.path.join(indir, rel_infilenameB)
                    
                    rel_outfileName = '_'.join(['crossing',
                                                thisRunA['charge_filename'],
                                                thisRunB['charge_filename'],
                                                '.npy'])
                    outfileName = os.path.join(outdir, rel_outfileName)

                    sbatch_cmd = " ".join(["sbatch ./closeness_batch.sh",
                                           infileNameA,
                                           infileNameB,
                                           outfileName])
                    os.system(sbatch_cmd)
                    
                    n_launched += 1

                    print ("input A: " + infileNameA)
                    print ("input B: " + infileNameB)
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
    parser.add_argument('--nJobs','-n',
			default = -1,
			type=int)
    parser.add_argument('--outdir','-o',
			required=True,
			type=str)

    args = parser.parse_args()

    main(**vars(args))
