import numpy as np
import matplotlib.pyplot as plt
import h5py

vd = 1.6

def pairWiseDist(posPairArray):
    return np.sqrt(np.sum(np.power(posPairArray[:,1,:] - posPairArray[:,0,:], 2), axis = -1))

def closeness(trackA, trackB, h5File, metric = 'hit'):
    if metric == 'hit':
        hitsA = h5File['hits'][trackA['hit_ref']]
        hitsB = h5File['hits'][trackB['hit_ref']]

        t0A = trackA['t0']
        t0B = trackB['t0']

        posA = np.array([hitsA['px'],
                         hitsA['py'],
                         vd*(hitsA['ts'] - t0A)])
        posB = np.array([hitsB['px'],
                         hitsB['py'],
                         vd*(hitsB['ts'] - t0B)])

        posPairs = np.array([[thisPosA, thisPosB]
                             for thisPosA in posA.T
                             for thisPosB in posB.T])
        return np.min(pairWiseDist(posPairs))

    elif metric == 'PCA':
        # TODO implement segment-based distance
        return 0
    
    else:
        raise ValueError ("this metric is not defined!")

def main(args):
    f = h5py.File(args.infile)

    rawTracks = np.array(f['tracks'])

    mask = rawTracks['length'] > 100.

    tracks = rawTracks[mask]
    
    hitDists = []
    Aid = []
    Bid = []

    if args.n > 0:
        Ntracks = args.n
    else:
        Ntracks = tracks.shape[0]

    for trackA in tracks[:Ntracks]:
        for trackB in tracks:
            if trackA['track_id'] > trackB['track_id']:
                print ("measuring track closeness between tracks:")
                print (" ".join(["trackA:",
                                 str(trackA['track_id']),
                                 "nhits:",
                                 str(trackA['nhit'])]))
                print (" ".join(["trackB:",
                                 str(trackB['track_id']),
                                 "nhits:",
                                 str(trackB['nhit'])]))

                Ai = trackA['track_id']
                Bi = trackB['track_id']
                d = closeness(trackA, trackB, f)

                Aid.append(Ai)
                Bid.append(Bi)
                hitDists.append(d)
                print ("closeness: " + str(d))
                

    if args.o:
        np.savetxt(args.o, np.array([Aid, Bid, hitDists]))
    else:
        plt.figure()
        plt.hist(hitDists)

        plt.show()
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Find the track-to-track distance between pairs of tracks in a given dataset.')
    parser.add_argument('infile',
                        help='input larpix data')
    parser.add_argument('-n',
                        default = -1,
                        type = int,
                        help='examine only the first n tracks.')
    parser.add_argument('-o',
                        default = "",
                        type = str,
                        help='Output file to save the results.')

    args = parser.parse_args()

    main(args)
