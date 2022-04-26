import numpy as np
import matplotlib.pyplot as plt
import h5py

vd = 1.62 # mm/us
drift_dist = 310. # mm
clock_interval = 0.1 # us

class POCA:
    def __init__(self, d, pointA, pointB):
        """
        This class contains data about a point of closest approach (POCA)
        d: distance of closest approach
        pointA: point from one track for which distance is minimized
        pointB: point from other track for which distance is minimized
        """

        self.d = d
        self.pointA = pointA
        self.pointB = pointB

def pairWiseDist(posPairArray):
    return np.sqrt(np.sum(np.power(posPairArray[:,1,:] - posPairArray[:,0,:], 2), axis = -1))

def hit_to_3d(hits, event):
    t0 = event['ts_start']

    x = hits['px']
    y = hits['py']
    t = clock_interval*(hits['ts'] - t0)
    grp = hits['iogroup']
    parity = np.power(-1, grp)
    z = parity*(drift_dist - t*vd)

    pos3d = np.array([x, y, z])
    return pos3d

def closeness(trackA, trackB, h5File, metric = 'hit'):
    if metric == 'hit':
        hitsA = h5File['hits'][trackA['hit_ref']]
        hitsB = h5File['hits'][trackB['hit_ref']]

        evA = h5File['events'][trackA['event_ref']]
        evB = h5File['events'][trackB['event_ref']]

        posA = hit_to_3d(hitsA, evA)
        posB = hit_to_3d(hitsB, evB)

        posPairs = np.array([[thisPosA, thisPosB]
                             for thisPosA in posA.T
                             for thisPosB in posB.T])
        pwd = pairWiseDist(posPairs)
        d = np.min(pwd)
        posPairMin = np.argmin(pwd)
        posAMin, posBMin = posPairs[posPairMin]
        return (d)

    elif metric == 'PCA':
        # TODO implement segment-based distance
        As = np.array([trackA['start'][0],
                       trackA['start'][1],
                       trackA['start'][2]])
        Ae = np.array([trackA['end'][0],
                       trackA['end'][1],
                       trackA['end'][2]])
        Ahat = As - Ae
        
        Bs = np.array([trackB['start'][0],
                       trackB['start'][1],
                       trackB['start'][2]])
        Be = np.array([trackB['end'][0],
                       trackB['end'][1],
                       trackB['end'][2]])
        Bhat = Bs - Be
                
        det = np.dot(Ahat, Ahat)*np.dot(Bhat, Bhat) - np.power(np.dot(Ahat, Bhat), 2)

        # length along the 'A' track.  Within the track is l in [0, 1]
        l = 1./det*(np.dot(Bhat, Bhat)*(np.dot(Be, Ahat) - np.dot(Ae, Ahat)) + np.dot(Ahat, Bhat)*(np.dot(Ae, Bhat) - np.dot(Be, Bhat)))

        # length along the 'B' track.  Within the track is s in [0, 1]
        s = 1./det*(np.dot(Ahat, Ahat)*(np.dot(Ae, Bhat) - np.dot(Be, Bhat)) + np.dot(Ahat, Bhat)*(np.dot(Be, Ahat) - np.dot(Ae, Ahat)))

        if l < 0:
            pocaA = Ae
        elif l > 1:
            pocaA = As
        else:
            pocaA = Ae + l*Ahat
            
        if s < 0:
            pocaB = Be
        elif s > 1:
            pocaB = Bs
        else:
            pocaB = Be + s*Bhat

        return np.power(np.dot(pocaA-pocaB, pocaA-pocaB), 0.5)
            
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
                dPCA = closeness(trackA, trackB, f, metric = 'PCA')

                Aid.append(Ai)
                Bid.append(Bi)
                hitDists.append(d)
                print ("closeness: " + str(d))
                print ("closeness (PCA): " + str(dPCA))
                

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
