import numpy as np
import matplotlib.pyplot as plt
import h5py

from utils import *

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

        # if the point is inside of both segments, that's it!
        # if the point is outside of at least one segment,
        # take the segment which ends farthest from its point
        if l < 0:
            dOutsideA = np.abs(l)
            pocaA = Ae
        elif l > 1:
            dOutsideA = l - 1
            pocaA = As
        else:
            dOutsideA = 0
            pocaA = Ae + l*Ahat

        if s < 0:
            dOutsideB = np.abs(s)
            pocaB = Be
        elif s > 1:
            dOutsideB = s - 1
            pocaB = Bs
        else:
            dOutsideB = 0
            pocaB = Be + s*Bhat
       
        if dOutsideA > dOutsideB:
            # then, using that as the new point, re-calculate for the other
            s = np.dot(Bhat, pocaA - Be)/np.dot(Bhat, Bhat)
            # if the point of the other is still outside of the segment,
            # then take the end which is closest
            if s < 0:
                pocaB = Be
            elif s > 1:
                pocaB = Bs
            else:
                pocaB = Be + s*Bhat

        elif dOutsideB > dOutsideA:
            l = np.dot(Ahat, pocaB - Ae)/np.dot(Ahat, Ahat)
            if l < 0:
                pocaA = Ae
            elif l > 1:
                pocaA = As
            else:
                pocaA = Ae + l*Ahat

        return np.power(np.dot(pocaA-pocaB, pocaA-pocaB), 0.5)
            
    else:
        raise ValueError ("this metric is not defined!")

def main(args):
    f = h5py.File(args.infile)

    rawTracks = np.array(f['tracks'])

    mask = rawTracks['length'] > 500.

    tracks = rawTracks[mask]
    
    hitDists = []
    PCAdists = []
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
                PCAdists.append(dPCA)
                print ("closeness: " + str(d))
                print ("closeness (PCA): " + str(dPCA))
                    
    if args.o:
        np.savetxt(args.o, np.array([Aid, Bid, hitDists]))
    else:
        plt.figure()
        bins = np.concatenate([np.linspace(0, 100, 10)[:-1],
                               np.linspace(100, 1200, 20)])
        plt.hist(hitDists,
                 histtype = 'step',
                 label = 'hit-to-hit',
                 density = True,
                 bins = bins)
        plt.hist(PCAdists,
                 histtype = 'step',
                 label = 'axis-to-axis',
                 density = True,
                 bins = bins)
        plt.xlabel(r'track-to-track distance (mm per bin width)')
        plt.legend()

        plt.figure()
        diff = np.abs(np.array(hitDists) - np.array(PCAdists))
        print ('mean diff: ', np.mean(diff))
        plt.hist(diff,
                 histtype = 'step',
                 density = True)
        plt.xlabel(r'distance metric disagreement (mm)')
        plt.semilogy()
        
        plt.figure()
        plt.scatter(hitDists, PCAdists)
        plt.xlabel(r'hit-to-hit distance (mm)')
        plt.ylabel(r'axis-to-axis distance (mm)')

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
