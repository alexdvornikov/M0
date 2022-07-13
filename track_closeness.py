import numpy as np
import matplotlib.pyplot as plt
import h5py

from utils import *

def pairWiseDist(posPairArray):
    return np.sqrt(np.sum(np.power(np.diff(posPairArray, axis = 1), 2), axis = -1))

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

def closeness(trackA, trackB, metric = 'hit'):
    if metric == 'hit':
        hitsA = fA['hits'][trackA['hit_ref']]
        hitsB = fB['hits'][trackB['hit_ref']]

        evA = fA['events'][trackA['event_ref']]
        evB = fB['events'][trackB['event_ref']]

        posA = hit_to_3d(my_geometry, hitsA, trackA['t0'])
        posB = hit_to_3d(my_geometry, hitsB, trackB['t0'])

        posPairs = np.array([[thisPosA, thisPosB]
                             for thisPosA in np.array(posA).T
                             for thisPosB in np.array(posB).T])
        pwd = pairWiseDist(posPairs)
        d = np.min(pwd)
        posPairMin = np.argmin(pwd)
        posAMin, posBMin = posPairs[posPairMin]

        midpoint = np.mean([posAMin, posBMin], axis = 0)
        return d, midpoint

    elif metric == 'PCA':
        As = np.array([trackA['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                       trackA['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                       trackA['start'][2] + my_geometry.tpc_offsets[0][2]*10])
        Ae = np.array([trackA['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                       trackA['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                       trackA['end'][2] + my_geometry.tpc_offsets[0][2]*10])
        Ahat = As - Ae
        
        Bs = np.array([trackB['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                       trackB['start'][1] + my_geometry.tpc_offsets[0][1]*10,
                       trackB['start'][2] + my_geometry.tpc_offsets[0][2]*10])
        Be = np.array([trackB['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                       trackB['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                       trackB['end'][2] + my_geometry.tpc_offsets[0][2]*10])
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

        midpoint = np.mean([pocaA, pocaB], axis = 0)
        return np.power(np.dot(pocaA-pocaB, pocaA-pocaB), 0.5), midpoint
            
    else:
        raise ValueError ("this metric is not defined!")

def main(args):
    global fA
    fA = h5py.File(args.infileA)
    infileAtag = hash(args.infileA)
    
    global fB
    fB = h5py.File(args.infileB)
    infileBtag = hash(args.infileB)

    global my_geometry
    my_geometry = DetectorGeometry(args.detector, args.geometry)

    minLength = 100. 
    # minLength = 500.
   
    rawTracksA = np.array(fA['tracks'])
    # maskA = rawTracksA['length'] > 500.
    maskA = rawTracksA['length'] > minLength
    tracksA = rawTracksA[maskA]

    tracksAsubset = tracksA[args.b::args.N]

    rawTracksB = np.array(fB['tracks'])
    # maskB = rawTracksB['length'] > 500.
    maskB = rawTracksB['length'] > minLength
    tracksB = rawTracksB[maskB]

    # hitDists = []
    PCAdists = []
    # hitMid = []
    PCAmid = []
    Aid = []
    Bid = []
    Afile = []
    Bfile = []

    for trackA in tracksAsubset:
        evA = fA['events'][trackA['event_ref']]

        if evA['n_ext_trigs'] >= 2:
            for trackB in tracksB:
                evB = fB['events'][trackB['event_ref']]

                if evB['n_ext_trigs'] >= 2:
                    if (trackA['track_id'] > trackB['track_id']) or (fA != fB):
                        
                        Ai = trackA['track_id']
                        Bi = trackB['track_id']
                        # dHit, mpHit = closeness(trackA, trackB)
                        dPCA, mpPCA = closeness(trackA, trackB, metric = 'PCA')

                        Aid.append(Ai)
                        Bid.append(Bi)
                        
                        Afile.append(infileAtag)
                        Bfile.append(infileBtag)

                        # hitDists.append(dHit)
                        PCAdists.append(dPCA)

                        # hitMid.append(mpHit)
                        PCAmid.append(mpPCA)

                        if args.verbose:
                            print ("measuring track closeness between tracks:")
                            print (" ".join(["trackA:",
                                             str(trackA['track_id']),
                                             "nhits:",
                                             str(trackA['nhit'])]))
                            print (" ".join(["trackB:",
                                             str(trackB['track_id']),
                                             "nhits:",
                                             str(trackB['nhit'])]))
                            # print ("closeness (Hit): " + str(dHit))
                            print ("closeness (PCA): " + str(dPCA))
                    
    if args.o:
        # hitMid = np.array(hitMid)
        PCAmid = np.array(PCAmid)
        
        outArray = np.array([Afile, Bfile,
                             Aid, Bid,
                             # hitDists,
                             PCAdists,
                             # hitMid[:,0],
                             # hitMid[:,1],
                             # hitMid[:,2],
                             PCAmid[:,0],
                             PCAmid[:,1],
                             PCAmid[:,2]])
        np.save(args.o, outArray)
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
    parser.add_argument('infileA',
                        help='input larpix data')
    parser.add_argument('infileB',
                        help='input larpix data')
    parser.add_argument('-b',
                        default = 0,
                        type = int,
                        help='batch number (from 0 to N-1)')
    parser.add_argument('-N',
                        default = 1,
                        type = int,
                        help='number of batches to divide the job into')
    parser.add_argument('-g', '--geometry',
                        default = './pixel_layouts/multi_tile_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d', '--detector',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-o',
                        default = "",
                        type = str,
                        help='Output file to save the results.')
    parser.add_argument('-v',
                        '--verbose',
                        action = 'store_true',
                        help = 'Verbose')

    args = parser.parse_args()

    main(args)
