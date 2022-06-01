import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import *
from plotting import *

def main(args):
    global my_geometry
    global TPC_bounds, anode_z, cathode_z, top, bottom, upstream, downstream
    global passing_counter
    global length_cut
    global epsilon

    my_geometry = DetectorGeometry(args.d, args.g)
    TPC_bounds = get_TPC_bounds()
    anode_z = TPC_bounds[1][2][0]
    cathode_z = TPC_bounds[1][2][1]
    top = TPC_bounds[0][1][1]
    bottom = TPC_bounds[0][1][0]
    upstream = TPC_bounds[0][0][1]
    downstream = TPC_bounds[0][0][0]
    length_cut = 100 #mm
    # epsilon = 50 #mm
    epsilon = 10 #mm

    endPoints = np.concatenate([np.load(infileName)
                                for infileName in args.infileList])

    crossingMask = approx_equals(endPoints[:,0], upstream, epsilon)
    crossingPoints = endPoints[crossingMask]

    dx = crossingPoints[:,0] - upstream
    # dz = crossingPoints[:,2] - cathode_z
    
    print (TPC_bounds)
    xmin = -320.
    xmax = 320.
    ymin = -600.
    ymax = 600.
    zmin = -302.768
    zmax = 302.768

    n_xbin = 65
    n_ybin = 121
    n_zbin = 61

    xgap = (xmax - xmin) / (n_xbin-1)
    ygap = (ymax - ymin) / (n_ybin-1)
    zgap = (zmax - zmin) / (n_zbin-1)

    x_grid = np.linspace(xmin, xmax, n_xbin)
    y_grid = np.linspace(ymin, ymax, n_ybin)
    z_grid = np.linspace(zmin, zmax, n_zbin)

    fig = plt.figure()
    totalD, bins, boop = np.histogram2d(crossingPoints[:,1] + my_geometry.tpc_offsets[0][1]*10,
                                        crossingPoints[:,2] + my_geometry.tpc_offsets[0][2]*10,
                                        bins = (y_grid, z_grid),
                                        weights = dz)
    counts, bins, bop = np.histogram2d(crossingPoints[:,1] + my_geometry.tpc_offsets[0][1]*10,
                                       crossingPoints[:,2] + my_geometry.tpc_offsets[0][2]*10,
                                       bins = (y_grid, z_grid))

    print (totalD/counts)

    plt.imshow((totalD/counts).T, origin='lower')
    plt.colorbar()
    
    plt.xlabel(r'y (vertical) [mm]')
    plt.ylabel(r'z (drift) [mm]')

    plt.show()
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='plotting track endpoints')
    parser.add_argument('infileList',
                        help = 'input numpy arrays containing track endpoint data',
                        nargs = '+')
    parser.add_argument('-g',
                        default = './pixel_layouts/multi_tile_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')


    args = parser.parse_args()

    main(args)
