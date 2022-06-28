import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage

from utils import *
from plotting import *

# Pretty fonts for figures 
# mpl.rc('text', usetex = True)
# mpl.rc('font', family='SignPainter')

# Ignore divide by zero warnings
np.seterr(divide='ignore', invalid='ignore')

def main(args):
    global my_geometry
    global TPC_bounds, anode_z, cathode_z, top, bottom, upstream, downstream
    global epsilon

    my_geometry = DetectorGeometry(args.d, args.g)
    TPC_bounds = get_TPC_bounds()
    anode_z = TPC_bounds[1][2][0]
    cathode_z = TPC_bounds[1][2][1]
    top = TPC_bounds[0][1][1]
    bottom = TPC_bounds[0][1][0]
    upstream = TPC_bounds[0][0][1]
    downstream = TPC_bounds[0][0][0]
    epsilon = 50 #mm

    endPoints = np.concatenate([np.load(infileName)
                                for infileName in args.infileList])

    if args.f == 'top':
        crossingMask = approx_equals(endPoints[:,1], top, epsilon)
        crossingPoints = endPoints[crossingMask]
        dy = top - crossingPoints[:,1]
    if args.f == 'bottom':
        crossingMask = approx_equals(endPoints[:,1], bottom, epsilon)
        crossingPoints = endPoints[crossingMask]
        dy = bottom - crossingPoints[:,1]
    if args.f == 'upstream':
        crossingMask = approx_equals(endPoints[:,0], upstream, epsilon)
        crossingPoints = endPoints[crossingMask]
        dx = upstream - crossingPoints[:,0] 
    if args.f == 'downstream':
        crossingMask = approx_equals(endPoints[:,0], downstream, epsilon)
        crossingPoints = endPoints[crossingMask]
        dx = downstream - crossingPoints[:,0] 

    print( 'Number of selected points:' + str( crossingPoints.shape ) )

    c = 0.443 #$coarseness (if 1, then every cm)

    n_xbin = int( round( (upstream + abs(downstream) )/(c*10) ) )
    n_ybin = int( round( ( top + abs(bottom) )/(c*10) ) )
    n_zbin = int( round( 2*anode_z/(c*10) ) )

    x_grid = np.linspace(downstream, upstream, n_xbin)
    y_grid = np.linspace(bottom, top, n_ybin)
    z_grid = np.linspace(-anode_z, anode_z, n_zbin)

    fig = plt.figure()

    if args.f == 'upstream' or args.f == 'downstream':
        totalD, bins, boop = np.histogram2d(crossingPoints[:,1],
                                            crossingPoints[:,2],
                                            bins = (y_grid, z_grid),
                                            weights = dx)
        counts, bins, bop = np.histogram2d(crossingPoints[:,1],
                                        crossingPoints[:,2],
                                        bins = (y_grid, z_grid))

    if args.f == 'top' or args.f == 'bottom':
        totalD, bins, boop = np.histogram2d(crossingPoints[:,2],
                                            crossingPoints[:,0],
                                            bins = (z_grid,x_grid),
                                            weights = dy)
        counts, bins, bop = np.histogram2d(crossingPoints[:,2],
                                        crossingPoints[:,0],
                                        bins = (z_grid,x_grid))


    to_m = 1000 #convert mm to m
    to_cm = 10 #convert mm to cm
    H = (totalD/counts)/to_cm
    im = plt.imshow(H, origin='lower', interpolation='none',cmap='winter',
        extent=[bop[0]/to_m, bop[-1]/to_m, bins[0]/to_m, bins[-1]/to_m]) #Gaussian interp if want smearing

    # plt.hist2d(crossingPoints[:,2],
    #         crossingPoints[:,1],
    #         bins = (z_grid, y_grid),cmap=plt.cm.jet)
    plt.tick_params(axis='both', which='both', labelsize = 15, direction = 'in')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=15)

    if args.f == 'upstream' or args.f == 'downstream':
        cbar.set_label(r'$\Delta x$ [cm]',fontsize = 15)
        plt.xlabel(r'z (drift) [cm]',fontsize = 15)
        plt.ylabel(r'y (vertical) [cm]',fontsize = 15)
        plt.xlim(-anode_z/to_m,anode_z/to_m)
        plt.ylim(bottom/to_m,top/to_m)

    if args.f == 'top' or args.f == 'bottom':
        cbar.set_label(r'$\Delta y$ [cm]',fontsize = 15)
        plt.xlabel(r'x (horizontal) [m]',fontsize = 15)
        plt.ylabel(r'z (drift) [m]',fontsize = 15)
        plt.xlim(downstream/to_m,upstream/to_m)
        plt.ylim(-anode_z/to_m,anode_z/to_m)

    # plt.show()
    plt.savefig( str(args.f) + '.png' )
        
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='plotting track endpoints')
    parser.add_argument('infileList',
                        help = 'input numpy arrays containing track endpoint data',
                        nargs = '+')
    parser.add_argument('-g',
                        default = './pixel_layouts/module1_layout-2.3.16.yaml',
                        # default = './pixel_layouts/multi_tile_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-f',
                        type = str,
                        help = 'face: top, bottom, upstream, downstream')


    args = parser.parse_args()

    main(args)


# python faceAna.py /global/project/projectdirs/dune/users/olexiy/M0/merged_array.npy -f upstream