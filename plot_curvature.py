# Example usage
# python3 plot_curvature.py distortions.npy


# If interested in timiming the code...
# from datetime import datetime
# startTime = datetime.now()
# Also need to uncomment a datetime line above __main__ at the bottom of the code. 

from utils_m1 import *
    
def main(args):
    global my_geometry
    global TPC_bounds, anode_z, cathode_z, top, bottom, upstream, downstream
    global length_cut
    global epsilon

    my_geometry = DetectorGeometry(args.detector, args.geometry)
    TPC_bounds = get_TPC_bounds()
    anode_z = TPC_bounds[1][2][0]
    cathode_z = TPC_bounds[1][2][1]
    top = TPC_bounds[0][1][1]
    bottom = TPC_bounds[0][1][0]
    upstream = TPC_bounds[0][0][1]
    downstream = TPC_bounds[0][0][0]


    hb_counts = np.load('hb_counts.npy', allow_pickle=True)
    hb_pos = np.load('hb_pos.npy', allow_pickle=True)

    counts_dx, counts_dy = hb_counts[0], hb_counts[1]
    x_dx,y_dx = hb_pos[0], hb_pos[1]
    x_dy,y_dy = hb_pos[2], hb_pos[3]

    #----------------------------------------------------------------------------#
    # Plot deviation vs drift coordinate
    #----------------------------------------------------------------------------#
    global plt
    fig, axes = plt.subplots(2,1,figsize=(8, 8),sharex=True,sharey=True)

    major_ticks_x = np.arange(-35, 5, 5)
    minor_ticks_x = np.arange(-35, 5, 1)
    major_ticks_y = np.arange(-21, 21, 5)
    minor_ticks_y = np.arange(-21, 21, 1)

    offset_range = 20

    ax = axes[0]
    # ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis')
    # counts, xedges, yedges, im = ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis', norm=LogNorm())
    # fig.colorbar(im, ax=ax)

    hb = ax.hexbin(x_dx, y_dx, C=counts_dx, gridsize = 50, cmap='viridis', bins='log')
    ax.axis([downstream/10, 0, -offset_range, offset_range])
    cb = fig.colorbar(hb, ax=ax)

    ax.set_ylabel('$\mathrm{\Delta x}$ [cm]')

    ax.set_xticks(major_ticks_x)
    ax.set_xticks(minor_ticks_x, minor=True)
    ax.set_yticks(major_ticks_y)
    ax.set_yticks(minor_ticks_y, minor=True)

    ax.grid(which='both', color = 'grey', linestyle = '--', linewidth = 0.5)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)




    ax = axes[1]
    hb = ax.hexbin(x_dy, y_dy, C=counts_dy, gridsize = 50, cmap='viridis', bins='log')
    # hb = ax.hexbin(z/10, dy/10, extent = (downstream/10, 0, -offset_range, offset_range), mincnt=1, gridsize = 50, cmap='viridis', bins='log')
    ax.axis([downstream/10, 0, -offset_range, offset_range])
    cb = fig.colorbar(hb, ax=ax)
    ax.set_ylabel('$\mathrm{\Delta y}$ [cm]')
    ax.set_xlabel('z [cm]')

    ax.set_xticks(major_ticks_x)
    ax.set_xticks(minor_ticks_x, minor=True)
    ax.set_yticks(major_ticks_y)
    ax.set_yticks(minor_ticks_y, minor=True)

    ax.grid(which='both', color = 'grey', linestyle = '--', linewidth = 0.5)
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)

    ax.set_xlim(downstream/10,0)
    ax.set_ylim(-offset_range, offset_range)

    plt.show()
    # plt.savefig('res.eps')
    #----------------------------------------------------------------------------#

    #----------------------------------------------------------------------------#
    # Plot face projection
    #----------------------------------------------------------------------------#
    # global plt
    # fig, axes = plt.subplots(1,2, figsize=(6, 6), sharex=True, sharey=True)

    # vmax = 3
    # bins = [5,10]
    # kwargs = dict(
    #     vmin=-vmax, vmax=vmax,
    #     cmap='winter',
    #     origin='lower',
    #     extent=[-32,0,-30,30],
    #     interpolation='none',
    # )
        

    # ax = axes[0]
    # pf = binned_statistic_2d(z/10, x/10, dx/10., bins=bins, range=[(-30,0), (-30,30)])
    # ax.imshow(pf.statistic.T, **kwargs)
    # ax.set_ylabel('$\mathrm{x_{r}}$ [cm]')
    # ax.set_xlabel('$\mathrm{z_{r}}$ [cm]')
    # ax.set_title('$\mathrm{\Delta x}$ [cm]')

    # ax = axes[1]
    # pf = binned_statistic_2d(z/10, x/10, dy/10., bins=bins, range=[(-30,0), (-30,30)])
    # im = ax.imshow(pf.statistic.T, **kwargs)
    # ax.set_title('$\mathrm{\Delta y}$ [cm]')

    # fig.tight_layout()

    # cax,kw = mpl.colorbar.make_axes(axes.flatten())
    # plt.colorbar(im, cax=cax, label='[cm]')
    # plt.show()
    #----------------------------------------------------------------------------#



    if args.plot:

        # fig.update_scenes(
        #     aspectmode='data', #cube
        #     xaxis_range=(downstream,upstream),
        #     yaxis_range=(bottom,top),
        #     zaxis_range=(-anode_z,anode_z),
        #     xaxis_title='x [mm]',
        #     yaxis_title='y [mm]',
        #     zaxis_title='z [mm]',
        # )

        fig.show()
        # plotly image will be rendered in web browser 
        # Can screen grab or save the image with the kaleido package (pip install)
        # fig.write_image("figure.pdf", engine="kaleido")

    # f.close() # Don't need this if using a numpy array

    # Show time elapsed for running the code
    # print(datetime.now() - startTime)

if __name__ == '__main__': 
    import argparse

    parser = argparse.ArgumentParser(description='Plot the first N tracks from a given file')
    # parser.add_argument('infile',
    #                     help = 'input larpix data with track reconstruction')
    parser.add_argument('-g', '--geometry',
                        default = './pixel_layouts/module1_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d', '--detector',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    parser.add_argument('-o', '--output',
                        default = '',
                        type = str,
                        help = 'save the data which passes the selection to a file')
    parser.add_argument('-p', '--plot',
                        default = True,
                        type = bool,
                        help = 'show plots')

    args = parser.parse_args()

    main(args)
