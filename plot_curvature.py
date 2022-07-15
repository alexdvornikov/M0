# Example usage
# python3 plot_curvature.py


# Put hb_counts.npy, hb_pos.npy
# and hist2d_zx.npy, hist2d_zy.npy
# in working directory or add appropriate paths to np.load() in this script. 

from utils_m1 import *
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')
    
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
    cm = 0.1 #Conversion from mm to cm


    #----------------------------------------------------------------------------#
    # Visuals for binned_statistic_2d plot
    #----------------------------------------------------------------------------#
    def h2d_config( ax, extent, fs ):

        major_ticks_x = np.arange(-40, 40, 10)
        # major_ticks_x = np.arange(-40, 10, 10)
        minor_ticks_x = np.arange(-40, 40, 5)
        major_ticks_y = np.arange(-90, 90, 10)
        minor_ticks_y = np.arange(-90, 90, 5)

        ax.set_xticklabels( major_ticks_x )
        ax.set_yticklabels( major_ticks_y ) 

        ax.set_xticks(major_ticks_x)
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.set_yticks(major_ticks_y)
        ax.set_yticks(minor_ticks_y, minor=True)
        ax.tick_params(axis='both', which='minor', labelsize=fs)

        ax.grid(which='both', color = 'white', linestyle = '--', linewidth = 1)
        ax.grid(which='minor', alpha=1)
        ax.grid(which='major', alpha=1)

        ax.set_xlim(extent[0],extent[1])
        ax.set_ylim(extent[2], extent[3])

        ax.axis([extent[0], extent[1], extent[2],extent[3]])

        ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    #----------------------------------------------------------------------------#



    #----------------------------------------------------------------------------#
    # Visuals for hexbin plot
    #----------------------------------------------------------------------------#
    def hexbin_config( ax, hb, fs ):

        offset_range = 2
        # offset_range = 20

        major_ticks_x = np.arange(-35, 35, 5)
        minor_ticks_x = np.arange(-35, 35, 1)
        major_ticks_y = np.arange(-25, 25, 5)
        minor_ticks_y = np.arange(-21, 21, 1)
        

        # ax.set_xticklabels( major_ticks_x )
        # ax.set_yticklabels( major_ticks_y ) 

        ax.hlines(y=minor_ticks_y , xmin = -anode_z, xmax = anode_z, colors='white', linestyles='--', lw=1)
        ax.vlines(x=major_ticks_x , ymin = -offset_range, ymax = offset_range, colors='white', linestyles='--', lw=1)

        ax.set_xticks(major_ticks_x)
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.set_yticks(major_ticks_y)
        ax.set_yticks(minor_ticks_y, minor=True)
        ax.tick_params(axis='both', which='minor', labelsize=fs)

        ax.grid(which='both', color = 'white', linestyle = '--', linewidth = 1)
        ax.grid(which='minor', alpha=1)
        ax.grid(which='major', alpha=1)

        # ax.set_xlim(-anode_z*cm,anode_z*cm)
        # ax.set_xlim(-anode_z*cm,cathode_z*cm)
        # ax.set_ylim(-offset_range, offset_range)

        # ax.axis([-anode_z*cm, anode_z*cm, -1, 1])
        ax.axis([-anode_z*cm, anode_z*cm, -offset_range, offset_range])
        # ax.axis([-anode_z*cm, cathode_z*cm, -offset_range, offset_range])
        cb = fig.colorbar(hb, ax=ax)
        cb.ax.tick_params(labelsize=fs)

        # Two lines to guide the eye at plus/minus 5 [cm]
        ax.hlines(y=[-5,5], xmin = -anode_z, xmax = anode_z, colors='r', linestyles='-', lw=2)
    #----------------------------------------------------------------------------#


    #----------------------------------------------------------------------------#
    # Load hexbin data
    #----------------------------------------------------------------------------#
    hb_counts = np.load('hb_counts_2anodes.npy', allow_pickle=True)
    hb_pos = np.load('hb_pos_2anodes.npy', allow_pickle=True)
    # hb_counts = np.load('hb_counts_merged.npy', allow_pickle=True)
    # hb_pos = np.load('hb_pos.npy', allow_pickle=True)

    counts_dx, counts_dy = hb_counts[0], hb_counts[1]
    x_dx,y_dx = hb_pos[0], hb_pos[1]
    x_dy,y_dy = hb_pos[2], hb_pos[3]

    # # Sanity check for array shapes
    # # print(counts_dx.shape)
    # # print(x_dx.shape)
    # # print(counts_dy.shape)
    # # print(x_dy.shape)

    #----------------------------------------------------------------------------#
    # Load 2d histogram data  
    #----------------------------------------------------------------------------#
    def load_hist(file):
        h2d_vals = np.load(file, allow_pickle=True)
        return h2d_vals[0].T, h2d_vals[1].T

    dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_2anodes.npy')
    dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_2anodes.npy')

    # dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_merged.npy')
    # dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_merged.npy')

    #----------------------------------------------------------------------------#
    # Plot hexbin deviation vs drift coordinate (using [cm] units)
    #----------------------------------------------------------------------------#
    global plt # Only need to do this once per file
    fs = 15 # Font size
    fig, axes = plt.subplots(2,1,figsize=(16, 4),sharex=True,sharey=True)
    ax = axes[0]
    hb = ax.hexbin(x_dx, y_dx, C=counts_dx, gridsize = 50, cmap='viridis') #bins ='log'
    ax.set_ylabel('$\mathrm{\Delta x}$ [cm]', fontsize = fs)
    hexbin_config(ax,hb,fs)


    ax = axes[1]
    hb = ax.hexbin(x_dy, y_dy, C=counts_dy, gridsize = 50, cmap='viridis') #bins ='log'
    ax.set_ylabel('$\mathrm{\Delta y}$ [cm]', fontsize = fs)
    ax.set_xlabel('z [cm]', fontsize = fs)
    hexbin_config(ax,hb,fs)

    # plt.show()

    # In want to try hist2d instead of hexbin...
    # ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis')
    # counts, xedges, yedges, im = ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis', norm=LogNorm())
    # fig.colorbar(im, ax=ax)
    #----------------------------------------------------------------------------#

    #----------------------------------------------------------------------------#
    # Plot zx face projection
    #----------------------------------------------------------------------------#
    # global plt
    fig, axes = plt.subplots(1,2, figsize=(4, 4), sharex=True, sharey=True)

    fs = 15 # Font size
    vmax = 2 # Range on colorbar [cm]
    extent = [-anode_z*cm,anode_z*cm,downstream*cm,upstream*cm]
    # extent = [-anode_z*cm,cathode_z*cm,downstream*cm,upstream*cm]
    kwargs = dict(
        vmin=-vmax, vmax=vmax,
        cmap='winter',
        origin='lower',
        extent=extent,
        interpolation='none',
    )
        
    ax = axes[0]
    im = ax.imshow(dx_hist_zx, **kwargs)
    ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    # ax.tick_params(axis='both', which='minor', labelsize=fs)
    h2d_config( ax, extent, fs )

    ax = axes[1]
    im = ax.imshow(dy_hist_zx, **kwargs)
    # ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    # ax.tick_params(axis='both', which='minor', labelsize=fs)
    h2d_config( ax, extent, fs )

    cax,kw = mpl.colorbar.make_axes(axes.flatten())
    cbar = plt.colorbar(im, cax=cax, label='[cm]')
    cbar.ax.tick_params(labelsize=fs)
    # plt.show()
    #----------------------------------------------------------------------------#


    #----------------------------------------------------------------------------#
    # Plot zy face projection
    #----------------------------------------------------------------------------#
    # global plt
    fig, axes = plt.subplots(1,2, figsize=(3, 6), sharex=True, sharey=True)

    fs = 15 # Font size
    vmax = 2 # Range on colorbar [cm]
    extent = [-anode_z*cm,anode_z*cm, bottom*cm,top*cm]
    # extent = [-anode_z*cm,cathode_z*cm, bottom*cm,top*cm]
    kwargs = dict(
        vmin=-vmax, vmax=vmax,
        cmap='winter',
        origin='lower',
        extent = extent,
        interpolation='none',
    )
        
    ax = axes[0]
    im = ax.imshow(dx_hist_zy, **kwargs)
    ax.set_ylabel('$\mathrm{y_{reco}}$ [cm]',fontsize=fs)
    # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    # ax.tick_params(axis='both', which='minor', labelsize=fs)
    h2d_config( ax, extent, fs )

    ax = axes[1]
    im = ax.imshow(dy_hist_zy, **kwargs)
    # ax.set_ylabel('$\mathrm{y_{reco}}$ [cm]',fontsize=fs)
    # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    # ax.tick_params(axis='both', which='minor', labelsize=fs)
    h2d_config( ax, extent, fs )


    cax,kw = mpl.colorbar.make_axes(axes.flatten())
    cbar = plt.colorbar(im, cax=cax, label='[cm]')
    cbar.ax.tick_params(labelsize=fs)

    plt.show()
    # Can save individual figures by...
    # plt.savefig('blah.eps')
    #----------------------------------------------------------------------------#

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
