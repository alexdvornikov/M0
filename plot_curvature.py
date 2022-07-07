# Example usage
# python3 plot_curvature.py distortions.npy
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
    cm = 0.1 #Conversion from mm to cm

    #----------------------------------------------------------------------------#
    # Visuals for hexbin plot
    #----------------------------------------------------------------------------#
    def hexbin_config( ax, hb, fs ):

        offset_range = 20

        major_ticks_x = np.arange(-35, 5, 5)
        minor_ticks_x = np.arange(-35, 5, 1)
        major_ticks_y = np.arange(-25, 25, 5)
        minor_ticks_y = np.arange(-21, 21, 1)

        ax.set_xticks(major_ticks_x)
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.set_yticks(major_ticks_y)
        ax.set_yticks(minor_ticks_y, minor=True)
        ax.tick_params(axis='both', which='minor', labelsize=fs)

        ax.grid(which='both', color = 'grey', linestyle = '--', linewidth = 0.5)
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        ax.set_xlim(downstream*cm,0)
        ax.set_ylim(-offset_range, offset_range)

        ax.axis([-anode_z*cm, cathode_z*cm, -offset_range, offset_range])
        cb = fig.colorbar(hb, ax=ax)
        cb.ax.tick_params(labelsize=fs)

        # Two lines to guide the eye at plus/minus 5 [cm]
        ax.hlines(y=[-5,5], xmin = -anode_z, xmax = cathode_z, colors='r', linestyles='--', lw=1)
    #----------------------------------------------------------------------------#


    #----------------------------------------------------------------------------#
    # Load hexbin data
    #----------------------------------------------------------------------------#
    hb_counts = np.load('hb_counts.npy', allow_pickle=True)
    hb_pos = np.load('hb_pos.npy', allow_pickle=True)

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

    dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx.npy')
    dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy.npy')

    #----------------------------------------------------------------------------#
    # Plot hexbin deviation vs drift coordinate (using [cm] units)
    #----------------------------------------------------------------------------#
    global plt # Only need to do this once per file
    fs = 15 # Font size
    fig, axes = plt.subplots(2,1,figsize=(10, 8),sharex=True,sharey=True)
    ax = axes[0]
    hb = ax.hexbin(x_dx, y_dx, C=counts_dx, gridsize = 50, cmap='viridis', bins='log')
    ax.set_ylabel('$\mathrm{\Delta x}$ [cm]', fontsize = fs)
    hexbin_config(ax,hb,fs)


    ax = axes[1]
    hb = ax.hexbin(x_dy, y_dy, C=counts_dy, gridsize = 50, cmap='viridis', bins='log')
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
    fig, axes = plt.subplots(1,2, figsize=(8, 8), sharex=True, sharey=True)

    fs = 15 # Font size
    vmax = 10 # Range on colorbar [cm]
    kwargs = dict(
        vmin=-vmax, vmax=vmax,
        cmap='winter',
        origin='lower',
        extent=[-anode_z*cm,cathode_z*cm,downstream*cm,upstream*cm],
        interpolation='none',
    )
        
    ax = axes[0]
    im = ax.imshow(dx_hist_zx, **kwargs)
    ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)

    ax = axes[1]
    im = ax.imshow(dy_hist_zx, **kwargs)
    ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)
    cax,kw = mpl.colorbar.make_axes(axes.flatten())
    cbar = plt.colorbar(im, cax=cax, label='[cm]')
    cbar.ax.tick_params(labelsize=fs)
    # plt.show()
    #----------------------------------------------------------------------------#


    #----------------------------------------------------------------------------#
    # Plot zy face projection
    #----------------------------------------------------------------------------#
    # global plt
    fig, axes = plt.subplots(1,2, figsize=(6, 8), sharex=True, sharey=True)

    fs = 15 # Font size
    vmax = 10 # Range on colorbar [cm]
    kwargs = dict(
        vmin=-vmax, vmax=vmax,
        cmap='winter',
        origin='lower',
        extent = [-anode_z*cm,cathode_z*cm, bottom*cm,top*cm],
        interpolation='none',
    )
        
    ax = axes[0]
    im = ax.imshow(dx_hist_zy, **kwargs)
    ax.set_ylabel('$\mathrm{y_{reco}}$ [cm]',fontsize=fs)
    ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)

    ax = axes[1]
    im = ax.imshow(dy_hist_zy, **kwargs)
    ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)
    # fig.tight_layout()
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
