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
    global fs 
    fs = 10 # Font size for plots

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

        ax.vlines(x=0 , ymin = -90, ymax = 90, colors='grey', linestyles='-', lw=2)

        ax.set_xticklabels( major_ticks_x , fontsize = fs)
        ax.set_yticklabels( major_ticks_y, fontsize=fs ) 

        ax.set_xticks(major_ticks_x)
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.set_yticks(major_ticks_y)
        ax.set_yticks(minor_ticks_y, minor=True)
        ax.tick_params(axis='both', which='minor', labelsize=fs)

        ax.grid(which='both', color = 'white', linestyle = '--', linewidth = 0)
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

        offset_range = 5
        # offset_range = 20

        major_ticks_x = np.arange(-35, 35, 5)
        minor_ticks_x = np.arange(-35, 35, 1)
        major_ticks_y = np.arange(-5, 5, 1)
        minor_ticks_y = np.arange(-5, 5, 1)

        ax.set_xticklabels( major_ticks_x, fontsize = fs )
        ax.set_yticklabels( major_ticks_y, fontsize = fs ) 
        
        # ax.set_xticklabels( major_ticks_x )
        # ax.set_yticklabels( major_ticks_y ) 

        # ax.hlines(y=minor_ticks_y , xmin = -anode_z, xmax = anode_z, colors='white', linestyles='--', lw=1)
        # ax.vlines(x=major_ticks_x , ymin = -offset_range, ymax = offset_range, colors='white', linestyles='--', lw=1)
        ax.vlines(x=0 , ymin = -offset_range, ymax = offset_range, colors='grey', linestyles='-', lw=2)

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
    hb_counts = np.load('hb_counts_2anodes_merged.npy', allow_pickle=True)
    # hb_counts = np.load('hb_counts_2anodes.npy', allow_pickle=True)
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
        # print(h2d_vals.shape)
        return h2d_vals[0].T, h2d_vals[1].T

    # dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_2anodes_merged.npy')
    # dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_2anodes_merged.npy')

    dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_2anodes.npy')
    dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_2anodes.npy')

    # dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_merged.npy')
    # dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_merged.npy')

    #----------------------------------------------------------------------------#
    def load_3dhist(file):
        h3d_vals = np.load(file, allow_pickle=True)
        print(h3d_vals.shape)
        # print(h3d_vals)

        print('Non-empty bins: ')
        print( np.count_nonzero( ~np.isnan( h3d_vals ) )/3 ) #Number of non NaNs 
        print('Empty bins: ')
        print( np.count_nonzero(np.isnan( h3d_vals ))/3 ) #Number of NaNs 

        h3d_vals[np.isnan(h3d_vals)] = 0 #Set NaNs to zeros
        # h3d_vals = h3d_vals[ ~np.isnan( h3d_vals ) ]
        # print(h3d_vals[0].shape)
        # print(h3d_vals)
        # print(h3d_vals[0])
        # print(h3d_vals[30])
        return h3d_vals[0], h3d_vals[1], h3d_vals[2]
        # return h3d_vals[0].T, h3d_vals[1].T, h3d_vals[2].T
    
    dx_3dhist, dy_3dhist, dz_3dhist = load_3dhist('hist3d.npy')


    # np.count_nonzero(~np.isnan( h3d_vals ))

    # #----------------------------------------------------------------------------#
    # # Plot hexbin deviation vs drift coordinate (using [cm] units)
    # #----------------------------------------------------------------------------#
    # global plt # Only need to do this once per file
    # # fs = 10 # Font size
    # grdsize = 2*round(2*anode_z*cm)
    # fig, axes = plt.subplots(2,1,figsize=(20, 4),sharex=True,sharey=True)
    # ax = axes[0]
    # hb = ax.hexbin(x_dx, y_dx, C=counts_dx, gridsize = grdsize, cmap='viridis', bins ='log') #bins ='log'
    # ax.set_ylabel('$\mathrm{\Delta x}$ [cm]', fontsize = fs)
    # hexbin_config(ax,hb,fs)


    # ax = axes[1]
    # hb = ax.hexbin(x_dy, y_dy, C=counts_dy, gridsize = grdsize, cmap='viridis', bins ='log') #bins ='log'
    # ax.set_ylabel('$\mathrm{\Delta y}$ [cm]', fontsize = fs)
    # ax.set_xlabel('z [cm]', fontsize = fs)
    # hexbin_config(ax,hb,fs)

    # # plt.show()

    # # In want to try hist2d instead of hexbin...
    # # ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis')
    # # counts, xedges, yedges, im = ax.hist2d(z/10, dx/10, bins = [20,20], range=[(-32,0), (-10,10)], cmap='viridis', norm=LogNorm())
    # # fig.colorbar(im, ax=ax)
    # #----------------------------------------------------------------------------#

    # print(dx_hist_zx.shape)
    # print(dy_hist_zx.shape)
    # print(dx_hist_zy.shape)
    # print(dy_hist_zy.shape)
    # # (62,61)
    # # #----------------------------------------------------------------------------#
    # # # Vector plot  zx (this is not functional, just messing around)
    # # #----------------------------------------------------------------------------#
    # zbins = round( (2*anode_z)*cm)
    # xbins = round( (upstream + abs(downstream))*cm )
    # ybins = round( (top + abs(bottom))*cm )
    # # # bns = [zbins,xbins]
    # # # extent = [(-anode_z*cm,anode_z*cm), (downstream*cm,upstream*cm)]
    # z_grid = np.linspace(-anode_z*cm,anode_z*cm,zbins)
    # x_grid = np.linspace(downstream*cm,upstream*cm,xbins)
    # y_grid = np.linspace(bottom*cm,top*cm, ybins)
    # print(z_grid.shape) #61
    # print(x_grid.shape) #62
    # print(y_grid.shape) #124
    # xx,yy,zz = np.meshgrid(x_grid, y_grid,z_grid)
    # print(xx.shape)
    # print(yy.shape)
    # print(zz.shape)

    # zi = len(z_grid)
    # # print(zi)
    # new_array_dx = []
    # for i in range(zi):
    #     b1, b2 = np.meshgrid( dx_hist_zx[:,i],dx_hist_zy[:,i] )
    #     print(b1.shape)
    #     print(b2.shape)
    #     new_array_dx.append( np.meshgrid( dx_hist_zx[:,i],dx_hist_zy[:,i] ) )
    # new_array_dx = np.array(new_array_dx)
    # print(new_array_dx.shape)
    # # print( dx_hist_zx[:,1].shape )
    # # dz = np.full(dx_hist_zx.shape, 0) # Just messing around for now 

    # # # fig, ax = plt.subplots(figsize=(4, 4))
    # # ax = plt.figure().add_subplot(projection='3d')
    # # ax.quiver(dx_hist_zx, dy_hist_zx, np.full(dx_hist_zx.shape, 0))
    # # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # # ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # # ax.xaxis.set_ticks([])
    # # ax.yaxis.set_ticks([])
    # # ax.set_aspect('equal')
    # # plt.show()
    # # # ax.quiver(X, Y, u, v)

    zbins = round( (2*anode_z)*cm)
    xbins = round( (upstream + abs(downstream))*cm )
    ybins = round( (top + abs(bottom))*cm )

    xbins, ybins, zbins = round(xbins), round(ybins), round(zbins)

    z_grid = np.linspace(-anode_z*cm,anode_z*cm,zbins)
    x_grid = np.linspace(downstream*cm,upstream*cm,xbins)
    y_grid = np.linspace(bottom*cm,top*cm, ybins)
    # print(z_grid.shape) #61
    # print(x_grid.shape) #62
    # print(y_grid.shape) #124
    xx,yy,zz = np.meshgrid(x_grid, y_grid,z_grid, indexing='ij')



    import matplotlib.cm as cm
    from matplotlib.colors import Normalize

    # fig, ax = plt.subplots(figsize=(4, 4))
    # ax = plt.figure().add_subplot(projection='3d')
    u, v, w = dx_3dhist, dy_3dhist, dz_3dhist
    # w = 0*w #Set z to zero

    # colors = np.sqrt( u*u + v*v + w*w)
    # norm = Normalize()
    # norm.autoscale(colors)

    # colormap = cm.inferno
    # print(xx.shape)
    # print(yy.shape)
    # print(zz.shape)
    # print(u.shape)
    # print(v.shape)
    # print(w.shape)
    # print(u)
    # print(v)
    # print(w)
    # blah = w.flatten()
    # print( np.count_nonzero(np.isnan(blah)) )
    # print( blah[0] )
    # print( np.mean( u.flatten() ) )
    # print( np.mean( v.flatten() ) )
    # print( np.mean( w.flatten() ) )

    # print(u.shape) 
    # print(v.shape)
    # print(w.shape)

    # print(xx.shape)
    # print(yy.shape)
    # print(zz.shape)
    # color=colormap(norm(colors))
    # q = ax.quiver(xx, yy, zz, u/colors, v/colors, w/colors)

    data_length = len( xx.flatten() )
    x=xx.flatten()
    y=yy.flatten()
    z=zz.flatten()
    u=dx_3dhist.flatten()
    v=dy_3dhist.flatten()
    w=dz_3dhist.flatten()

    # mags = np.sqrt( u*u + v*v + w*w)
    # u = u/mags
    # u = u/mags
    # v = v/mags

    X = []
    Y = []
    Z = []
    U = []
    V = []
    W = []

    for i in range(data_length):
        if u[i] != 0:
            X.append(x[i])
            Y.append(y[i])
            Z.append(z[i])

            U.append(u[i])
            V.append(v[i])
            W.append(w[i])

    print(len(U))
    print( max(U) )
    print( min(U) )

    fig = go.Figure(data = go.Cone(
        x=X,y=Y,z=Z, u=U,v=V,w=W,
        colorscale='Blues',
        sizemode="scaled", #sizemode="absolute", sizeref=40
        sizeref=1))

    # fig = go.Figure(data = go.Cone(
    #     x=xx.flatten(),
    #     y=yy.flatten(),
    #     z=zz.flatten(),
    #     u=dx_3dhist.flatten(),
    #     v=dy_3dhist.flatten(),
    #     w=dz_3dhist.flatten(),
    #     colorscale='Blues',
    #     sizemode="absolute",
    #     sizeref=40),
    #     render_mode='webgl' )

    fig.show()


    # ax.quiver(dx_hist_zx, dy_hist_zx, np.full(dx_hist_zx.shape, 0))
    # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # ax.xaxis.set_ticks([])
    # ax.yaxis.set_ticks([])
    # ax.set_aspect('equal')
    # plt.show()
    # ax.quiver(X, Y, u, v)
    #----------------------------------------------------------------------------#
    # Plot zx face projection
    #----------------------------------------------------------------------------#
    # global plt
    fig, axes = plt.subplots(1,2, figsize=(4, 4), sharex=True, sharey=True)

    # fs = 10 # Font size
    vmax = 0.5 # Range on colorbar [cm]
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
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(label='[cm]', size=fs)
    cbar.ax.tick_params(labelsize=fs)
    # plt.show()
    #----------------------------------------------------------------------------#


    #----------------------------------------------------------------------------#
    # Plot zy face projection
    #----------------------------------------------------------------------------#
    # global plt
    fig, axes = plt.subplots(1,2, figsize=(3, 6), sharex=True, sharey=True)

    # fs = 10 # Font size
    vmax = 0.5 # Range on colorbar [cm]
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
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label(label='[cm]', size=fs)
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
