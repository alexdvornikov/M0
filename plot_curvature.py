# Example usage
# python3 plot_curvature.py

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
    print(anode_z) #304.31 [mm]
    print(cathode_z) #1.5875 [mm] or 1/16''
    print(top) #402.524
    print(bottom) #839
    print(upstream) #310.38
    print(downstream) #-310.38
    # The cathode is 1/8'' (so 1/16'' into each TPC)?


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


    # #----------------------------------------------------------------------------#
    # # Load hexbin data
    # #----------------------------------------------------------------------------#
    # hb_counts = np.load('hb_counts_2anodes_merged.npy', allow_pickle=True)
    # # hb_counts = np.load('hb_counts_2anodes.npy', allow_pickle=True)
    # hb_pos = np.load('hb_pos_2anodes.npy', allow_pickle=True)
    # # hb_counts = np.load('hb_counts_merged.npy', allow_pickle=True)
    # # hb_pos = np.load('hb_pos.npy', allow_pickle=True)

    # counts_dx, counts_dy = hb_counts[0], hb_counts[1]
    # x_dx,y_dx = hb_pos[0], hb_pos[1]
    # x_dy,y_dy = hb_pos[2], hb_pos[3]

    # #----------------------------------------------------------------------------#
    # # Load 2d histogram data  
    # #----------------------------------------------------------------------------#
    # def load_hist(file):
    #     h2d_vals = np.load(file, allow_pickle=True)
    #     return h2d_vals[0].T, h2d_vals[1].T

    # dx_hist_zx, dy_hist_zx = load_hist('hist2d_zx_2anodes_merged.npy')
    # dx_hist_zy, dy_hist_zy = load_hist('hist2d_zy_2anodes_merged.npy')

    #----------------------------------------------------------------------------#
    def load_3dhist(file):
        h3d_vals = np.load(file, allow_pickle=True)
        # notNaN = np.count_nonzero( ~np.isnan( h3d_vals ) )
        # isNaN = np.count_nonzero(np.isnan( h3d_vals ))
        # print( 'Percent occupied: ', 100*(notNaN/isNaN) )
        h3d_vals[np.isnan(h3d_vals)] = 0 #Set NaNs to zeros
        return h3d_vals[0], h3d_vals[1], h3d_vals[2]
    
    dx_3dhist, dy_3dhist, dz_3dhist = load_3dhist('hist3d_merged.npy')

    # print( (dx_3dhist.T).shape )
    # print( ( (dx_3dhist.T)[0] ).T.shape )

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


    # #----------------------------------------------------------------------------#
    # # Set up the quiver plot. 
    # #----------------------------------------------------------------------------#
    # zbins = round( (2*anode_z)*cm)
    # xbins = round( (upstream + abs(downstream))*cm )
    # ybins = round( (top + abs(bottom))*cm )
    # xbins, ybins, zbins = round(xbins), round(ybins), round(zbins)

    # z_grid = np.linspace(-anode_z*cm,anode_z*cm,zbins)
    # x_grid = np.linspace(downstream*cm,upstream*cm,xbins)
    # y_grid = np.linspace(bottom*cm,top*cm, ybins)
    # xx,yy,zz = np.meshgrid(x_grid, y_grid,z_grid, indexing='ij')


    # data_length = len( xx.flatten() )
    # x=xx.flatten()
    # y=yy.flatten()
    # z=zz.flatten()
    # u=dx_3dhist.flatten()
    # v=dy_3dhist.flatten()
    # w=dz_3dhist.flatten()
    # w = 0*w #Set z to zero (we're at some z slices)

    # X, Y, Z = [], [], []
    # U, V, W = [], [], []

    # for i in range(data_length):
    #     if u[i] != 0:
    #         X.append(x[i])
    #         Y.append(y[i])
    #         Z.append(z[i])

    #         U.append(u[i])
    #         V.append(v[i])
    #         W.append(w[i])


    # fig = go.Figure( layout=solarized.layout3d_dark )
    # # Surface for cathode
    # fig.add_trace( go.Surface( x = x_grid, y = y_grid, z= 0*np.ones((ybins,xbins)), 
    #                                 colorscale='Gray', showscale=False,
    #                                 opacity = 0.5) )

    # # Use plotly's Cone plot for the vector field (aka quiver plot)
    # fig.add_trace( go.Cone(
    #     x=X,y=Y,z=Z, u=U,v=V,w=W,
    #     # cmax = 0.5, cmin = 0, 
    #     colorscale='viridis', 
    #     sizemode="absolute", #sizemode="scaled",
    #     sizeref = 6,
    #     opacity=1) )


    # fig.update_scenes(
    #     aspectmode='data', #cube
    #     xaxis_range=(downstream*cm,upstream*cm),
    #     yaxis_range=(bottom*cm,top*cm),
    #     zaxis_range=(-anode_z*cm,anode_z*cm),
    #     xaxis_title='x [cm]',
    #     yaxis_title='y [cm]',
    #     zaxis_title='z [cm]',
    # )

    # fig.write_html("quiver_test.html")
    # fig.show()






    #----------------------------------------------------------------------------#
    # Set up a 2d quiver plot. 
    #----------------------------------------------------------------------------#
    import plotly.figure_factory as ff
    z_slice = 30
    # U = ( (dx_3dhist.T)[z_slice] ).T #get slice from 3d vector field and flip it accordingly 
    # V = ( (dy_3dhist.T)[z_slice] ).T
    # zbins = round( (2*anode_z)*cm)
    xbins = round( (upstream + abs(downstream))*cm )
    ybins = round( (top + abs(bottom))*cm )
    xbins, ybins = round(xbins), round(ybins)
    # xbins, ybins, zbins = round(xbins), round(ybins), round(zbins)

    # z_grid = np.linspace(-anode_z*cm,anode_z*cm,zbins)
    x_grid = np.linspace(downstream*cm,upstream*cm,xbins)
    y_grid = np.linspace(bottom*cm,top*cm, ybins)
    # X,Y = np.meshgrid(x_grid, y_grid, indexing='ij')
    xx,yy = np.meshgrid(x_grid, y_grid, indexing='ij')
    # xx,yy,zz = np.meshgrid(x_grid, y_grid,z_grid, indexing='ij')

    # print(xx.shape)
    # print(yy.shape)


    data_length = len( xx.flatten() )
    x=xx.flatten()
    y=yy.flatten()
    # z=zz.flatten()
    print(  (dx_3dhist.T).shape )
    print( ( ( (dx_3dhist.T)[z_slice] ).T ).shape )
    u = ( ( (dx_3dhist.T)[z_slice] ).T ).flatten()
    v = ( ( (dy_3dhist.T)[z_slice] ).T ).flatten()
    # w=dz_3dhist.flatten()
    # w = 0*w #Set z to zero (we're at some z slices)

    X,Y = [],[]
    U,V = [],[]

    # X, Y, Z = [], [], []
    # U, V, W = [], [], []

    for i in range(data_length):
        if u[i] != 0:
            X.append(x[i])
            Y.append(y[i])
            # Z.append(z[i])

            U.append(u[i])
            V.append(v[i])
            # W.append(w[i])


    fig, ax = plt.subplots(figsize =(14, 8))
    ax.quiver(X, Y, U, V, units='xy', scale=0.1)
    # plt.figure(figsize=(10, 10))
    # plt.streamplot(X,Y,U,V,linewidth=None)
    
    # ax.xaxis.set_ticks([])
    # ax.yaxis.set_ticks([])
    # ax.axis([-0.3, 2.3, -0.3, 2.3])
    # ax.set_aspect('equal')
    
    # show plot
    plt.show()



    # fig = go.Figure( layout=solarized.layout3d_dark )
    # fig = ff.create_quiver(x=X, y=Y, u=U, v=V, 
    #                             scale = 5, arrow_scale = 0.5,line_width = 2)
    # Surface for cathode
    # fig.add_trace( go.Surface( x = x_grid, y = y_grid, z= 0*np.ones((ybins,xbins)), 
    #                                 colorscale='Gray', showscale=False,
    #                                 opacity = 0.5) )

    # # Use plotly's Cone plot for the vector field (aka quiver plot)
    # fig.add_trace( go.Cone(
    #     x=X,y=Y,z=Z, u=U,v=V,w=W,
    #     # cmax = 0.5, cmin = 0, 
    #     colorscale='viridis', 
    #     sizemode="absolute", #sizemode="scaled",
    #     sizeref = 6,
    #     opacity=1) )


    # fig.update_scenes(
    #     aspectmode='data', #cube
    #     xaxis_range=(downstream*cm,upstream*cm),
    #     yaxis_range=(bottom*cm,top*cm),
    #     # zaxis_range=(-anode_z*cm,anode_z*cm),
    #     xaxis_title='x [cm]',
    #     yaxis_title='y [cm]',
    #     # zaxis_title='z [cm]',
    # )

    # # fig.write_html("quiver_test_2d.html")
    # fig.show()

    # #----------------------------------------------------------------------------#
    # # Plot zx face projection
    # #----------------------------------------------------------------------------#
    # # global plt
    # fig, axes = plt.subplots(1,2, figsize=(4, 4), sharex=True, sharey=True)

    # # fs = 10 # Font size
    # vmax = 0.5 # Range on colorbar [cm]
    # extent = [-anode_z*cm,anode_z*cm,downstream*cm,upstream*cm]
    # # extent = [-anode_z*cm,cathode_z*cm,downstream*cm,upstream*cm]
    # kwargs = dict(
    #     vmin=-vmax, vmax=vmax,
    #     cmap='winter',
    #     origin='lower',
    #     extent=extent,
    #     interpolation='none',
    # )
        
    # ax = axes[0]
    # im = ax.imshow(dx_hist_zx, **kwargs)
    # ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    # # ax.tick_params(axis='both', which='minor', labelsize=fs)
    # h2d_config( ax, extent, fs )

    # ax = axes[1]
    # im = ax.imshow(dy_hist_zx, **kwargs)
    # # ax.set_ylabel('$\mathrm{x_{reco}}$ [cm]',fontsize=fs)
    # # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    # # ax.tick_params(axis='both', which='minor', labelsize=fs)
    # h2d_config( ax, extent, fs )

    # cax,kw = mpl.colorbar.make_axes(axes.flatten())
    # cbar = plt.colorbar(im, cax=cax)
    # cbar.set_label(label='[cm]', size=fs)
    # cbar.ax.tick_params(labelsize=fs)
    # # plt.show()
    # #----------------------------------------------------------------------------#


    # #----------------------------------------------------------------------------#
    # # Plot zy face projection
    # #----------------------------------------------------------------------------#
    # # global plt
    # fig, axes = plt.subplots(1,2, figsize=(3, 6), sharex=True, sharey=True)

    # # fs = 10 # Font size
    # vmax = 0.5 # Range on colorbar [cm]
    # extent = [-anode_z*cm,anode_z*cm, bottom*cm,top*cm]
    # # extent = [-anode_z*cm,cathode_z*cm, bottom*cm,top*cm]
    # kwargs = dict(
    #     vmin=-vmax, vmax=vmax,
    #     cmap='winter',
    #     origin='lower',
    #     extent = extent,
    #     interpolation='none',
    # )
        
    # ax = axes[0]
    # im = ax.imshow(dx_hist_zy, **kwargs)
    # ax.set_ylabel('$\mathrm{y_{reco}}$ [cm]',fontsize=fs)
    # # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # ax.set_title('$\mathrm{\Delta x}$ [cm]',fontsize=fs)
    # # ax.tick_params(axis='both', which='minor', labelsize=fs)
    # h2d_config( ax, extent, fs )

    # ax = axes[1]
    # im = ax.imshow(dy_hist_zy, **kwargs)
    # # ax.set_ylabel('$\mathrm{y_{reco}}$ [cm]',fontsize=fs)
    # # ax.set_xlabel('$\mathrm{z_{reco}}$ [cm]',fontsize=fs)
    # ax.set_title('$\mathrm{\Delta y}$ [cm]',fontsize=fs)
    # # ax.tick_params(axis='both', which='minor', labelsize=fs)
    # h2d_config( ax, extent, fs )


    # cax,kw = mpl.colorbar.make_axes(axes.flatten())
    # cbar = plt.colorbar(im, cax=cax)
    # cbar.set_label(label='[cm]', size=fs)
    # cbar.ax.tick_params(labelsize=fs)

    plt.show()
    # # Can save individual figures by...
    # # plt.savefig('blah.eps')
    # #----------------------------------------------------------------------------#

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
