# Example usage
# python3 curvatureAna.py events_2022_02_08_07_36_25_CET.gz.h5


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
    length_cut = 100 #mm
    epsilon = 10 #mm


    global f
    f = h5py.File(args.infile, 'r')
    data = Data(f)
    #----------------------------------------------------------------------------#

    # # If plotting residuals
    # fig, axes = plt.subplots(2,1,figsize=(8, 8),sharex=True,sharey=True)
    
    # If plotting hits and PCAs
    if args.plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits import mplot3d as plt3d
        from plotting import draw_boundaries

        fig = go.Figure(layout=plot_theme.layout3d_white)
        colors = cycle(plot_theme.DEFAULT_SEQUENCE)
        
    trackMask = (data.rawTracks['length'] > length_cut) & (data.rawTracks['nhit'] > 0)
    track_idx = np.where(trackMask)[0]

    if args.n >= 0:
        N = args.n
    else:
        N = len(track_idx)

    max_trks_show = args.nshow
    n_trks = 0

    output = []

    print(N)

    for thisTrack_idx in track_idx[:N]:
        thisTrack = data.rawTracks[thisTrack_idx]
        thisEvent = dereference(thisTrack_idx, data.track_ref, data.rawEvents, region=data.track_reg, ref_direction=(1,0))
        evt_idx = thisEvent['id']
        t0_dtype = dereference(evt_idx, data.t0_ref, data.rawT0, region=data.t0_reg, ref_direction=(1,0))
        t0_unflattened = t0_dtype['ts'] #Light trigger
        t0 = t0_unflattened[0][0]

        if thisEvent['n_ext_trigs'] < 2:
            continue
        #----------------------------------------------------------------------------#
        # Get hits
        #----------------------------------------------------------------------------#
        theseHits = dereference(thisTrack_idx, data.hit_ref, data.rawHits, region=data.hit_reg, ref_direction=(0,1))
        pos3d = hit_to_3d(my_geometry, theseHits[0], t0)
        pos3d = np.array([np.array(entry) for entry in pos3d]) #Nested list to nested array
        #----------------------------------------------------------------------------#
        # Get endpoints
        #----------------------------------------------------------------------------#
        corrected_endpoints = []
        startHitPos, endHitPos = get_extreme_hit_pos(t0,thisTrack,theseHits[0],my_geometry)
        corrected_endpoints.append(startHitPos)
        corrected_endpoints.append(endHitPos)

        # # Focus on a single TPC for now (avoid positive drift coordinates)
        # if startHitPos[2] > 0 or endHitPos[2] > 0:
        #     continue
        #----------------------------------------------------------------------------#
        # Check anode crossing and plot
        # #----------------------------------------------------------------------------#
        # print(startHitPos[2])
        # print(endHitPos[2])
        if is_both_anodes_piercer(startHitPos, endHitPos, anode_z, epsilon):
        # if is_anode_piercer(startHitPos, endHitPos, anode_z, epsilon):

            # For more on get_pca_endpts see utils_m1.py
            ds = distortions(t0, my_geometry, theseHits[0], pos3d, near_anode = True, nhit = 10)
            output.append(ds)

            if args.plot:
                endpts = ds['endpts']

                # Plot projected PCA endpoints
                fig.add_trace(
                    go.Scatter3d(
                        x=endpts[:,0],
                        y=endpts[:,1],
                        z=endpts[:,2],
                        mode='lines',
                        line_width=5,
                        line_color=plot_theme.TURQ,
                    )
                )

                # Plot hits
                fig.add_trace(
                    go.Scatter3d(
                        x=pos3d[0], y=pos3d[1], z=pos3d[2],
                        mode='markers',
                        marker_size=3,
                        marker_color=next(colors),
                    )
                )


                n_trks += 1
                if n_trks >= max_trks_show:
                    break 


    # if args.output:
    #     np.save(args.output, corrected_endpoints)

    #----------------------------------------------------------------------------#
    # Get reco and true coords
    #----------------------------------------------------------------------------#
    print( len(output) )
    concat = lambda key : np.concatenate([out[key] for out in output])
    reco_coords = concat('reco')
    true_coords = concat('true')

    x, y, z = reco_coords.T
    dx, dy, dz = (reco_coords - true_coords).T
    # #----------------------------------------------------------------------------#
    # # Plot deviation vs drift coordinate
    # #----------------------------------------------------------------------------#
    # global plt
    # fig, axes = plt.subplots(2,1,figsize=(8, 8),sharex=True,sharey=True)

    # major_ticks_x = np.arange(-35, 5, 5)
    # minor_ticks_x = np.arange(-35, 5, 1)
    # major_ticks_y = np.arange(-15, 15, 5)
    # minor_ticks_y = np.arange(-15, 15, 1)

    # ax = axes[0]
    # # ax.scatter(z/10, dx/10)
    # ax.hexbin(z/10, dx/10, mincnt=1, cmap='viridis', bins='log')
    # ax.set_ylabel('$\mathrm{\Delta x}$ [cm]')
    # # ax.set_xticks(np.arange(0, 35, 5))
    # # ax.set_yticks(np.arange(-10, 10, 5))
    # # ax.grid(color = 'grey', linestyle = '--', linewidth = 0.5)

    # ax.set_xticks(major_ticks_x)
    # ax.set_xticks(minor_ticks_x, minor=True)
    # ax.set_yticks(major_ticks_y)
    # ax.set_yticks(minor_ticks_y, minor=True)

    # ax.grid(which='both', color = 'grey', linestyle = '--', linewidth = 0.5)
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', alpha=0.5)

    # ax = axes[1]
    # # ax.scatter(z/10, dy/10)
    # ax.hexbin(z/10, dy/10, mincnt=1, cmap='viridis', bins='log')
    # ax.set_ylabel('$\mathrm{\Delta y}$ [cm]')
    # ax.set_xlabel('z [cm]')
    # # ax.set_xticks(np.arange(0, 35, 5))
    # # ax.set_yticks(np.arange(-10, 10, 5))

    # ax.set_xticks(major_ticks_x)
    # ax.set_xticks(minor_ticks_x, minor=True)
    # ax.set_yticks(major_ticks_y)
    # ax.set_yticks(minor_ticks_y, minor=True)

    # ax.grid(which='both', color = 'grey', linestyle = '--', linewidth = 0.5)
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', alpha=0.5)

    # ax.set_xlim(-32,0)
    # ax.set_ylim(-10, 10)

    # plt.show()
    # # plt.savefig('res.eps')
    # #----------------------------------------------------------------------------#

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

        fig.update_scenes(
            aspectmode='data', #cube
            xaxis_range=(downstream,upstream),
            yaxis_range=(bottom,top),
            zaxis_range=(-anode_z,anode_z),
            xaxis_title='x [mm]',
            yaxis_title='y [mm]',
            zaxis_title='z [mm]',
        )

        fig.show()
        # plotly image will be rendered in web browser 
        # Can screen grab or save the image with the kaleido package (pip install)
        # fig.write_image("figure.pdf", engine="kaleido")

    f.close()

    # Show time elapsed for running the code
    # print(datetime.now() - startTime)

if __name__ == '__main__': 
    import argparse

    parser = argparse.ArgumentParser(description='Plot the first N tracks from a given file')
    parser.add_argument('infile',
                        help = 'input larpix data with track reconstruction')
    parser.add_argument('-n',
                        default = -1,
                        type = int,
                        help = 'evaluate the selection criteria over the first N tracks')
    parser.add_argument('-nshow',
                        default = 10,
                        type = int,
                        help = 'number of tracks to show')
    parser.add_argument('-g', '--geometry',
                        default = './pixel_layouts/module1_layout-2.3.16.yaml',
                        type = str,
                        help = 'path to the pixel layout YAML')
    parser.add_argument('-d', '--detector',
                        default = './detector_properties/module0.yaml',
                        type = str,
                        help = 'path to the detector properties YAML')
    # parser.add_argument('-o', '--output',
    #                     default = '',
    #                     type = str,
    #                     help = 'save the data which passes the selection to a file')
    parser.add_argument('-p', '--plot',
                        default = False,
                        type = bool,
                        help = 'show a plot of the selected tracks')

    args = parser.parse_args()

    main(args)
