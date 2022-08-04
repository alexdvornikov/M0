import numpy as np
import os
from tqdm import tqdm

bins = (np.linspace(-300, 300, 31),
        np.linspace(-800, 380, 61),
        np.linspace(-300, 300, 31),
        )

denom_hist = np.zeros(shape = [len(bins[0])-1,
                               len(bins[1])-1,
                               len(bins[2])-1])
dx_hist = np.zeros(shape = [len(bins[0])-1,
                            len(bins[1])-1,
                            len(bins[2])-1])

inputDir = './cathode_crossers_volumetric'
for infile in tqdm(os.listdir(inputDir)):
    path = os.path.join(inputDir, infile)

    true = np.load(path)[0]
    reco = np.load(path)[1]
    pca_dir = np.load(path)[2]
    if np.any(reco):
        ds = np.diff(np.array([true, reco]) , axis = 0)[0]

        denom_hist += np.histogramdd(reco,
                                     bins = bins)[0]
        dx_hist += np.histogramdd(reco,
                                  weights = ds[:,0],
                                  bins = bins)[0]

dx_mean = dx_hist/denom_hist

from matplotlib import cm
import matplotlib.pyplot as plt

nonNaN = dx_mean[~np.isnan(dx_mean)]
print(nonNaN)
vext = np.max([np.max(nonNaN), -np.min(nonNaN)])
zBinCenters = 0.5*(bins[2][1:] + bins[2][:-1])
for i, binCenter in enumerate(zBinCenters):
    fig = plt.figure()
    plt.imshow(dx_mean[:,:,i].T,
               origin = 'lower',
               extent = (np.min(bins[0]), np.max(bins[0]),
                         np.min(bins[1]), np.max(bins[1])),
               cmap = 'RdYlBu',
               vmin = vext,
               vmax = -vext)
    plt.xlabel(r'x [mm]')
    plt.ylabel(r'y [mm]')
    plt.title(r'z = '+str(binCenter)+ r' mm')
    plt.gca().set_aspect('equal')
    cb = plt.colorbar()
    cb.set_label(r'$\Delta x$ [mm]')
    plt.subplots_adjust(left=-0.9, bottom=0.125, right=1, top=0.9)

    # plt.savefig('cathode_crosser_volume_plots/slice_z_'+str(binCenter)+'.png',
    #             dpi=300)
    plt.show()
