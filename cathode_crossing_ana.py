import numpy as np
import os

TPC1_true = []
TPC1_reco = []

TPC2_true = []
TPC2_reco = []

inputDir = './cathode_crossers'

TPC1_true = np.concatenate([np.load(os.path.join(inputDir, infile))[0,:,:]
                            for infile in os.listdir(inputDir)])
TPC1_reco = np.concatenate([np.load(os.path.join(inputDir, infile))[1,:,:]
                            for infile in os.listdir(inputDir)])
TPC2_true = np.concatenate([np.load(os.path.join(inputDir, infile))[2,:,:]
                            for infile in os.listdir(inputDir)])
TPC2_reco = np.concatenate([np.load(os.path.join(inputDir, infile))[3,:,:]
                            for infile in os.listdir(inputDir)])
v_dir = np.concatenate([np.load(os.path.join(inputDir, infile))[4,:,:]
                        for infile in os.listdir(inputDir)])

import matplotlib.pyplot as plt
import matplotlib

print (TPC1_true.shape)

print(np.min(TPC1_true[:,1]))

bins = (np.linspace(-300, 300, 31),
        np.linspace(-800, 380, 61))

TPC1_ds = np.diff(np.array([TPC1_true, TPC1_reco]), axis = 0)[0]
TPC2_ds = np.diff(np.array([TPC2_true, TPC2_reco]), axis = 0)[0]

TPC1_dxHist = np.histogram2d(TPC1_reco[:,0],
                             TPC1_reco[:,1],
                             weights = TPC1_ds[:,0],
                             bins = bins)[0]
TPC1_dyHist = np.histogram2d(TPC1_reco[:,0],
                             TPC1_reco[:,1],
                             weights = TPC1_ds[:,1],
                             bins = bins)[0]
TPC1_dzHist = np.histogram2d(TPC1_reco[:,0],
                             TPC1_reco[:,1],
                             weights = TPC1_ds[:,2],
                             bins = bins)[0]
TPC1_denomHist = np.histogram2d(TPC1_reco[:,0],
                                TPC1_reco[:,1],
                                bins = bins)[0]

TPC1_dx_mean = TPC1_dxHist.T/TPC1_denomHist.T
TPC1_dy_mean = TPC1_dyHist.T/TPC1_denomHist.T
TPC1_dz_mean = TPC1_dzHist.T/TPC1_denomHist.T

TPC2_dxHist = np.histogram2d(TPC2_reco[:,0],
                             TPC2_reco[:,1],
                             weights = TPC2_ds[:,0],
                             bins = bins)[0]
TPC2_dyHist = np.histogram2d(TPC2_reco[:,0],
                             TPC2_reco[:,1],
                             weights = TPC2_ds[:,1],
                             bins = bins)[0]
TPC2_dzHist = np.histogram2d(TPC2_reco[:,0],
                             TPC2_reco[:,1],
                             weights = TPC2_ds[:,2],
                             bins = bins)[0]
TPC2_denomHist = np.histogram2d(TPC2_reco[:,0],
                                TPC2_reco[:,1],
                                bins = bins)[0]

TPC2_dx_mean = TPC2_dxHist.T/TPC2_denomHist.T
TPC2_dy_mean = TPC2_dyHist.T/TPC2_denomHist.T
TPC2_dz_mean = TPC2_dzHist.T/TPC2_denomHist.T

# print (hist[0])

from matplotlib import cm

# nonNaN = TPC1_dx_mean[~np.isnan(TPC1_dx_mean)]
# vext = np.max([np.max(nonNaN), -np.min(nonNaN)])
# plt.imshow(TPC1_dx_mean,
#            origin = 'lower',
#            extent = (np.min(bins[0]), np.max(bins[0]),
#                      np.min(bins[1]), np.max(bins[1])),
#            cmap = 'RdYlBu',
#            vmin = vext,
#            vmax = -vext)
# plt.imshow(TPC1_dy_mean,
#            origin = 'lower')
# plt.imshow(TPC1_dz_mean,
#            origin = 'lower')
nonNaN = TPC2_dx_mean[~np.isnan(TPC2_dx_mean)]
vext = np.max([np.max(nonNaN), -np.min(nonNaN)])
print ("vext", vext)
plt.imshow(TPC2_dx_mean,
           origin = 'lower',
           extent = (np.min(bins[0]), np.max(bins[0]),
                     np.min(bins[1]), np.max(bins[1])),
           cmap = 'RdYlBu',
           vmin = vext,
           vmax = -vext)

# plt.imshow(TPC2_dy_mean,
#            origin = 'lower')
# plt.imshow(TPC2_dz_mean,
#            origin = 'lower')
plt.xlabel(r'x [mm]')
plt.ylabel(r'y [mm]')
plt.title(r'TPC 2')
plt.gca().set_aspect('equal')
cb = plt.colorbar()
cb.set_label(r'$\Delta x$ [mm]')
plt.show()
