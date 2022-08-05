import numpy as np
import os
from tqdm import tqdm

projection = 'x'
scanDir = 'y'

if scanDir == 'x':
    plotAx = ('z', 'y')
elif scanDir == 'y':
    plotAx = ('x', 'z')
elif scanDir == 'z':
    plotAx = ('x', 'y')

sigma_axial = 1.
# sigma_axial = 1.
sigma_radial = 1.

unit_vectors = {'x': np.array([1, 0, 0]),
                'y': np.array([0, 1, 0]),
                'z': np.array([0, 0, 1]),
                }

bins = {'x': np.linspace(-300, 300, 31),
        'y': np.linspace(-800, 380, 61),
        'z': np.linspace(-300, 300, 31),
        }

histkwargs = {'bins': (bins['x'], bins['y'], bins['z'])}

denom_hist = np.zeros(shape = [len(bins['x'])-1,
                               len(bins['y'])-1,
                               len(bins['z'])-1])
disp_hist = np.zeros(shape = [len(bins['x'])-1,
                              len(bins['y'])-1,
                              len(bins['z'])-1])

inputDir = './cathode_crossers_volumetric'
for infile in tqdm(os.listdir(inputDir)[:10]):
    path = os.path.join(inputDir, infile)

    true = np.load(path)[0]
    reco = np.load(path)[1]
    pca_dir = np.load(path)[2]
    if np.any(reco):
        ds = np.diff(np.array([true, reco]) , axis = 0)[0]
        sigma = np.sqrt(np.power(np.dot(pca_dir,
                                        unit_vectors[projection])
                                 * sigma_axial, 2)
                        + np.power(np.linalg.norm(np.cross(pca_dir,
                                                           unit_vectors[projection]),
                                                  axis = -1)
                                   * sigma_radial, 2))

        displacement_proj = np.dot(ds, unit_vectors[projection])

        denom_hist += np.histogramdd(reco,
                                     weights = 1./sigma,
                                     **histkwargs)[0]
        disp_hist += np.histogramdd(reco,
                                    weights = displacement_proj/sigma,
                                    **histkwargs)[0]

disp_mean = disp_hist/denom_hist

from matplotlib import cm
import matplotlib.pyplot as plt

nonNaN = disp_mean[~np.isnan(disp_mean)]
vext = np.max([np.max(nonNaN), -np.min(nonNaN)])

scanBinCenters = 0.5*(bins[scanDir][1:] + bins[scanDir][:-1])
for i, binCenter in enumerate(scanBinCenters):
    # fig = plt.figure()
    if scanDir == 'x':
        thisSlice = disp_mean[i,:,:]
    elif scanDir == 'y':
        thisSlice = disp_mean[:,i,:].T
    elif scanDir == 'z':
        thisSlice = disp_mean[:,:,i].T

    plt.imshow(thisSlice,
               origin = 'lower',
               extent = (np.min(bins[plotAx[0]]),
                         np.max(bins[plotAx[0]]),
                         np.min(bins[plotAx[1]]),
                         np.max(bins[plotAx[1]])),
               cmap = 'RdYlBu',
               vmin = vext,
               vmax = -vext)
    plt.xlabel(plotAx[0]+r' [mm]')
    plt.ylabel(plotAx[1]+r' [mm]')
    plt.title(scanDir+r' = '+str(binCenter)+ r' mm')
    plt.gca().set_aspect('equal')
    cb = plt.colorbar()
    cb.set_label(r'$\Delta '+projection+r'$ [mm]')
    plt.subplots_adjust(left=-0.9, bottom=0.125, right=1, top=0.9)

    plt.savefig('cathode_crosser_volume_plots/delta_'+projection+'_slice_'+scanDir+'_'+str(binCenter)+'.png',
                dpi=300)
    # plt.show()
    plt.clf()
