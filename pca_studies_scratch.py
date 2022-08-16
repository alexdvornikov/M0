import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from utils_m1 import *
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')


def SciForm(number):
    return "{:e}".format(number)

def get_coord(coord_idx):
    # coord_idx is eitehr 0, 1, or 2 (for x,y,z)

    q = []
    for a in pts:
        for b in a:
            q.append( b[coord_idx] )
    q = np.array(q, dtype=object)

    blah = []
    for e in q:
        blah.append( e.tolist() )

    Q = [x for xs in blah for x in xs] #Flatten list of lists

    return Q



def distortions_2anodes(pts,endpts):

    pos3d = np.column_stack( [ pts[0], pts[1], pts[2] ] )
    endpts3d = np.column_stack( [ endpts[0], endpts[1], endpts[2] ] )
   
    pca = PCA(1)
    pca.fit(pos3d)
    
    v_dir = pca.components_[0]
    # want v_dir to have well-defined handedness
    # to make a consistent sense of TPC ordering
    if v_dir[2] > 0:
        v_dir *= -1

    pca = PCA(1)
    pca.fit(endpts3d)
    AA_dir = pca.components_[0]
    if AA_dir[2] > 0:
        AA_dir *= -1

    dot_product = np.dot(v_dir, AA_dir)
    rad = np.arccos(dot_product)
    deg = np.rad2deg(rad)

    return deg

#--------------------------------------------------------------------------------------------------------------#

endpts = []

dir = 'path/to/dir'
for filename in os.listdir(dir):
    file = os.path.join(dir, filename)
    d = ( ( np.load(file, allow_pickle=True) )['data'] ).item()
    endpts.append(  np.array( d['start'] )  )
    endpts.append( np.array( d['end'] )  )

endpts = np.array(endpts, dtype=object)

zs = []
for pts in endpts:
    zs.append( pts[:,2] )
zs = np.array(zs, dtype=object)


zs_lists = []
for e in zs:
    zs_lists.append( e.tolist() )

z_flat = [x for xs in zs_lists for x in xs] #Flatten list of lists

#--------------------------------------------------------------------------------------------------------------#
pts = []
endpts = []

dir = '/Users/alex/Desktop/m1_AA_crossers'
for filename in os.listdir(dir):
    file = os.path.join(dir, filename)
    d = ( ( np.load(file, allow_pickle=True) )['data'] ).item()

    d['pos']
    d['start']
    d['end']



    pts.append( d['pos']  )
    endpts.append(  np.array( d['start'] )  )
    endpts.append( np.array( d['end'] )  )

pts = np.array(pts, dtype=object)
endpts = np.array(endpts, dtype=object)






X = get_coord(0)
Y = get_coord(1)
Z = get_coord(2)