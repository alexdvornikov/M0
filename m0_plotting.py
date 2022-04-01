import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as plt3d
import h5py


filename = 'datalog_2021_04_02_01_30_31_CESTevd.h5'
f = h5py.File(filename)
f['events']

for event in f['events']:
    print (event.dtype)
    print (event)
    eventHits = f['hits'][event['hit_ref']]
    
    t0 = event['ts_start']
    
    xs = eventHits['px']
    ys = eventHits['py']
    ts = eventHits['ts'] - t0

    q = eventHits['q']
 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    
    # Time is in 0.1 microseconds ticks (10 MHz clock). 
    # Divide by 10 to convert the tick numbers to microseconds.  
    ax.scatter(xs, ys, ts/10, c = q, cmap = 'plasma', alpha = 1) 

    #Set background to white
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax.set_xlim3d(-300, 300)
    ax.set_ylim3d(-600, 600)
    ax.set_zlim3d(0, 500)

    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('t [us]')
    # ax.set_zlabel('t [$\mu$ s]') #If have LaTex

    plt.show()
    plt.close()