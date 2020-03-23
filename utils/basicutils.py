import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "/home/redo/Sources/GridDataFormats/")

from gridData import OpenDX
from ase.io.cube import read_cube

def cubetodx(infname, outfname, valperl = 3, scientfn=False):

     OpenDX.global_values_per_line.append(valperl)
     if scientfn:
         OpenDX.global_float_format.append("e")
         OpenDX.global_float_precision.append(8)

     conv = 0.529177

     fp = open(infname)
     cubef = read_cube(fp)
     fp.close()

     origin = cubef['origin']
     grid_shape = cubef['data'].shape
     data = cubef['data']

     with open(infname) as fp:    
         next(fp)
         next(fp)
         array = []
         for i in range(4):
             head = next(fp)
             array.append([float(x) for x in head.split()])
     header = np.array(array)
     delta = [float(header[m,m])*conv for m in range(1,4)]
     delta = np.diagflat(delta)
     fp.close()

     new = OpenDX.field('density')
     new.add('positions', OpenDX.gridpositions(1, grid_shape, origin, delta))
     new.add('connections', OpenDX.gridconnections(2, grid_shape))
     new.add('data',OpenDX.array(3, data, type='double'))
     new.write(outfname)

def exportheatmap(occups, virtua, valuesin, xlabel = "Occ", \
            ylabel = "Virt", title = "Analysis of ...", \
            outfname = "test.jpg", norm = 1.0):


    values = np.divide(valuesin, norm)

    df = pd.DataFrame(data=values, \
            index=virtua, \
            columns=occups)
    
    labels = np.asarray(["%4.2f"%v \
            for v in list(values.flatten())]).reshape(values.shape)
    
    fig, ax = plt.subplots(figsize=(values.shape[0]*2, values.shape[1]+2))
    
    plt.title(title, fontsize=18)
    
    ttl = ax.title
    ttl.set_position([0.5, 1.05])
    
    #ax.set_xticks(occups)
    #ax.set_yticks(virtua)
    #ax.axis('off')
    
    sns.heatmap(df, annot=labels, fmt="", cmap="Blues", \
            linewidths=1, linecolor='black',ax=ax)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    #plt.show()
    fig.savefig('test.jpg')
