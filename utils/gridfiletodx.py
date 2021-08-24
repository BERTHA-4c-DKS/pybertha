import sys
import numpy as np
from gridData import Grid

data = np.genfromtxt(sys.argv[1], delimiter=" ")  

xmin = min(data[:,0])
ymin = min(data[:,1])
zmin = min(data[:,2])
density = data[:,4]

dstep = data[0,2] - data[1,2]
density.shape = [50, 50, 50]

g = Grid(density, origin=[xmin, ymin, zmin], \
        delta=[dstep, dstep, dstep])

g.export(sys.argv[1].replace(".txt", ".dx"))
