import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata as gd
from scipy.interpolate import interpn as itp

from gridData import Grid

N = 40
WEI = False

gridfilename = ""

if len(sys.argv) != 2:
    print("usage: ", sys.argv[0] , " gridfilename.txt")
    exit(1)
else:
    gridfilename = sys.argv[1]

fp = open(gridfilename)
w = []

xs = []
ys = []
zs = []

for line in fp:
    sline = line.split()
    if len(sline) == 4:
        x = np.float64(sline[0])
        y = np.float64(sline[1])
        z = np.float64(sline[2])
        v = np.float64(sline[3])
        w.append(v)
        xs.append(x)
        ys.append(y)
        zs.append(z)
    else:
        print("ERROR at line ", line)

fp.close()

#print(max(xs), max(ys), max(zs))

dx = (abs(min(xs)))/int(N/2)
dy = (abs(min(ys)))/int(N/2)
dz = (abs(min(zs)))/int(N/2)

#print("Regular_grid")
print((N+1)**3)
x = min(xs)
for ix in range(N+1):
    y = min(ys)
    for iy in range(N+1):
        z = min(zs)
        for iz in range(N+1):
            w = dx*dy*dz
            print("%.15e %.15e %.15e %.5e"%(x, y, z, w))
            z += dz
        y += dy
    x += dx
