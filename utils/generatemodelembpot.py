import math 
import sys

import numpy as np
from numpy.lib.polynomial import _raise_power

gridfilename = ""
avgspherefname = ""

if len(sys.argv) == 3:
    gridfilename = sys.argv[1]
    avgspherefname = sys.argv[2]
else:
    print("usage: ", sys.argv[0], " gridfilename.txt phereavgfilename.txt")
    exit(1)

fp = open(gridfilename, "r")

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

rpotvalues = {}

fp = open(avgspherefname, "r")

for l in fp:
    sline = l.split()
    if len(sline) == 2:
        r = float(sline[0])
        v = float(sline[1])

        rpotvalues[r] = v

    else:
        print("Error at line ", sline)        

fp.close()

s = []
for i, v in enumerate(w):
    x = xs[i]
    y = ys[i]
    z = zs[i]

    r = math.sqrt(x**2 + y**2 + z**2)

    inserted = False
    rmin = 0
    for ra in rpotvalues:
        potval = rpotvalues[ra]
        if r >= rmin and r < ra:
            s.append(potval)
            inserted = True
        rmin = ra

    if not inserted:
         s.append(0.0)

for v in s:
    print(v)
