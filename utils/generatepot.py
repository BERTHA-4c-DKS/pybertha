import math 
import sys

import numpy as np

gridfilename = ""

if len(sys.argv) == 2:
    gridfilename = sys.argv[1]
else:
    print("usage: ", sys.argv[0], " gridfilename.txt")
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

rmin = 5.8
d = 1.9 

s = []
for i, v in enumerate(w):
    x = xs[i]
    y = ys[i]
    z = zs[i]

    r = math.sqrt(x**2 + y**2 + z**2)

    if r >= rmin and r <= rmin+d:
        s.append(0.302 / w[i])
    else:
        s.append(0.000)

for v in s:
    print(v)