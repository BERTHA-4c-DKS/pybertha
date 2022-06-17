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

fp.close() 

#print("X limits: ", min(xs), max(xs))
#print("Y limits: ", min(ys), max(zs))
#print("Z limits: ", min(zs), max(zs))

#rmin = 5.8
#d = 1.9 
#potval = -0.30134
#
#s = []
#for i, v in enumerate(w):
#    x = xs[i]
#    y = ys[i]
#    z = zs[i]
#
#    r = math.sqrt(x**2 + y**2 + z**2)
#
#    if r >= rmin and r <= rmin+d:
#        #s.append(0.302 / w[i])
#        s.append(potval)
#    else:
#        s.append(0.000)

#rminp = 4.5
#rmaxp = 5.8
#rminn = 5.8
#rmaxn = 7.8
#npotval = -0.2
#ppotval = 0.1
#
#s = []
#for i, v in enumerate(w):
#    val = 0.0
#
#    x = xs[i]
#    y = ys[i]
#    z = zs[i]
#
#    r = math.sqrt(x**2 + y**2 + z**2)
#
#    if (r >= rminp) and (r < (rmaxp)):
#        val = ppotval
#
#    if (r >= rminn) and (r <= (rmaxn)):
#        val = npotval
#
#    s.append(val)

rminn_1 = 2.6
rmaxn_1 = 4.5
rminp = 4.5
rmaxp = 5.8
rminn = 5.8
rmaxn = 7.8
npotval_1 = -0.005
ppotval = 0.1
npotval = -0.2

s = []
for i, v in enumerate(w):
    val = 0.0

    x = xs[i]
    y = ys[i]
    z = zs[i]

    r = math.sqrt(x**2 + y**2 + z**2)

    if (r >= rminn_1) and (x < (rmaxn_1)):
        vsl = npotval_1

    if (x >= rminp) and (x < (rmaxp)):
        val = ppotval

    if (x >= rminn) and (x <= (rmaxn)):
        val = npotval

    s.append(val)

for v in s:
    print(v)
