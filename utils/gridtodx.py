import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata as gd

from gridData import Grid


N = 50
WEI = False

gridfilename = ""
fieldfilename = ""

if len(sys.argv) != 4:
    print("usage: ", sys.argv[0] , " gridfilename.txtn fieldfilename.txt WEI=[0,1]")
    exit(1)
else:
    gridfilename = sys.argv[1]
    fieldfilename = sys.argv[2]

    if (int(sys.argv[3]) == 1):
        WEI = True
    else:
        WEI = False


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

fp = open(fieldfilename)
s = []

rmin = np.float64("inf")
rmax = np.float64("-inf")

for i, line in enumerate(fp):
    v = np.float64(line)
    if v < rmin :
        rmin = v
    if v > rmax :
        rmax = v
    if WEI :
        s.append(v*w[i])
    else:
        s.append(v)

fp.close()

if len(s) != len(xs):
    print("Error in file sizes")
    exit(1)

#generate new grid
xi, yi, zi = np.ogrid[min(xs):max(xs):N*1j, \
    min(ys):max(ys):N*1j, \
        min(zs):max(zs):N*1j]

X1 = xi.reshape(xi.shape[0],)
Y1 = yi.reshape(yi.shape[1],)
Z1 = zi.reshape(zi.shape[2],)

ar_len = len(X1)*len(Y1)*len(Z1)
X = np.arange(ar_len,dtype=float)
Y = np.arange(ar_len,dtype=float)
Z = np.arange(ar_len,dtype=float)

l=0
for i in range(0,len(X1)):
    for j in range(0,len(Y1)):
        for k in range(0,len(Z1)):
            X[l]=X1[i]
            Y[l]=Y1[j]
            Z[l]=Z1[k]
            l=l+1

dxstep = (max(xs)-min(xs))/N
dzstep = (max(ys)-min(ys))/N
dystep = (max(zs)-min(zs))/N

xold = X[0]
for xv in X:
    if xv != xold:
       dxstep = abs(xold - xv)
       break 

yold = Y[0]
for yv in Y:
    if yv != yold:
       dystep = abs(yold - yv)
       break 

zold = Z[0]
for zv in Z:
    if zv != zold:
       dzstep = abs(zold - zv)
       break 

print("Steps: ", dxstep, dzstep, dystep)
#print("        ", X[1]-X[0], Y[1]-Y[0], Z[1]-Z[0])
#print(X)
#print(Y)
#print(Z)

print("Interpolate...")
S = gd((xs,ys,zs), s, (X,Y,Z), method='linear')
print("")

#s = np.sin(x*y*z)/(x*y*z)
S = S.reshape(N, N, N)

#print(S.shape, type(S))

g = Grid(S, origin=[min(xs), min(ys), min(zs)], \
        delta=[dxstep, dystep, dzstep])

g.export(fieldfilename.replace(".txt", ".dx"))

"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

def f(x, y, z):
  return 2 * x**3 + 3 * y**2 - z

x = np.linspace(1, 4, 11)
y = np.linspace(4, 7, 22)

z = np.linspace(7, 9, 33)
xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)

print(xg.shape)
print(yg.shape)
print(zg.shape)

data = f(xg, yg, zg)

print(data.shape)

my_interpolating_function = RegularGridInterpolator((x, y, z), data)
pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])

my_interpolating_function(pts)

"""
