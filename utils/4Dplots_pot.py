from os import path
import numpy as np
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata as gd


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

TYPES = 4
N = 50
CUTTHR = -0.1

gridfilename = ""
fieldfilename = ""

if len(sys.argv) != 4:
    print("usage: ", sys.argv[0] , " gridfilename.txtn fieldfilename.txt TYPES")
    exit(1)
else:
    gridfilename = sys.argv[1]
    fieldfilename = sys.argv[2]
    TYPES =int(sys.argv[3])

fp = open(gridfilename)
grid = []
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
        grid.append((x, y, z, v))
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

for line in fp:
    v = np.float64(line)
    if v < rmin :
        rmin = v
    if v > rmax :
        rmax = v
    if v < CUTTHR:
        v = 0.0
    s.append(v)

fp.close()

if len(s) != len(xs):
    print("Error in file sizes")
    exit(1)

#print(len(set(xs)), len(xs))
#print(len(set(ys)), len(ys))
#print(len(set(zs)), len(zs))

vvmin = np.percentile(np.array(s), 20)
vvmax = np.percentile(np.array(s), 80)
print("20 percetile  and 80 percentile: ", vvmin, vvmax)
print("Real Min and Max: ", rmin, rmax)
print("After cut Min and Max: ", min(s), max(s))

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

print("Interpolate...")
S = gd((xs,ys,zs), s, (X,Y,Z), method='linear')
print("")

#s = np.sin(x*y*z)/(x*y*z)
#print(S.shape, type(S))
#print(min(xs), max(xs))
#mlab.contour3d(s)


if TYPES == 0:
    n_bins = 100

    plt.hist(s, bins = n_bins)

    plt.show()

elif TYPES == 1:
    #Plot original values
    fig1 = plt.figure()
    ax1=fig1.gca(projection='3d')
    sc1=ax1.scatter(xs, ys, zs, c=s, cmap=plt.hot())
    plt.colorbar(sc1)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')

    #Plot interpolated values
    fig2 = plt.figure()
    ax2=fig2.gca(projection='3d')
    sc2=ax2.scatter(X, Y, Z, c=S, cmap=plt.hot())
    plt.colorbar(sc2)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')

    #Show plots
    plt.show()

elif TYPES == 2:
    S = S.reshape(N, N, N)
    mlab.pipeline.volume(mlab.pipeline.scalar_field(S))
    mlab.show()
elif TYPES == 3:
    S = S.reshape(N, N, N)
    mlab.pipeline.volume(mlab.pipeline.scalar_field(S), vmin=vvmin, vmax=vvmax)
    mlab.show()
elif TYPES == 4:
    S = S.reshape(N, N, N)

    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(S),
                            plane_orientation='x_axes',
                            slice_index=int(N/2),
                        )
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(S),
                            plane_orientation='y_axes',
                            slice_index=int(N/2),
                        )
    mlab.outline()
    mlab.volume_slice(S, plane_orientation='x_axes', slice_index=30)
    mlab.show()
elif TYPES == 5:
    S = S.reshape(N, N, N)
    mlab.contour3d(S)
    mlab.show()

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
