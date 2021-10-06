from os import path
import numpy as np
from mayavi import mlab
from scipy.interpolate import RegularGridInterpolator

import sys

filename = ""

if len(sys.argv) != 2:
    print("usage: ", sys.argv[0] , " filename.txt")
    exit(1)
else:
    filename = sys.argv[1]

fp = open(filename)
grid = []
s = []

xs = []
ys = []
zs = []

for line in fp:
    sline = line.split()
    if len(sline) == 4:
        x = float(sline[0])
        y = float(sline[1])
        z = float(sline[2])
        v = float(sline[3])
        grid.append((x, y, z, v))
        s.append(v)
        xs.append(x)
        ys.append(y)
        zs.append(z)
    else:
        print("ERROR at line ", line)

print(len(set(xs)), len(xs))

print(len(set(ys)), len(ys))

print(len(set(zs)), len(zs))

x, y, z = np.ogrid[min(xs):max(xs):100j, \
    min(ys):max(ys):100j, \
        min(zs):max(zs):100j]

s = np.sin(x*y*z)/(x*y*z)

print(s.shape)

print(min(xs), max(xs))

#mlab.contour3d(s)


#mlab.pipeline.volume(mlab.pipeline.scalar_field(s), vmin=0, vmax=0.8)


mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='x_axes',
                            slice_index=50,
                        )
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='y_axes',
                            slice_index=50,
                        )
mlab.outline()
mlab.volume_slice(s, plane_orientation='x_axes', slice_index=30)

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