import numpy as np
from mayavi import mlab

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


x, y, z = np.ogrid[-10:10:20j, -10:10:20j, -10:10:20j]
s = np.sin(x*y*z)/(x*y*z)

print(s.shape)

print(min(xs), max(xs))

#mlab.contour3d(s)


#mlab.pipeline.volume(mlab.pipeline.scalar_field(s), vmin=0, vmax=0.8)


mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='x_axes',
                            slice_index=10,
                        )
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='y_axes',
                            slice_index=10,
                        )
mlab.outline()
mlab.volume_slice(s, plane_orientation='x_axes', slice_index=10)

mlab.show()