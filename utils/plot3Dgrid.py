import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mayavi import mlab
from scipy.interpolate import griddata

def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)

data = np.genfromtxt(sys.argv[1], delimiter=" ")  

print(data)
x = data[:,0]
y = data[:,1]
z = data[:,2]
density = data[:,4]

density = [0.0 if abs(a_) < 1.0e-4 else a_ for a_ in density]
xtp = []
ytp = []
ztp = []
densitytp = []

for idx, a in enumerate(density):
    if a != 0.0:
        xtp.append(x[idx])
        ytp.append(y[idx])
        ztp.append(z[idx])
        densitytp.append(a)

#dx, pts = 2, 100j
#R = np.column_stack((x, y, z))
#X,Y,Z = np.mgrid[-dx:dx:pts, -dx:dx:pts, -dx:dx:pts]
#F = griddata(R, density, (X,Y,Z))
F = np.column_stack((x, y, z, density))
#print(F)
#print(density.shape[0]**(1/3))
#density.shape = (50, 50, 50)
#source = mlab.pipeline.scalar_field(density)
#min = density.min()
#max = density.max()
#vol = mlab.pipeline.volume(source) 

#figure = mlab.figure('DensityPlot')
pts = mlab.points3d(xtp, ytp, ztp, densitytp, opacity=0.2, scale_mode='none', scale_factor=0.07)
#pts = mlab.contour3d(F,contours=8, opacity=.2)
#pts = mlab.contour_surf(F)
#pts = mlab.plot3d(x, y, z, density)
#mesh = mlab.pipeline.delaunay2d(pts)
#surf = mlab.pipeline.surface(mesh)

mlab.axes()
mlab.show()

"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

limit = 1e-5

# draw a sphere for each data point
for (xi,yi,zi,ri) in zip(x,y,z,r):
    rp = ri/10.0

    if abs(rp) > limit:
        if rp < 0.0:
            (xs,ys,zs) = drawSphere(xi,yi,zi,-1.0*rp)
            ax.plot_wireframe(xs, ys, zs, color="r")
        else:
            (xs,ys,zs) = drawSphere(xi,yi,zi,-1.0*rp)
            ax.plot_wireframe(xs, ys, zs, color="b")

plt.show()
"""