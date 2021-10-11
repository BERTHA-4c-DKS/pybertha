import math 
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

r = 1

N = 1000
ncount = 0
a = (4.0*math.pi*r**2)/N
d = math.sqrt(a)
Mt = int(math.pi/d)
dt = math.pi/Mt
dg = a/dt

xs = []
ys = []
zs = []

#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')

for m in range(Mt):
    t = math.pi*(m+0.5)/Mt
    Mg = int(2.0 * math.pi * math.sin(t) / dg)
    for n in range(Mg):
        g = 2.0*math.pi*n / Mg
        x = r * math.sin(t) * math.cos(g)
        y = r * math.sin(t) * math.sin(g)
        z = r * math.cos(t)

        #ax.scatter(x, y, z)

        #print (x, y, z)
        xs.append(x)
        ys.append(y)
        zs.append(z)
        ncount += 1

#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

print(min(xs), max(xs))
print(min(ys), max(ys))
print(min(zs), max(zs))

#plt.show()
