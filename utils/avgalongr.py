import sys
import math

from gridData import Grid

filename = ""

if len(sys.argv) != 2:
    print("usage: ", sys.argv[0] , " gridfilename.txt ")
    exit(1)
else:
    filename = sys.argv[1]

gin = Grid(filename)
g = gin.resample_factor(4)

nx = g.grid.shape[0]
ny = g.grid.shape[1]
nz = g.grid.shape[2]

deltax = g.delta[0]
deltay = g.delta[1]
deltaz = g.delta[2]

dr = 3.0 * max(deltax, deltay, deltaz)

xmin = g.origin[0]
ymin = g.origin[1]
zmin = g.origin[2]

xmax = xmin + (nx-1)*deltax
ymax = ymin + (ny-1)*deltay
zmax = zmin + (nz-1)*deltaz

Rmax = min(xmax-xmin, ymax-ymin, zmax-zmin)/2.0
nr = int(Rmax/dr)-1

sumofpoints = []
numofpoints = []
rminramx = []

r = 0.0
for i in range(nr):
    sumofpoints.append(0.0)
    numofpoints.append(0)
    rminramx.append((r, r+dr))
    r += dr

print(dr, Rmax, nr, file=sys.stderr)

for ix in range(nx):
    print(ix , " of ", nx, file=sys.stderr)
    x = xmin + ix*deltax
    for iy in range(ny):
        y = ymin + iy*deltay
        for iz in range(nz):
            z = zmin + iz*deltaz

            dist = math.sqrt(x**2 + y**2 + z**2)

            for i in range(nr):
                if dist >= rminramx[i][0] and dist < rminramx[i][1]:
                    sumofpoints[i] += g.grid[ix, iy, iz]
                    numofpoints[i] += 1
                    break;

r = dr
for i in range(nr):
    print(r, sumofpoints[i]/numofpoints[i])
    r += dr