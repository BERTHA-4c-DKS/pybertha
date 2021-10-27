import sys

from gridData import Grid

filename = ""
iadd = 0

if len(sys.argv) != 3:
    print("usage: ", sys.argv[0] , " gridfilename.txtb extrai")
    exit(1)
else:
    filename = sys.argv[1]
    iadd = int(sys.argv[2])

#g = Grid(filename)
gin = Grid(filename)
g = gin.resample_factor(4)

nx = g.grid.shape[0]
ny = g.grid.shape[1]
nz = g.grid.shape[2]

deltax = g.delta[0]
deltay = g.delta[1]
deltaz = g.delta[2]

#print(nx, ny, nz)
i0x = 0
x0 = 0.0
x = g.origin[0]
diff = float("inf")
for ix in range(nx):
    #print(ix, x)
    if abs(x) < diff:
        diff = abs(x)
        i0x = ix
        x0 = x
    x = x + deltax

i0y = 0
y0 = 0.0
y = g.origin[1]
diff = float("inf")
for iy in range(ny):
    #print(ix, x)
    if abs(y) < diff:
        diff = abs(y)
        i0y = iy
        y0 = y
    y = y + deltay

i0z = 0
z0 = 0.0
z = g.origin[2]
diff = float("inf")
for iz in range(nz):
    #print(ix, x)
    if abs(z) < diff:
        diff = abs(z)
        i0z = iz
        z0 = z
    z = z + deltaz

#print(i0x, i0y, i0z)
#print(x0, y0, z0)

"""
rmin = 5.8
d = 1.9 
potval = -0.30134

x = g.origin[0]
for i in range(nx):
    if abs(x) >= rmin and abs(x) <= (rmin + d):
        v = potval 
    else :
        v = 0.0 
    print(x, v)
    x = x + deltax

"""

i0z += iadd

x = g.origin[0]
for i in range(nx):
    print(x, g.grid[i, i0y, i0z])
    x = x + deltax

"""
y = g.origin[1]
for i in range(ny):
    print(y, g.grid[i0x, i, i0z])
    y = y + deltay

z = g.origin[2]
for i in range(nz):
    print(z, g.grid[i0x, i0y, i])
    z = z + deltaz
"""