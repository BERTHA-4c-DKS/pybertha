from pymol.cgo import *
from pymol import cmd
import glob
import os
import sys

lista = []
xyzfile = ""
#lista = ["diff_tot.cube", "diff_tot_ortho.cube", "nocv-1.cube", "nocv+1.cube", "pair1.cube"]
#cmd.load('AuCn+.xyz')

print(sys.argv)

if len(sys.argv) != 5:
    print("usage: ", sys.argv[0], " cube xyzfile")
    exit(1)
else:
    lista.append(sys.argv[3])
    xyzfile = sys.argv[4]

cmd.load(xyzfile)

cmd.set_bond ('stick_radius', 0.1, 'all', 'all')
cmd.set ('sphere_scale', 0.15, 'all')
cmd.show ('sticks', 'all')
cmd.show ('spheres', 'all')
cmd.set ('stick_ball_color', 'atomic')

cmd.color('gray20', '(name C*)')
cmd.set ('transparency_mode', 1)

w = 0.03 # cylinder width
l = 1.65 # cylinder length
h = 0.2 # cone hight
d = w * 1.618 # cone base diameter

obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   3*l, w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
       CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
       CONE, 0.0, 0.0,   3*l, 0.0, 0.0, h+3*l, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]

cmd.load_cgo(obj, 'axes')

cmd.set ('label_font_id', 16)
cmd.set ('label_size', 20)
cmd.pseudoatom('xatom', pos=[2.5,0,0], label="x")
cmd.pseudoatom('yatom', pos=[0,2.5,0], label="y")
cmd.pseudoatom('zatom', pos=[0,0,5.5], label="z")


cmd.bg_color('white')
cmd.set('ray_opaque_background', 'off')
cmd.set('transparency', 0.50)

for idx, name in enumerate(lista):
    cmd.load(name)
    basename = os.path.splitext(name)[0]
    cmd.isosurface(basename+'_sp', basename, +0.0014)
    cmd.color('blue', basename+'_sp')
    cmd.isosurface(basename+'_sm', basename, -0.0014)
    cmd.color('red', basename+'_sm')


cmd.set_view ('\
    -0.086753346,    0.836276770,    0.541441798,\
     0.670052588,   -0.353231400,    0.652936995,\
     0.737265468,    0.419426382,   -0.529700875,\
     0.000000000,    0.000000000,  -34.167503357,\
     0.108778238,    0.131587982,    2.662951946,\
    26.937919617,   41.397087097,  -20.000000000' )

cmd.ray(600,600)
cmd.png('this.png', dpi=150)

cmd.quit()

