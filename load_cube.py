from pymol.cgo import *
from pymol import cmd
import glob
import os

lista = ["diff_tot.cube", "diff_tot_ortho.cube", "nocv-1.cube", "nocv+1.cube", "pair1.cube"]

cmd.load('aucn.xyz')

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
