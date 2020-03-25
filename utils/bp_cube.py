# axes.py
from pymol.cgo import *
from pymol import cmd

#cmd.set('connect_cutoff', 1.7)
cmd.load('_xyzfile.xyz')
#cmd.bond('id 1', 'id 2')
#cmd.bond('id 1', 'id 3')
#cmd.bond('id 1', 'id 4')

#cmd.set_bond ('stick_color', 'white', 'all', 'all')
#cmd.set_bond ('stick_radius', -0.14, 'all', 'all')
#cmd.set ('stick_ball', 1)
#cmd.set ('stick_ball_ratio', -1)
#cmd.set ('stick_ball_color', 'atomic')
#cmd.show ('sticks', 'all')

#cmd.color('black', 'id 1')
#cmd.color('gray', '(name Au*)')

cmd.set_bond ('stick_radius', 0.1, 'all', 'all')
cmd.set ('sphere_scale', 0.15, 'all')
cmd.show ('sticks', 'all')
cmd.show ('spheres', 'all')
cmd.set ('stick_ball_color', 'atomic')


cmd.color('gray20', '(name C*)')

cmd.set ('transparency_mode', 1)
 
#w = 0.01 # cylinder width 
#l = 0.5 # cylinder length
#h = 0.15 # cone hight
#d = w * 1.618 # cone base diameter
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

cmd.bg_color('white')


#cmd.set_bond ('stick_color', 'white', 'clauc2h2', 'clauc2h2')
#cmd.set_bond ('stick_radius', -0.14, 'clauc2h2', 'clauc2h2')
#cmd.set ('stick_ball', 1)
#cmd.set ('stick_ball_ratio', -0.5)
#cmd.set ('stick_ball_color', 'atomic')
#cmd.show ('sticks', 'clauc2h2')
##cmd.set ('sphere_scale', 0.25, 'clauc2h2')
##cmd.show ('spheres', 'clauc2h2')

cmd.load('this.dx')
cmd.isosurface('sp', 'this', _iso)
cmd.color('_colorp', 'sp')
cmd.isosurface('sm', 'this', -_iso)
cmd.color('_colorm', 'sm')

#cmd.set ('label_font_id', 16)
#cmd.set ('label_size', 24)
#cmd.pseudoatom('xatom', pos=[1,0,0], label="x")
#cmd.pseudoatom('yatom', pos=[0,1,0], label="y")
#cmd.pseudoatom('zatom', pos=[0,0,2], label="z")
cmd.set ('label_font_id', 16)
cmd.set ('label_size', 20)
cmd.pseudoatom('xatom', pos=[2.5,0,0], label="x")
cmd.pseudoatom('yatom', pos=[0,2.5,0], label="y")
cmd.pseudoatom('zatom', pos=[0,0,5.5], label="z")

cmd.set('ray_opaque_background', 'on')
cmd.set('transparency', 0.50)

cmd.set_view ('\
    -0.086753346,    0.836276770,    0.541441798,\
     0.670052588,   -0.353231400,    0.652936995,\
     0.737265468,    0.419426382,   -0.529700875,\
     0.000000000,    0.000000000,  -34.167503357,\
     0.108778238,    0.131587982,    2.662951946,\
    26.937919617,   41.397087097,  -20.000000000' )

cmd.pseudoatom('latom', pos=[3.5,-2.5,0], label="_label")

cmd.ray(1024,1024)
cmd.png('this.png', dpi=600)

cmd.quit()
