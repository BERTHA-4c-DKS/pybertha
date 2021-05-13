# axes.py
from pymol.cgo import *
from pymol import cmd

cmd.set('connect_cutoff', 2.5)
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

cmd.set_bond ('stick_radius', 0.05, 'all', 'all')
cmd.set ('sphere_scale', 0.15, 'All')
cmd.show ('sticks', 'all')
cmd.show ('spheres', 'all')
cmd.set ('stick_ball_color', 'atomic')
cmd.color('gray20', '(name C*)')
cmd.color('blue', '(name L*)')
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
       CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   2*l, w, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 
       CONE, 0.0, 0.0,   2*l, 0.0, 0.0, h+2*l, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]

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

cmd.load('this.cube')
cmd.isosurface('sp', 'this', _iso)
cmd.color('green', 'sp')
cmd.isosurface('sm', 'this', -_iso)
cmd.color('_colorm', 'sm')


#cmd.load('dif_B1.cube')
#cmd.isosurface('spp', 'dif_B1', 0.00008)
#cmd.color('blue', 'spp')
#cmd.isosurface('smp', 'dif_B1', -0.00008)
#cmd.color('red', 'smp')


#cmd.set ('label_font_id', 16)
#cmd.set ('label_size', 24)
#cmd.pseudoatom('xatom', pos=[1,0,0], label="x")
#cmd.pseudoatom('yatom', pos=[0,1,0], label="y")
#cmd.pseudoatom('zatom', pos=[0,0,2], label="z")
cmd.set ('label_font_id', 16)
cmd.set ('label_size', 20)
cmd.pseudoatom('xatom', pos=[2.5,0,0], label="x")
cmd.pseudoatom('yatom', pos=[0,2.5,0], label="y")
cmd.pseudoatom('zatom', pos=[0,0,4.0], label="z")

cmd.set('ray_opaque_background', 'off')
cmd.set('transparency', 0.30)

#cmd.set_view ('\
#    -0.086753346,    0.836276770,    0.541441798,\
#     0.670052588,   -0.353231400,    0.652936995,\
#     0.737265468,    0.419426382,   -0.529700875,\
#     0.000000000,    0.000000000,  -34.167503357,\
#     0.108778238,    0.131587982,    1.662951946,\
#     26.937919617,   41.397087097,  -20.000000000' )
#

cmd.set_view ('\
    -0.036295272,    0.997309446,    0.063655443,\
    -0.180431277,   -0.069193237,    0.981151760,\
     0.982915938,    0.024127077,    0.182457983,\
     0.000000000,    0.000000000,  -42.977798462,\
     0.235994816,    0.235994816,    0.235994816,\
    33.884025574,   52.071571350,  -20.000000000' )




#cmd.pseudoatom('latom', pos=[3.5,-2.5,0], label="_label")
#cmd.pseudoatom('latom', pos=[-3,0,0], label=u"\u0394\u03C1'\u2081".encode('utf-8'))
cmd.pseudoatom('latom', pos=[-3,0,0], label=u"_label".encode('utf-8'))
cmd.ray(1200,1200)
cmd.png('this.png', dpi=600)

#cmd.quit()
