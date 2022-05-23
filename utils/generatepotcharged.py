import math 
import sys

import numpy as np

rmin = 5.8
d = 1.9 
charge = 1 

s = []
for r in np.arange(0.1, (rmin+d), 0.1, dtype=float):

    if r >= rmin and r <= rmin+d:
        #s.append(0.302 / w[i])
        s.append(charge / (rmin + d))
    else:
        s.append(charge/r)

for v in s:
    print(v)
