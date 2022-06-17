import math 
import sys


max = 35.0
dx = 0.1

n = int(max/dx)

rmin = 5.8
d = 1.9 
potval = -0.30134


#x = 0.0
#for i in range(n):
#    v = 0.0
#
#    if (x >= rmin) and (x <= (rmin+d)):
#        v = potval
#    
#    print("%8.3f %8.3f"%(x, v))
#    x += dx 

#rminp = 4.5
#rmaxp = 5.8
#rminn = 5.8
#rmaxn = 7.8
#npotval = -0.2
#ppotval = 0.1
#
#x = 0.0
#for i in range(n):
#    v = 0.0
#
#    if (x >= rminp) and (x < (rmaxp)):
#        v = ppotval
#
#    if (x >= rminn) and (x <= (rmaxn)):
#        v = npotval
#    
#    print("%8.3f %8.3f"%(x, v))
#    x += dx 

rminn_1 = 2.6
rmaxn_1 = 4.5
rminp = 4.5
rmaxp = 5.8
rminn = 5.8
rmaxn = 7.8
npotval_1 = -0.005
ppotval = 0.1
npotval = -0.2

x = 0.0
for i in range(n):
    v = 0.0

    if (x >= rminn_1) and (x < (rmaxn_1)):
        v = npotval_1

    if (x >= rminp) and (x < (rmaxp)):
        v = ppotval

    if (x >= rminn) and (x <= (rmaxn)):
        v = npotval
    
    print("%8.3f %8.3f"%(x, v))
    x += dx 
