import sys
import numpy

fp = open(sys.argv[1], "r")

mx = 0.0e0 

for l in fp:
    sl = l.split()

    v1 = numpy.float64(sl[0]) 
    v2 = numpy.float64(sl[1]) 

    print(v2-v1)

    if (v2-v1) > mx:
        mx = v2-v1


print("Maxval: ", mx)