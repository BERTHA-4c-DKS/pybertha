import sys
import math
import numpy as np

import matplotlib.pyplot as plt
import sklearn

from scipy.optimize import curve_fit, leastsq

T = 1000.0
k = 0.00001

###########################################################

def func(x, f, w):

    x = x[...,np.newaxis] 

    return np.sum (f * x/w * \
            (np.sin((x - w)*T) / (math.pi * (x - w))),  \
            axis=-1)

###########################################################

def sfunc(x, f, w0):

    return (f * (x/w0) * \
            (np.sin((x - w0)*T) / (math.pi * (x - w0)) ))

###########################################################


filename = ""

xcutmin = 0.26
xcutmax = 0.36

if len(sys.argv) == 4:
    filenname = sys.argv[1]
    xcutmin = np.float64(sys.argv[2])
    xcutmax = np.float64(sys.argv[3])
else:
    print("usage: ", sys.argv[0], " filename xcutmin xcutmax")
    exit(1)

fp = open(filenname)

x_data = []
y_data = []

ymax = float("-inf")
xmax = xcutmin 

for line in fp:
    sline = line.split()

    x = np.float64(sline[0])
    y = np.float64(sline[1])

    y = y * -2.0/3.0/math.pi/k

    if x <= xcutmax and x >= xcutmin:
        x_data.append(x)
        y_data.append(y)

        if (y> ymax):
            ymax = y
            xmax = x

fp.close()

X = np.asarray(x_data)
Y = np.asarray(y_data)
single = True
N = 1

if single:

    ymax = ymax/300.0
    xmax -= 0.0001

    p0=[ymax, xmax]

    print("Initial %.8e %.8e"%(xmax, ymax))
    print("")

    for idx in range(0,5):
        params, params_covariance = curve_fit(sfunc, x_data, y_data,
                p0)
    
        print(idx+1, " --> %.8e %.8e"%(params[0], params[1]))
        p0=params

    print("")
    print("Final value:")
    print("%.8e %.8e"%(params[0], params[1]))

    plt.scatter(x_data, y_data, label='Data', s=1)
    plt.plot(x_data, sfunc(X, params[0], params[1]),
                 label='Fitted function')
    plt.show()
else:
    sortedvalues = np.sort(y_data)

    R0 = np.asarray(np.random.rand(N))
    S0 = np.asarray(np.random.rand(N))
    
    for i in range(N):
        ymax = sortedvalues[-(i+1)]
        R0[i] = ymax / 100.0
        index = np.where(y_data == ymax)
        xmax = x_data[index[0][0]]
        S0[i] = x_data[index[0][0]]-0.00001
        print(i, " S0 %.8e "%S0[i], " selected ")
        print(i, " R0 %.8e "%R0[i], " selected ")
    
    print(" ")

    for idx in range(0,5):
 
       popt,pcov = curve_fit(
               lambda x,*params: func(x,params[:N],params[N:]),
                     X,Y, np.r_[R0,S0],
               )
       
       R = popt[:N]
       S = popt[N:]
       
       # fit using leastsq
       popt,ier = leastsq(
               lambda params: func(X,params[:N],params[N:]) - Y,
               np.r_[R0,S0],
               )
       
       R = popt[:N]
       S = popt[N:]
       
       for i in range(N):
           R0[i] = R[i]
           S0[i] = S[i]

       print (idx+1, " --> %.8e %.8e"%(R, S))


    
    plt.scatter(x_data, y_data, label='Data', s=3)
    plt.plot(x_data, func(X, R, S),
                     label='Fitted function')
    
    plt.show()
