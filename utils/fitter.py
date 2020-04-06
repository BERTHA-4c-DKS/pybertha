import sys
import math
import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit, leastsq

T = 1000.0

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

if len(sys.argv) == 2:
    filenname = sys.argv[1]
else:
    exit(1)

fp = open(filenname)

x_data = []
y_data = []

xcutmin = 0.2
xcutmax = 0.40

for line in fp:
    sline = line.split()

    x = np.float64(sline[0])
    y = np.float64(sline[1])

    if x <= xcutmax and x >= xcutmin:
        x_data.append(x)
        y_data.append(y)

fp.close()

X = np.asarray(x_data)
Y = np.asarray(y_data)
single = False

if single:
    params, params_covariance = curve_fit(sfunc, x_data, y_data,
        p0=[0.01, 0.26])
    
    print(params)
    plt.scatter(x_data, y_data, label='Data', s=1)
    plt.plot(x_data, sfunc(X, params[0], params[1]),
                 label='Fitted function')
    plt.show()
else:
    N = 2
    
    R0 = np.asarray(np.random.rand(N))
    S0 = np.asarray(np.random.rand(N))

    S0[0] = 0.26
    S0[1] = 0.33

    #R0[0] = -7.0e-06
    #R0[1] = -7.0e-06
    
    print(R0, S0)
    
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
    
    print (R, S)
    
    plt.scatter(x_data, y_data, label='Data', s=3)
    plt.plot(x_data, func(X, R, S),
                     label='Fitted function')
    
    plt.show()
