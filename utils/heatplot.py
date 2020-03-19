import numpy as np

import basicutils

if __name__ == "__main__":

    occups = np.array([1,2,3,4])
    virtua = np.array([12,11,10,13])
    values = np.array([[0.0,  1.0, 0.0, 0.0],\
                       [0.0, 10.0, 0.0, 0.0],\
                       [0.0,  0.0, 8.0, 0.0],\
                       [0.0,  3.0, 0.0, 5.6]])

    xlabel = "Occupied"
    ylabel = "Virtual"
    title =  "Analysis of ..."

    outfname = "test.jpg"
    norm = 10.0
    
    basicutils.exportheatmap(occups, virtua, values, xlabel, \
            ylabel, title, outfname, norm)

