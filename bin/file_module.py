# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 10:01:25 2014

@author: petr
"""
import numpy as np

def read_data(file):
    
    with open(file) as f:
        # read first the dimension of the grid
        d = np.fromfile(f, np.int32, 1)
        
        # read grid size
        gr_dim = np.fromfile(f, np.int32, 3)
    
        # read resolution
        dx = np.fromfile(f, np.float64, 3)
    
        # bottom left corner
        x_lo = np.fromfile(f, np.float64, 3)
    
        # upper right corner
        x_hi = np.fromfile(f, np.float64, 3)
    
        nCells = np.prod(gr_dim)
    
        
    
        # read the data
        data = np.fromfile(f, np.float64)
    


    grid = {'dim' : d, 'size' : gr_dim, 'x_lo' : x_lo, 'x_hi' : x_hi}

    data = data.reshape(gr_dim, order='F')
    data = data.transpose((1,0,2))
        
    
    return grid, data