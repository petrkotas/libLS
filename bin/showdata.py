import numpy as np
from mayavi import mlab

if __name__ == "__main__":
    with open("sphere_100.data") as f:
        # read first the dimension of the grid
        d = np.fromfile(f, np.int32, 1)
        print d
        # read grid size
        gr_dim = np.fromfile(f, np.int32, 3)
        print gr_dim
        # read resolution
        dx = np.fromfile(f, np.float64, 3)
        print dx
        # bottom left corner
        x_lo = np.fromfile(f, np.float64, 3)

        # upper right corner
        x_hi = np.fromfile(f, np.float64, 3)

        nCells = np.prod(gr_dim)
        print nCells
        # read the data
        data = np.fromfile(f, np.float64)
        print data.size
        print data
        #data = data.reshape(gr_dim, order='F')
        #data = data.transpose((1,0,2))

        # show the data
        #mlab.contour3d(data, contours=[0,]) 
        #mlab.show()
        
        
        
        
