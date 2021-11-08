import datetime
import pickle

import numpy as np
from numpy.random import RandomState
import numpy as np
import graphutils
import hdfutils
import windows_directories, ascii_info
from functools import partial


# Get data
data = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                          ascii_info.fullpanda)

# Set up L via vec_L
vec_L = data['vec_L'].tolist()
vec_L = np.array(vec_L)

# Get covariance matrices
cov_L = data['covtrix'].tolist()
cov_L = np.array(cov_L)


# Define rng
np.random.seed(11) # make sure to set the seed.
rng=np.random.default_rng()



"""
# Return probabilistic selection for bounding box as True or False (a generator object.)
def getCuboid(L_data):
    # Specify the bounds for the cuboid
    bounds = np.array([[-8e3,6e3],
                       [-10e3,10e3],
                       [-7e3,7e3]])
    coords = L_data
    box_limits = bounds
    # Get probs: 1 or 0.
    return (coords[:,0] > box_limits[0,0]) & (coords[:,0] < box_limits[0,1]) \
         & (coords[:,1] > box_limits[1,0]) & (coords[:,1] < box_limits[1,1]) \
         & (coords[:,2] > box_limits[2,0]) & (coords[:,2] < box_limits[2,1])

# Apply on data
cov_L = cov_L[getCuboid(vec_L)]
vec_L = vec_L[getCuboid(vec_L)]




mus = np.array([[426, -4950, 1436],
                [-0.7, 11, 3],
                [144, -278, 122],
                [-3905, -2322, -4664]])
covtrix_1 = [[655**2, 0, 0],
             [0, 1255**2, 0],
             [0, 0, 659**2]]
covtrix_2 = [[11**2, 0, 0],
             [0, 12**2, 0],
             [0, 0, 11**2]]
covtrix_3 = [[1755**2, 0, 0],
             [0, 1926**2, 0],
             [0, 0, 1733**2]]
covtrix_4 = [[13**2, 0, 0],
             [0, 12**2, 0],
             [0, 0, 13**2]]
covtrices = np.array([covtrix_1, covtrix_2, covtrix_3, covtrix_4])


graphutils.threed_graph().kmeans_L_multinormal(vec_L, mus, covtrices, 1000)


"""