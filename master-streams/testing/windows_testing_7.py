import datetime

import numpy as np
import pygmmis
from numpy.random import RandomState

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
rng=np.random.default_rng()

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

# Run HDBSCAN on dataset
