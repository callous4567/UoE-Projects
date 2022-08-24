import multiprocessing
import numpy as np

import ascii_info_new
import windows_multiprocessing_new

# Generate aics/bics for a bunch of clusterings.

# Set up coordinate grid
ntheta = 20
nphi = 20
thetas = np.linspace(0,180,ntheta)
phis = np.linspace(0,360,nphi)

# Specify dtheta/radmin for gcc cut
dtheta = 5 #in degrees
radmin = 20 # in kpc

# The maximum number of clusters for bic
k_max = 10

# Specify the table you're working this on: for now full table. Check later.
tableinfo = [ascii_info_new.fullgroup, ascii_info_new.fullset]

# Create the variable grid for pool
mapped_variables = []
for theta in thetas:
    for phi in phis:
        gccpar = [theta, phi, dtheta, radmin]
        savedex = ("{0:.2f}_{1:.2f}").format(theta,phi)
        clustpar = k_max, savedex
        full_variable = [tableinfo, gccpar, clustpar]
        mapped_variables.append(full_variable)

"""
tableinfo, gcc_params, clustering_params
table, theta, phi, dtheta, radmin
table, k_max, savedex):
"""
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    pool.map(windows_multiprocessing_new.do_gaussaicsbics, mapped_variables)