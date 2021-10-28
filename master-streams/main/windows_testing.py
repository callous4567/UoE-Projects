import math
import multiprocessing
import numpy as np
from astropy.table import Table
import time
import ascii_info
import galcentricutils
import graphutils
import hdfutils
import windows_directories
import windows_multiprocessing


# Set up coordinate grid
n_points = 1.2e4
coordgrid = galcentricutils.greatcount().fullgridgen(n_points) # IN RADIANS!
#graphutils.threed_graph().unitsphere(coordgrid[0], coordgrid[1], coordgrid[2])

# Specify dtheta/radmin for gcc cut: n_phi too for number of phi regions. dtheta is half-width of cell.
dtheta = 2.5 #in degrees
radmin = 20 # in kpc
n_phi = 24 # number of phi regions

# Calculate the probability to get a single star within the segment- see notes.
P = (2*(2*np.pi/n_phi)*np.sin(math.radians(dtheta))) / (4*np.pi)

# Specify the table you're working this on: for now full table. Check later.
tableinfo = [ascii_info.fullgroup, ascii_info.fullset]

# Create the variable grid for pool
mapped_variables = []
for theta, phi, index in zip(*coordgrid):
    theta, phi = math.degrees(theta), math.degrees(phi)
    gccpar = [theta, phi, index, dtheta, radmin, n_phi]
    savedex = ("{0:.2f}_{1:.2f}").format(theta,phi)
    full_variable = [tableinfo, gccpar, P, savedex]
    mapped_variables.append(full_variable)

# Run the function...
def main():
    if __name__ == '__main__':
        pool = multiprocessing.Pool(4)
        result = pool.map(func=windows_multiprocessing.do_gccsplit_confidence, iterable=mapped_variables)
        pool.close()
        return result

results = main()
if type(results) == list:
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    results = np.array(results).T
    thetas, phis, indices, sigmins = results
    table = Table()
    table['theta'],table['phi'],table['index'],table['sigmin'] = thetas, phis, indices, sigmins
    writer.write_table(tableinfo[0], "GCC_TEST_CONFIDENCE", table)






