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

# Note that this tool is deprecated and not useful for our dataset, at least not that useful.
# Only Sagittarius is discoverable as a GCC (since other targets are hazed out by the data/noise/etc.)
"""
This tool will generate confidence maps for the sky, following the works of 
IDENTIFYING STAR STREAMS IN THE MILKY WAY HALO
Charles King III et al.
dtheta is the half-width of the cell (different from kings paper.)  
nphi is the number of phi-regions per GCC.
This does not account for the bonferroni correction- to do this, see original paper.
"""

# Set up the grid of nphi/dtheta just as they were in windows_testing
dthetarange = np.arange(0.25, 0.76, 0.25)
nphiarange = np.arange(13, 30, 4)

# Iterate for a range of dtheta/n_phi to generate confidence maps
for dtheta in dthetarange:
    for n_phi in nphiarange:
        # Set up coordinate grid
        n_points = 2.4e4
        coordgrid = galcentricutils.greatcount().fullgridgen(n_points) # IN RADIANS!
        #graphutils.threed_graph().unitsphere(coordgrid[0], coordgrid[1], coordgrid[2])

        # Specify dtheta/radmin for gcc cut: n_phi too for number of phi regions. dtheta is half-width of cell.
        #dtheta = 2.5 #in degrees
        radmin = 20 # in kpc
        #n_phi = 1 # number of phi regions

        # For saving within the original group
        string_format = ("{0:.2f}_{1:.0f}").format(dtheta, n_phi)

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
                pool = multiprocessing.Pool(8)
                result = pool.map(func=windows_multiprocessing.do_gccsplit_confidence, iterable=mapped_variables)
                pool.close()
                return result

        # Get results and write the "confidence map" to a table
        results = main()
        if type(results) == list:
            writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
            results = np.array(results).T
            thetas, phis, indices, nsplitsigmins, sigmins = results
            table = Table()
            table['theta'],table['phi'], table['index'], table['phiregion_n'], table['sigmin'] = thetas, phis, \
                                                                                                 indices, \
                                                                                                 nsplitsigmins, \
                                                                                                 sigmins
            while True:
                try:
                    print("mapped and now writing, ", dtheta, n_phi)
                    writer.write_table(tableinfo[0], string_format, table)
                    break
                except:
                    continue


