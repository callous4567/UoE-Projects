import multiprocessing
import os
import time

from astropy.table import Table, vstack

import ascii_info
import asciiutils
import galcentricutils
import hdfutils
import windows_directories
import windows_multiprocessing

"""
Multiprocessed Monte-Carlo. Also sorts out Angular Momenta. 
Multiprocess the monte-carlo thing.
The process needs to be ENTIRELY SELF CONTAINED for EACH TABLE: no calling externals.
"""

# NOTE THAT WE'RE DEALING WITH COVARIANCE MATRICES! THE RESULTS FROM HEREON ARE IN PANDAS DATAFRAMES!!!

# Grab groups/sets necessary
groups = ascii_info.all_groups
sets = [ascii_info.set_raw for group in groups]
savesets = [ascii_info.panda_raw for group in groups]
zipped = list(zip(groups, sets, savesets))
windows_multiprocessing.n_monte = 200
windows_multiprocessing.sourcecoord = "solar_info.dat"
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(4)
    pool.map(windows_multiprocessing.do_covmonte_table, zipped)


