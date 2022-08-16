import multiprocessing
import os
import time
import numpy as np
from astropy.table import Table, vstack
import ascii_info
import asciiutils
import galcentricutils
import hdfutils
import windows_directories
import windows_multiprocessing
import windows_stack
from windows_stack import dataframes_to_tables

"""
Multiprocessed Monte-Carlo. Also sorts out Angular Momenta and 4D Clustering parameters we wanted. 
Multiprocess the monte-carlo thing.
The process needs to be ENTIRELY SELF CONTAINED for EACH TABLE: no calling externals.
Does covmonte for all data. 
"""


# NOTE THAT WE'RE DEALING WITH COVARIANCE MATRICES! THE RESULTS FROM HEREON ARE IN PANDAS DATAFRAMES!!!
# Grab groups/sets necessary
groups = ascii_info.all_groups
astropysets = [ascii_info.set_raw for group in groups]
pandasets = [ascii_info.panda_raw for group in groups]
zipped = list(zip(groups, astropysets, pandasets))
# Number of datapoints to generate errors over (higher=slower.) Final run has 500.
windows_multiprocessing.n_monte = 200
windows_multiprocessing.sourcecoord = "solar_info.dat"
results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    results = pool.map(windows_multiprocessing.do_covmonte_table, zipped)
    pool.close()

# Once we've covered all sets, this will be true.
if len(results) > 2:
    # quick cleanup of all dataframes into tables (removing covtrix/vec_L nparray columns)
    dataframes_to_tables()
    # stack 'em up :)
    windows_stack.stacking().tables()
    windows_stack.stacking().dataframes()

