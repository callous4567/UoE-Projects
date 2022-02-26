import itertools
import multiprocessing
import numpy as np
import windows_multiprocess_functions
import time
from SimAndVis.C2 import twod_gol

# Just to keep track of how long this sim takes.
time_start = time.time()

# Number of runs
nrun = 10000
zipped = np.arange(0, nrun, 1)

results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocess_functions.twod_run, zipped)
    pool.close()

if len(results) == len(zipped):
    time_end = time.time()
    print("Total Multirun took ",time_end - time_start)
    uwu = twod_gol.twod_gol().multigraph_absorbing()
