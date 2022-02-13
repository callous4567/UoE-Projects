import multiprocessing
import numpy as np
from astropy.table import Table
import os
import hdfutils
import windows_xy_functions
import time

from SimAndVis.C1_2.xy_ising import xy_ising

time_start = time.time()
temperatures = np.linspace(0.5, 2, 30)

# This is to generate all the simulation runs for each temperature
results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(10)
    results = pool.map(windows_xy_functions.model_run, temperatures)
    pool.close()

if len(results) == len(temperatures):
    time_end = time.time()
    print(time_end - time_start)
    uwu = xy_ising()
    uwu.multigraph()