import multiprocessing
import numpy as np
from astropy.table import Table
import os
import hdfutils
import windows_multiprocess_functions
import time

time_start = time.time()
dynamics = 'g'
if dynamics == 'g':
    temperatures = np.linspace(1, 3, 60)
    dynamics = [dynamics for T in temperatures]
if dynamics == 'k':
    temperatures = np.linspace(2, 3, 30)
    #temps_2 = np.array([10**x for x in temps_2])
    #temperatures = np.linspace(0.1, 2, 60)
    #temperatures = np.concatenate([temps_1, temps_2]) # , temps_2, temps_4, temps_5, temps_6])
    dynamics = [dynamics for T in temperatures]
zipped = list(zip(temperatures, dynamics))

# This is to generate all the simulation runs for each temperature

results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(10)
    results = pool.map(windows_multiprocess_functions.twod_regenerate_averages, zipped)
    pool.close()

if len(results) == len(temperatures):
    time_end = time.time()
    print(time_end - time_start)