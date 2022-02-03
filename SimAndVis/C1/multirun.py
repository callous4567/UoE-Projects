import multiprocessing
import numpy as np
from astropy.table import Table
import os
import hdfutils
import windows_multiprocess_functions
import time

dynamics = 'g'
if dynamics == 'g':
    temps_1 = np.linspace(1, 2, 10)
    temps_2 = np.linspace(2, 2.2, 20)
    temps_4 = np.linspace(2.2, 2.4, 40)
    temps_5 = np.linspace(2.4, 2.5, 10)
    temps_6 = np.linspace(2.5, 3.5, 10)
    temperatures = np.concatenate([temps_1, temps_2, temps_4, temps_5, temps_6])
    dynamics = [dynamics for T in temperatures]
if dynamics == 'k':
    temps_1 = np.linspace(0.0001, 0.1, 30)
    temps_2 = np.linspace(0.1, 2, 5)
    temps_4 = np.linspace(2.2, 2.4, 10)
    temps_5 = np.linspace(2.4, 2.5, 10)
    temps_6 = np.linspace(2.5, 3.5, 10)
    temperatures = np.concatenate([temps_1, temps_2, temps_4, temps_5, temps_6])
    dynamics = [dynamics for T in temperatures]
zipped = list(zip(temperatures, dynamics))

# This is to generate all the simulation runs for each temperature

results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocess_functions.twod_run, zipped)
    pool.close()

# This is to regenerate all the averages
#for tempdyn in zipped:
#    windows_multiprocess_functions.twod_regenerate_averages(tempdyn)
