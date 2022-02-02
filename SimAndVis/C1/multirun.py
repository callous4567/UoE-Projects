import multiprocessing
import numpy as np
from astropy.table import Table
import os
import hdfutils
import windows_multiprocess_functions


# Set range for simulation
#temperatures = np.linspace(0.5, 5, 100)
#dynamics = ['k' for d in temperatures]

# Second run for Glauber! Given the speed-up from removing final-point energy-calculations :)
# Set range for simulation
temperatures = np.linspace(2.5, 5, 1000)
dynamics = ['g' for d in temperatures]

zipped = list(zip(temperatures, dynamics))
results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocess_functions.twod_run, zipped)
    pool.close()

"""
# Once we've covered all sets, this will be true.
# quick fix: see pickling/autosaving. This threw an error due to writing to disk and I cba dealing with that. 
# HDF parallel write did me in. 
if len(results) > 2:
    # Save the temperatures we've simulated under a group "dynamic" and dataset "temperatures"
    labels = ['avg_M', 'avg_E', 'avg_M_err', 'avg_E_err', 'chi_true', 'chi_error', 'c_true', 'c_error', 'T']
    results = list(np.array(results).T)
    results.append(temperatures)
    results = np.array(results).T
    table = Table(results, names=labels)
    writer = hdfutils.hdf5_writer(os.getcwd(), "multirun_data.hdf5")
    writer.write_table(group=dynamics[0], dataset="multirun_data", astropytable=table) """