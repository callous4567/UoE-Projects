import itertools
import multiprocessing
import numpy as np
import windows_multiprocess_functions
import time

# Just to keep track of how long this sim takes.
time_start = time.time()

# Generate the parameter space for the runs.
one_range = np.linspace(0, 1, 50)
two_range = np.linspace(0, 1, 50)
tre_range = np.linspace(0, 1, 50)
zipped = list(itertools.product(one_range, two_range, tre_range))

results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocess_functions.twod_run, zipped)
    pool.close()

if len(results) == len(zipped):
    time_end = time.time()
    print(time_end - time_start)