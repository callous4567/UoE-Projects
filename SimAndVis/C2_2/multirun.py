import itertools
import multiprocessing
import numpy as np
import windows_multiprocess_functions
import time

"""
Question 2: 

# Generate the parameter space for the runs (for probability.)
one_range = np.linspace(0, 1, 10)
two_range = np.linspace(0, 1, 10)
tre_range = np.linspace(0, 1, 10)
zipped = list(itertools.product(one_range, two_range, tre_range))

With GRAPH = TRUE! (Generate lots and lots of them with graphs to examine.

"""

"""
Question 3: 

# p2 is fixed, p1-p3 plane. 

# Generate the parameter space for the runs (for probability.)
one_range = np.linspace(0, 1, 20)
two_range = np.array([0.5]) # set p2 = 0.5 
tre_range = np.linspace(0, 1, 20)
zipped = list(itertools.product(one_range, two_range, tre_range))
Set graph to False.  

"""

"""
Question 4: 

# 1D Slice with varied p1: p2 and p3 set to 0.5. 

# Generate the parameter space for the runs (for probability.)
one_range = np.linspace(0.2, 0.5, 20)
two_range = np.array([0.5]) # set p2 = 0.5
tre_range = np.array([0.5]) # set p3 = 0.5, too
zipped = list(itertools.product(one_range, two_range, tre_range))

"""

# Just to keep track of how long this sim takes.
time_start = time.time()

# Generate the parameter space for the runs (for probability.)
one_range = np.array([0.5])
two_range = np.array([0.5]) # set p2 = 0.5
tre_range = np.array([0.5]) # set p3 = 0.5, too
immu_range = np.linspace(0, 1, 20)
zipped = list(itertools.product(one_range, two_range, tre_range, immu_range))

results = ["null"]
# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(10)
    results = pool.map(windows_multiprocess_functions.twod_run, zipped)
    pool.close()

if len(results) == len(zipped):
    time_end = time.time()
    print(time_end - time_start)