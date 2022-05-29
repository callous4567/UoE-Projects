# Set the variables you desire
import time

import numpy as np
from astropy.stats import sigma_clipped_stats
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from SimAndVis.Exam import exam

# Parameters for the simulation
dx = 1
dt = 0.1
nx = 50
D = 1
q = 1
p = 0.5
max_sweeps = int(1000/dt)
varray = [dx, dt, nx, D, q, p, max_sweeps]

# Number of runs to average for the absorbing state
nruns = 200

# Run it over (and save.)
times = []
time_start = time.time()
for i in range(nruns):
    print("Running ", i)
    model = exam.model(*varray)
    TT = model.run_until_absorption()
    if TT != False:
        times.append(TT)
time_end = time.time()
total_time = time_end - time_start
print("Took ", total_time)
times=np.array(times)
print(times)

# Mean and std
mean, median, std = sigma_clipped_stats(times)

# Generate a residual chart
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
axs.hist(times, bins=50, density=True)
axs.axvspan(mean-std, mean+std, color='lime', alpha=0.25)
axs.axvline(mean, 0, 1, color='green')
axs.set(xlabel="T", ylabel=r'$\rho$')
plt.title(("{0}-clipped mean of {1:.3f} with {2} of {3:.3f}").format(r'$3\sigma$', mean, r'$\sigma$', std))
plt.savefig("part_c.png", dpi=300)
plt.show()
