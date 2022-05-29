import time

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from scipy.fft import fft, fftfreq
from scipy.optimize import curve_fit

from SimAndVis.Exam import exam

dx = 1
dt = 0.01
nx = 50
D = 0.5
q = 1
p = 2.5
max_sweeps = int(1000/dt)
viewing_binrate = int(10)
varray = [dx, dt, nx, D, q, p, max_sweeps, viewing_binrate]

# Run
uwu = exam.model(*varray)
times = np.arange(0, uwu.max_sweeps, 1)*dt
types_1, types_2 = uwu.a_over_time_twopoints([10,10], [20,20])
types = [types_1, types_2]
types = [np.array(d) for d in types]
types = [d[int(max_sweeps/5):] for d in types]

# Plots
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
colors=['blue', 'red']
labels = ["10,10", "20,20"]
for num,type in enumerate(types):

    plt.plot(times[int(max_sweeps/5):], type - np.mean(type), color=colors[num], label=labels[num])
axs.grid(which='major', color='pink')
axs.legend(loc='upper right')
axs.set(xlabel="T", ylabel="a")
plt.savefig("fraction_time_q5.png", dpi=300)
plt.show()
plt.close()

# Fit sines to both sets of data.
"""
general sine
asin(bt+c) + d 
"""

type = types[0]

# sample spacing
T = dt
x = times[int(max_sweeps/5):][0:80000]
y = type[0:80000]
freq_mags = fft(y)
frequencies = fftfreq(len(x)) * 1/dt

plt.stem(frequencies, np.abs(freq_mags))
plt.xlim([-0.25,0.25])
plt.grid()
plt.show()