

# Set the variables you desire
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from SimAndVis.Exam import exam

dx = 1
dt = 0.01
nx = 50
D = 1
q = 1
p = 0.5
max_sweeps = int(20000)
times = np.arange(0, max_sweeps, 1)*dt
varray = [dx, dt, nx, D, q, p, max_sweeps]
model = exam.model(*varray)
type_over_time = model.type_over_time()
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
fracs = ["a", "b", "c"]
colours = ['r', 'g', 'b']
patches = [Patch(facecolor=c, edgecolor='black', label=l) for c,l in zip(colours,fracs)]
for num,type in enumerate(type_over_time):
    axs.plot(times, type, color=colours[num])
axs.grid(which='major', color='pink')
axs.set(xlabel="T", ylabel="Fraction")
axs.legend(handles=patches, loc='upper right')
plt.savefig("fraction_over_time_Q2.png", dpi=300)
plt.show()