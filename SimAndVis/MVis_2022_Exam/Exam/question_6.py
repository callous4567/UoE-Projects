from matplotlib import pyplot as plt

from SimAndVis.Exam import exam
from SimAndVis.Exam.fast_exam import typefield, correlation_func

dx = 1
dt = 0.01
nx = 50
D = 0.3
q = 1
p = 2.5
max_sweeps = int(1000/dt)
viewing_binrate = int(10)
varray = [dx, dt, nx, D, q, p, max_sweeps, viewing_binrate]

# Run
uwu = exam.model(*varray)
while uwu.sweep <= int(max_sweeps/5):
    uwu.fast_model()

# Grab the correlation function for this
r_range, probs = correlation_func(typefield(uwu.a, uwu.b, uwu.c, uwu.d))

# Create plot
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
axs.plot(r_range, probs, color='red')
axs.grid(which='major', color='pink')
axs.set(xlabel="Euclidean Distance", ylabel="Probability (Not Normalised)")
plt.savefig("correlation_plot_0.3.png", dpi=300)
plt.show()