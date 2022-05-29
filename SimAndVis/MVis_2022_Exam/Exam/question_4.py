# Set the variables you desire
from SimAndVis.Exam import exam

dx = 1
dt = 0.3
nx = 256
D = 0.5
q = 4
p = 2
max_sweeps = int(1e4)
viewing_binrate = int(10)
varray = [dx, dt, nx, D, q, p, max_sweeps, viewing_binrate]

# Run
uwu = exam.model(*varray)
uwu.run()