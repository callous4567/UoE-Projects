import exam

# Set the variables you desire
dx = 1
dt = 0.1
nx = 50
D = 1
q = 1
p = 0.5
max_sweeps = int(1e4)
viewing_binrate = int(10)
varray = [dx, dt, nx, D, q, p, max_sweeps, viewing_binrate]

# Run
uwu = exam.model(*varray)
uwu.run()