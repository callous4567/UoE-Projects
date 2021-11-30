import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import G
import astropy.constants.iau2015 as iau
import astropy.units as u

# Define properties of M and R. We're going to use Earth, just because. Assume average smooth density.
m = iau.M_earth.to(u.kg).value
r = iau.R_earth.to(u.m).value

# Define r'$\phi(r)$' and r'$|g(r)|$'
phi = lambda m, r: -1*G*m/r
grav = lambda m, r: G*m/(r**2)

# Set up m1m2/r1r2
m1 = 3*m
m2 = m
r1 = ((m1/m)**(1/3)) * r
r2 = ((m2/m)**(1/3)) * r
D = 10*r

# Define x-axis in units of r
x = np.linspace(0, D, 1000)

# Get potential as function of x
phi_of_x = phi(m1, np.abs(x)) + phi(m2, np.abs(x-D))

# Also get gravity. Note that grav is not signed and this must be considered.
grav_of_x = -1*grav(m1,np.abs(x))*(x-0)/np.abs(x-0) + -1*grav(m2,np.abs(x-D))*(x-D)/np.abs(x-D)

# Plot it: we'll introduce asymptotes at r1/r2 to indicate where our approx breaks down
fig, ax = plt.subplots(nrows=2,ncols=1, sharex='all')
plt.subplots_adjust(wspace=0, hspace=0)

# Set up axes
ax[0].plot(x/r, phi_of_x)
ax[1].plot(x/r, grav_of_x)
ax[0].set(ylabel=r'$\phi(x)$',
          xlabel=r'$\frac{r}{r_{earth}}$')
ax[1].set(ylabel=r'$g(x)$', ylim=[-40,40],
          xlabel=r'$\frac{r}{r_{earth}}$')
ax[0].axvline(r1/r,ymin=0,ymax=1, linewidth=1, linestyle='dashed', color='red')
ax[0].axvline((D-r2)/r,ymin=0,ymax=1, linewidth=1, linestyle='dashed', color='red')
ax[1].axvline(r1/r,ymin=0,ymax=1, linewidth=1, linestyle='dashed', color='red')
ax[1].axvline((D-r2)/r,ymin=0,ymax=1, linewidth=1, linestyle='dashed', color='red')
plt.savefig("hello.png", dpi=300)
plt.show()

