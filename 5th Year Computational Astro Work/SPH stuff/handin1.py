# Handin 1, Sebastian Straszak (S1728659), Part 3 :D
####################################################
import time

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from numpy import random, ravel

# Define all various parameters and get the sample.
"""
In this case I'm using a fixed smoothing length alongside equal masses
The mass distribution is a Gaussian in regular cartesian space, drawn from a multivariate normal
"""
nrows = 50 # number of samples
mass = 0.1 #example mass
smoothing_length = 0.5
mean, stdev = 0, 3
total_mass = nrows*mass
particles = random.default_rng().normal(loc=mean, scale=stdev, size=(nrows, 3))

# Example plot of points.
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(*particles.T, color='black',s=1)
ax.set(xlabel="x",
       ylabel="y",
       zlabel="z")
plt.show()

# Define the smoothing kernel, Wab(h). Using Gaussian for simplicity: 3D.
"""
W(a,b,h) = $\frac{1}{h^3\pi^\frac{3}{2}}\exp{\left(-\frac{|\vec{r_a}-\vec{r_b}|^2}{h^2}\right)}$
Taken from Sadegh's notes- I'm assuming it's all fine and normalized/etc in this form.
"""
def kernel(posa, posb, h):
    return (1/((h**3)*(np.pi**(3/2))))*np.exp(-(np.linalg.norm(posa-posb)**2)/h**2)

# Calculate the density at some point [r]
"""
Standard form for a variable A given a fixed smoothing length:
A_i = \displaystyle \sum_{j=1}^N \frac{A_jm_j}{\rho_j}W_ij(h)
"""
def density(x, y, z,smoothing_length=smoothing_length):
    # Create vector for the position you're calculating density for (useful for map)
    r = np.array([x,y,z])

    rho = float(0)
    for particle in particles:
        rho += mass*kernel(r, particle, smoothing_length)
    return rho

"""
Get the total mass: untested. 
"""
# Build domain in each axis
oned_domain = np.array([mean-(3*stdev), mean+(3*stdev)])
oned_points = np.linspace(*oned_domain, nrows)
volume = ((max(oned_domain)-min(oned_domain))/(nrows-1))**3
# Define grid
x,y,z = oned_points,oned_points,oned_points
# Find density over grid elements
densities = density(x[:,None,None], y[None,:,None], z[None,None,:]) # https://stackoverflow.com/questions/22774726/numpy-evaluate-function-on-a-grid-of-points
# Get mass for each element + sum for total
mass_elements = densities*volume*mass
intmass = np.sum(mass_elements)
print(intmass)

