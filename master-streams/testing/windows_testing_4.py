import math
import numpy as np
import scipy
import ascii_info
import galcentricutils
import matplotlib.pyplot as plt
import mpsci
from scipy.stats import poisson


# Set up coordinate grid
n_points = 2.4e2
coordgrid = galcentricutils.greatcount().fullgridgen(n_points) # IN RADIANS!
#graphutils.threed_graph().unitsphere(coordgrid[0], coordgrid[1], coordgrid[2])

# Specify dtheta/radmin for gcc cut: n_phi too for number of phi regions. dtheta is half-width of cell.
dtheta = 2.5 #in degrees
radmin = 20 # in kpc
n_phi = 14 # number of phi regions

# Calculate the probability to get a single star within the segment- see notes.
P = (2*(2*np.pi/n_phi)*np.sin(math.radians(dtheta))) / (4*np.pi)
radnumber, n, sig = 10864, 100, 0.0
n_range = np.arange(0, 100, 1)
cdf_range = [mpsci.distributions.binomial.cdf(n, radnumber, P, method='sumpmf') for n in n_range]
print(cdf_range)
plt.scatter(n_range, cdf_range)
plt.xlim(0,100)
plt.show()