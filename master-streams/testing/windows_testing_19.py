import pickle
import time
import astropy.units as u
from galpy.util import bovy_coords
import galpy.orbit
import numpy as np
from matplotlib import pyplot as plt
import ascii_info
import hdfutils
import windows_directories
from energistics import orbigistics

"""
Galpy Orbit Fitting using orbit.from_fit, for the data.
"""

# Load member table (of the members who survived the initial monte-carlo'ing)
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
membership_table = writer.read_table(ascii_info.fullgroup, "percent_table")

# Grab the data table that you are fitting
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)

# Decide a cluster to fit an orbit to, just
clust_to_fit = 1
clustering = membership_table['probable_clust']
clustselec = [True if d == clust_to_fit else False for d in clustering]

# Clip data
data_to_fit = table[clustselec]

# Set up orbigistics (for potential)
orbigist = orbigistics()

# Use the solarmotion='schoenrich' local solar velocities
# (as noted by Drimmel and Poggio as being decent- these are the defaults.)
solarmotion='schoenrich'

# Set up the position vectors for galpy, which will use galactic coordinates straight-up, and also use errors
vxvv = np.empty(shape=(len(data_to_fit), 6))
vxvv[:,0], vxvv[:,1], vxvv[:,2], vxvv[:,3],vxvv[:,4], vxvv[:,5] = data_to_fit['l'],\
                                                                  data_to_fit['b'],\
                                                                  data_to_fit['dist'],\
                                                                  data_to_fit['dmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                                  data_to_fit['dmu_b'],\
                                                                  data_to_fit['vlos']


# Also set up the errors. Don't forget: edmu is multiplied by cos(dec) in Galpy.
evxvv = np.empty(shape=(len(data_to_fit), 6))
evxvv[:,0], evxvv[:,1], evxvv[:,2], evxvv[:,3],evxvv[:,4], evxvv[:,5] = np.zeros(len(data_to_fit)),\
                                                                        np.zeros(len(data_to_fit)),\
                                                                        data_to_fit['edist'],\
                                                                        data_to_fit['edmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                                        data_to_fit['edmu_b'],\
                                                                        data_to_fit['evlost']

# Set up the initial orbit guess
init_vxvv = np.mean(vxvv, axis=0)

load_fit = True
if load_fit == False:
    # Run and save the fit
    fit = galpy.orbit.Orbit.from_fit(init_vxvv,
                                     vxvv,
                                     evxvv,
                                     orbigist.pot,
                                     ro=orbigist.rovo[0]*u.kpc,
                                     vo=orbigist.rovo[1]*u.km/u.s,
                                     zo=orbigist.zo*u.kpc,
                                     solarmotion=solarmotion,
                                     lb=True)
    with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt", 'wb') as f:
        pickle.dump(obj=fit, file=f)
else:
    # Load the fit
    with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt", 'rb') as f:
        fit = pickle.load(file=f)

# Get a "flip" for backward orbit
fit_backward = fit.flip()

# Get integrals
fit.integrate(np.linspace(0, 0.1e9, 1000)*u.yr, orbigist.pot, 'rk4_c')
fit_backward.integrate(np.linspace(0, 0.1e9, 1000)*u.yr, orbigist.pot, 'rk4_c')
# Get orbits
orbits = fit.getOrbit()
orbits_backward = fit_backward.getOrbit()
orbits_full = np.concatenate((orbits, orbits_backward), axis=0)
R, vR, vT, z, vz, phi = orbits_full.T  # the fitted parameters, in units of rovo/etc.
X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
#r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
plt.scatter(R/orbigist.rovo[0], l, color='red', marker='x', s=0.1)

# Clip data and get galpy cylindrical coordinates, alongside getting data as an array.
R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
#r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
plt.scatter(R, l, color='black', marker='x')
plt.savefig(windows_directories.imgdir + "\\" + str(clust_to_fit) + "_all_orbits_test.png", dpi=300)
plt.xlabel('l')
plt.ylabel('b')
plt.show()


