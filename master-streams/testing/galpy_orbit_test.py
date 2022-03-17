import astropy.units

from galpy.orbit import Orbit
import numpy as np
import ascii_info
import hdfutils
import windows_directories
from energistics import orbigistics
import astropy.units as u
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
                                                                  data_to_fit['dist']/orbigist.rovo[0],\
                                                                  data_to_fit['dmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                                  data_to_fit['dmu_b'],\
                                                                  data_to_fit['vlos']


# Also set up the errors. Don't forget: edmu is multiplied by cos(dec) in Galpy.
evxvv = np.empty(shape=(len(data_to_fit), 6))
evxvv[:,0], evxvv[:,1], evxvv[:,2], evxvv[:,3],evxvv[:,4], evxvv[:,5] = np.zeros(len(data_to_fit)),\
                                                                        np.zeros(len(data_to_fit)),\
                                                                        data_to_fit['edist']/orbigist.rovo[0],\
                                                                        data_to_fit['edmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                                        data_to_fit['edmu_b'],\
                                                                        data_to_fit['evlost']

# Set up the initial orbit guess
init_vxvv = np.mean(vxvv, axis=0)

# Get LMC from actual catalogue
o= Orbit.from_name('LMC')

# Define it from Vizier manually
l, b, dist, dmul, dmub, vlos = 280.4652, -32.8884, 49.97, 1.910, 0.229*np.cos(np.radians(-32.8884)), 262.2
newo = Orbit([l, b, dist, dmul, dmub, vlos], ro=orbigist.rovo[0]*u.kpc, vo=orbigist.rovo[1]*u.km/u.s, zo=orbigist.zo*u.kpc, lb=True)
print(o.ll(), o.bb(), o.R())
print(newo.ll(), newo.bb(), newo.R())