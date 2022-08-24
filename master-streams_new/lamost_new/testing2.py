# Load in the original fits
import numpy as np
from astropy.io import fits
from astropy.table import Table

import windows_directories_new

with fits.open(windows_directories_new.lamodir + "\\dr7_2MASS_aidists.fits", mode='readonly', memmap=True) as hdul:
    # Data (which is an astropy table when loaded.)
    # "gaia_source_id" is dr2_source_id inside this.
    data = Table(hdul[1].data)

data.rename_columns(['rv','erv', 'd_kpc_', 'e_d_kpc_'],['vlos','evlost','dist','edist'])


# Get the galactic coordinates from these equatorial ones
from energistics_new import orbigistics
from galcentricutils_new import angular
orbigist = orbigistics()
data = orbigist.converter.nowrite_ICRS_to_GAL(data, has_cosfactor=True)
data = orbigist.converter.nowrite_GAL_to_GALCENT(data)
data = angular().get_momentum(data)

# Remove nuisance stars (outside min_radius, but within the max_radius of interest.)
min_radius = 12
max_radius = 100
data = data[[True if max_radius > r > min_radius else False for r in data['r']]]

# Grab the table of globular clusters
with fits.open(windows_directories_new.baumgardt_fits, mode='readonly', memmap=False) as hdul:
    gc_table = Table(hdul[1].data)

# Add dummy values for proper motions
gc_table.add_column(0., name='pmra')
gc_table.add_column(0., name='pmdec')
gc_table.add_column(0., name='pmra_error')
gc_table.add_column(0., name='pmdec_error')
gc_table.add_column(0., name='vlos')
gc_table.add_column(0., name='evlost')

# Get Galcent
gc_table = orbigist.converter.nowrite_ICRS_to_GAL(gc_table, has_cosfactor=True)
gc_table = orbigist.converter.nowrite_GAL_to_GALCENT(gc_table)

# Convert tidal radius to kpc
gc_table['rt']/=1000

# Remove stars within gcs
import gcc_utils
gc_to_remove = gcc_utils.remove_gcs(np.array(data['x'], float), np.array(data['y'], float), np.array(data['z'], float),
                                    np.array(gc_table['x'], float), np.array(gc_table['y'], float), np.array(gc_table['z'], float),
                                    np.array(gc_table['rt'], float))
data = data[gc_to_remove]


from galcentricutils_new import cluster3d
import graphutils_new
from energistics_new import fast_energistics_new

data['x'], data['y'], data['z'], \
data['vx'], data['vy'], data['vz'] = np.array(data['x'], float), np.array(data['y'], float), np.array(data['z'], float), \
                                                          np.array(data['vx'], float), np.array(data['vy'], float), np.array(data['vz'], float)
master_fits = fast_energistics_new().default_E_c(data)
table = fast_energistics_new().default_E_c(data)
import matplotlib.pyplot as plt
#plt.scatter(master_fits['Lz'], master_fits['E'], color='blue', s=1)
#plt.scatter(table['Lz'], table['E'], color='red', s=0.5)

#plt.hist(table['vlos'], bins=100)
plt.hist(master_fits['vlos'], bins=100)
plt.show()
data = np.array([master_fits['Lx'], master_fits['Ly'], master_fits['Lz']]).T
clustered = cluster3d().listhdbs(data, [1000, 20])
data = np.array(data)
graphutils_new.threed_graph().kmeans_L_array(data, clustered, False, browser=True, outliers=True)





