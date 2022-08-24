import os

import numpy as np
from astropy.io import fits
from astropy.table import Table

import windows_directories_new


import graphutils_new
from energistics_new import fast_energistics_new
import hdfutils
import matplotlib.pyplot as plt
import graphutils_new

finalstack_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_final.fits")

import graphutils_new
from energistics_new import fast_energistics_new
import hdfutils
import matplotlib.pyplot as plt
import graphutils_new

with fits.open(finalstack_dir, mode='readonly', memmap=True) as hdul:
    # Data (which is an astropy table when loaded.)
    # "gaia_source_id" is dr2_source_id inside this.
    table = Table(hdul[1].data)
    table = table[[True if source == "LAMOST_sstraszak" else False for source in table['source']]]
    table = table[[True if feh < -0.5 else False for feh in table['feh']]]


# data = data[[True if 1 - abs(circ) > 0.9 else False for circ in data['circ']]]

# data = data[[True if abs(d) < 300 else False for d in data['vx']]]
# data = data[[True if abs(d) < 300 else False for d in data['vy']]]
# data = data[[True if abs(d) < 300 else False for d in data['vz']]]

# data = fast_energistics_new().default_E_c(data)

table['x'], table['y'], table['z'], \
table['vx'], table['vy'], table['vz'] = np.array(table['x'], float), \
                                     np.array(table['y'], float), \
                                     np.array(table['z'], float), \
                                     np.array(table['vx'], float), \
                                     np.array(table['vy'], float), \
                                     np.array(table['vz'], float)

# plt.scatter(master_fits['Lz'], master_fits['E'], color='blue', s=1)
# plt.scatter(table['Lz'], table['E'], color='red', s=0.5)

# plt.hist(table['vlos'], bins=100)
# plt.hist(data['vlos'], bins=100)
# plt.show()
datatab = np.array([table['Lx'], table['Ly'], table['Lz']]).T
# clustered = cluster3d().listhdbs(data, [1000, 20])
graphutils_new.threed_graph().kmeans_L_array(datatab, [1 for d in datatab[:, 0]], False, browser=True, outliers=True)

# writer = hdfutils.hdf5_writer(windows_directories_new.datadir, "stardata.hdf5")
# data = writer.read_table("LAMOST_K_FULL_edr3", "astrotable")
# data = data[[True if r > 15 else False for r in data['r']]]

plt.hist(table['vx'], bins=1000, label='vx')
plt.hist(table['vy'], bins=1000, label='vy')
plt.hist(table['vz'], bins=1000, label='vz')
plt.legend()
plt.show()
clustered = [1 for d in table['source']]
data = fast_energistics_new().default_E_c(table)
graphutils_new.spec_graph().energy_plot(data,
                                    clustered,
                                    list(set(clustered)),
                                    list(set(clustered)),
                                    [-10000, 10000],
                                    [-200000, 10000],
                                    "lamost_test_finale.png",
                                    False,
                                    False)