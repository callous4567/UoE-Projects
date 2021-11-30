import numpy as np
import yt
from yt.units import *
import astropy.units as u
from dm_overdensity import _overdensity
import matplotlib.pyplot as plt
import yt.units as u
from yt.analysis_modules.halo_analysis.api import *

"""
Note that results is quite hefty- don't spam too many results up...
KEEP ALL STEPS IN ORDER and for the love of god almighty, Label. Each. Step. WELL. 
"""

# This has the working directory for all (directly) related Enzo stuff.
resdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm/results"
supdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm"

# WSdir for this workshop (i.e. save images here.)
wsdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/ws2"

# TODO: Steps for Q1
"""
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(("gas", "dm_overdensity"), function=_overdensity, sampling_type="cell")
ad0 = ds0.all_data()
den0 = ad0['all_cic']
oden0 = ad0['dm_overdensity']

# in gcm^-3 CGS UNITS!!!
min = den0.min()
max = den0.max()
print("density minmax", min, max)
print(max - min)

# in gcm^-3 CGS UNITS!!!
min = oden0.min()
max = oden0.max()
print("overdensity minmax", min, max)
print(max - min)

# Obtain the grid zone size (Q4)
redshift = 60
boxsize_init = 20 # Mpc/h
h = 0.7 # dimensionless hub
boxsize_init = boxsize_init / h # Mpc
new_a = (1 + redshift)**(-1) # via redshift relation
new_boxsize = boxsize_init*new_a
per_cell = new_boxsize/128 # in Mpc
per_cell_mpc = per_cell*u.Mpc
per_cell_cm = per_cell_mpc.to(u.cm)
per_cell_volume = per_cell_cm**3

# Mean particle mass
partimass = ad0['particle_mass'].mean() # in grams already.

# Density minimum non-zero assuming a mean particle
denminnonzero = partimass / per_cell_volume
print("denminnonzero " , denminnonzero) """


# TODO: Steps for Q2
"""
# Z = 60
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 4
ds1 = yt.load(resdir + "/RD0002/RD0002")
ds1.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 1
ds2 = yt.load(resdir + "/RD0005/RD0005")
ds2.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 0
ds3 = yt.load(resdir + "/RD0009/RD0009")
ds3.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")

# ZLIST
zlist = [60, 4, 1, 0]
dslis = [ds0, ds1, ds2, ds3]
odenvars = []
aass = []
for ds, z in zip(dslis, zlist):
    ad = ds.all_data()
    den = ad['all_cic']
    oden = ad['dm_overdensity']
    odenmin, odenmax = oden.min(), oden.max()
    odenvar = oden.var(axis=0)
    odenvars.append(odenvar)
    a = (1 + z)**(-1)
    aass.append(a)
    print(odenmin, odenmax, odenvar, "min, max, var")

plt.figure()
plt.plot(aass, odenvars)
plt.show() """

# TODO: Steps for Q3
"""
# Z = 60
ds0 = yt.load(resdir + "/RD0000/RD0000")
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 4
ds1 = yt.load(resdir + "/RD0002/RD0002")
ds1.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 1
ds2 = yt.load(resdir + "/RD0005/RD0005")
ds2.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 0
ds3 = yt.load(resdir + "/RD0009/RD0009")
ds3.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")

# ZLIST
zlist = [60, 4, 1, 0]
ymaxes = np.array([145000, 300000, 200000, 150000])/2
dslis = [ds0, ds1, ds2, ds3]
odenvars = []
aass = []
for ds, z, ymax in zip(dslis, zlist, ymaxes):
    ad = ds.all_data()
    den = ad['all_cic']
    oden = ad['dm_overdensity']
    oden += 1e-10
    oden = np.log10(oden)
    histoden, binoden = np.histogram(oden, bins=100)
    plt.hist(oden, bins=200)
    plt.xlim([-6, 4])
    plt.ylim([0, ymax])
    plt.xlabel(r'$\log_{10}\delta_{DM}$')
    plt.ylabel(r'$n$')
    plt.suptitle("Overdensity histogam for z=" + str(z))
    plt.savefig("histoden" + str(z) + ".png", dpi=300)
    plt.show()
"""


# TODO: Steps for Q4
"""
# Z = 4
ds1 = yt.load(resdir + "/RD0002/RD0002")
ds1.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 0
ds3 = yt.load(resdir + "/RD0009/RD0009")
ds3.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")

# ZLIST
zlist = [4, 0]
dslis = [ds1, ds3]
odenvars = []
aass = []
for ds, z in zip(dslis, zlist):
    # Obtain the grid zone size (Q4)
    redshift = z
    boxsize_init = 20  # Mpc/h
    h = 0.7  # dimensionless hub
    boxsize_init = boxsize_init / h  # Mpc
    new_a = (1 + redshift) ** (-1)  # via redshift relation
    new_boxsize = boxsize_init * new_a
    print(new_boxsize)

    ad = ds.all_data()
    dense_ad = ad.cut_region(["obj['dm_overdensity'] > 200"])
    proj_dense = yt.ProjectionPlot(ds, 'z', "dm_overdensity",
                                   weight_field="dm_overdensity",
                                   data_source=dense_ad,
                                   width=1,
                                   center=[0.5, 0.5, 0.5])
    proj_dense.save("projected_UNcut" + str(z) + ".png")"""

# TODO: Steps for Q5. Downgrade YT to 3.0.4 so that you can do the assessment.
# Z = 4
ds1 = yt.load(resdir + "/RD0002/RD0002")
ds1.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 0
ds3 = yt.load(resdir + "/RD0009/RD0009")
ds3.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
zlist = [4, 0]
dslis = [ds1, ds3]
odenvars = []
aass = []
for ds, z in zip(dslis, zlist):
    # Obtain the grid zone size (Q4)
    redshift = z
    boxsize_init = 20  # Mpc/h
    h = 0.7  # dimensionless hub
    boxsize_init = boxsize_init / h  # Mpc
    new_a = (1 + redshift) ** (-1)  # via redshift relation
    new_boxsize = boxsize_init * new_a
    ad = ds.all_data()
    hc = HaloCatalog(data_ds=ds, finder_method='hop')
    hc.load()
