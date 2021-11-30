import numpy as np
import yt
from yt.units import *
import astropy.units as u
from dm_overdensity import _overdensity
import matplotlib.pyplot as plt
import yt.units as u
from yt.analysis_modules.halo_analysis.api import *

# This has the working directory for all (directly) related Enzo stuff.
resdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm/results"
supdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/sup-pm"
z60 = resdir + "/z60"
z10 = resdir + "/z10"

# WSdir for this workshop (i.e. save images here.)
wsdir = "/home/callicious/Documents/pycharm/CompAstro/pmesh/ws3"

# Code copied from WS2 TODO: Q2
# Z = 60
ds0 = yt.load(z60 + '/RD0000/RD0000')
ds0.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")
# Z = 10
ds2 = yt.load(z10 + '/RD0000/RD0000')
ds2.add_field(name=("gas",'dm_overdensity'),
              function=_overdensity,
              sampling_type="cell")

# ZLIST
zlist = [60, 10]
mass_list = [1e11, 3e11, 5e11, 7e11, 9e11, 11e11, 13e11, 15e11, 17e11, 19e11]
dslis = [ds0, ds2]
odenvars = []
aass = []
for ds, z in zip(dslis, zlist):
    """
    ad = ds.all_data()
    den = ad['all_cic']
    oden = ad['dm_overdensity']
    odenmin, odenmax = oden.min(), oden.max()
    odenvar = oden.var(axis=0)
    odenvars.append(odenvar)
    a = (1 + z)**(-1)
    aass.append(a)
    mean = oden.mean()
    print(odenmin, odenmax, odenvar, mean, "min, max, var, mean")
    """
    # For the entire simulation volume
    slc0 = yt.SlicePlot(ds, 'z', 'all_cic', center=[0.5, 0.5, 0], width=(1))
    slc0.save('densiproj_' + str(z) + '.png')

    # Find halos
    hc = HaloCatalog(data_ds=ds, finder_method='hop')
    hc.create()
    nhalos = []
    for mass in mass_list:
        hc.add_filter("quantity_value", "particle_mass", ">", 1e8, "Msun")
        number = len(hc.halo_list) # get the number above this mass. Cumulative dist...
        nhalos.append(number)

    # Plot the mass distribution function (unnormalized.)
    # NONE OF THIS IS WORKING! No time to fix it.
    plt.figure()
    plt.plot(mass_list, nhalos)
    plt.xlabel("Mass")
    plt.ylabel("Halo Count")
    plt.savefig("halo_mass_dist_funct_unnormalized.png")

