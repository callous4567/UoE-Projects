import multiprocessing
import time

import matplotlib
import matplotlib as matplotlib
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import timeit

import graphutils
import windows_directories
import windows_multiprocessing
plt.rcParams["font.family"] = "serif"
from astropy.coordinates import Galactic
import astropy.units as u
import galcentricutils
from windows_directories import sourcedir, imgdir, datadir
import ascii_info
from galcentricutils import monte_angular, galconversion
import hdfutils
import math
import timeit


"""
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
bhb_data = writer.read_table(ascii_info.bhb, ascii_info.set_raw)
# Grab just the ra/dec/etc for one row
row = bhb_data[1000]
vlos, dist = row['vlos', 'dist']
l,b,dmul,dmub = row['l'],row['b'],row['dmu_l'],row['dmu_b']
edist, edmul, demub, evlos = row['edist'],row['edmu_l'],row['edmu_b'],row['evlost']
vec = [l,b,dist,dmul,dmub,vlos, edist, edmul, demub, evlos]

monte = galcentricutils.monte_angular()
monte.galdefine(sourcedir, "solar_info.dat")
test_monte=monte.vec_covmonte(vec,100)
for test in test_monte:
    print(test)
"""

"""
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
full_data = writer.read_table(ascii_info.bhb, ascii_info.set_raw)
converter = galcentricutils.galconversion()
converter.solinfo_grab(sourcedir, "solar_info.dat")
converter.solgal_set()
grapher = graphutils.twod_graph()
grapher.lbplot(full_data)
grapher.haitoff(full_data)
full_data_icrs = converter.nowrite_GAL_to_ICRS(full_data)
grapher.radecplot(full_data_icrs) """


theta,phi = galcentricutils.angular().vec_polar([3000,2000,-2500])[1:]
#theta,phi = [272,103]
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
full_data = writer.read_table(ascii_info.fullgroup, ascii_info.fullset)
greatcircled = galcentricutils.greatcount().gcc_table(full_data, theta, phi, 90, 20)
#clustered_table = galcentricutils.cluster().gaussmix(greatcircled,k,savedex="test",browser=True,graph=True)
#dbsclusterd = galcentricutils.cluster().dbs(greatcircled, eps=1.1e3, min_samples=15, browser=True)
hdbscluster = galcentricutils.cluster().hdbs(greatcircled, browser=True)
#afproped = galcentricutils.cluster().afprop(greatcircled)

k = 10
for i in range(k):
    table_by_k = hdbscluster.group_by("k_index")
    k_mask = table_by_k.groups.keys["k_index"] == i
    masked_by_k = table_by_k.groups[k_mask]
    L = np.array([masked_by_k['Lx'],masked_by_k['Ly'],masked_by_k['Lz']])
    print(np.mean(L, axis=1))





