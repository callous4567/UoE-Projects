import multiprocessing
import time

import matplotlib
import matplotlib as matplotlib
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import timeit

import graphutils_new
import windows_directories_new
import windows_multiprocessing_new
plt.rcParams["font.family"] = "serif"
from astropy.coordinates import Galactic
import astropy.units as u
import galcentricutils_new
from windows_directories_new import sourcedir, imgdir, datadir
import ascii_info_new
from galcentricutils_new import monte_angular, galconversion
import hdfutils
import math
import timeit



writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
bhb_data = writer.read_table(ascii_info_new.bhb, ascii_info_new.set_raw)
test_bhb_rows = bhb_data[0:5]
#monte_dataframe = galcentricutils_new.monte_angular().table_covmonte(test_bhb_rows, 10)
#writer.write_df("test_dataframe","test_monte",monte_dataframe)
#monte_read = writer.read_df("test_dataframe","test_monte")
#print(monte_read['covtrix'][0])
"""
# Grab just the ra/dec/etc for one row
row = bhb_data[1000]
vlos, dist = row['vlos', 'dist']
l,b,dmul,dmub = row['l'],row['b'],row['dmu_l'],row['dmu_b']
edist, edmul, demub, evlos = row['edist'],row['edmu_l'],row['edmu_b'],row['evlost']
vec = [l,b,dist,dmul,dmub,vlos, edist, edmul, demub, evlos]

monte = galcentricutils_new.monte_angular()
monte.galdefine(sourcedir, "solar_info.dat")
test_monte=monte.vec_covmonte(vec,100)
for test in test_monte:
    print(test)
"""

"""
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
full_data = writer.read_table(ascii_info_new.bhb, ascii_info_new.set_raw)
converter = galcentricutils_new.galconversion()
converter.solinfo_grab(sourcedir, "solar_info.dat")
converter.solgal_set()
grapher = graphutils_new.twod_graph()
grapher.lbplot(full_data)
grapher.haitoff(full_data)
full_data_icrs = converter.nowrite_GAL_to_ICRS(full_data)
grapher.radecplot(full_data_icrs) """

"""
theta,phi = galcentricutils_new.angular().vec_polar([3000,2000,-2500])[1:]
#theta,phi = [272,103]
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)
full_data = writer.read_table(ascii_info_new.fullgroup, ascii_info_new.fullset)
greatcircled = galcentricutils_new.greatcount().gcc_table(full_data, theta, phi, 90, 20)
#clustered_table = galcentricutils_new.cluster().gaussmix(greatcircled,k,savedex="test",browser=True,graph=True)
#dbsclusterd = galcentricutils_new.cluster().dbs(greatcircled, eps=1.1e3, min_samples=15, browser=True)
hdbscluster = galcentricutils_new.cluster().hdbs(greatcircled, browser=True)
#afproped = galcentricutils_new.cluster().afprop(greatcircled)

k = 10
for i in range(k):
    table_by_k = hdbscluster.group_by("k_index")
    k_mask = table_by_k.groups.keys["k_index"] == i
    masked_by_k = table_by_k.groups[k_mask]
    L = np.array([masked_by_k['Lx'],masked_by_k['Ly'],masked_by_k['Lz']])
    print(np.mean(L, axis=1)) """





