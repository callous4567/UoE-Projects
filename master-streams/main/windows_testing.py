import matplotlib
import matplotlib as matplotlib
import numpy as np
from matplotlib import rcParams
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
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


uwu = [1,2]
owo = [3,4]
gamma = [-3,-5]
array = np.array([uwu, owo, gamma])
print(np.average(array, axis=0))

#writer = hdfutils.hdf5_writer(datadir, ascii_info.asciiname)
#bhb_data = writer.read_table(ascii_info.bhb, ascii_info.set_raw)


"""
# Grab just the ra/dec/etc for one row
row = bhb_data[1000]
vlos, dist = row['vlos', 'dist']
l,b,dmul,dmub = row['l'],row['b'],row['dmu_l'],row['dmu_b']
x,y,z,vx,vy,vz = row['x','y','z','vx','vy','vz']
edist, edmul, edmub, evlos = row['edist'],row['edmu_l'],row['edmu_b'], row['evlost']
original_galactic_coords = [l,b,dist,dmul,dmub,vlos,edist,edmul,edmub,evlos, edist, edmul, edmub, evlos]

mh = monte_angular()
mh.galdefine(sourcedir, "solar_info.dat")

# This is just some testing to see how convergence of the Monte-Carlo behaves.
n_range = np.append(np.arange(1,100,1), np.arange(100,1000,10))
meanstd_list = []
for n in n_range:
    meanstd_list.append(mh.astropy_vec_monte(original_galactic_coords, n))
meanstd_list = np.array(meanstd_list).T
plt.plot(n_range, meanstd_list[0], color="red", marker="x", lw=1)
plt.plot(n_range, meanstd_list[0+3], color="green", marker="x", lw=1)
legend_elements = [Patch(edgecolor='black', facecolor="red", label=r'$\mu$'),
                   Patch(edgecolor='black', facecolor="green", label=r'$\sigma$')]
plt.legend(handles=legend_elements, loc="upper right")
plt.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
plt.savefig(imgdir + "\\monte_x_test.png")
plt.show()
plt.plot(n_range, meanstd_list[1], color="red", marker="x", lw=1)
plt.plot(n_range, meanstd_list[1+3], color="green", marker="x", lw=1)
legend_elements = [Patch(edgecolor='black', facecolor="red", label=r'$\mu$'),
                   Patch(edgecolor='black', facecolor="green", label=r'$\sigma$')]
plt.legend(handles=legend_elements, loc="upper right")
plt.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
plt.savefig(imgdir + "\\monte_y_test.png")
plt.show()
plt.plot(n_range, meanstd_list[2], color="red", marker="x", lw=1)
plt.plot(n_range, meanstd_list[2+3], color="green", marker="x", lw=1)
legend_elements = [Patch(edgecolor='black', facecolor="red", label=r'$\mu$'),
                   Patch(edgecolor='black', facecolor="green", label=r'$\sigma$')]
plt.legend(handles=legend_elements, loc="upper right")
plt.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)
plt.savefig(imgdir + "\\monte_z_test.png")
plt.show()
"""
