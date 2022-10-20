
# Obtain the monte carlo realisation membership
import pickle

import numpy as np
from matplotlib import pyplot as plt, rc
from matplotlib.patches import Patch

import ascii_info_new
import energistics_new_constants
import hdfutils
import windows_directories_new

# Load writer
from energistics_new import orbigistics
from galcentricutils_new import angular

# Etc etc
# Enable TeX
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Writer
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.asciiname)

# Grab final membership
membership_table = writer.read_table(ascii_info_new.fullgroup, "percent_table_greatfitted")
clustering = membership_table['probable_clust']

with open(windows_directories_new.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
    clustering_prelim = pickle.load(file=f)

# Grab entropy table (final)
entropy_table = writer.read_table("entropy_fit_witherrors", "per_cluster_haloonly")


# Grab radii/etc
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.fullset)
table = angular().get_polar(table)

# Set up orbigist
orbigist = orbigistics()

# Radii and Masses
radii, masses = [],[]
rerr, merr = [],[]
all_rs, all_ms = [],[]

# Scatter plot
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

# For each cluster in the entropy table, estimate mass interior of the cluster (roughly.)
colors = []
for cluster, M_d, M_nfw, c_nfw in zip(entropy_table['cluster'], entropy_table['M_d'], entropy_table['M_nfw'], entropy_table['c_nfw']):

    if cluster == 4:
        pass
    elif cluster == 11:
        pass
    else:
        # Truefalse
        truefalse = [True if d == cluster else False for d in clustering]

        # Isolate the r
        r_vals = table[truefalse]['r']
        mean_r = np.mean(r_vals)
        all_rs += list(r_vals)

        # Obtain mass estimates for this
        m_vals = [orbigist.ret_halo_disk_manudisk(r, M_d*energistics_new_constants.M_d, M_nfw*energistics_new_constants.M_nfw, c_nfw*energistics_new_constants.c_nfw) for r in r_vals]
        mean_mass = np.mean(m_vals)
        all_ms += m_vals

        # Append
        radii.append(mean_r)
        masses.append(mean_mass)
        rerr.append(np.std(r_vals))
        merr.append(np.std(m_vals))
        axs.scatter(all_rs, all_ms, s=1, zorder=0, color="red")


axs.scatter(x=radii, y=masses, color="blue", marker='x', s=40, zorder=1)

# Also grab the Cautun potential/etc
radii = np.linspace(10, 55, 50)
masses = [orbigist.ret_cautun(d) for d in radii]
masses_2 = [orbigist.ret_bovy(d) for d in radii]
axs.set(xlim=[10,55],
        ylim=[0.6e11, 3.9e11])

axs.plot(radii, masses, color='black')
axs.plot(radii, masses_2, color='gray')

# Legend
legend_elements = [Patch(edgecolor='black',
                         facecolor='black',
                         label="Cautun 2020"),
                   Patch(edgecolor='gray',
                         facecolor='gray',
                         label="Bovy 2014")]
plt.legend(handles=legend_elements, loc='lower right')

axs.grid(which='major', color='pink')
axs.set(xlabel="r / kpc")
axs.set(ylabel=r'$\frac{M}{M_\odot}(r)$')
plt.savefig(windows_directories_new.imgdir + "\\entropy_oned_halo", dpi=300)
plt.show()
plt.close()



"""

# Load the entropy fit tables
entropy_table = entropy_table.to_pandas(index=False)

dataframe = DataFrame(data_array, columns=["Cluster", "Preliminary", "Monte", "Monte GCCs"])
dataframe.to_latex(caption="Cluster memberships determined for the dataset: for the preliminary clustering, quasi-soft clustering, and quasi-soft refined clustering. The cluster index -1 corresponds to the stars labelled as ``noise'' by \verb|hdbscan|.",
                   buf=windows_directories_new.imgdir + "\\memberships.tex",
                   label="tab:memberships",
                   float_format="%.0d",
                   index=False) """