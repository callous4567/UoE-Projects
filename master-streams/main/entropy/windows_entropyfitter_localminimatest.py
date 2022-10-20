import multiprocessing
import pickle
import time
import warnings
import numpy as np
from matplotlib import pyplot as plt, rc, cm
import ascii_info
import energistics_constants
import graphutils
import hdfutils
import windows_directories
import energistics
from energistics_constants import M_nfw, M_d, a_b, a_d, M_b, a_nfw, b_d, c_nfw
import mpl_scatter_density
# Etc etc
# Enable TeX
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

warnings.filterwarnings("ignore",
                        message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")

if __name__ == "__main__":

    # Load in the "mean" percenttable and map
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")

    # Get clustering data and select for relevant data.
    clustering_final = membership_table['greatcircle_probable_clust']

    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
        clustering_prelim = pickle.load(file=f)

    # Set
    clustering = clustering_final

    # Grab the data you're going to use
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.set_raw)

    # Decide monte points
    nmonte = 10000

    # Decide group and the clusters to cluster
    group = ascii_info.fullgroup
    clust_to_run = [1,3,4,8]
    ranges = np.array([[1,1],
                       [0.5,2],
                       [0.5,2]]) # Md Mnfw Cnfw

    # Grab and produce required data lists
    x, y, z, vx, vy, vz = [[] for d in range(6)]
    for clust_to_fit in clust_to_run:
        data_to_fit = table[[True if d == clust_to_fit else False for d in clustering]]
        xx, yy, zz = data_to_fit['x'], data_to_fit['y'], data_to_fit['z']
        vxx, vyy, vzz = data_to_fit['vx'], data_to_fit['vy'], data_to_fit['vz']
        x.append(xx)
        y.append(yy)
        z.append(zz)
        vx.append(vxx)
        vy.append(vyy)
        vz.append(vzz)

    # Set up entrofitter
    entrofit = energistics.entro_fitter(x, y, z,
                                        vx, vy, vz,
                                        ranges[0]*M_d, ranges[1]*M_nfw, ranges[2]*c_nfw,
                                        [M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw],
                                        nmonte,
                                        0.2, 100)

    # Orbi
    orbigist = energistics.orbigistics()

    # Set up points
    points = entrofit.genpoints()
    M_ds, M_Nfws, c_nfws = points.T

    # Bin the points into the number of cores we'll use
    ncores = 10
    multipoints = np.split(points, ncores)

    # Generate the entropies for the points
    pool = multiprocessing.Pool(ncores)
    print("Running Pool")
    entropies = pool.map(entrofit.list_entropy, multipoints)
    print("Pool Done")
    pool.close()

    # Concatenate them all
    entropies = np.concatenate(entropies)
    masses = []

    mean_r = 15

    # Obtain mass estimates for this in all points
    for num,point in enumerate(points):
        print(num)
        mass = orbigist.ret_halo_disk_manudisk(mean_r, *point)
        masses.append(mass)

    # Get the best points (minimum locus.)
    best_vecs, bestropies = entrofit.lowest_monte_X_list(200, points, entropies)
    best_vecs = np.array(best_vecs).T
    bestmds, bestmnfs, bestcnfws = best_vecs

    # Divide masses by the highest mass possible
    masses = np.array(masses)

    # set up plot
    cmap = cm.get_cmap("hot")
    cmap.set_under('w')
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(6,6))
    axis_1 = axs.scatter(M_Nfws, c_nfws, cmap=cmap, c=masses, zorder=0)  # , norm=norm)
    axs.scatter(bestmnfs, bestcnfws, zorder=2, color='red', s=5)
    axs.grid(which='major', color='pink')
    plt.colorbar(axis_1, label='Enclosed mass / ' + r'$M_\odot$')
    axs.set(xlabel=r'$\frac{M_H}{M_\odot}$',
            ylabel=r'$c_H$',
            xlim=ranges[1]*energistics_constants.M_nfw,
            ylim=ranges[2]*energistics_constants.c_nfw)
    plt.savefig(windows_directories.imgdir + "\\entropyorphan_massmap.png", dpi=300)
    plt.show(dpi=72)







