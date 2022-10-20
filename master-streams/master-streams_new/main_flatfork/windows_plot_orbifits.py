import copy
import os
import numpy as np
from matplotlib import pyplot as plt
import astropy.units as u
import ascii_info_new
import energistics_new
import graphutils_new
import hdfutils
import windows_directories_new

# Take data greattable and obtain clusters to orbifit (length =/= 0)
clusters_to_orbifit = []
total_percent_table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                              "total_percent_table_greatfitted")

for num,clust in enumerate(total_percent_table['cluster']):

    if clust not in clusters_to_orbifit:

        if clust != -1:

            # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
            if total_percent_table[num]['membership_fraction'] < 0.1:

                clusters_to_orbifit.append(clust)

# Print it
print("Orbifitting for ", clusters_to_orbifit)
import time
time.sleep(10)

# Load in the table, greattable, and so forth.
writer = hdfutils.hdf5_writer(windows_directories_new.datadir, ascii_info_new.flatfork_asciiname)
membership_table = writer.read_table(ascii_info_new.fullgroup, "percent_table_greatfitted")

# Get clustering data and select for relevant data.
clustering_final = membership_table['greatcircle_probable_clust']

# Grab the data you're going to use
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                              ascii_info_new.set_raw)

# Prelim clust
prelim_clust = table['prelim_clust']

# Set up spec-grapher
specgrapher = graphutils_new.spec_graph()

do_clustplots = True
do_lbplots = True
do_energyplots = True

# Iterate over clusters
if do_lbplots == True:

    savedir = os.path.join(windows_directories_new.imgdir, "flatfork_lb_plots")

    try:
        os.mkdir(savedir)
    except:
        pass
    for clust in clusters_to_orbifit:

        # Grab an orbit fit
        orbifit = energistics_new.orbifitter().flatfork_galpy_final_fitting(table,
                                                                        clust,
                                                                        2000,
                                                                        0.3e9,
                                                                        1000,
                                                                        True,
                                                                        False,
                                                                        False,
                                                                        False)

        # Forward Integral
        forward = copy.deepcopy(orbifit)
        backward = copy.deepcopy(orbifit)
        forward.integrate((np.linspace(0, 1e9, 2000)*u.yr), energistics_new.orbigistics().pot)
        backward.integrate((np.linspace(0, -1e9, 2000)*u.yr), energistics_new.orbigistics().pot)
        llsf, bbsf, llsb, bbsb = forward.ll((np.linspace(0, 0.7e9, 2000)*u.yr)).value, \
                                 forward.bb((np.linspace(0, 0.7e9, 2000)*u.yr)).value, \
                                 backward.ll((np.linspace(0, -0.15e9, 2000)*u.yr)).value, \
                                 backward.bb((np.linspace(0, -0.15e9, 2000)*u.yr)).value
        lls, bbs = np.concatenate([llsf, llsb]), np.concatenate([bbsf, bbsb])
        lls = [d - 360 if d > 180 else d for d in lls]

        specgrapher = graphutils_new.spec_graph()
        fig, axs = specgrapher.lb_orbits(table[[True if d == clust
                                                else False
                                                for d in clustering_final]],
                                         0.3e9, [-180, 180],
                                         [-90, 90],None,
                                         line=False,points=4000)
        axs.scatter(lls, bbs, color='lime', marker='o', s=3)
        # Misc axis things
        axs.set_facecolor("k")
        axs.grid(True, which='major', alpha=1, linewidth=0.25, color='white')
        axs.grid(True, which='minor', alpha=1, linewidth=0.25, color='white')
        axs.grid(color="white")

        # Savedir
        savepath = os.path.join(savedir,str(clust) + ".png")
        plt.savefig(savepath, dpi=300, transparent=False)
        plt.close()

if do_energyplots == True:

    savedir = os.path.join(windows_directories_new.imgdir, "flatfork_energy_plots")
    try:
        os.mkdir(savedir)
    except:
        pass
    prelim_path = os.path.join(savedir, "flatfork_streamplot_prelim.png")
    final_path = os.path.join(savedir, "flatfork_streamplot_final.png")

    specgrapher.energy_plot(table, prelim_clust, clusters_to_orbifit, [str(d) for d in clusters_to_orbifit],
                            [-7000, 7000], [-150000, 0], prelim_path,
                            mcmillan=False, kpcmyr=False)
    plt.close()
    specgrapher.energy_plot(table, clustering_final, clusters_to_orbifit, [str(d) for d in clusters_to_orbifit],
                            [-7000,7000], [-150000,0], final_path,
                            mcmillan=False, kpcmyr=False)
    plt.close()

if do_clustplots == True:

    savedir = os.path.join(windows_directories_new.imgdir, "flatfork_clustplots")
    try:
        os.mkdir(savedir)
    except:
        pass

    # Trim noise. TODO: Make it so that we can feed in the clustercount, nice and ordered (to parity colour between two.)
    clustering_final = np.array(clustering_final)
    truefalse = [False if d == -1 else True for d in clustering_final]
    tablee, clustering_final = table[truefalse], clustering_final[truefalse]
    graphutils_new.twod_graph().tripL_colour(tablee, 15, clustering_final, os.path.join(savedir, "tripL_final.png"))

    plt.close()

    clustering_prelim = np.array(prelim_clust)
    truefalse = [False if d == -1 else True for d in clustering_prelim]
    tablee, clustering_prelim = table[truefalse], clustering_prelim[truefalse]
    graphutils_new.twod_graph().tripL_colour(tablee, 15, clustering_prelim, os.path.join(savedir, "tripL_prelim.png"))

    plt.close()