import copy
import os
import numpy as np
from matplotlib import pyplot as plt
import astropy.units as u
import ascii_info
import energistics
import graphutils
import hdfutils
import windows_directories

# Take data greattable and obtain clusters to orbifit (length =/= 0)
clusters_to_orbifit = []

# Grab the data you're going to use
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.set_raw)

# Prelim clust
prelim_clust = table['prelim_clust']

# Decide ones we are fitting
set = list(set(prelim_clust))
for clust in set:

    if clust not in clusters_to_orbifit:

        if clust != -1:

            size_fraction = len(np.where(prelim_clust==clust)[0])/len(prelim_clust)

            # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
            if size_fraction < 0.1:

                clusters_to_orbifit.append(clust)

# Print it
print("Orbifitting for ", clusters_to_orbifit)
import time
time.sleep(10)

# Set up the data.

# Load in the table, greattable, and so forth.
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)

# Grab the data you're going to use
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.set_raw)

# Prelim clust
prelim_clust = table['prelim_clust']

# Set up spec-grapher
specgrapher = graphutils.spec_graph()

do_lbplots = True

# Iterate over clusters
if do_lbplots == True:
    savedir = os.path.join(windows_directories.imgdir, "lb_plots_maindata")
    try:
        os.mkdir(savedir)
    except:
        pass
    for clust in clusters_to_orbifit:

        # Grab an orbit fit
        orbifit = energistics.orbifitter().galpy_fitting_nomemb(table,
                                                                         prelim_clust,
                                                                         clust,
                                                                         2000,
                                                                         0.3e9,
                                                                         1000,
                                                                         True,
                                                                         False,
                                                                         False,
                                                                         False,
                                                                         extra_text="for_lbplots")

        # Forward Integral
        forward = copy.deepcopy(orbifit)
        backward = copy.deepcopy(orbifit)
        forward.integrate((np.linspace(0, 1e9, 2000)*u.yr), energistics.orbigistics().pot)
        backward.integrate((np.linspace(0, -1e9, 2000)*u.yr), energistics.orbigistics().pot)
        llsf, bbsf, llsb, bbsb = forward.ll((np.linspace(0, 0.7e9, 2000)*u.yr)).value, \
                                 forward.bb((np.linspace(0, 0.7e9, 2000)*u.yr)).value, \
                                 backward.ll((np.linspace(0, -0.15e9, 2000)*u.yr)).value, \
                                 backward.bb((np.linspace(0, -0.15e9, 2000)*u.yr)).value
        lls, bbs = np.concatenate([llsf, llsb]), np.concatenate([bbsf, bbsb])
        lls = [d - 360 if d > 180 else d for d in lls]

        specgrapher = graphutils.spec_graph()
        fig, axs = specgrapher.lb_orbits(table[[True if d == clust
                                                else False
                                                for d in prelim_clust]],
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
