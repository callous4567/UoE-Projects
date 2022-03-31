import multiprocessing
import pickle
import time

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

import ascii_info
import energistics
import graphutils
import hdfutils
import windows_directories
from energistics_constants import M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw


# Set up the data.
if __name__ == "__main__":

    # Load in the "mean" percenttable and map
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")

    # Get clustering data and select for relevant data.
    clustering_final = membership_table['greatcircle_probable_clust']

    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
        clustering_prelim = pickle.load(file=f)

    # Grab the data you're going to use
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.set_raw)
    table = table[[True if d == 4 else False for d in clustering_final]]
    print(np.max(table['r']))
    time.sleep(5000)
    #truefalse = [True if d == 5 else False for d in clustering_prelim]
    #orbigist = energistics.orbigistics()
    #table = table[truefalse]
    #table = orbigist.converter.nowrite_GAL_to_ICRS(table)
    #table.write(windows_directories.datadir + "\\zero.xml", table_id="zeroth_table", format='votable', overwrite='true')
    #print(table[['ra','dec','vlos']])
    # Set clusters
    #clusts_to_fitTEST = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    #clusts_to_energyplot = [0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17]
    #streams_to_energyplot = [11, 15, 10, 3, 8, 12, 0, 7, 14, 1, 9, 2, 4]
    rest_to_energyplot = [16,13,0,1,2,3,4,5,6,7,8,9,10,11,12,14,15]
    rest_to_skyplot = [0,1,2,3,4,5,6,7,8,9,10,11,12,14,15]
    final_rest_to_energyplot = [16,13,1,3,4,5,6,8,11]
    final_to_skyplot = [1,3,4,5,6,8,11]
    greatpercent_streams = [13]
    #streams_surviving_greatpercent = [1] # ,3,5,6,4,11,13,8,16

    #specgrapher = graphutils.spec_graph()

    #graphutils.twod_graph().tripL_colour(table, 15, clustering_prelim, windows_directories.imgdir + "\\tripL_preliminary.png")

    # Trim noise
    #clustering_final = np.array(clustering_final)
    #truefalse = [False if d == -1 else True for d in clustering_final]
    #table, clustering_final = table[truefalse], clustering_final[truefalse]
    #graphutils.twod_graph().tripL_colour(table, 15, clustering_final, windows_directories.imgdir + "\\tripL_final.png")

    clustering_prelim = np.array(clustering_prelim)
    truefalse = [False if d == -1 else True for d in clustering_prelim]
    table, clustering_prelim = table[truefalse], clustering_prelim[truefalse]
    graphutils.twod_graph().tripL_colour(table, 15, clustering_prelim, windows_directories.imgdir + "\\tripL_prelim.png")



    #for stream in rest_to_skyplot:
    #    specgrapher.lb_orbits(table[[True if d == stream else False for d in clustering_prelim]], 0.3e9, [-180,180], [-90,90],
    #                          windows_directories.imgdir + "\\" + "maindata" + str(stream) + ".png", line=False,
    #                          points=4000)
    #for stream in final_to_skyplot:
    #    specgrapher.lb_orbits(table[[True if d == stream else False for d in clustering_final]], 0.3e9, [-180,180], [-90,90],
    #                          windows_directories.imgdir + "\\" + "finalgreatcount" + str(stream) + ".png", line=False,
    #                          points=4000)
    """
    specgrapher.energy_plot(table, clustering_prelim, rest_to_energyplot, [str(d) for d in rest_to_energyplot],
                            [-7000, 7000], [-150000, 0], windows_directories.imgdir + "\\streamplot_preliminary.png",
                            mcmillan=False, kpcmyr=False)
    plt.close()
    specgrapher.energy_plot(table, clustering_final, final_rest_to_energyplot, [str(d) for d in final_rest_to_energyplot],
                            [-7000,7000], [-150000,0], windows_directories.imgdir + "\\streamplot_final.png", mcmillan=False, kpcmyr=False)
    plt.close()
    specgrapher.clust_lb_colors(table, clustering_prelim, rest_to_skyplot,
                             windows_directories.imgdir + "\\streamskytest_prelim.png", False)
    plt.close()
    specgrapher.clust_lb_colors(table, clustering_final, final_to_skyplot,
                             windows_directories.imgdir + "\\streamskytest_final.png", False)
    plt.close()
    specgrapher.clust_lb_colors(table, clustering_final, [13],
                             windows_directories.imgdir + "\\streamskytest_sagittarius.png", True) """

    """
    pointtable = Table(names=["clust", "M_d", "M_nfw", "c_nfw"], data=np.array([streams_surviving_greatpercent,
                                                                                np.zeros_like(streams_surviving_greatpercent),
                                                                                np.zeros_like(streams_surviving_greatpercent),
                                                                                np.zeros_like(streams_surviving_greatpercent)]).T)
    pointtable['clust'] = streams_surviving_greatpercent
    points_saved = []
    for clusts_to_fit in streams_surviving_greatpercent:
        clusts_to_fit = [clusts_to_fit]
        # Grab and produce required data lists
        x,y,z,vx,vy,vz = [[] for d in range(6)]
        for clust_to_fit in clusts_to_fit:
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
        entrofit = energistics.entro_fitter(x,y,z,
                                            vx,vy,vz,
                                            [M_d, M_d], [0.1*M_nfw, 2*M_nfw], [0.5*c_nfw, 2*c_nfw],
                                            [M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw],
                                            1000000,
                                            0.2, 100)

        # Set up points
        points = entrofit.genpoints()
        M_ds, M_Nfws, c_nfws = points.T

        # Bin the points into the number of cores we'll use
        ncores = 10
        multipoints = np.split(points, ncores)

        # Generate the entropies for the points
        pool = multiprocessing.Pool(ncores)
        entropies = pool.map(entrofit.list_entropy, multipoints)
        pool.close()

        # Concatenate them all
        entropies = np.concatenate(entropies)

        # Get the minimum point
        best_point = points[np.argmin(entropies)]/np.array([M_d, M_nfw, c_nfw])
        points_saved.append(best_point)

    points_saved = np.array(points_saved).T
    pointtable['M_d'], pointtable['M_nfw'], pointtable['c_nfw'] = points_saved
    writer.write_table("entropy_fit", "orphan_greatpercent_test_onlyhalo_andconc", pointtable)

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))
    plt.subplots_adjust(bottom=0.15)
    axs.scatter(M_Nfws / M_nfw, entropies, marker='.', s=0.1, color='red')
    plt.grid(which='major', color='pink')
    plt.xlabel(r'$\frac{M}{M_H}$')
    plt.ylabel(r'$\tilde{H}$')
    #plt.savefig(windows_directories.imgdir + "\\entropy_test_MNFW_orphan.png", dpi=300)
    plt.savefig(windows_directories.imgdir + "\\entropy_test_MNFW_orphan.pgf")
    plt.show()

    """

