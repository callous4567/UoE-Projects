import multiprocessing
import pickle
import time

import numpy as np
from astropy.table import Table

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
    clustering = membership_table['greatcircle_probable_clust']

    with open(windows_directories.clusterdir + "\\" + "fullgroup.cluster.txt", 'rb') as f:
        clustering = pickle.load(file=f)

    # Grab the data you're going to use
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.set_raw)

    # Set clusters
    clusts_to_fitTEST = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    clusts_to_energyplot = [0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17]
    streams_to_energyplot = [11, 15, 10, 3, 8, 12, 0, 7, 14, 1, 9, 2, 4]
    rest_to_energyplot = [13,0,1,2,3,4,5,7,8,9,10,11,12,14,17]
    streams_to_check = [1, 3, 8, 10]
    streams_surviving_greatpercent = [1,3,5,6,4,11,13,8,16]

    #specgrapher = graphutils.spec_graph()
    #for stream in streams_surviving_greatpercent:
    #    specgrapher.lb_orbits(table[[True if d == stream else False for d in clustering]], 0.3e9, [-180,180], [-90,90],
    #                          windows_directories.imgdir + "\\" + "greatmemberpercent_unlimited" + str(stream) + ".png", line=False,
    #                          points=4000)
    #specgrapher.energy_plot(table, clustering, rest_to_energyplot, [str(d) for d in rest_to_energyplot],
    #                        [-7000,7000], [-150000,0], windows_directories.imgdir + "\\streamplot.png", mcmillan=False, kpcmyr=False)
    #specgrapher.clust_lb_colors(table, clustering, streams_to_check,
    #                         windows_directories.imgdir + "\\streamskytest.png")

    pointtable = Table(names=["clust", "M_d", "M_nfw", "c_nfw"], data=np.array([clusts_to_fitTEST,
                                                                                np.zeros_like(clusts_to_fitTEST),
                                                                                np.zeros_like(clusts_to_fitTEST),
                                                                                np.zeros_like(clusts_to_fitTEST)]).T)
    pointtable['clust'] = clusts_to_fitTEST
    points_saved = []
    for clusts_to_fit in clusts_to_fitTEST:
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
                                            [0.2*M_d, 2*M_d], [0.1*M_nfw, 2*M_nfw], [0.5*c_nfw, 2*c_nfw],
                                            [M_b, a_b, M_d, a_d, b_d, M_nfw, a_nfw, c_nfw],
                                            1000000,
                                            0.3, 250)

        # Set up points
        points = entrofit.genpoints()

        # Bin the points into the number of cores we'll use
        ncores = 8
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
    writer.write_table("entropy_fit", "maindata", pointtable)
