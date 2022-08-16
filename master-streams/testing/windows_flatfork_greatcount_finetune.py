# Automated process proposal.
"""

For each cluster found in maindata
Isolate the greatcircle
Optimize for the width of the greatcircle
Fine-tune the parameters for HDBSCAN such that the cluster is isolated perfectly
Train the clusterer on this greatfit just as we did with the flat clustering previously
Run the entire process though
If it improves/adds any stars, adjust, if not, pass (just like we did.)

"""
"""

To optimize for the width
Use 5x the standard deviation in distance to greatcircle
Or use maximum distance
Winners choice. 
"""
"""

To optimize for the cluster being isolated perfectly by HDBSCAN
Use the min cluster size that the usual clusterer uses (just as we did in cluster_maindata clusterer object)
Adjust min_samples iteratively until the original cluster is recovered- the *entire* cluster- at a minimum.  
If this fails, manual intervention is required and a flat clusterer manually adjusted to provide results is necessary. 

"""

# Imports
import os
owd = os.getcwd()
from galcentricutils import greatcount, compclust, cluster3d
import numpy as np
from hdbscan import HDBSCAN
import hdfutils
import windows_directories
import ascii_info

# set up dependent classes
gc = greatcount()
cc = compclust()
c3d = cluster3d()
c3d.sizediflim = 2

# Preliminary Clustering for Greatcount
for i in range(25):
    clust_to_finetune = i
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.flatfork_asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)

    greattable = hdfutils.hdf5_writer(windows_directories.datadir,
                                             ascii_info.flatfork_asciiname).read_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles")
    greattable['max_dist'] += 1
    row = np.where(greattable['clust_to_try'] == clust_to_finetune)[0][0]
    thetaphiwidth = greattable[row][['theta', 'phi', 'max_dist']]
    table = gc.gcc_table_retain(table, *thetaphiwidth)
    cut_table = table[table['greatcount']]
    cut_clustered = cut_table['prelim_clust']
    is_gse = np.array([True if d == 23 else False for d in cut_clustered])
    is_gse = np.sum(is_gse)*(1.1)
    cut_array = np.array([cut_table['Lx'], cut_table['Ly'], cut_table['Lz']]).T
    max_clust_size = is_gse
    best_samples, passfail, \
    numclusts, scores, sizedifs = c3d.minsamples_finetune(cut_array, cut_clustered,
                        clust_to_finetune, [4,20], ascii_info.fulldata_minpar[0],
                        20, 0.8, 0.8, max_clust_size)

    min_clust_size = ascii_info.fulldata_minpar[0]
    min_samples = best_samples

    # Cluster it accordingly
    """
    trial_clustered = flat.HDBSCAN_flat(
        X=cut_array,
        n_clusters=numclusts,
        min_cluster_size=int(min_clust_size),
        min_samples=int(min_samples),
        metric='l2',
        algorithm='best',
        prediction_data=True
    ).labels_
    """
    trial_clustered = HDBSCAN(min_cluster_size=int(min_clust_size),
                              max_cluster_size=is_gse,
                              min_samples=int(min_samples),
                              algorithm='best',
                              metric='l2').fit_predict(cut_array)

    import graphutils
    dd = "test_" + str(clust_to_finetune)
    graphutils.threed_graph().kmeans_L_array(cut_array, cut_clustered, dd, browser=False, outliers=True)
    # graphutils.threed_graph().kmeans_L_array(cut_array, trial_clustered, False, browser=True, outliers=True)
    """
    import pickle
    with open(os.path.join(owd,"dump1.txt"),"wb") as f:
        pickle.dump(obj=cut_clustered, file=f)
    with open(os.path.join(owd,"dump2.txt"),"wb") as f:
        pickle.dump(obj=trial_clustered, file=f)
    """
    trial_clustered = cc.compclust_multilabel_julia(cut_clustered, trial_clustered, 30)
    dd = "test_done_" + str(clust_to_finetune)
    graphutils.threed_graph().kmeans_L_array(cut_array, trial_clustered, dd, browser=False, outliers=True)

    sum = 0
    for d in cut_clustered:
        if d == clust_to_finetune:
            sum += 1
    print("original sum", sum)
    sum = 0
    for d in trial_clustered:
        if d == clust_to_finetune:
            sum += 1
    print("trial sum", sum)
