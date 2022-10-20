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
from galcentricutils_new import greatcount, compclust, cluster3d
import numpy as np
from hdbscan import HDBSCAN, prediction, flat
import hdfutils
import windows_directories_new
import ascii_info_new


# Suppress Loky Warnings
import warnings
warnings.filterwarnings("ignore", message="Loky-backed parallel loops cannot be called in a multiprocessing, setting n_jobs=1")



if __name__ == "__main__":


    # Parameters for finetune (max_clust_size is adaptive.)
    minsamples_range, \
    min_clust_size, \
    trials_per_samples, \
    minimum_trial_score, \
    minimum_trial_analytically = [4,14], \
                                 ascii_info_new.fulldata_minpar[0], \
                                 20, \
                                 4/5, \
                                 4/5
    finetune_args = minsamples_range, \
                    min_clust_size, \
                    trials_per_samples, \
                    minimum_trial_score, \
                    minimum_trial_analytically


    # set up dependent classes
    gc = greatcount()
    cc = compclust()
    c3d = cluster3d()
    c3d.sizediflim = 2

    # Grab table + greattable
    table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                 ascii_info_new.flatfork_asciiname).read_table(ascii_info_new.fullgroup,
                                                                  ascii_info_new.fullset)

    greattable = hdfutils.hdf5_writer(windows_directories_new.datadir,
                                             ascii_info_new.flatfork_asciiname).read_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles")
    greattable['max_dist'] += 1
    max_clust = np.max(table['prelim_clust'])

    # Get GCC tables/indices/etc
    gcc_tables = []
    clusts_to_finetune = np.arange(0, max_clust+1, 1)
    finetuneargs_lists = []

    # Set up the finetune dir
    try:
        os.mkdir(os.path.join(windows_directories_new.imgdir, "flatfork_finetune"))
    except:
        pass

    # Iterate
    for clust_to_finetune in clusts_to_finetune:

        row = np.where(greattable['clust_to_try'] == clust_to_finetune)[0][0]
        thetaphiwidth = greattable[row][['theta', 'phi', 'max_dist']]
        cut_table = gc.gcc_table_retain(table, *thetaphiwidth)
        cut_table = table[table['greatcount']]
        gcc_tables.append(cut_table)
        finetuneargs_lists.append(finetune_args)

    # Prepare the pool args
    zipped = list(zip(gcc_tables, clusts_to_finetune, finetuneargs_lists))

    # Regularly map/pool :)
    import multiprocessing
    import windows_multiprocessing_new
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing_new.flatfork_do_finetune, zipped)
    pool.close()

    # Array up
    results = np.array(results, int).T
    best_samples, inisums, finsums = results
    try:
        greattable.add_column(0, name="best_samples")
        greattable.add_column(0, name="inisum")
        greattable.add_column(0, name="finsum")
    except:
        pass

    # Modify rows (note that the greattable is not ordered- I am an idiot, I know.)
    for num, clust_to_finetune in enumerate(clusts_to_finetune):
        row_index = np.where(greattable['clust_to_try']==clust_to_finetune)[0][0]
        greattable[row_index]["best_samples"], \
        greattable[row_index]["inisum"], \
        greattable[row_index]["finsum"] = best_samples[num],inisums[num],finsums[num]

    # Write
    greattable['max_dist'] -= 1
    hdfutils.hdf5_writer(windows_directories_new.datadir,
                                             ascii_info_new.flatfork_asciiname).write_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles",
        greattable)
