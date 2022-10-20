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
import hdfutils
import windows_directories
import ascii_info


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
                                 ascii_info.fulldata_minpar[0], \
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
    writer = hdfutils.hdf5_writer(windows_directories.datadir,
                                  ascii_info.finetune_asciiname)
    table = writer.read_table(ascii_info.fullgroup,ascii_info.fullset)

    greattable = writer.read_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles"
    )
    greattable['max_dist'] += 1
    clustered = table['prelim_clust']

    # Choose clusts to finetune
    clusts_to_finetune = []
    for clust in list(set(clustered)):

        if clust not in clusts_to_finetune:

            if clust != -1:

                # Only bother if smaller than 10% fraction (i.e. Sgr/etc)
                if len(np.where(clustered==clust)[0])/len(clustered) < 0.1:

                    clusts_to_finetune.append(clust)

    # Get GCC tables/indices/etc
    gcc_tables = []
    finetuneargs_lists = []

    # Set up the finetune dir
    try:
        os.mkdir(os.path.join(windows_directories.imgdir, "finetune_finetune"))
    except:
        pass

    # Iterate
    for clust_to_finetune in clusts_to_finetune:

        row = np.where(greattable['clust_to_try'] == clust_to_finetune)[0][0]
        thetaphiwidth = greattable[row][['theta', 'phi', 'max_dist']]
        cut_table = gc.gcc_table_retain(table, *thetaphiwidth)
        cut_table = cut_table[cut_table['greatcount']]
        gcc_tables.append(cut_table)
        finetuneargs_lists.append(finetune_args)

    # Prepare the pool args
    zipped = list(zip(gcc_tables, clusts_to_finetune, finetuneargs_lists))

    # Regularly map/pool :)
    import multiprocessing
    import windows_multiprocessing
    pool = multiprocessing.Pool(8)
    results = pool.map(windows_multiprocessing.finetune_do_finetune, zipped)
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
    writer.write_table(
        "greatcircles_nonmemberpercent_preliminary",
        "greatcircles",
        greattable
    )
