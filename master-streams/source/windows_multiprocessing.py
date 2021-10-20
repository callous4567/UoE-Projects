import time
import ascii_info
import galcentricutils
import hdfutils
import windows_directories

"""
Monte-Carlo for a table [group,set]: for multiprocessed Monte-Carlo simulation. Default 4 workers.
"""
# Default monte-carlo parameters for the multiprocessed table. Change externally if needed.
n_monte = 5
sourcecoord = "solar_info.dat"
def do_table(groupset):
    group,set = groupset
    # Grab table
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(group, set)
    # Set up Monte
    sourcedir = windows_directories.sourcedir
    monte = galcentricutils.monte_angular()
    monte.galdefine(sourcedir, sourcecoord)
    # Run Monte on table
    table = monte.table_monte(table,n_monte)
    # Save the table: watch out for if file is already being written to (retry if it fails.)
    while True:
        try:
            writer.write_table(group, set, table)
            break
        except:
            time.sleep(5)
            continue


"""
Example clustering processing to try and guess good DBS parameters for the circle fitting.
"""
def do_dbstest(paramarray):
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    full_data = writer.read_table(ascii_info.fullgroup, ascii_info.fullset)
    cluster = galcentricutils.cluster()
    cluster.dbs(full_data, *paramarray)

"""
Basic multi-processed KMeans routine for great circle counts, 
- generate grid in theta/phi 
- great circle count each grid
- kmeans each great circle table 
- get a score for how good the kmeans is 
- save html for each great circle plot alongside returning the GCC parameter and GCC score 
- decide a score to use: there's sklearn score default, etc.
"""
def do_clustering(vararray):
    tableinfo, gcc_params, clustering_params = vararray
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(*tableinfo)
    # Great Circle the table
    table = galcentricutils.greatcount().gcc_table(table, *gcc_params)
    # Cluster the table
    clustered_table = galcentricutils.cluster().gaussmix(table,*clustering_params)
def do_gaussaicsbics(vararray):
    tableinfo, gcc_params, clustering_params = vararray
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    table = writer.read_table(*tableinfo)
    table = galcentricutils.greatcount().gcc_table(table, *gcc_params)
    galcentricutils.cluster().gaussaicsbics(table,*clustering_params)
