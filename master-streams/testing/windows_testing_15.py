import pickle
import numpy as np
import ascii_info
import hdfutils
import graphutils
import windows_directories


panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.bhb,
                                                           ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.bhb,
                                                             ascii_info.set_raw)

# Load in the lists as a list of lists!
arrays = [windows_directories.duplimontedir + "\\" + ascii_info.bhb + "\\" + d + ".cluster.txt" \
          for d in ascii_info.duplimonte_saveids]
arrtest = None
with open(arrays[0], 'rb') as f:
    arrtest = pickle.load(f)
graphutils.spec_graph().clust_radec(table, arrtest, cluster_id=0, vasiliev=True)
#graphutils.spec_graph().clust_radec(table, arrtest, cluster_id=False, vasiliev=True)
