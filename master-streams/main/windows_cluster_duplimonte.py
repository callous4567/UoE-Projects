import multiprocessing

import ascii_info, windows_directories
import pickle
import numpy as np
import graphutils
import hdfutils
import galcentricutils

"""
This file will run clustering for all the individual duplimonte'd angulars we have and produce clustered lists (pickled)
"""

"""
# Generate arrayinfo, minpar = arrayinfominpar
# Arrayinfo should be [group, saveid] for the duplimonte: minpar as in cluster3d()
import windows_multiprocessing

arrayinfominpars = []
for group, minpar in zip(ascii_info.all_groups, ascii_info.minpars_allgroups):
    for saveid in ascii_info.duplimonte_saveids:
        arrayinfominpars.append([[group,saveid],minpar])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    results = pool.map(windows_multiprocessing.do_hdbscan, arrayinfominpars)
    pool.close() """

# The above, but for the full group/data. More convenient than dealing with stacking/etc.
import windows_multiprocessing

arrayinfominpars = []
group = ascii_info.fullgroup
minpar = ascii_info.fulldata_minpar
for saveid in ascii_info.duplimonte_saveids:
    arrayinfominpars.append([[group,saveid],minpar])

# Regularly map/pool :)
if __name__ == '__main__':
    pool = multiprocessing.Pool(6)
    results = pool.map(windows_multiprocessing.do_hdbscan, arrayinfominpars)
    pool.close()















# For parameter selection.
"""
data = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                          ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                           ascii_info.set_raw)
#L_ex = [1179, -4450, -1640]
#r, theta, phi = galcentricutils.angular().asSpherical(L_ex)
#indices = galcentricutils.greatcount().gcc_table_indices(t ble, theta, phi, 10, 20)

vec_L = data['vec_L'].tolist()
vec_L = np.array(vec_L)#[indices]
vec_L = vec_L[galcentricutils.cluster3d().getCuboid(vec_L)]
minclustsize = 8
minsamples = 12
minpar = [minclustsize, minsamples]
clustered = galcentricutils.cluster3d().listhdbs(vec_L, minpar)
graphutils.threed_graph().kmeans_L_array(vec_L, clustered, False, True) """

