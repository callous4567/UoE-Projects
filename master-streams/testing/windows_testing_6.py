import pickle
import galcentricutils
import hdfutils
import windows_directories, ascii_info
from astropy.table import Table
import numpy as np
import pygmmis


group = ascii_info.lamostk
panda = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname).read_df(group, ascii_info.panda_raw)


with open(windows_directories.duplimontedir + "\\" + group + "\\L_0.txt", 'rb') as f:
    exdata = pickle.load(f)

exdata_values = np.array(exdata).T
table = Table()
table['Lx'],table['Ly'],table['Lz'] = exdata_values

galcentricutils.cluster3d().hdbs(table, True)

