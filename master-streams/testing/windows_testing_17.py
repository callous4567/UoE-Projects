import ascii_info
import galcentricutils
import hdfutils
import windows_directories
from pandas import DataFrame

table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)
panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.fullgroup,
                                                              ascii_info.fullpanda)
galcentricutils.cluster3d().L_clean(panda, 3)
