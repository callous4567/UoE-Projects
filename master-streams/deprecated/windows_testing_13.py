import galcentricutils
import ascii_info
import hdfutils
import windows_directories

panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.bhb,
                                                           ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.bhb,
                                                             ascii_info.set_raw)
galcentricutils.cluster3d().hdbs(table, True, ascii_info.bhb_minpar)