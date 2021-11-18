import ascii_info
import energistics
import hdfutils
import windows_directories

energistics = energistics.energistics()

data = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_df(ascii_info.gcs,
                                                          ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.gcs,
                                                             ascii_info.set_raw)

# Evaluate energy/etc
table = energistics.pot_eval(table)
