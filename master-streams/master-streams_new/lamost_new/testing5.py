import os
from astropy.io import fits
from astropy.table import Table
from pandas import DataFrame
import hdfutils
import windows_directories_new

writer = hdfutils.hdf5_writer(windows_directories_new.rootdir, "test_file.hdf5")
writer.create()
dataframe = DataFrame()
writer.write_df("test", "test", dataframe)
