import os
import pickle

from matplotlib import pyplot as plt, rc
from pandas import DataFrame

import ascii_info
import hdfutils
import windows_directories


# Writer
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)

# Grab entropy table (final)
entropy_table = writer.read_table("entropy_fit_witherrors", "per_cluster_haloonly")[['cluster','M_nfw','M_nfw_err']]

entropy_pandas = entropy_table.to_pandas(index=False)
entropy_pandas.to_latex(caption="placeholder",
                        buf=windows_directories.imgdir + "\\per_cluster_haloonly.tex",
                        label="tab:per_cluster_haloonly",
                        float_format="%.2f",
                        index=False)
