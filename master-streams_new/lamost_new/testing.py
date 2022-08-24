import multiprocessing
import os
import numba
import numpy as np
from astropy.table import Table
from astroquery.gaia import Gaia
from numba import njit
from astropy.io import votable, fits

import ascii_info
import hdfutils
import windows_directories
import windows_directories_new
from energistics_new import orbigistics

table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)

import galcentricutils_quality_cuts
indices = galcentricutils_quality_cuts.quality_cuts().fancy_LSR_cut(table, 210)
print(len(indices), len(table))