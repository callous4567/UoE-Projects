import os
import pickle
import time

import numpy as np
from matplotlib import pyplot as plt, rc
from pandas import DataFrame

import ascii_info
import graphutils
import hdfutils
import windows_directories

writer = hdfutils.hdf5_writer(windows_directories.lamodir, "lamodata.hdf5")
final_table = writer.read_table("lamodata", "astrotable")
print([not d for d in final_table['dr2_radial_velocity'].mask][0:4])