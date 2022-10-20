import os
import pickle
import time

import numpy as np
from matplotlib import pyplot as plt, rc
from pandas import DataFrame

import ascii_info_new
import graphutils_new
import hdfutils
import windows_directories_new

writer = hdfutils.hdf5_writer(windows_directories_new.lamodir, "lamodata.hdf5")
final_table = writer.read_table("lamodata", "astrotable")
print([not d for d in final_table['dr2_radial_velocity'].mask][0:4])