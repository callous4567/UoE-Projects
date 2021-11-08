import math
import mpmath
import numpy as np
import scipy.stats

# Create the variable grid for pool
from past.builtins import execfile

import ascii_info
import graphutils
import hdfutils
import windows_multiprocessing
import windows_testing
from windows_directories import datadir
import os
import windows_testing





"""
mapped_variables = []
phi, theta_latipolar = 23,175 # latipolar
theta_polar = 90 - theta_latipolar
dtheta = 2.5
radmin = 20
n_phi = 15
index = 1
P = (2*(2*np.pi/n_phi)*np.sin(math.radians(dtheta))) / (4*np.pi)
tableinfo = [ascii_info.fullgroup, ascii_info.fullset]
gccpar = [theta_polar, phi, index, dtheta, radmin, n_phi]
savedex = ("{0:.2f}_{1:.2f}").format(theta_polar,phi)
full_variable = [tableinfo, gccpar, P, savedex]
done = windows_multiprocessing.do_gccsplit_confidence(full_variable)
significance = done[3]
print(significance)
print(np.log10(1/significance))

data = hdfutils.hdf5_writer(datadir, ascii_info.asciiname).read_table(ascii_info.fullgroup, ascii_info.fullset)
graphutils.twod_graph().latipolar(data, theta_polar, phi) """