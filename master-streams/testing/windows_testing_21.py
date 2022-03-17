import pickle
import time
import timeit

import astropy.units as u
import galpy.orbit
import numpy as np
from astropy.table import Table

import graphutils
from galpy import orbit
from matplotlib import pyplot as plt
import ascii_info
import galcentricutils
import hdfutils
import windows_directories
from energistics import orbigistics
from galcentricutils import compclust

"""
Method for orbit fitting:
- Establish the range of orbits that fit the data
- Set up parameter space (min-max) of this
- Monte-carlo generate a set of orbits within this space (n-iterations to test.) 
- Get least-squares estimate of best fit
This will get you the orbital parameters of the best-fitting orbit for that set. 
Galpy orbit generation is the bottleneck: make sure to save the monte-carlo'd orbits for re-use. 
Once a good orbit has been generated, a markov-chain can be used to search near it to see if an improved orbit exists.
The same basic monte-carlo'd orbital elements should realistically be usable for all sets of the data, saving comp time.
"""

# Load member table (of the members who survived the initial monte-carlo'ing)
writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
membership_table = writer.read_table(ascii_info.fullgroup, "percent_table")

# Grab the data table that you are fitting
table = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                              ascii_info.fullset)

# Test
grapher = graphutils.twod_graph()
grapher.hist_fancy(table['vlos'], 100, r'$\textrm{vlos} / \textrm{kms}^{-1}$', r'$\rho$', savepath=None, dencurve=False)