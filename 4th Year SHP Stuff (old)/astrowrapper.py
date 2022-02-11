# This file must be included with main_runetime.py for the latter to be viable, and kept in the same directory.
# Specify rootdir for the operating directory under which the file structure resides.
# Import Preamble includes packages that may/may not be used, or were used for testing purposes. Do not remove any T_T.
import nexusformat.nexus as nx
import IPython
import astropy
from MeanStars import MeanStars
import random
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import from_table, writeto
from astroquery.gaia import Gaia
import copy
from threading import Thread
import functools
import gc
import yaml
import matplotlib as mpl
from adjustText import adjust_text
from astropy.io.misc import hdf5
from astropy.io import fits
from astropy import wcs as wcs
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table, Column
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import pandas
import h5py
import shutil
from astroquery.astrometry_net import AstrometryNet
from photutils import psf
from astropy.modeling import fitting
from photutils import background
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import SkyCircularAnnulus
from photutils import SkyCircularAperture
from photutils import datasets # Sample data to use for testing
from photutils import IRAFStarFinder
from photutils import DAOStarFinder
from photutils import aperture_photometry
import numpy as np
import matplotlib.pyplot as plt
import sys
import io
import os
import glob
import time
import scipy
import itertools
import multiprocessing
from matplotlib.patches import Patch

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from skimage import transform
from scipy.ndimage.interpolation import shift
# Set AstrometryNet() API key.
ast = AstrometryNet()
ast.api_key = 'rkwvuedpwhdrllyp'
# Set GAIA login
gaia_username = "sstrasza"
gaia_password = "CENT@rino4657898"
# Set the root directory in WINDOWS format. This file is made for usage on PyCharm on a Windows OS.
rootdir = 'D:\Prog\pycharm\SHP2021\Code'
# This is a bunch of info on the file structure, for folders + directories, alongside the h5py file.
"""

DIRECTORIES
===========
Sci -> Cluster Folders + 'STD' Folder -> U B V -> CALIBRATED -> WCSDONE (which has the final files)
Cali -> U B V DARKS BIAS -> masterflatUBV / masterbias / masterdark


h5py/hdf5 "data.hdf5"
=====================
Groups -> Datasets 
------------------
111773/2039 -> depth_u,b,v depth_err_u,b,v alt_dif_u,b,v 

"""

# DISCLAIMER FOR REUSE.
"""
If you intend to reuse this code in part or in full, make reference to the original author and/or source.
Callicious/plasmolian@gmail.com. Proudly coded by Callicious! <3 
"""


# Holds various utilities for photometric analysis for TGP2020 Project.
# INCLUDES CALIBRATION
# INCLUDES OPTICAL DEPTH STUFF
# INCLUDES ASTROMETRY
# INCLUDES FILE UTILITIES (SO FAR)

# Because I'm sick of dealing with a huge "utils" class, I'm writing a new class that draws from utils.
# This is exclusively for alignment of FITS files. Takes "fwhm" and "threshold" as initialization parameter, for DAOPHOT starfinder.
# Threshold is now deprecated (uses static values based on filter band)
# INCLUDES FITS ALIGNMENT PROCESSES


# h5py/hdf5 file-writing utilities for the 'data.hdf5' file in the specified directory.
# Directory is explicit (i.e. D:\TGP2020\etc). Only singular arrays to each.
# INCLUDES HDF5 FILE HANDLING UTILITIES
class hdf5_writer(object):
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename
    # Creates the file if it doesn't exist, does nothing otherwise.
    def create(self):
        os.chdir(self.directory)
        file_created = h5py.File(self.filename, 'a')
        file_created.close()
        os.chdir(self.directory)
    # Returns the list of groups for the file.
    def info(self):
        os.chdir(self.directory)
        f = nx.nxload(self.filename)
        print(f.tree)
    # Writes ARRAY to the DATASET inside GROUP of FILENAME. Can also write SINGLE STRINGS.
    def write(self, group, dataset, array):
        os.chdir(self.directory)
        # First try to create group and add to dataset.
        # Except: The group exists. Access the group instead.
        # Further Except: The group exists, AND the dataset exists. Delete old dataset and retry.

        # In each case, consider whether the data is in the form of a single string or not.
        try:
            with h5py.File(self.filename, 'a') as f:
                grouper = f.create_group(group, 'a')
                if type(array) == str:
                    stringer = array.encode("ascii", "ignore")
                    datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                else:
                    datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)
        except:
            try:
                with h5py.File(self.filename, 'a') as f:
                    grouper = f[group]
                    if type(array) == str:
                        stringer = array.encode("ascii", "ignore")
                        datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                    else:
                        datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)
            except:
                with h5py.File(self.filename, 'a') as f:
                    del f[group][dataset]
                    grouper = f[group]
                    if type(array) == str:
                        stringer = array.encode("ascii", "ignore")
                        datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                    else:
                        datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)


    # Write an astropy table
    def write_table(self, group, dataset, astropytable):
        os.chdir(self.directory)
        astropy.io.misc.hdf5.write_table_hdf5(astropytable, self.filename, path=(group + "/" + dataset), append=True,overwrite=True, serialize_meta=True)

    # Read an astropy table
    def read_table(self, group, tablename):
        os.chdir(self.directory)
        tab = astropy.io.misc.hdf5.read_table_hdf5(input=self.filename, path=(group + "/" + tablename), character_as_bytes=False)
        return tab

    # Reads ARRAY from DATASET inside GROUP of FILENAME
    def read(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a np array.
            data = np.array(f[group][dataset])
        return data

    # Read a string instead, singular, not tested on arrays.
    def read_string(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a string (well, bytes.)
            data = str(np.array(f[group][dataset]), "ascii")
        return data

    # Converts a table to votable xml format + exports. Exportname automatically gets .xml export key.
    def write_votable(self, group, dataset, exportname):
        table = self.read_table(group, dataset)
        votable = from_table(table)
        writeto(votable,exportname + ".xml")

