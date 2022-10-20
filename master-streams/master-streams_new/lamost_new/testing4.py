import os

import numpy as np
from astropy.io import fits
from astropy.table import Table

import windows_directories_new

finalstack_dir = os.path.join(windows_directories_new.lamodir, "LAMOST_final.fits")

import graphutils_new
from energistics_new import fast_energistics_new
import hdfutils
import matplotlib.pyplot as plt
import graphutils_new

with fits.open(finalstack_dir, mode='readonly', memmap=True) as hdul:
    # Data (which is an astropy table when loaded.)
    # "gaia_source_id" is dr2_source_id inside this.
    data = Table(hdul[1].data)
    data = data[[True if source == "LAMOST_sstraszak" else False for source in data['source']]]

print(*data[0][['ra','dec']])