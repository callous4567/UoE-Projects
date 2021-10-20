import numpy as np
from astropy.coordinates import Galactocentric

import ascii_info
import asciiutils
import os

import galcentricutils
import windows_directories
from windows_directories import datadir, asciidir, sourcedir
from ascii_info import asciiname, fullgroup, fullset, set_raw

# Future notes for using columns
"""
To grab the data:
go under group
grab raw via hdf table reader (this is the astropy table)

To grab the unit type for each column:
go under group
use column id as dataset
read string using this dataset
returns original unit string from ascii file
"""
# ASCII FILE IMPORT!

# Navigate to asciidir, grab filenames with/without .txt
os.chdir(asciidir)
ascii_list = os.listdir()

# Set up ascii-port objects for all of them.
asciis = [asciiutils.ascii_port(asciidir, d) for d in ascii_list]

# Grab the name of the ascii without the .txt extension
stringnames = [d.replace(".txt","") for d in ascii_list]

# Get ascii tables and save them (individually for debug)
for n,d in enumerate(asciis):
    d.astrotabify()
    d.table['dmu_l'] = d.table['dmu_l']/np.cos(np.deg2rad(d.table['b']))
    d.save(datadir, asciiname, stringnames[n], set_raw)


# All this is for conversion stuff.
# Set up handler for conversion
galcent = galcentricutils.galconversion()

# First remove Jorge's conversion from the data, moving back to Galactic
galcent.solinfo_grab(sourcedir, "solar_info_jorge.dat")
galcent.solgal_set()

# For each original ascii
groups = ascii_info.all_groups
for group in groups:
    galcent.GALCENT_to_GAL(windows_directories.datadir,
                           ascii_info.asciiname,
                           group,
                           ascii_info.set_raw)

# For the "full" ascii
galcent.GALCENT_to_GAL(windows_directories.datadir,
                       ascii_info.asciiname,
                       ascii_info.fullgroup,
                       ascii_info.fullset)

# This gives us back the original Galactic data. Next: update Galactocentric using our own conversion (updated values.)
galcent.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
galcent.solgal_set()
# For each original ascii
groups = ascii_info.all_groups
for group in groups:
    galcent.GAL_to_GALCENT(windows_directories.datadir,
                           ascii_info.asciiname,
                           group,
                           ascii_info.set_raw)