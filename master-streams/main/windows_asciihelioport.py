import numpy as np
import asciiutils
import os
import galcentricutils
import windows_stack
from windows_builddir import dirbuilder
from windows_directories import datadir, asciidir, sourcedir
from ascii_info import asciiname, fullgroup, fullset, set_raw

"""
This file will handle import and pre-processing of data from the provided ASCII's we've got.
"""

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

# Try to make dir
dirbuilder()

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
    d.table['dmu_l'] = d.table['dmu_l']/np.cos(np.radians(d.table['b'])) # because they were stored multiplied by cosdec

    # Set up handler for conversion (using our coordinates)
    galcent = galcentricutils.galconversion()

    # Get galcentric transformation
    galcent.solinfo_grab(sourcedir, "solar_info.dat")
    galcent.solgal_set()

    # Convert to Galcentric + ICRS + and Get Momenta (angular.)
    d.table = galcent.nowrite_GAL_to_GALCENT(d.table)
    d.table = galcent.nowrite_GAL_to_ICRS(d.table)
    d.table = galcentricutils.angular().get_momentum(d.table)

    # Rename the awkward FeH with feh (to match LAMOST_2MASS_GAIA.py)
    try:
        d.table.rename_column("FeH", "feh")
    except:
        pass

    # Also get the Galcentric Polars
    d.table = galcentricutils.angular().get_polar(d.table)
    minradius = 15

    # Remove stars within a certain radius. Originally used 20. Try with 15 next.
    d.table = galcentricutils.cluster3d().r_clean(d.table, minimum_radius=minradius)

    # Save Tables
    print(datadir, asciiname, stringnames[n], set_raw)
    d.save(datadir, asciiname, stringnames[n], set_raw)

# Stack all tables (temporarily disabled for MSc publication, given that we manually do this with our custom LAMOST)
windows_stack.stacking().tables()