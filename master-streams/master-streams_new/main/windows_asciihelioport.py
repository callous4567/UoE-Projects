import numpy as np
import asciiutils
import os
import galcentricutils_new
import windows_stack
from windows_builddir_new import dirbuilder
from windows_directories_new import datadir, asciidir, sourcedir
from ascii_info_new import asciiname, set_raw

"""
Version for the new master_streams datasets.
"""

def default():

    # Try to make dir
    dirbuilder()

    # Navigate to asciidir, grab filenames with/without .txt
    os.chdir(asciidir)
    ascii_list = os.listdir()

    # Set up ascii-port objects for all of them.
    asciis = [asciiutils.ascii_port(asciidir, d) for d in ascii_list]

    # Grab the name of the ascii without the .txt extension
    stringnames = [d.replace(".txt", "") for d in ascii_list]

    # Get ascii tables and save them (individually for debug)
    for n, d in enumerate(asciis):
        d.astrotabify()
        d.table['dmu_l'] = d.table['dmu_l'] / np.cos(
            np.radians(d.table['b']))  # because they were stored multiplied by cosdec

        # Set up handler for conversion (using our coordinates)
        galcent = galcentricutils_new.galconversion()

        # Get galcentric transformation
        galcent.solinfo_grab(sourcedir, "solar_info.dat")
        galcent.solgal_set()

        # Convert to Galcentric + ICRS + and Get Momenta (angular.)
        d.table = galcent.nowrite_GAL_to_GALCENT(d.table)
        d.table = galcent.nowrite_GAL_to_ICRS(d.table)
        d.table = galcentricutils_new.angular().get_momentum(d.table)

        # Rename the awkward FeH with feh (to match LAMOST_2MASS_GAIA.py)
        try:
            d.table.rename_column("FeH", "feh")
        except:
            pass

        # Also get the Galcentric Polars
        d.table = galcentricutils_new.angular().get_polar(d.table)
        minradius = 15

        # Remove stars within a certain radius. Originally used 20. Try with 15 next.
        d.table = galcentricutils_new.cluster3d().r_clean(d.table, minimum_radius=minradius)

        # Save Tables
        d.save(datadir, asciiname, stringnames[n], set_raw)

    # Stack all tables (temporarily disable for MSc publication, given that we manually do this with our custom LAMOST)
    windows_stack.stacking().tables()

def GAIA_GAIA_stack():

    # Navigate to asciidir, grab filenames with/without .txt
    os.chdir(asciidir)
    ascii_list = os.listdir()

    # Set up ascii-port objects for all of them.
    asciis = [asciiutils.ascii_port(asciidir, d) for d in ascii_list]

    # Get ascii tables
    for n, d in enumerate(asciis):
        d.astrotabify()
        d.table['dmu_l'] = d.table['dmu_l'] / np.cos(
            np.radians(d.table['b']))  # because they were stored multiplied by cosdec

        # Set up handler for conversion (using our coordinates)
        galcent = galcentricutils_new.galconversion()

        # Get galcentric transformation
        galcent.solinfo_grab(sourcedir, "solar_info.dat")
        galcent.solgal_set()

        # Convert to ICRS. Do not get the momenta. Record source as new text column.
        d.table = galcent.nowrite_GAL_to_ICRS(d.table)
        d.table['source'] = ascii_list[n].split(".")[0]

        # Rename the awkward FeH with feh (to match LAMOST_GAIA-GAIA.py)
        try:
            d.table.rename_column("FeH", "feh")
        except:
            pass

    # Get just the tables
    asciis = [d.table for d in asciis]
    from astropy.table import vstack
    asciis = vstack(asciis)
    return asciis

if __name__ == "__main__":
    default()
