import numpy as np
import asciiutils
import os
import galcentricutils_new
from windows_directories_new import asciidir, sourcedir

"""
Version for the new master_streams datasets.
"""

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

