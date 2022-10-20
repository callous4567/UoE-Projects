import galpy.util.bovy_coords
import ascii_info_new
import hdfutils
import windows_directories_new
from galcentricutils_new import galconversion
import numpy as np


"""
panda = hdfutils.hdf5_writer(windows_directories_new.datadir,
                             ascii_info_new.asciiname).read_df(ascii_info_new.fullgroup,
                                                           ascii_info_new.fullpanda)
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                            ascii_info_new.asciiname).read_table(ascii_info_new.fullgroup,
                                                             ascii_info_new.fullset)

# Evaluate energy/etc
hey = energistics_new().pot_eval(table)
"""

"""
We want to test that our conversion of coordinate system via astropy is the same as galpy.
Take the GCs catalogue (small and easy)
Convert that to ICRS 
Convert this to GAL using GALPY
Convert to Galcent using both Galpy and our own work 
"""


# Testing to see that our coord transformation is alright
"""
table = hdfutils.hdf5_writer(windows_directories_new.datadir,
                            ascii_info_new.asciiname).read_table(ascii_info_new.gcs,
                                                             ascii_info_new.set_raw)

# Convert to ICRS. Success.
galcon = galconversion()
galcon.solinfo_grab(windows_directories_new.sourcedir, "solar_info.dat")
galcon.solgal_set()
table2 = galcon.nowrite_GAL_to_ICRS(table)
# Next, take the table and convert it to XYZ (radec->XYZ) using Galpy
XYZ = galpy.util.bovy_coords.lbd_to_XYZ(table2['l'],table2['b'],table2['dist'],degree=True)
# Note about sol_params
"""
#0 is pos, 1 is err, 2 is vel, 3 is err
#x y z respectively.
"""
X_sun = np.abs(galcon.sol_params[0][0])
Z_sun = np.abs(galcon.sol_params[0][2])
xyz = galpy.util.bovy_coords.XYZ_to_galcenrect(*(XYZ.T), X_sun, Z_sun, True)
for i in range(len(table)):
    print(table[i]['x'], xyz[i][0])

# Alright. Our coordinate transformations are on point- aside from the different coordinate system definition
# Astropy uses right-handed. See https://galaxiesbook.org/chapters/A.-Coordinate-systems.html

"""