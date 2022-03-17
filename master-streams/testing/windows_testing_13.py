import ascii_info
import hdfutils
import windows_directories
from energistics import orbigistics
from galpy.util import bovy_coords
import numpy as np

"""
Verify that coordinate transformations we have are identical to Galpy. 
"""

panda = hdfutils.hdf5_writer(windows_directories.datadir,
                             ascii_info.asciiname).read_df(ascii_info.bhb,
                                                           ascii_info.panda_raw)
table = hdfutils.hdf5_writer(windows_directories.datadir,
                            ascii_info.asciiname).read_table(ascii_info.bhb,
                                                             ascii_info.set_raw)

# Attempt our own conversion
orbigist = orbigistics()
R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(table)

# Attempt to use Galpy to get it, too
newcoords = bovy_coords.lbd_to_XYZ(table['l'], table['b'], table['dist'], degree=True)
RR, phiphi, zz = bovy_coords.XYZ_to_galcencyl(*newcoords.T, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T # galcencyl R, phi, z

print(np.max(RR - R))
print(np.max(phiphi - phi))
print(np.max(zz - z))