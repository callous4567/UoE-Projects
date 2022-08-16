import astropy.units

from galpy.orbit import Orbit
import numpy as np
import ascii_info
import hdfutils
import windows_directories
from energistics import orbigistics
import astropy.units as u
"""
Galpy Orbit Fitting using orbit.from_fit, for the data.
"""



# Get LMC from actual catalogue
o= Orbit.from_name('Sgr')
print(o.x(), o.y(), o.z())
print(o.vx(), o.vy(), o.vz(), o.vphi())