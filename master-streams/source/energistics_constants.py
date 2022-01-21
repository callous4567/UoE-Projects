"""
This is the set of constants for energistics to model the bulge, disk, and halo.
Note that these have all been taken from Helmer-Helmi 2019 paper.
https://www.aanda.org/articles/aa/abs/2019/05/aa34769-18/aa34769-18.html
Energistics will use these in calculating amplitudes and etc for using GALPY.

You should note that energistics introduces astropy units/quantities on top of these-
Provide the values here RAW and WITHOUT UNITS.
"""

# Hernquist Bulge. Mass in solmas. Scale in kpc.
M_b = 3e10 # bulge mass
a_b = 0.7 # radial scale

# Miyamoto-Nagai Constants.
M_d = 9.3e10 # disk mass
a_d, b_d = 6.5, 0.26 # thick and thin scale

# NFW Profile Constants.
M_nfw = 1e12 # virial halo mass
c_nfw = 12 # concentration ratio (virial radius divided by scale radius)
a_nfw = 21.5 # radial scale
