import galpy.potential as potential
import astropy.units as u
import astropy.constants.iau2015 as iau
import numpy as np

# Hernquist Bulge. Mass in solmas. Scale in kpc.
M_b = 3e10 # bulge mass
a_b = 0.7 # radial scale

# Specify custom ro/vo
ro = 8.178
vo = 245.6

# Example coordinate to test
x, y, z = 1, 2, 3
xx, yy ,zz = x*u.kpc, y*u.kpc, z*u.kpc

# Cylindrically- galpy doesn't seem to want quantities for R and z.
R = np.sqrt(x**2 + y**2)*u.kpc

# Hernquist Potential. Takes quantities.
def hernquist(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return (-iau.GM_sun*M_b/(r + (a_b*u.kpc))).to(u.km**2/u.s**2)

# Galpy Potential.
amp_hernquist = 2*M_b*iau.GM_sun
bulge = potential.HernquistPotential(amp=amp_hernquist,
                                     a=a_b*u.kpc,
                                     normalize=False,
                                     ro=ro*u.kpc,
                                     vo=vo*u.km/u.s)

# Evaluate using manual/Galpy
manual = hernquist(xx, yy, zz)
galpy = potential.evaluatePotentials(bulge, R, zz)
print(manual, galpy)


