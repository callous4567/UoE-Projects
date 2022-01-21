import time

import numpy as np
from astropy.table import Table, QTable
from galpy import potential, orbit
from galpy.potential import plotRotcurve, MWPotential2014, PowerSphericalPotentialwCutoff, MiyamotoNagaiPotential, \
    NFWPotential
from matplotlib import pyplot as plt
import astropy.units as u
import astropy.constants.iau2015 as iau
import galcentricutils
from energistics_constants import a_b, a_d, b_d, a_nfw, M_b, M_d, M_nfw, c_nfw
import windows_directories
import scipy.interpolate as interpolate
import scipy.spatial as spatial

# This file handles various energy-orbit related things-
# separated because galpy takes a while to load on Windows.
# Note that we're doing our best to avoid units in our tables-
# If you intend to save a table, while hdfutils supports it,
# SAVE AS A BASIC TABLE WITHOUT AN ASTROPY UNITS INSIDE!!!!!!
# Earlier code was made without considering units inside tables
# Do not go against this.


"""
# A note on Galpy Config:
# It's inside "Users\\Callicious".
# Enable astropy-units == True for this to work at all.
# Another note: galpy functions do NOT take keywords as arguments- input must be in order, without specifying R=R.
# See our GitHub Issue for this- https://github.com/jobovy/galpy/issues/462
# TODO: See about leveraging QTables more- they're basically the new version of Table and better in most ways.
"""


# Generate energy and other associated potential properties for data, including circularity.
# Built for Astropy Tables- the table should have no quantities attached, but otherwise be in standard units.
class energistics(object):
    def __init__(self):
        # Grab the parameters of the solar system in the galactocentric frame- use galcentricutils.
        converter = galcentricutils.galconversion()
        converter.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
        """
        Galpy converts everything into internal units (radii and velocities.) 
        The scale is ro and vo- radius of sol galcentrically, and velocity.
        We're going to use ro as the x-radius in galcentric (cylindrical radius)
        We're going to use vo as the y-velocity in galcentric (cylindrical velocity)
        Note:
        GALPY REALLY REALLY LIKES CYLINDRICAL COORDINATES! :D 
        """
        ro = np.abs(converter.sol_params[0][0])
        vo = np.abs(converter.sol_params[2][1])
        self.rovo = [ro,vo]
        self.zo = np.abs(converter.sol_params[0][2]) # baggage for orbigistics

        # Next up is to set up the potential amplitudes for galpy. Instantiate all this in natural units.
        amp_hernquist = 2*M_b*iau.GM_sun
        amp_miyanagai = M_d*iau.GM_sun
        A_nfw = np.log(1+c_nfw) - (c_nfw/(1+c_nfw)) # see material.
        amp_nfw = M_nfw*iau.GM_sun/A_nfw

        # Set up the potentials...

        # Hernquist Bulge, Miyamoto-Nagai Disk, NFW Halo. 
        self.bulge = potential.HernquistPotential(amp=amp_hernquist,
                                                  a=a_b*u.kpc,
                                                  normalize=False,
                                                  ro=ro*u.kpc,
                                                  vo=vo*u.km/u.s)
        self.disk = potential.MiyamotoNagaiPotential(amp=amp_miyanagai,
                                                     a=a_d*u.kpc,
                                                     b=b_d*u.kpc,
                                                     normalize=False,
                                                     ro=ro*u.kpc,
                                                     vo=vo*u.km/u.s)
        self.halo = potential.NFWPotential(amp=amp_nfw,
                                           a=a_nfw*u.kpc,
                                           normalize=False,
                                           ro=ro*u.kpc,
                                           vo=vo*u.km/u.s)
        self.pot = self.bulge+self.disk+self.halo

        # This is the Bovy 2015 potential, but with our own ro/vo
        # https://docs.galpy.org/en/v1.7.0/reference/potential.html
        bp = PowerSphericalPotentialwCutoff(alpha=1.8, rc=1.9 / ro, normalize=0.05, ro=ro, vo=vo)
        dp = MiyamotoNagaiPotential(a=3. / ro, b=0.28 / ro, normalize=.6, ro=ro, vo=vo)
        hp = NFWPotential(a=16 / ro, normalize=.35, ro=ro,vo=vo)
        self.pot2014 = bp + dp + hp

    # DEPRECATED
    # This is purely for debugging- plot the rotcurve for our potential vs. Bovy 2015 potential.
    def rotcurve(self):
        plotRotcurve(MWPotential2014, label=r'$\mathrm{MWPotential2014 Bovy}$', ro=self.rovo[0],vo=self.rovo[1])
        plotRotcurve(self.pot, overplot=True, label=r'$Model$')
        plotRotcurve(self.pot2014, overplot=True, label=r'Bovy 2014, with mod')
        plt.legend()
        plt.show()

    # DEPRECATED since Orbigistics can do this directly.
    # Evaluate: Potential, Kinetic, Circularity.
    """
    v^2/2 + phi = total 
    """
    def pot_eval(self, table):
        # QTable it- instead of the entire column label being "unit" each value is labelled "unit."
        qtable = QTable(table)
        # Get cylindrical R and Z in units
        qtable['R'] = np.sqrt(qtable['x']**2 + qtable['y']**2) * u.kpc
        qtable['z'] = qtable['z'] * u.kpc
        # Evaluate potential
        table['pot'] = potential.evaluatePotentials(self.pot,
                                                    qtable['R'], #R
                                                    qtable['z'], #Z
                                                    )
        # Evaluate kinetic
        table['kinetic'] = (1/2)*np.sqrt(qtable['vx']**2 +
                                         qtable['vy']**2 +
                                         qtable['vz']**2) * u.km**2/u.s**2
        # Set total
        table['energy'] = table['pot'] + table['kinetic']
        # Return.
        return table

"""
Energistics but manually defined- see notes.
All given positions must be as astropy quantities.
This was just for debugging to ensure that the galpy cooperates- see https://github.com/jobovy/galpy/issues/462
"""
class energistics_manual(object):
    def __init__(self):
        # Grab the parameters of the solar system in the galactocentric frame- use galcentricutils.
        converter = galcentricutils.galconversion()
        converter.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
        """
        Galpy converts everything into internal units (radii and velocities.) 
        The scale is ro and vo- radius of sol galcentrically, and velocity.
        We're going to use ro as the x-radius in galcentric (cylindrical radius)
        We're going to use vo as the y-velocity in galcentric (cylindrical velocity)
        Note:
        GALPY REALLY REALLY LIKES CYLINDRICAL COORDINATES! :D 
        """
        ro = np.abs(converter.sol_params[0][0])
        vo = np.abs(converter.sol_params[2][1])
        self.rovo = [ro,vo]

        # Calculate the value of A_NFW (saves processing time- it's a chonker.)
        A_nfw = np.log(1+c_nfw) - (c_nfw/(1+c_nfw))
        self.ampnfw = -iau.GM_sun*M_nfw/A_nfw # -1*GMvir/A

    # Hernquist Potential.
    def hernquist(self, x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)
        return (-iau.GM_sun*M_b/(r + (a_b*u.kpc))).to(u.km**2/u.s**2)

    # Miyamoto-Nagai
    def nagai(self, x, y, z):
        R_cisq = x**2 + y**2
        Z_circ = np.sqrt(z**2 + (u.kpc * b_d)**2)
        full_R = np.sqrt(R_cisq + (Z_circ + (u.kpc * a_d))**2)
        return (-iau.GM_sun*M_d/full_R).to(u.km**2/u.s**2)

    # Navarro-Frenk-White
    def nfw(self, x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2)
        return ((self.ampnfw*np.log(1 + r/(a_nfw*u.kpc)))/r).to(u.km**2/u.s**2)

    # Full Potential
    def pot(self, x, y, z):
        return self.hernquist(x,y,z) + self.nagai(x,y,z) + self.nfw(x,y,z)

    # Evaluate potentials for Table alongside Total Energy
    """
    v^2/2 + phi = total 
    """
    def pot_eval(self, table):
        qtable = QTable(table)
        qtable['x'],qtable['y'],qtable['z'] = qtable['x']*u.kpc,\
                                              qtable['y']*u.kpc,\
                                              qtable['z']*u.kpc
        # Evaluate potential
        qtable['pot'] = self.pot(qtable['x'],qtable['y'],qtable['z'])
        # Evaluate for the solar system, too.debug
        """
        R_sol = self.rovo[0]
        z_sol = 0 # approx.
        r_sol = np.sqrt(R_sol**2 + z_sol**2)
        solunitstab = Table([[R_sol], [z_sol], [r_sol]], names=['R','z','r'])
        solunitstab['pot'] = potential.evaluatePotentials(self.pot,
                                                          R=solunitstab['R'],
                                                          z=solunitstab['z'],
                                                          t=0) """

        # Evaluate kinetic
        qtable['kinetic'] = (1/2)*np.sqrt(qtable['vx']**2 +
                                          qtable['vy']**2 +
                                          qtable['vz']**2) * u.km**2/u.s**2

        # Set total
        qtable['energy'] = qtable['pot'] + qtable['kinetic']
        fig = plt.figure(figsize=(20,20))
        plt.scatter(qtable['Lz'],qtable['energy'],marker="x",s=2, color='green')
        plt.xlim([-6000,6000])
        plt.ylim([-175000,0])
        plt.xlabel(r'$L_z$')
        plt.ylabel(r'$E$')

        # Als evaluate and overplot using automatic energistics
        table_orig = energistics().pot_eval(table)
        plt.scatter(table_orig['Lz'],table_orig['energy'],marker="x",s=0.5,color='red')

        # Save and show.
        plt.savefig(windows_directories.imgdir + "\\test_energy.png",dpi=600)
        plt.show(dpi=300)

        # Return.
        return table

# Class for dealing with galpy orbits. Inherits from energistics.
# See https://docs.galpy.org/en/v1.7.1/reference/orbit.html for information on handling orbits.
class orbigistics(energistics):
    def __init__(self):
        energistics.__init__(self)

    # Get list of orbit objects for table. Returns list of orbit objects.
    def orbits(self, table):
        # Cylindrify and qtable-up
        table = galcentricutils.angular().get_cylindrical(table)
        qtable = QTable(table)
        qtable['R'],qtable['vR'],qtable['vT'],qtable['z'],qtable['vZ'],qtable['phi'] = qtable['R']*u.kpc,\
                                                                                       qtable['vR']*u.km/u.s,\
                                                                                       qtable['vT']*u.km/u.s,\
                                                                                       qtable['z']*u.kpc,\
                                                                                       qtable['vz']*u.km/u.s,\
                                                                                       qtable['phi']*u.deg
        # Set up orbit objects for each row.
        orbits = orbit.Orbit(vxvv=[qtable['R'],qtable['vR'],qtable['vT'],qtable['z'],qtable['vZ'],qtable['phi']],
                             ro=self.rovo[0]*u.kpc,
                             vo=self.rovo[1]*u.km/u.s,
                             zo=self.zo*u.kpc)

        return orbits

    # Evaluate E, Phi, kin, and circ and dump inside the table, provided the orbits and interpolation range for circ.
    def circularity(self, orbits, table, interpar):
        # Evaluate orbit energies
        table['E'] = orbits.E(pot=self.pot).to(u.km**2 / u.s**2).value \
                     + 0.5*(table['vx']**2 + table['vy']**2 + table['vz']**2)
        table['bound'] = [True if E <= 0 else False for E in table['E']]

        # Set up an interp1d object for potential and vcirc over our range.
        R_vals = np.linspace(interpar[0][0], interpar[0][1], interpar[1])*u.kpc
        E_vals = potential.evaluatePotentials(self.pot, R_vals, 0) + (1/2)*potential.vcirc(self.pot, R_vals)**2
        R_func = interpolate.interp1d(x=E_vals, y=R_vals, kind='cubic', fill_value='extrapolate')


        # Get the R for all our energies (corresponding to circular orbit)
        R_E = R_func(table['E'])*u.kpc
        vcirc_E = potential.vcirc(self.pot, R_E)
        L_E = R_E*vcirc_E

        # Finally, circularity
        circ = (table['Lz']*u.kpc*u.km/u.s)/L_E
        table['circ'] = circ.value
        return table

    # Orbits and Circularity
    def orbilarity(self, table):
        interpar = [[0.1,50], 4000]
        orbits = self.orbits(table)
        circus = self.circularity(orbits, table, interpar)
        return circus


