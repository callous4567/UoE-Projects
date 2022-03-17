import copy
import os
import pickle

from matplotlib.patches import Patch

import ascii_info
import hdfutils
from galpy.util import bovy_coords

import galcentricutils
import time
import numpy as np
from astropy.table import Table, QTable
from numba import njit, types
from numba.experimental import jitclass
import galpy.potential.Potential
from galpy import potential, orbit
from galpy.potential import plotRotcurve, MWPotential2014, PowerSphericalPotentialwCutoff, MiyamotoNagaiPotential, \
    NFWPotential
from matplotlib import pyplot as plt
import astropy.units as u
import astropy.constants.iau2015 as iau
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

# Important note regarding Galpy coordinate system TODO: IMPORTANT!

"""
The x coordinate galpy takes is negative the x coordinate that astropy takes
The y coordinate is unchanged
The z coordinate is unchanged
(I.e. a left-right clip: galpy left with index toward us, astropy right with index away.)

Consequently:

The theta coordinate is unchanged
The phi coordinate from galpy needs to be subtracted from 180 degrees to be the astropy equivalent



"Specifically, we will define Galactocentric coordinates as the 
left-handed coordinate system that has the Galactic center at 
its center, with x increasing in the disk midplane towards the 
location of the Sun, y increasing in the disk midplane 
perpendicular to x in the direction of Galactic rotation at 
the Sun’s position, and z increasing towards the direction of 
the North Galactic Pole."
"""
# Generate energy and other associated potential properties for data, including circularity.
# Built for Astropy Tables- the table should have no quantities attached, but otherwise be in standard units.
class energistics(object):
    def __init__(self):
        # Grab the parameters of the solar system in the galactocentric frame- use galcentricutils.
        self.converter = galcentricutils.galconversion()
        self.converter.solinfo_grab(windows_directories.sourcedir, "solar_info.dat")
        """
        Galpy converts everything into internal units (radii and velocities.) 
        The scale is ro and vo- radius of sol galcentrically, and velocity.
        We're going to use ro as the x-radius in galcentric (cylindrical radius)
        We're going to use vo as the y-velocity in galcentric (cylindrical velocity)
        Note:
        GALPY REALLY REALLY LIKES CYLINDRICAL COORDINATES! :D 
        Note that ro is negative inside our solar_info.dat file: this is due to astropy being right handed
        """
        ro = np.abs(self.converter.sol_params[0][0])
        vo = np.abs(self.converter.sol_params[2][1])
        self.rovo = [ro,vo] # in units of kpc and kms^-1 cylindrically.
        self.zo = np.abs(self.converter.sol_params[0][2]) # baggage for orbigistics

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
    def pot_eval(self, table, savedexdir):
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
        plt.savefig(windows_directories.imgdir + savedexdir + ".png",dpi=600)
        plt.show(dpi=300)

        # Return.
        return table

# Class for dealing with galpy orbits. Inherits from energistics.
# See https://docs.galpy.org/en/v1.7.1/reference/orbit.html for information on handling orbits.
class orbigistics(energistics):
    def __init__(self):
        energistics.__init__(self)

    # Just get a regular array, non-vectorized, [R, vR...etc]. Phi here are in radians. [-pi,pi.]
    def get_leftgalpy(self, table):
        # Get left-handed cylindrical coordinates for galpy.
        x, y, z, vx, vy, vz = table['x'], \
                              table['y'], \
                              table['z'], \
                              table['vx'], \
                              table['vy'], \
                              table['vz']
        data_array = np.empty((6, len(x)))
        data_array[0, :], data_array[1, :], data_array[2, :], data_array[3, :], data_array[4, :], data_array[5, :] = \
            x, y, z, \
            vx, vy, vz
        R, \
        vR, \
        vT, \
        z, \
        vz, \
        phi = galcentricutils.angular().left_numba_cylindrical(data_array)  # phi in radians.
        #phi = [d + 2 * np.pi if d < 0 else d for d in phi]  # make sure all phi are [0,360]
        #phi = np.array(phi)
        return R, vR, vT, z, vz, phi

    # Just get a qtable with R, vR, vT, z, vz, phi, in the Galpy system. Phi are in degrees and with units. [-pi,pi]
    def get_orbitqtable(self, table):
        # Get left-handed cylindrical coordinates for galpy.
        R, vR, vT, z, vz, phi = self.get_leftgalpy(table)
        phi = np.degrees(phi)

        # Get qtable of results
        qtable = QTable(table)
        qtable['R'], qtable['vR'], qtable['vT'], qtable['z'], qtable['vZ'], qtable['phi'] = R * u.kpc, \
                                                                                            vR * u.km / u.s, \
                                                                                            vT * u.km / u.s, \
                                                                                            z * u.kpc, \
                                                                                            vz * u.km / u.s, \
                                                                                            phi * u.deg
        return qtable

    # Get list of orbit objects for table. Returns list of orbit objects.
    def orbits(self, table):
        # Get left-handed cylindrical coordinates for galpy.
        qtable = self.get_orbitqtable(table)

        # Set up orbit objects for each row.
        orbits = orbit.Orbit(vxvv=[qtable['R'],qtable['vR'],qtable['vT'],qtable['z'],qtable['vZ'],qtable['phi']],
                             ro=self.rovo[0]*u.kpc,
                             vo=self.rovo[1]*u.km/u.s,
                             zo=self.zo*u.kpc)
        # Option to return the qtable, also
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

# Class for fitting Galpy related things. Packed with static methods! Inherits from Angular for cylindrical coords.
class orbifitter():
    def __init__(self):
        self.ang = galcentricutils.angular()

    # Least-squares of 6D data vs. model R, vR, vT, z, vz, phi (phi in radians.)
    @staticmethod
    @njit(fastmath=True)
    def least_squares(data, model):
        # Get all metrics to generate least squares
        length = len(data)
        x, y, z = data[:, 0] * np.cos(data[:, 5]), data[:, 0] * np.sin(data[:, 5]), data[:, 3]
        x_model, y_model, z_model = model[:, 0] * np.cos(model[:, 5]), model[:, 0] * np.sin(model[:, 5]), model[:, 3]
        vR, vT, vz = data[:, 1], data[:, 2], data[:, 4]
        vR_model, vT_model, vz_model = model[:, 1], model[:, 2], model[:, 4]
        # Get leastsq
        leastsq = 0
        for xx, yy, zz, vrr, vtt, vzz in zip(x, y, z, vR, vT, vz):
            spatial_square = (xx - x_model) ** 2 + (yy - y_model) ** 2 + (zz - z_model) ** 2
            velocity_squar = (vrr - vR_model) ** 2 + (vtt - vT_model) ** 2 + (vzz - vz_model) ** 2
            spatial_square /= np.sum(spatial_square) / length
            velocity_squar /= np.sum(velocity_squar) / length
            total_square = spatial_square + velocity_squar
            leastsq += np.min(total_square)
        return leastsq

    # Spatial + Velocity (normalized sum) fitting, weighted by the clustering probabilities.
    @staticmethod
    @njit(fastmath=True)
    def least_multisquares(the_data, the_models, clustering_probability):
        # Define within, since you can't use external (but we want fastmath.)
        def least_squares(data, model, probabilities):
            # Get all metrics to generate least squares
            x, y, z = data[:, 0] * np.cos(data[:, 5]), data[:, 0] * np.sin(data[:, 5]), data[:, 3]
            x_model, y_model, z_model = model[:, 0] * np.cos(model[:, 5]), model[:, 0] * np.sin(model[:, 5]), model[:,3]
            vR, vT, vz = data[:, 1], data[:, 2], data[:, 4]
            vR_model, vT_model, vz_model = model[:, 1], model[:, 2], model[:, 4]
            # Get leastsq
            leastsq = 0
            for xx, yy, zz, vrr, vtt, vzz, probability in zip(x, y, z, vR, vT, vz, probabilities):
                # array of deviations from model for all rr
                spatial_square = (xx - x_model) ** 2 + (yy - y_model) ** 2 + (zz - z_model) ** 2
                # array of deviations for all vv
                velocity_squar = (vrr - vR_model) ** 2 + (vtt - vT_model) ** 2 + (vzz - vz_model) ** 2
                total_square = spatial_square + velocity_squar
                leastsq += probability*np.min(total_square)
            return leastsq
        # Generate least-square list
        least_list = np.empty(len(the_models), dtype=types.float32)
        for i in range(len(the_models)):
            least_list[i] = least_squares(the_data, the_models[i], clustering_probability)
        # Returns.
        argmin = np.argmin(least_list)
        return argmin, least_list

    # The above, but just with velocity- seems to work better: spatial component is letting us down, for some reason.
    # Spatial + Velocity (normalized sum) fitting.
    @staticmethod
    @njit(fastmath=True)
    def least_multisquares_velocityonly(the_data, the_models):
        # Define within, since you can't use external (but we want fastmath.)
        def least_squares(data, model):
            # Get all metrics to generate least squares
            vR, vT, vz = data[:, 1], data[:, 2], data[:, 4]
            vR_model, vT_model, vz_model = model[:, 1], model[:, 2], model[:, 4]
            # Get leastsq
            leastsq = 0
            for vrr, vtt, vzz in zip(vR, vT, vz):
                velocity_squar = (vrr - vR_model) ** 2 + (vtt - vT_model) ** 2 + (vzz - vz_model) ** 2
                leastsq += np.min(velocity_squar)
            return leastsq
        # Generate least-square list
        least_list = np.empty(len(the_models), dtype=types.float32)
        for i in range(len(the_models)):
            least_list[i] = least_squares(the_data, the_models[i])
        # Returns.
        argmin = np.argmin(least_list)
        return argmin, least_list

    # The above, but fitting with position on the sky (i.e. phi, and theta.)
    @staticmethod
    @njit(fastmath=True)
    def least_multisquares_sky(the_data, the_models):
        # Define within, since you can't use external (but we want fastmath.) R, vR, vT, z, vz, phi
        def least_squares(data, model):
            # Rtan(latitude) = z
            latitude = np.arctan(data[:,3]/data[:,0])
            model_latitude = np.arctan(model[:,3]/model[:,0])
            phis = data[:,5]
            model_phi = model[:,5]
            # Convert latitude to spherical polar
            thetas = np.pi/2 - latitude
            model_theta = np.pi/2 - model_latitude
            # Generate spherical unit vectors
            model_univects = np.empty((len(model_latitude), 3))
            model_univects[:,0], model_univects[:,1], model_univects[:,2] = np.cos(model_phi) * np.sin(model_theta), np.sin(model_phi) * np.sin(model_theta), np.cos(model_theta)
            # Get least-sq
            leastsq = 0
            for theta, phi in zip(thetas, phis):
                x,y,z = np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)
                x_dev,y_dev,z_dev = x-model_univects[:,0], y-model_univects[:,1], z-model_univects[:,2]
                deviation = x_dev**2 + y_dev**2 + z_dev**2
                leastsq += np.min(deviation)
            return leastsq
        # Generate least-square list
        least_list = np.empty(len(the_models), dtype=types.float32)
        for i in range(len(the_models)):
            least_list[i] = least_squares(the_data, the_models[i])
        # Returns.
        argmin = np.argmin(least_list)
        return argmin, least_list

    """
    Least Squares Preliminary Fitting.
    - iteration defines number of orbtis to generate. Orphan is good for 4,000
    - time_to_integrate is number of years forward or backward for the orbit. Orphan is good for 0.1e9
    - number_of_steps defines the number of points to generate for the integration. Orphan is good for 1000
    - clust_to_fit defines which cluster you want to fit, from the data. 
    - try_load specifies whether we should try reloading old orbits to save on generation time: saves a lot of time for MC.
    """

    def least_squares_preliminary_orbit_fitting_thing_with_a_long_name(self,
                                                                       table,
                                                                       clust_to_fit,
                                                                       iterations,
                                                                       time_to_integrate,
                                                                       number_of_steps,
                                                                       membership_table=None,
                                                                       clustering=None,
                                                                       try_load=True,
                                                                       graph=False,
                                                                       try_save=True,
                                                                       extra_text=None):
        """
        Method for orbit fitting:
        - Establish the range of orbits that fit the data
        - Set up parameter space (assume gaussians, use covariance matrix and means, generate set.)
        - Monte-carlo generate a set of orbits within this space (n-iterations to test.)
        - Get least-squares estimate of best fit
        This will get you the orbital parameters of the best-fitting orbit for that set.
        Galpy orbit generation is the bottleneck: make sure to save the monte-carlo'd orbits for re-use.
        Once a good orbit has been generated, a markov-chain can be used to search near it to see if an improved orbit exists.
        The same basic monte-carlo'd orbital elements should realistically be usable for all sets of the data, saving comp time.


        :param table: data table
        :param clust_to_fit: which cluster to fit
        :param iterations: number of monte-carlo generations
        :param time_to_integrate: time to integrate for
        :param number_of_steps: resolution of orbit generations
        :param membership_table: the membership table (if applicable: default None)
        :param clustering: the clustering (if applicable: default None)
        :param try_load: should I try to load orbits?
        :param graph: should I produce and save graphs?
        :param try_save: attempt saving the orbits or not? Default True.
        :return: best fitting galpy.orbit.Orbit object.
        """

        # Galactocentric Coordinate Frame Least-squares Orbit Fitting for the data.
        """
        TODO:
        Change fit to fit inside galactic coordinates, instead of ra/dec space. 
        """

        # Set up integration time
        integrate_time = np.linspace(0, time_to_integrate, number_of_steps) * u.yr  # in years.

        # If a membership table exists...
        if membership_table != None:
            # Grab the clustering you're fitting from the membership table and probabilities for weighting.
            clustering = membership_table['probable_clust']
            clustering_probability = membership_table['probability']
        # If the membership table doesn't exist, but a clustering does... assume flat probabilities of unity
        else:
            clustering_probability = np.ones_like(table['l'])

        clustselec = [True if d == clust_to_fit else False for d in clustering]

        # Try to craft specific directory
        try:
            os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables_guess")
        except:
            pass

        if extra_text == None:
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_guess" + "\\" + str(clust_to_fit)
        else:
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_guess" + "\\" + str(clust_to_fit) + "_" + extra_text

        try:
            os.mkdir(savedir)
        except:
            pass

        # Set up orbigistics (for potential) and cylindrification, and also the orbit fitter
        orbigist = orbigistics()
        # Set up the orbit fitter
        orbifitt = orbifitter()

        # Clip data
        data_to_fit = table[clustselec]

        # Attempt to get the data in left-handed galpy, for the fitting.
        try:
            R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
        # Didn't work: the data may not be in cartesian form. Try this.
        except:
            try:
                data_to_fit = orbigist.converter.nowrite_GAL_to_GALCENT(data_to_fit)
                R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
            # Still didn't work. Fudge
            except:
                print("Sorry- broken.")

        # Renormalize our data to Galpy ROVO, since that's what orbits are outputted in.
        R /= orbigist.rovo[0]
        vR /= orbigist.rovo[1]
        vT /= orbigist.rovo[1]
        z /= orbigist.rovo[0]
        vz /= orbigist.rovo[1]
        # Set up table values of these re-normalized values.
        data_to_fit['R_galpy'], \
        data_to_fit['vR_galpy'], \
        data_to_fit['vT_galpy'], \
        data_to_fit['z_galpy'], \
        data_to_fit['vz_galpy'], \
        data_to_fit['phi_galpy'] = R, \
                                   vR, \
                                   vT, \
                                   z, \
                                   vz, \
                                   phi
        data_array = np.empty((6, len(R)))
        data_array[0, :], \
        data_array[1, :], \
        data_array[2, :], \
        data_array[3, :], \
        data_array[4, :], \
        data_array[5, :] = R, vR, vT, z, vz, phi  # replace xyzvxvyvz with R,vR,vT,z,vz,phi (in radians)
        data_array = data_array.T  # convert to vectors. note that we're now in a left-handed system.

        # TODO: Important Note for the Report!!! Defines local solar.
        """
        Use the solarmotion='schoenrich' local solar velocities. 
        (as noted by Drimmel and Poggio as being decent- these are the defaults.)
        solarmotion='schoenrich'
        """

        # Try to load orbits if try_load is true.
        try:
            if try_load == True:
                with open(
                        windows_directories.orbitsdir + "\\" + (
                        "preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_{2:.0f}_timesteps").format(clust_to_fit,
                                                                                                         iterations,
                                                                                                         number_of_steps) + ".txt",
                        'rb') as f:
                    orbit_list = pickle.load(file=f)
                with open(
                        windows_directories.orbitsdir + "\\" + (
                        "preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_elements_{2:.0f}_timesteps").format(
                            clust_to_fit, iterations, number_of_steps)
                                .format(clust_to_fit, iterations) + ".txt", 'rb') as f:
                    orbit_elements = pickle.load(file=f)
                print("Orbits successfully loaded for " + str(clust_to_fit) + " with an iteration number of " + str(
                    iterations))
            else:
                raise ValueError("Not loading orbits, Generating them.")
        # Failed to load them. Generate them.
        except Exception as e:
            if try_load != False:
                print("Orbits failed to load. Generating preliminary Monte orbits.")
            else:
                print(e)

            # Grab the data in galactic coords
            l, b, dist, dmu_l_cosdec, dmu_b, vlos = data_to_fit['l'], \
                                                    data_to_fit['b'], \
                                                    data_to_fit['dist'], \
                                                    data_to_fit['dmu_l'] * np.cos(np.radians(data_to_fit['b'])), \
                                                    data_to_fit['dmu_b'], \
                                                    data_to_fit['vlos']  # deg deg kpc mas/yr mas/yr km/s

            # Set up the means and deviations
            lm, bm, distm, dmulcosdecm, dmubm, vlosm = np.mean(l), np.mean(b), np.mean(dist), \
                                                       np.mean(dmu_l_cosdec), np.mean(dmu_b), np.mean(vlos)
            ldev, bdev, distdev, dmulcosdecdev, dmubdev, vlosdev = l - lm, b - bm, dist - distm, \
                                                                   dmu_l_cosdec - dmulcosdecm, dmu_b - dmubm, vlos - vlosm
            devlist = [ldev, bdev, distdev, dmulcosdecdev, dmubdev, vlosdev]

            # Set up covariance matrix and mean vector, then generate points
            mean_vector = np.array([lm, bm, distm, dmulcosdecm, dmubm, vlosm])
            covtrix = np.empty(shape=(6, 6))
            for i in range(6):
                for j in range(6):
                    covtrix[i, j] = np.mean(devlist[i] * devlist[j])
            orbit_elements = np.random.default_rng().multivariate_normal(mean=mean_vector, cov=covtrix, size=iterations)

            # For all this orbit elements, generate orbit data (list of arrays of orbit vectors: each array is [vec1, vec2...]
            start = time.time()
            orbit_list = []
            for element in orbit_elements:
                element = list(element)
                element[0] *= u.deg
                element[1] *= u.deg
                element[2] *= u.kpc
                element[3] *= u.mas / u.yr
                element[4] *= u.mas / u.yr
                element[5] *= u.km / u.s
                orbit_forward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0] * u.kpc, vo=orbigist.rovo[1] * u.km / u.s,
                                            zo=orbigist.zo * u.kpc, lb=True)
                orbit_backward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0] * u.kpc,
                                             vo=orbigist.rovo[1] * u.km / u.s, zo=orbigist.zo * u.kpc, lb=True)
                orbit_forward.integrate(integrate_time, orbigist.pot, 'rk4_c')
                orbit_backward.integrate(-integrate_time, orbigist.pot, 'rk4_c')
                orbits = np.concatenate([orbit_forward.getOrbit(), orbit_backward.getOrbit()], axis=0)
                orbit_list.append(orbits)
                # Galpy returns are in units of RO and VO from orbigist, here at least.
                # Phi is also [-pi,pi] instead of [0,2pi.] Be aware of this.
                # To rectify this, re-normalize our data to rovo: we took care of this earlier for this reason.

                """
                # Try Plotting (for test.) DEBUG DEBUG DEBUG 
                R, vR, vT, z, vz, phi = orbits.T  # the fitted parameters, which now need to be converter.
                phi = [np.pi - d for d in phi]  # convert galpy left to astropy right
                phi = np.array(phi)
                x, y = R * np.cos(phi), R * np.sin(phi)
                table = Table()
                table['x'], table['y'], table['z'] = x, y, z
                table = galcentricutils.angular().get_polar(table)  # get the fit in galactocentric polars
                plt.scatter(table['theta'], table['phi'], color='green', marker='x', s=0.05)  # plot the fit. """
            end = time.time()
            print(
                "Integration of initial orbits too " + str(end - start) + " seconds, for cluster " + str(clust_to_fit))
            orbit_list = np.array(orbit_list)

            # Try to save the orbits that have been generated
            if try_save == True:
                try:
                    with open(windows_directories.orbitsdir + "\\" + (
                    "preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_{2:.0f}_timesteps").format(clust_to_fit,
                                                                                                     iterations,
                                                                                                     number_of_steps) + ".txt",
                              'wb') as f:
                        pickle.dump(obj=orbit_list, file=f)
                except:
                    pass

                # Also try to save the orbital elements we've generated.
                try:
                    with open(windows_directories.orbitsdir + "\\" + (
                    "preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_elements_{2:.0f}_timesteps").format(clust_to_fit,
                                                                                                              iterations,
                                                                                                              number_of_steps) + ".txt",
                              'wb') as f:
                        pickle.dump(obj=orbit_elements, file=f)
                except:
                    pass

        # Generate preliminary MC least squares fits
        argmin, least_list = orbifitt.least_multisquares(data_array, orbit_list, clustering_probability)

        # Grab the correct orbit element and integrate/etc, get the fit, plot it, etc, in galactocentric polars.
        argmin_element = orbit_elements[argmin]  # l, b, dist, dmu_l_cosdec, dmu_b, vlos
        argmin_element = list(argmin_element)
        argmin_element[0] *= u.deg
        argmin_element[1] *= u.deg
        argmin_element[2] *= u.kpc
        argmin_element[3] *= u.mas / u.yr
        argmin_element[4] *= u.mas / u.yr
        argmin_element[5] *= u.km / u.s
        best_fit = orbit.Orbit(vxvv=argmin_element, ro=orbigist.rovo[0] * u.kpc, vo=orbigist.rovo[1] * u.km / u.s,
                               zo=orbigist.zo * u.kpc, lb=True)

        # Produce graphs if true.
        if graph == True:
            best_fit_reverse = orbit.Orbit(vxvv=argmin_element, ro=orbigist.rovo[0] * u.kpc,
                                           vo=orbigist.rovo[1] * u.km / u.s, zo=orbigist.zo * u.kpc, lb=True)
            best_fit.integrate(integrate_time, orbigist.pot, method='rk4_c')
            best_fit_reverse.integrate(-integrate_time, orbigist.pot, method='rk4_c')
            best_fits_full = np.concatenate([np.flipud(best_fit_reverse.getOrbit()), best_fit.getOrbit()], axis=0)
            R, vR, vT, z, vz, phi = best_fits_full.T  # the fitted parameters, which now need to be converted.
            R *= orbigist.rovo[0]
            vR *= orbigist.rovo[1]
            vT *= orbigist.rovo[1]
            z *= orbigist.rovo[0]
            vz *= orbigist.rovo[1]
            modelfits = [R, vR, vT, z, vz, phi]
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T

            # Main example sky plot in galactic longitude/latitude
            data_array = data_array.T
            Rdata, vRdata, vTdata, zdata, vzdata, phidata = data_array[0, :], \
                                                            data_array[1, :], \
                                                            data_array[2, :], \
                                                            data_array[3, :], \
                                                            data_array[4, :], \
                                                            data_array[5, :]
            Rdata *= orbigist.rovo[0]
            vRdata *= orbigist.rovo[1]
            vTdata *= orbigist.rovo[1]
            zdata *= orbigist.rovo[0]
            vzdata *= orbigist.rovo[1]
            datafits = [Rdata, vRdata, vTdata, zdata, vzdata, phidata]
            Xd, Yd, Zd = bovy_coords.galcencyl_to_XYZ(Rdata, phidata, zdata, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            ld, bd, dd = bovy_coords.XYZ_to_lbd(Xd, Yd, Zd, degree=True).T
            plt.scatter(ld, bd, color='black', marker='x')
            plt.scatter(l, b, color='red', marker='x', s=0.5)  # plot the fit.
            plt.xlabel('l')
            plt.ylabel('b')
            plt.savefig(savedir + "\\" + str(clust_to_fit) + "_all_orbits_testgalactic.png", dpi=300)
            plt.close()

            # Set up permutations for comparative plots
            permutations = [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5], [1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [2, 4],
                            [2, 5], [3, 4], [3, 5], [4, 5]]
            titles = ["R vs. vR",
                      "R vs. vT",
                      "R vs. z",
                      "R vs. vz",
                      "R vs. phi",
                      "vR vs. vT",
                      "vR vs. z",
                      "vR vs, vz",
                      "vR vs. phi",
                      "vT vs. z",
                      "vT vs. vz",
                      "vT vs. phi",
                      "z vs. vz",
                      "z vs. phi",
                      "vz vs. phi"]

            # Create 1-1 plots to demonstrate error directions
            for permutation, title in zip(permutations, titles):
                x, y = permutation
                plt.scatter(datafits[x], datafits[y], color='black', marker='x', s=0.1)
                plt.scatter(modelfits[x], modelfits[y], color='red', marker='x', s=0.1)
                plt.title(title)
                plt.savefig(savedir + "\\" + str(clust_to_fit) + "_" + title.replace(" ", "_") + ".png", dpi=300)
                plt.close()

            # Return the best-fitting orbit object (to use as a "base" for the galpy orbit fitting suite.)
            return best_fit
        # Else, just return the best_fit.
        else:
            return best_fit

    """
    Galpy Final Fitting
    - Use least_squares as initial fit/init_vxvv
    - Generate lots of plots
    - profit.
    - Needs membership table.
    """
    def galpy_final_fitting(self,
                            table,
                            clust_to_fit,
                            iterations,
                            time_to_integrate,
                            number_of_steps,
                            try_load=True,
                            graph=False,
                            load_fit=False,
                            try_save=False,
                            debug_graph=None):

        """
        :param table: The data table, astropy or pandas.
        :param clust_to_fit: Which cluster do you want to fit from percent_table_greatfitted
        :param iterations: How many monte-carlo iterations to generate for least-squares guess
        :param time_to_integrate: In years, how long to forward/backward tail for least-squares (and graphs?)
        :param number_of_steps: Resolution of the least-squares fit
        :param try_load: Should I try loading the Monte-Carlo orbits that already exist to save time?
        :param graph: Should I produce graphing?
        :param load_fit: Should I try loading the previous orbit fit?
        :param try_save: Should I try saving the fit?
        :return: The orbit fit
        """
        # Load in the "mean" percenttable and map
        writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
        membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")

        # Get an initial guess for the orbit with out code (will also produce plots for the orbit- see folder.)
        guess = self.least_squares_preliminary_orbit_fitting_thing_with_a_long_name(table,
                                                                                    clust_to_fit,
                                                                                    iterations,
                                                                                    time_to_integrate,
                                                                                    number_of_steps,
                                                                                    membership_table=membership_table,
                                                                                    try_load=try_load,
                                                                                    graph=graph,
                                                                                    try_save=True)

        # Get clustering data and select for relevant data.
        clustering = membership_table['probable_clust']
        clustselec = [True if d == clust_to_fit else False for d in clustering]
        data_to_fit = table[clustselec]

        # Set up orbigistics (for potential and various parameters related to it/rovozo.)
        orbigist = orbigistics()

        # Use the solarmotion='schoenrich' local solar velocities
        # (as noted by Drimmel and Poggio as being decent- these are the defaults in Galpy anyway but we're specifying.)
        solarmotion = 'schoenrich'

        # Set up the position vectors for galpy, which will use galactic coordinates straight-up, and also use errors
        vxvv = np.empty(shape=(len(data_to_fit), 6))
        vxvv[:, 0], vxvv[:, 1], vxvv[:, 2], vxvv[:, 3], vxvv[:, 4], vxvv[:, 5] = data_to_fit['l'], \
                                                                                 data_to_fit['b'], \
                                                                                 data_to_fit['dist'], \
                                                                                 data_to_fit['dmu_l'] * np.cos(
                                                                                     np.radians(data_to_fit['b'])), \
                                                                                 data_to_fit['dmu_b'], \
                                                                                 data_to_fit['vlos']

        # Also set up the errors. Don't forget: edmu is multiplied by cos(dec) in Galpy.
        evxvv = np.empty(shape=(len(data_to_fit), 6))
        evxvv[:, 0], \
        evxvv[:, 1], \
        evxvv[:, 2], \
        evxvv[:, 3], \
        evxvv[:, 4], \
        evxvv[:, 5] = 0.05 * np.ones_like(data_to_fit['l']), \
                      0.05 * np.ones_like(data_to_fit['l']), \
                      data_to_fit['edist'], \
                      data_to_fit['edmu_l'] * np.cos(np.radians(data_to_fit['b'])), \
                      data_to_fit['edmu_b'], \
                      data_to_fit['evlost']  # Jo recommended using a few arcmins of error for l/b

        # Set up the initial orbit guess. Per Bovy, I'm just guessing that pmll here already has the cosdec factor... guessing.
        init_vxvv = [guess.ll(ro=orbigist.rovo[0]),
                     guess.bb(ro=orbigist.rovo[0]),
                     guess.dist(ro=orbigist.rovo[0]),
                     guess.pmll(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1]),
                     guess.pmbb(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1]),
                     guess.vlos(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1])]

        if load_fit == False:
            # Run and save the fit
            fit = galpy.orbit.Orbit.from_fit(init_vxvv,
                                             vxvv,
                                             evxvv,
                                             orbigist.pot,
                                             ro=orbigist.rovo[0] * u.kpc,
                                             vo=orbigist.rovo[1] * u.km / u.s,
                                             zo=orbigist.zo * u.kpc,
                                             solarmotion=solarmotion,
                                             lb=True)
            if try_save == True:
                with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt",
                          'wb') as f:
                    pickle.dump(obj=fit, file=f)
        else:
            # Load the fit
            with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt",
                      'rb') as f:
                fit = pickle.load(file=f)

        # If Graph is True, generate graphs
        if graph == True:

            # Get a "flip" for backward orbit
            fit_backward = copy.copy(fit)

            # Get integrals
            fit.integrate(np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')
            fit_backward.integrate(-np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')

            # Get orbits
            orbits = fit.getOrbit()
            orbits_backward = fit_backward.getOrbit()
            orbits_full = np.concatenate((orbits, orbits_backward), axis=0)
            R, vR, vT, z, vz, phi = orbits_full.T  # the fitted parameters, in units of rovo/etc.
            R *= orbigist.rovo[0]
            z *= orbigist.rovo[0]
            x, y = R * np.cos(phi), R * np.sin(phi)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T

            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T

            fig = plt.figure()
            plt.grid(which="major", color="pink")
            plt.scatter(l, b, color='red', marker='x', s=0.1)

            # Clip data and get galpy cylindrical coordinates, alongside getting data as an array.
            R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T

            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T

            plt.scatter(data_to_fit['l'], data_to_fit['b'], color='black', marker='x')

            # Generate the range for the plotting.
            xlength = np.max(data_to_fit['l']) - np.min(data_to_fit['l'])
            ylength = np.max(data_to_fit['b']) - np.min(data_to_fit['b'])
            xlim = [np.min(data_to_fit['l']) - 0.05 * xlength, np.max(data_to_fit['l']) + 0.05 * xlength]
            ylim = [np.min(data_to_fit['b']) - 0.05 * ylength, np.max(data_to_fit['b']) + 0.05 * ylength]
            plt.xlim(xlim)
            plt.ylim(ylim)
            legend_elements = [Patch(edgecolor='gray', facecolor='red', label='Fit'),
                               Patch(edgecolor='gray', facecolor='black', label='Data')]
            plt.legend(handles=legend_elements, loc='upper right')

            try:
                os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables")
            except:
                pass
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables" + "\\" + str(clust_to_fit)
            try:
                os.mkdir(savedir)
            except:
                pass
            plt.savefig(savedir + "\\" + str(clust_to_fit) + "_final_galpy_fit.png", dpi=300)
            plt.xlabel('l')
            plt.ylabel("b")
            plt.show()

        # If debug_graph isn't None, generate graphs for debug
        if debug_graph != None:

            # Get a "flip" for backward orbit
            fit_backward = copy.copy(fit)

            # Get integrals
            fit.integrate(np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')
            fit_backward.integrate(-np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')

            # Get orbits
            orbits = fit.getOrbit()
            orbits_backward = fit_backward.getOrbit()
            orbits_full = np.concatenate((orbits, orbits_backward), axis=0)
            R, vR, vT, z, vz, phi = orbits_full.T  # the fitted parameters, in units of rovo/etc.
            R *= orbigist.rovo[0]
            z *= orbigist.rovo[0]
            x, y = R * np.cos(phi), R * np.sin(phi)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
            RR, phiphi = copy.deepcopy(R), copy.deepcopy(phi)

            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T

            fig = plt.figure()
            plt.grid(which="major", color="pink")
            plt.scatter(l, b, color='red', marker='x', s=0.1)

            # Clip data and get galpy cylindrical coordinates, alongside getting data as an array.
            try:
                R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
            except:
                data_to_fit = orbigist.converter.nowrite_GAL_to_GALCENT(data_to_fit)
                R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T

            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T

            plt.scatter(data_to_fit['l'], data_to_fit['b'], color='black', marker='x')

            # Generate the range for the plotting.
            xlength = np.max(data_to_fit['l']) - np.min(data_to_fit['l'])
            ylength = np.max(data_to_fit['b']) - np.min(data_to_fit['b'])
            xlim = [np.min(data_to_fit['l']) - 0.05 * xlength, np.max(data_to_fit['l']) + 0.05 * xlength]
            ylim = [np.min(data_to_fit['b']) - 0.05 * ylength, np.max(data_to_fit['b']) + 0.05 * ylength]
            plt.xlim(xlim)
            plt.ylim(ylim)
            legend_elements = [Patch(edgecolor='gray', facecolor='red', label='Fit'),
                               Patch(edgecolor='gray', facecolor='black', label='Data')]
            plt.legend(handles=legend_elements, loc='upper right')

            try:
                os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables_debugging")
            except:
                pass
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_debugging" + "\\" + str(clust_to_fit)
            try:
                os.mkdir(savedir)
            except:
                pass

            plt.xlabel('l')
            plt.ylabel("b")
            plt.savefig(savedir + "\\" + str(debug_graph) + "_final_galpy_fit.png", dpi=300)

            # Also produce for R vs. phi (data)
            plt.close()
            fig = plt.figure()

            # Throw in the fit
            plt.scatter(RR, phiphi, color='red', marker='x')

            # And the data
            plt.scatter(R, phi, color='black', marker='x')

            # Generate the range for the plotting.
            xlength = np.max(R) - np.min(R)
            ylength = np.max(phi) - np.min(phi)
            xlim = [np.min(R) - 0.05 * xlength, np.max(R) + 0.05 * xlength]
            ylim = [np.min(phi) - 0.05 * ylength, np.max(phi) + 0.05 * ylength]
            plt.xlim(xlim)
            plt.ylim(ylim)
            legend_elements = [Patch(edgecolor='gray', facecolor='red', label='Fit'),
                               Patch(edgecolor='gray', facecolor='black', label='Data')]
            plt.legend(handles=legend_elements, loc='upper right')
            plt.xlabel('R')
            plt.ylabel(r'$\phi$')
            plt.grid(which='major', color='pink')
            plt.savefig(savedir + "\\" + str(debug_graph) + "_final_galpy_fit_RPHI.png", dpi=300)
        return fit

    """
    Galpy Final Fitting
    - Use least_squares as initial fit/init_vxvv
    - Generate lots of plots
    - profit.
    - Doesn't need membership table.
    - Fairly deprecated but whatever. 
    """
    def galpy_fitting_nomemb(self,
                             table,
                             clustering,
                             clust_to_fit,
                             iterations,
                             time_to_integrate,
                             number_of_steps,
                             try_load=True,
                             graph=False,
                             load_fit=False,
                             try_save=False,
                             extra_text="maindata_guess"):
        """
        :param table: the data table
        :param clustering: the clustering
        :param clust_to_fit: Which cluster do you want to fit from percent_table_greatfitted
        :param iterations: How many monte-carlo iterations to generate for least-squares guess
        :param time_to_integrate: In years, how long to forward/backward tail for least-squares (and graphs?)
        :param number_of_steps: Resolution of the least-squares fit
        :param try_load: Should I try loading the Monte-Carlo orbits that already exist to save time? Default True.
        :param graph: Should I produce graphing? Default False.
        :param load_fit: Should I load a previous fit? Default False.
        :param try_save: Should I try saving the fit?
        :return: The orbit fit
        """

        # Get an initial guess for the orbit with out code (will also produce plots for the orbit- see folder.)
        guess = self.least_squares_preliminary_orbit_fitting_thing_with_a_long_name(table,
                                                                                    clust_to_fit,
                                                                                    iterations,
                                                                                    time_to_integrate,
                                                                                    number_of_steps,
                                                                                    clustering=clustering,
                                                                                    try_load=try_load,
                                                                                    graph=graph,
                                                                                    try_save=True,
                                                                                    extra_text=extra_text)

        # Get clustering data and select for relevant data.
        clustselec = [True if d == clust_to_fit else False for d in clustering]
        data_to_fit = table[clustselec]

        # Set up orbigistics (for potential and various parameters related to it/rovozo.)
        orbigist = orbigistics()

        # Use the solarmotion='schoenrich' local solar velocities
        # (as noted by Drimmel and Poggio as being decent- these are the defaults in Galpy anyway but we're specifying.)
        solarmotion = 'schoenrich'

        # Set up the position vectors for galpy, which will use galactic coordinates straight-up, and also use errors
        vxvv = np.empty(shape=(len(data_to_fit), 6))
        vxvv[:, 0], vxvv[:, 1], vxvv[:, 2], vxvv[:, 3], vxvv[:, 4], vxvv[:, 5] = data_to_fit['l'], \
                                                                                 data_to_fit['b'], \
                                                                                 data_to_fit['dist'], \
                                                                                 data_to_fit['dmu_l'] * np.cos(
                                                                                     np.radians(data_to_fit['b'])), \
                                                                                 data_to_fit['dmu_b'], \
                                                                                 data_to_fit['vlos']

        # Also set up the errors. Don't forget: edmu is multiplied by cos(dec) in Galpy.
        evxvv = np.empty(shape=(len(data_to_fit), 6))
        evxvv[:, 0], \
        evxvv[:, 1], \
        evxvv[:, 2], \
        evxvv[:, 3], \
        evxvv[:, 4], \
        evxvv[:, 5] = 0.05 * np.ones_like(data_to_fit['l']), \
                      0.05 * np.ones_like(data_to_fit['l']), \
                      data_to_fit['edist'], \
                      data_to_fit['edmu_l'] * np.cos(np.radians(data_to_fit['b'])), \
                      data_to_fit['edmu_b'], \
                      data_to_fit['evlost']  # Jo recommended using a few arcmins of error for l/b

        # Set up the initial orbit guess. Per Bovy, I'm just guessing that pmll here already has the cosdec factor... guessing.
        init_vxvv = [guess.ll(ro=orbigist.rovo[0]),
                     guess.bb(ro=orbigist.rovo[0]),
                     guess.dist(ro=orbigist.rovo[0]),
                     guess.pmll(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1]),
                     guess.pmbb(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1]),
                     guess.vlos(ro=orbigist.rovo[0],
                                vo=orbigist.rovo[1])]

        if load_fit == False:
            # Run and save the fit
            fit = galpy.orbit.Orbit.from_fit(init_vxvv,
                                             vxvv,
                                             evxvv,
                                             orbigist.pot,
                                             ro=orbigist.rovo[0] * u.kpc,
                                             vo=orbigist.rovo[1] * u.km / u.s,
                                             zo=orbigist.zo * u.kpc,
                                             solarmotion=solarmotion,
                                             lb=True)
            if try_save == True:
                with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + "_maindata_galpy.orbit_fit.txt",
                          'wb') as f:
                    pickle.dump(obj=fit, file=f)
        else:
            # Load the fit
            with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + "_maindata_galpy.orbit_fit.txt",
                      'rb') as f:
                fit = pickle.load(file=f)

        # If Graph is true, generate graphs
        if graph == True:
            # Get a "flip" for backward orbit
            fit_backward = copy.copy(fit)

            # Get integrals
            fit.integrate(np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')
            fit_backward.integrate(-np.linspace(0, time_to_integrate, number_of_steps) * u.yr, orbigist.pot, 'rk4_c')
            # Get orbits
            orbits = fit.getOrbit()
            orbits_backward = fit_backward.getOrbit()
            orbits_full = np.concatenate((orbits, orbits_backward), axis=0)
            R, vR, vT, z, vz, phi = orbits_full.T  # the fitted parameters, in units of rovo/etc.
            R *= orbigist.rovo[0]
            z *= orbigist.rovo[0]
            x, y = R * np.cos(phi), R * np.sin(phi)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
            fig = plt.figure()
            plt.grid(which="major", color="pink")
            plt.scatter(l, b, color='red', marker='x', s=0.1)

            # Clip data and get galpy cylindrical coordinates, alongside getting data as an array.
            R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
            X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
            l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
            # r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
            plt.scatter(data_to_fit['l'], data_to_fit['b'], color='black', marker='x')

            # Generate the range for the plotting.
            xlength = np.max(data_to_fit['l']) - np.min(data_to_fit['l'])
            ylength = np.max(data_to_fit['b']) - np.min(data_to_fit['b'])
            xlim = [np.min(data_to_fit['l']) - 0.05 * xlength, np.max(data_to_fit['l']) + 0.05 * xlength]
            ylim = [np.min(data_to_fit['b']) - 0.05 * ylength, np.max(data_to_fit['b']) + 0.05 * ylength]
            plt.xlim(xlim)
            plt.ylim(ylim)
            legend_elements = [Patch(edgecolor='gray', facecolor='red', label='Fit'),
                               Patch(edgecolor='gray', facecolor='black', label='Data')]
            plt.legend(handles=legend_elements, loc='upper right')

            try:
                os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables_maindata")
            except:
                pass
            savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_maindata" + "\\" + str(clust_to_fit)
            try:
                os.mkdir(savedir)
            except:
                pass
            plt.savefig(savedir + "\\" + str(clust_to_fit) + "_maindata_galpy_fit.png", dpi=300)
            plt.xlabel('l')
            plt.ylabel("b")
            plt.show()

        return fit

    # Get the means and errors for a set of orbits provided to it.
    def orbistatistics(self, orbits, integration_time, number_of_steps, savedir, save_unique, method='dopr54_c'):

        # Set up orbigist
        orbigist = orbigistics()

        # Lists
        ees, periggs = [],[]
        EEs = []
        Lzs = []

        # For the orbits, integrate them
        for orbit in orbits:

            # Set up a copy for the forward and the backward
            forward = copy.deepcopy(orbit)
            backward = copy.deepcopy(orbit)

            # Do integrals and append these orbit objects to lists
            forward.integrate(t=np.linspace(0, integration_time, number_of_steps)*u.yr,
                              pot=orbigist.pot,
                              method=method)
            backward.integrate(t=-1*np.linspace(0, integration_time, number_of_steps)*u.yr,
                               pot=orbigist.pot,
                               method=method)

            # Grab the eccentricity
            e = forward.e(analytic=True, pot=orbigist.pot)
            ees.append(e)

            # Get the radii involved
            for_R, back_R = np.array([forward.R(t).to(u.kpc).value for t in np.linspace(0, integration_time, number_of_steps)*u.yr]),\
                            np.array([backward.R(t).to(u.kpc).value for t in -1*np.linspace(0, integration_time, number_of_steps)*u.yr])
            for_z, back_z = np.array([forward.z(t).to(u.kpc).value for t in np.linspace(0, integration_time, number_of_steps)*u.yr]),\
                            np.array([backward.z(t).to(u.kpc).value for t in -1*np.linspace(0, integration_time, number_of_steps)*u.yr])

            # Concatenate
            R_concat = np.concatenate([for_R, back_R])
            z_concat = np.concatenate([for_z, back_z])

            # Get the radius
            r = np.sqrt(R_concat**2 + z_concat**2)

            # Grab pergalacticons/closest approaches.
            perigalacticon = np.min(r)
            periggs.append(perigalacticon)

            # Also evaluate the energy (conserved after all) and angular momentum in Z
            E, Lz = orbit.E(0*u.yr, orbigist.pot, orbigist.rovo[1]*u.km/u.s), orbit.Lz(0*u.yr, orbigist.rovo[0],
                                                                                       orbigist.rovo[1])
            EEs.append(E.value)
            Lzs.append(Lz.value)

        # Get mean/stdev
        meanee, stdee = np.mean(ees), np.std(ees)
        meanpg, stdpg = np.mean(periggs), np.std(periggs)
        meanEE, stdEE = np.mean(EEs), np.std(EEs)
        meanLz, stdLz = np.mean(Lzs), np.std(Lzs)

        # Create basic plots
        return ees, meanee, stdee, periggs, meanpg, stdpg, EEs, meanEE, stdEE, Lzs, meanLz, stdLz
