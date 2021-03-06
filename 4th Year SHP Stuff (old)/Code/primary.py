# Code for SHP
# Sebastian Straszak, 2021

# Menagerie of imports of varying use.
import pickle

import nexusformat.nexus as nx
import IPython
import astropy
from MeanStars import MeanStars
import random
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import from_table, writeto
from astroquery.gaia import Gaia
import copy
from threading import Thread
import functools
import gc
import yaml
import matplotlib as mpl
from adjustText import adjust_text
from astropy.constants import iau2015 as astroconst
from astropy.constants import codata2018 as codaconst
from astropy.io.misc import hdf5
from astropy.io import fits
from astropy import wcs as wcs
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table, Column
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import pandas
import h5py
import shutil
from astroquery.astrometry_net import AstrometryNet
from photutils import psf
from astropy.modeling import fitting
from photutils import background
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import SkyCircularAnnulus
from photutils import SkyCircularAperture
from photutils import datasets # Sample data to use for testing
from photutils import IRAFStarFinder
from photutils import DAOStarFinder
from photutils import aperture_photometry
import numpy as np
import matplotlib.pyplot as plt
import sys
import io
import os
import glob
import time
import scipy
import scipy.constants
import itertools
import multiprocessing
from matplotlib.patches import Patch
from mpmath import mpf

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from skimage import transform
from scipy.ndimage.interpolation import shift

from SHP2021.Code import astrowrapper
font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}

mpl.rc('font', **font)

filer = astrowrapper.hdf5_writer("D:\Prog\pycharm\SHP2021\Code", "shp.hdf5")

# ALL THE PARAMETERS THAT ONE MIGHT NEED! Give in the form of "density parameters" since all surveys quote like that.
# Current values are random picked from book. If you give as raw densities, then alter the code to to account for this.
# 0.31, 9E-5, 0.69
# Grab scalar field stuff


# We're moving all definition/etc down below deck.
"""
dens_mat = 0.3
dens_rad = 0
dens_vac = 0.7
den_array = [dens_mat, dens_vac, dens_rad]
w_array = [0, -1, 1/3]
hub = 68 # Kms^-1Mpc^-1
current_age = 13.7e9 # Years
year = 365*86400
mpc = 3.0857*1e22
G = 6.6743e-11
w_inteconst = 4*np.pi*G * (year**2) / 3 # 4 * pi * G / 3 (in the form of year^-2) for integration routines."""
# Universe class. Follows LCDM model. Initialize with start parameters. Has functions for stepping forward/backward (in years).
# Rho should be in units of kgm^-3 and H0 in kms^-1Mpc^-1.
# Default R is 1
# X is currently in GIGAYEARS. Operational hubble in gigayear*s^-1
class Universe():
    def __init__(self, dV,dR ,dM, H0, age):
        self.Ov, self.Or, self.Om = dV, dR, dM
        self.H = H0
        self.H0 = H0
        self.addconst = 0 # Use const_calc to change if uni is curved
        self.R = 1
        self.age = age # Gy
    # Converts between the forms of h0 (year vs integrator units): does not alter self.value (needs to be done manually.)
    def converter(self, to_operator):
        if to_operator == True: # convert inputs into the operational format for integrator/etc
            h0 = self.H * (1e9) * year * 1e3 / mpc # In units of (gigayearyear * s^-1) for the purpose of timesteps by 1 year.
            return h0
        if to_operator == False: # take operational format of h0 and return to "visual" format
            h0 = self.H * mpc / (1e9 * year * 1e3)
            return h0
    # Calculate the addition constant (do at start) (only useful for curved universes)
    def const_calc(self):
        self.addconst = 1 - (self.Ov + self.Or + self.Om)
    # Step the system forward by X GIGA_years.
    def forstep(self, X):
        # Step R
        #print(self.R)
        self.R = self.R + X*(((self.R**2) * (self.H**2))**0.5)
        # Recalculate H
        Hsquared = (self.H0**2)*((self.Or/(self.R**4) + self.Om/(self.R**3) + self.Ov) + (self.addconst/(self.R**2)))
        self.H = np.sqrt(Hsquared)
        # Alter age GOOD
        self.age += X
    # Step the system forward by X GIGA_years.
    def backstep(self, X):
        # Step R
        self.R = self.R - X*(((self.R**2) * (self.H**2))**0.5)
        # Recalculate H
        Hsquared = (self.H0**2)*((self.Or/(self.R**4) + self.Om/(self.R**3) + self.Ov) + (self.addconst/(self.R**2)))
        self.H = np.sqrt(Hsquared)
        # Alter age GOOD
        self.age -= X
    # Dump age/R/H
    def dump(self):
        return [self.age, self.R, self.H]
# Universe class for the w-model. MATTER - VACUUM - RADATION = order.
# TODO: verify that changes to Gy instead of Yr actually work
class Universe_w():
    def __init__(self,w, rho, R0, RPrime0, h0, T0, X):
        self.w = w # EOS matrix!
        self.rho = rho # Density matrix (not density parameter matrix!!!)
        self.R = R0 # Initial scale factor
        self.Rprime = RPrime0 # R' / (dR/dt) in Gy*s^-1
        self.h = h0 # In the form of the integrator (pre-defined. see notes.)
        self.t = T0 # initial age in Gy
        self.step = X # timestep in Gy

    # Calculate R'' given parameters for current self.setup.
    def accel_calc(self, rho, w):
        summation = 0
        for num, par in enumerate(w):
            summation += rho[num] * (self.R**(-3*(1 + par))) * (1 + 3*par)
        accel = w_inteconst*(-1 * self.R)*summation
        return accel

    # Verlet forward_step
    def forstep(self):
        accel = self.accel_calc(self.rho, self.w)
        self.R += self.Rprime*self.step + 0.5*accel*(self.step**2)
        self.Rprime += 0.5*(accel + self.accel_calc(self.rho, self.w))*self.step
        self.t += self.step

    # Verlet back step
    def backstep(self):
        accel = self.accel_calc(self.rho, self.w)
        self.R -= self.Rprime*self.step + 0.5*accel*(self.step**2)
        self.Rprime -= 0.5*(accel + self.accel_calc(self.rho, self.w))*self.step
        self.t -= self.step

    # Dump current details
    def dump(self):
        return [self.R, self.Rprime, self.t]
# Throwaway generate the standard lambda-CDM model given parameters.
# TODO: densities over time, graphing, etcetera.
class lambda_cdm():
    def __init__(self, densities, hubble):
        self.densarray = densities
        self.hub = hubble # In the standard kms^-1kmpc^-1 format

    # Will return R(t) for given t, provided in gigayears
    def lambda_timecalc(self, t_list, t_shift):
        # Shift analytic model by amount t_shift
        #t_shift = 0.174
        t_list = [d + t_shift for d in t_list]

        # Do all standard shit
        lamfrac = (self.densarray[1] / self.densarray[0]) ** (1 / 3)
        hnought = self.hub * 1e3 / mpc
        u = ((1*10**9)*365 * 86400)
        t_list = [d*u for d in t_list]
        lamt = (2 / (3 * hnought)) / (self.densarray[1] ** 0.5)
        r_vals = [lamfrac * ((np.sinh(t / lamt)) ** (2 / 3)) for t in t_list]
        return r_vals
# main class that holds all stuff like integrators/etc.
# Group is save group identity. Step in years for integrators.
# This is for the finite step model that Andy said to use.
class primary():
    def __init__(self, group, set, age_min, age_max, minRmaxR, multistep, Rtransition):
        self.group, self.set = group, set
        self.agemin = age_min
        self.agemax = age_max
        self.minRmaxR = minRmaxR
        self.multistep = multistep
        self.Rtransition = Rtransition


    # Integrator. Define two Universes (one going back, one going forward.)
    def integrator(self):
        # Forward and backward going Universes
        forwd, backwd = Universe(dens_vac, dens_rad, dens_mat, hub, current_age), Universe(dens_vac, dens_rad, dens_mat, hub, current_age)
        # Calculate operator values
        forwd.H0, backwd.H0 = forwd.converter(True), backwd.converter(True)
        forwd.H, backwd.H = forwd.converter(True), backwd.converter(True)
        forwd.const_calc(), backwd.const_calc()
        # Dump lists. Initial dump is provided.
        Backdump, Fordump = [backwd.dump()],[forwd.dump()]
        # LOOPS FOR AGES!
        if self.minRmaxR == False:
            # Loop for FORWARD integrator
            while True:
                if forwd.R >= self.Rtransition:
                    # Greater than transition. Use large step.
                    forwd.forstep(self.multistep[0])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
                else:
                    # Less than transition. Use smol step.
                    forwd.forstep(self.multistep[1])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
            # Loop for BACKWARD integrator
            while True:
                # Less than transition. Use small step.
                if backwd.R < self.Rtransition:
                    backwd.backstep(self.multistep[1])
                    Backdump.append(backwd.dump())
                    if backwd.age <= self.agemin:
                        break
                # Greater than or equal to transition. Use large step.
                else:
                    backwd.backstep(self.multistep[0])
                    Backdump.append(backwd.dump())
                    if backwd.age <= self.agemin:
                        break
        if self.minRmaxR != False:
            # Loop for FORWARD integrator
            minR,maxR = self.minRmaxR
            while True:
                if forwd.R > maxR:
                    break
                # Above transition: large step
                if forwd.R >= self.Rtransition:
                    forwd.forstep(self.multistep[0])
                    Fordump.append(forwd.dump())
                # Below. Small step.
                else:
                    forwd.forstep(self.multistep[1])
                    Fordump.append(forwd.dump())

            # Loop for BACKWARD integrator
            while True:
                # Below. Use small step.
                if backwd.R < self.Rtransition:
                    backwd.backstep(self.multistep[1])
                    if backwd.R >= minR:
                        Backdump.append(backwd.dump())
                    else:
                        break
                # Above. Large step
                else:
                    backwd.backstep(self.multistep[0])
                    if backwd.R >= minR:
                        Backdump.append(backwd.dump())
                    else:
                        break
        # Format + save data. Both have same initial value (same age/etc) so clip the first element of backdump.
        Backdump = Backdump[1:]
        Backdump.reverse()
        Backdump += Fordump
        table = Table()
        table['t'], table['R'], table['H'] = np.array(Backdump).T

        # Also save H in a more user-friendly way.
        # h0 = self.H * (1e9) * year * 1e3 / mpc
        table['H_normal'] = table['H'] /( (1e9) * year * 1e3 / mpc)# kms^-1 Mpc^-1

        print("Sim done. Pickling.")

        # PICKLE THE BASTARD
        # Write to file
        f_myfile = open(self.set + '.pickle', 'wb')
        pickle.dump(table, f_myfile)
        f_myfile.close()

        print("Pickled!")

        # For the sake of the scalar field integrator, we require initial conditions. We'll just quickly grab those from this.
        # Segment a part of the table (say the first 100 elements) and save those to a separate file, HDF5.
        subtable = table[0:100]
        filer.write_table(self.group, self.set + "_subtable", subtable)

    # Visualize/Graph
    def R_grapher(self, tlims, rlims, t_shift):
        # Get data from the pickle!
        # Read from file
        print("Loading Pickle.")
        f_myfile = open(self.set + '.pickle', 'rb')
        data = pickle.load(f_myfile)  # variables come out in the order you put them in
        f_myfile.close()
        print("Loaded!")


        # Also generate analytic solution for this situation.
        lam_class = lambda_cdm(den_array, hub)
        analytic_y = lam_class.lambda_timecalc(data['t'], t_shift)
        # Generate rescale factor necessary for analytic_y
        scale_nought = lam_class.lambda_timecalc([current_age], t_shift)[0]
        analytic_y = np.array(analytic_y) / scale_nought  # assumes R0 = 1
        data['R_ana'] = analytic_y
        print("Got analytic data using time array.")

        x, y = data['t'], data['R']
        x, y = np.array(x), np.array(y)

        # Generate gradient for plots (same method as with integrator.) in log space.
        y_prime = []
        for num, y0 in enumerate(y):
            try:
                grad = (np.log(y[num + 1]) - np.log(y0)) / np.log(x[num + 1]/x[num])
                y_prime.append(grad)
            except:
                y_prime.append(y_prime[num - 1])

        # Some R might be a NaN
        for num, i in enumerate(y):
            if np.isnan(i) == True:
                y[num] = 0


        fig, axes = plt.subplots(2, figsize=(10, 10), dpi=72)
        fig.suptitle(("Evolution of R(t) for (t + shift) of {0:.3f}").format(t_shift))
        axes[0].plot(x,y, lw=1, color="black")
        axes[0].plot(x, analytic_y, lw=1, color="red")
        axes[0].set(xlabel="t / Gy",
                    ylabel="R(t)",
                    xlim=tlims,
                    ylim=rlims)
        legend_elements = [Patch(facecolor='black', edgecolor='black', label="Model"),
                           Patch(facecolor='red', edgecolor='black', label="Analytic")]
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[0].legend(handles=legend_elements, loc='upper right')

        # Calculate fractional
        analytic_fracs = []
        for num, val in enumerate(analytic_y):
            fractional = (analytic_y[num] - y[num]) / y[num]
            analytic_fracs.append(fractional)

        # Grab areas where fractional is to within 1% (alterable)
        percentage = 1
        analytic_satisfied, t_satisfied = [],[]
        for num,frac in enumerate(analytic_fracs):
            if np.abs(frac)*100 <= percentage:
                analytic_satisfied.append(frac*100), t_satisfied.append(x[num])


        # Error plotter
        axes[1].plot(x, np.array(analytic_fracs)*100, lw=1, color="red")
        axes[1].plot(t_satisfied, analytic_satisfied, lw=2, color="green")
        axes[1].set(xlabel="t / Gy",
                    ylabel="% Difference R(t)",
                    xlim=tlims,
                    ylim=[-10,10])
        legend_elements = [Patch(facecolor='red', edgecolor='black', label="Analytic - Model"),
                           Patch(facecolor='green', edgecolor='black', label="Err <= 1%")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[1].legend(handles=legend_elements, loc='upper right')
        plt.savefig(astrowrapper.rootdir + "\\Figures\\" + "roft_default.png", dpi=300)
        plt.show()
# This is for the w-model simulator.
# Takes various w models (for matter, vac, rad)
# Also takes current density parameters, current Hubble. Uses density parameters to calculate the actual densities involved.
# This one works with densities directly (will be adapted in future to deal with density parameters once adequate equation is given.
# w and den in order of "matter, vacuum, radiation"
# TODO: fix/etc like we did for other class
class primary_w():
    def __init__(self, w, den, H0, R0, T0, X, agerange, group, set):
        self.w = w
        self.den = den
        self.R0, self.T, self.H0 = R0, T0, H0
        self.Rprime = 0
        self.h0 = H0 * year * (1e3/mpc) # For integrator.
        self.OriH = H0
        self.rho = [0,0,0]
        self.X = X
        self.T0 = T0
        self.agerange = agerange
        self.group, self.set = group, set

    # Get initial densities via density parameters + a critical density
    def init_density(self):
        rho_crit = (3*(self.H0**2)/(8*np.pi*G))*(1e6)*(mpc**-2) # kgm^-3
        self.rho = [d*rho_crit for d in self.den]

    # Get the current value of R prime/gradient of R (per year.)
    def init_Rprime(self):
        self.Rprime = self.R0*self.h0

    # Integrator for the w-universe.
    def integrator(self):
        # Instantiate universes for forward/backward
        forwd, backwd = Universe_w(self.w, self.rho, self.R0, self.Rprime, self.h0, self.T0, self.X), Universe_w(self.w, self.rho, self.R0, self.Rprime, self.h0, self.T0, self.X)
        # Set up limits for ages
        minage, maxage = self.agerange
        # Set up variable tables
        forwd_dumps, backwd_dumps = [forwd.dump()],[backwd.dump()]
        # Run forward integrator
        while True:
            forwd.forstep()
            forwd_dumps.append(forwd.dump())
            if forwd.t >= maxage:
                break
            if forwd.R <= 0:
                break
        # Run backward integrator
        while True:
            backwd.backstep()
            backwd_dumps.append(backwd.dump())
            if backwd.t <= minage:
                break
            if forwd.R <= 0:
                break
        # Clip the first element of the forward integrator, reverse backwd, and get R/R'/t
        forwd_clipped = forwd_dumps[1:]
        backwd_dumps.reverse()

        R, Rprime, t = np.array(backwd_dumps + forwd_clipped).T
        logR, logRprime, logt = np.log(R), np.log(Rprime), np.log(t)

        # Create table
        table = Table()
        table['R'], table['Rprime'], table['t'], table['logR'], table['logRprime'], table['logt'] = R, Rprime, t, logR, logRprime, logt
        filer.write_table(self.group, self.set, table)

    # Visualize/Graph
    # Graphing/visualization
    def R_grapher(self, tlims):
        # Get data
        data = filer.read_table(self.group, self.set)
        x, y = data['t'], data['R']
        x, y = np.array(x), np.array(y)

        # Generate gradient for plots (same method as with integrator.) in log space.
        y_prime = []
        for num, y0 in enumerate(y):
            try:
                grad = (np.log(y[num + 1]) - np.log(y0)) / np.log(x[num + 1]/x[num])
                y_prime.append(grad)
            except:
                y_prime.append(y_prime[num - 1])

        # Some R might be a NaN
        for num, i in enumerate(y):
            if np.isnan(i) == True:
                y[num] = 0

        # Generate log plots w.r.t normalized t (against hubble t = 1/h0)
        #x_norm = year*x*hub*1000/(mpc)
        #logxnorm = np.log(x_norm)
        #logy = np.log(y)

        # Overplotting for the solution to lambda CDM, for FLAT cosmology with no radiation.
        fig, axes = plt.subplots(2, figsize=(10, 10), dpi=72)
        fig.suptitle("Evolution of R(t)")
        axes[0].plot(x, y, lw=1, color="black")
        axes[0].plot(x, self.lambda_Cdm(x), lw=1, color="red")
        axes[0].set(xlabel="t / Yr",
                    ylabel="R(t)",
                    xlims=tlims)
        legend_elements = [Patch(facecolor='black', edgecolor='black', label="Model"),
                           Patch(facecolor='red', edgecolor='black', label="Analytic")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[0].legend(handles=legend_elements, loc='upper right')

        axes[1].plot(x, self.lambda_Cdm(x) - y, lw=1, color="green")
        axes[1].set(xlabel="t / Yr",
                    ylabel="Dif. R(t)",
                    xlims=tlims)
        legend_elements = [Patch(facecolor='green', edgecolor='black', label="Analytic - Model")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[1].legend(handles=legend_elements, loc='upper right')



        """
        fig, axes = plt.subplots(3, figsize=(10, 10), dpi=72)
        fig.suptitle("Evolution of R(t)")
        axes[0].plot(x,y, lw=1, color="black")
        axes[0].set(xlabel="t / Yr",
                    ylabel="R(t)")
        axes[1].plot(logxnorm, logy, lw=1, color="black")
        axes[1].set(xlabel="ln(" + r'$H_0$' + "t)",
                    ylabel="ln(R(t))")
        axes[2].plot(x, y_prime, lw=1, color="black")
        axes[2].set(xlabel="t",
                    ylabel="R'")"""

        plt.show()



"""
FINAL CONSTANTS FOR INTEGRATION
t = -0.015299719978238067
R = 9.039240480615748E-7 
H = 4.402022513073902E10 
"""

# Legacy Parameters.
dens_mat, dens_vac, dens_rad = 0.3, 0.7, 0
w_array = [0, -1, 1/3] # Equation of state parameters, MATTER:VACUUM:RADIATION
multi_time_step_in_gigayears = [10**-3, 10**-4] # Test to be used with scalar field integrator for when R is very small (<= 10^-3) to improve sim accuracy.
transition_R = 0.1 # for transition from the low_step ([1]) to fast_step ([0]) regime.
age_R_ONE, hub_R_ONE = 13.7, 68  # 13.7971, 67.32 # leave as default from Planck 2018
scalar_init_R_ONE = 1
scalar_minmaxR = [0.010000000000001,15]
t_shift = 0.1737
scalar_filename = "alpha_zero"
scalar_minmaxage = [-0.3,30] # [-0.2,100]
scalar_minmaxage_nought = copy.deepcopy(scalar_minmaxage)

# SIMPLIFICATION.
# Scalar field model takes INITIAL MATTER DENSITY
# This is set by the initialization R value + modifies our present-day density.
# Scalar field model takes INITIAL SCALAR FIELD DENSITY
# This is set by some factor F, the ratio of scalar divided by matter.
# Scalar field model takes INITIAL FIELD VALUE AND DERIVATIVE
# Combined with initial density, defines the mass parameter.
# Initial age is given analogously to LCDM for the initial R/etc.
current_age, hub = -0.15518124978059802, hub_R_ONE # 37245.173264532576
scalar_init_R = 0.01
field_params = [2e-7, 0, [1, 0]] # Ratio, Alpha, [Psi, Psiprime]


# Misc formatting of parameters/etc for all the various integrators we have set up.
den_array = [dens_mat, dens_vac, dens_rad] # Array of the three, MATTER:VACUUM:RADIATION
year = 365*86400 # seconds
gigayear = year*(1e9)
mpc = 1e6*scipy.constants.parsec
G = scipy.constants.physical_constants['Newtonian constant of gravitation'][0] # standard units
planck_time = mpf(scipy.constants.physical_constants['Planck time'][0])
gyrplanckdiv = gigayear/planck_time
w_inteconst = 4*np.pi*G * ((1e9*year)**2) / 3 # 4 * pi * G / 3 (in the form of gigayear^-2) for integration routines.
scalar_minmaxage = np.array(scalar_minmaxage)/planck_time * 1e9 * year

# Calculate density of universe at current epoch (13.7 Gyr)
hub_full_current = hub_R_ONE*(1e3)*(1/mpc)
current_critical_density = (3 * (hub_full_current**2)) / (8 * np.pi * G)
current_matter_density_absolute = current_critical_density * dens_mat

# Also get this in terms of mp^2 (i.e the same units as the scalar density. Multi by tp^2)
current_matter_density_absolute_in_planck_mass_to_four = current_matter_density_absolute * planck_time**4
current_matter_density_absolute_planck_mass = current_matter_density_absolute_in_planck_mass_to_four / planck_time**2

# Calculate density of universe at start of scalar integrator
hub_full_init = hub*(1e3)*(1/mpc) # per second
integrator_critical_density_init = (3 * (hub_full_init**2)) / (8 * np.pi * G)
integrator_matter_density_absolute = current_matter_density_absolute * (scalar_init_R**-3)

# Scalar Constants.
scalar_hubble = hub*(1e3)*(1/mpc)*planck_time # units of tp^-1
multi_time_step_in_planck_time = np.array(multi_time_step_in_gigayears)*gigayear/planck_time
current_age_planck_time = current_age*gigayear/planck_time

# Scalar Integrator/Universe Parameter Arrays
universe_object_parameters = [current_age_planck_time, scalar_init_R, field_params] # the scalar value is INITIAL R!
primary_scalar_params = [scalar_filename, *scalar_minmaxage, universe_object_parameters, scalar_minmaxR, multi_time_step_in_planck_time, transition_R]



# Default Integrator
#main = primary("Final", "Finalset", *scalar_minmaxage, scalar_minmaxR, multi_time_step_in_gigayears, transition_R) # Finalset has 1e-7 steps
#main.integrator()
#main.R_grapher([-0.18,14], [0,1], t_shift)

# According to the primary model... we'll start near R = 10^-5 (a few timesteps before zero...)
"""
FINAL CONSTANTS FOR INTEGRATION
t = -0.015299719978238067
R = 9.039240480615748E-7 
H = 4.402022513073902E10 
"""


#According to the primary model... we'll start near R = 10^-5 (a few timesteps before zero...)
"""
t = -0.1726673736404837 Gyr 
R = 1.065708378629629e-5 
H = 1.0705623595312455e9 kms^-1 Mpc^-1 
hub, current_age = -0.1726673736404837, 1.0705623595312455e9
(for these simulation params)
scalar_init_R = 1.065708378629629e-5 
# LET'S GO FOR A LATER START INSTEAD!!! (for accuracy/fast-runs)
# ELEMENT 13
current_age, hub = -0.1726666736404837, 5.0292516379746115e8 
scalar_init_R = 1.7635080617517368e-5
"""


# More. 1E-2
"""
current_age, hub = -0.15518124978059802, 38.064743302015735
scalar_init_R = 0.010000000733600039

"""

# DEPRECATED STUFF FROM SECONDARY



"""


# Misc multiprocessing stuff to make things run fast. Masking is per sim.
def masking(scalarparams):
    # Run all sim stuff
    scalobject = primary_scalar(*scalarparams)
    scalobject.main_multi()
    # Remove the pickle
    try:
        os.remove(astrowrapper.rootdir + "\\" + scalarparams[0] + '.pickle')
    except:
        pass
    # Pickle the primary_scalar object (for the sake of later analyses)
    f_myfile = open(astrowrapper.rootdir + "\\Trials\\" + scalobject.set + '.pickle', 'wb')
    pickle.dump(scalobject, f_myfile)
    f_myfile.close()
    
# Multiprocess the sims. We've modified the code!!!
def miscellaneous_simrunner():
    # We need to set this up unique for each range.
    possible_alpha = [3, 4] # -2, 1, 2,
    fracexparange = np.linspace(-7, 0, 40)
    possible_fractions = [10**float(d) for d in fracexparange]
    psiexparange = np.linspace(-2, 1, 30)
    possible_psis = [10**float(d) for d in psiexparange]
    psiprimearange = np.linspace(-65, -59, 10)
    possible_primes = [10 ** float(d) for d in psiprimearange]
    possible_primes.append(0)
    possible_psiprimes_mintwo = [-1*d for d in possible_primes]
    parameter_arrays = []

    possible_psialphas = [possible_psis for d in possible_alpha]
    possible_psiprimesalphas = [possible_psiprimes_mintwo, possible_primes,possible_primes,possible_primes,possible_primes]

    #new_alphas = [2]
    #fracs = possible_fractions
    #psis = [possible_psis]
    #psiprimes = [possible_primes]

    for num, alpha in enumerate(possible_alpha):
        for psi in possible_psialphas[num]:
            for psiprime in possible_psiprimesalphas[num]:
                for denfrac in possible_fractions:
                    try:
                        filename = ("{0:.2e}_{1:.2e}_{2:.1e}_{3:.0f}").format(psi, psiprime, denfrac, alpha)
                        newfieldparams = [denfrac, alpha, [psi,psiprime]]
                        array = copy.deepcopy(primary.primary_scalar_params)
                        array[0] = filename
                        array[3][2] = newfieldparams
                        parameter_arrays.append(array)
                    except Exception as e:
                        print("fucked", e)
                        pass



    if __name__ == '__main__':
        pool = multiprocessing.Pool(10)
        pool.map(masking, parameter_arrays)
    else:
        print("Not main, mate.")
  
    
# Takes a list of scalar_field objects and does various plots. Modified to handle "FRAC".
def pickle_plotter(scalar_list, alpha):
    # Generic scatter plot for all points. First collect up the [psi, psiprime] for each point.
    fracs = [d.UNIVERSEPARAMS[2][0] for d in scalar_list]
    psis = [d.UNIVERSEPARAMS[2][2][0] for d in scalar_list]
    satisfieds = [d.satisfied for d in scalar_list]
    greens, reds = [],[]
    psigruns, psirots = [],[]
    for num, sat in enumerate(satisfieds):
        if sat == "True":
            greens.append(fracs[num])
            psigruns.append(psis[num])
        if sat == "False":
            reds.append(fracs[num])
            psirots.append(psis[num])



    # Next produce the plot
    fig, axs = plt.subplots(1)
    axs.scatter(psirots,reds, color="red", label="Invalid", s=16)
    axs.scatter(psigruns,greens, color="green",label="Valid", s=16)
    axs.set(xlabel="psi",
            ylabel="frac",
            title=alpha + ", psiprime is 0",
            xlim=[0,1])# ,ylim=[0.0675e-60, 0.09e-60]
    axs.legend()
    plt.show()
# This misc function will collect up all the pickles and produce graphs for validity using pickle_plotter.
def pickle_analysis(alphas):
    # Collect up all pickles
    os.chdir(astrowrapper.rootdir + "\\Trials")
    files = os.listdir()
    scalar_field_objects = []
    for file in files:
        f_myfile = open(file, 'rb')
        field_object = pickle.load(f_myfile)
        scalar_field_objects.append(field_object)
        f_myfile.close()

    # Set up alpha strings + relevant lists to hold tuples for each alpha (for plotting)
    alphastr = [("{0:.0f}").format(d) for d in alphas] # string it up
    alpha_lists = [[] for d in alphastr]

    # Split up the pickles into relevant alpha_list
    object_alphas = [("{0:.0f}").format(d.UNIVERSEPARAMS[2][1]) for d in scalar_field_objects]
    for num, alpha in enumerate(alphastr):
        for object_num, object_alpha in enumerate(object_alphas):
            if alpha == object_alpha:
                alpha_lists[num].append(scalar_field_objects[object_num])


    # Then run the plotting function for each relevant scalar object list.
    for num,alphalist in enumerate(alpha_lists):
        pickle_plotter(alphalist, alphastr[num])
#pickle_analysis([-2,1,2,3,4])


"""