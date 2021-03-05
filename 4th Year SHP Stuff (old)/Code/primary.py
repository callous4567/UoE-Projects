# Code for SHP
# Sebastian Straszak, 2021

# Menagerie of imports of varying use.
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
import itertools
import multiprocessing
from matplotlib.patches import Patch

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from skimage import transform
from scipy.ndimage.interpolation import shift

from SHP2021.Code import astrowrapper

# Set up data file in case it's needed.
filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "shp.hdf5")
#filer.create()

# ALL THE PARAMETERS THAT ONE MIGHT NEED! Give in the form of "density parameters" since all surveys quote like that.
# Current values are random picked from book. If you give as raw densities, then alter the code to to account for this.
# 0.31, 9E-5, 0.69
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
w_inteconst = 4*np.pi*G * (year**2) / 3 # 4 * pi * G / 3 (in the form of year^-2) for integration routines.

# Universe class. Follows LCDM model. Initialize with start parameters. Has functions for stepping forward/backward (in years).
# Rho should be in units of kgm^-3 and H0 in kms^-1Mpc^-1
class Universe():
    def __init__(self, dV,dR ,dM, H0, age):
        self.Ov, self.Or, self.Om = dV, dR, dM
        self.H = H0
        self.H0 = False
        self.addconst = False
        self.R = 1
        self.age = age
    # Converts between the forms of h0 (year vs natural units)
    def converter(self, to_operator):
        if to_operator == True: # convert inputs into the operational format for integrator/etc
            h0 = self.H * year * 1e3 / mpc # In units of (year * s^-1) for the purpose of timesteps by 1 year.
            return h0
        if to_operator == False: # take operational format of h0 and return to "visual" format
            h0 = self.H * mpc / (year * 1e3)
            return h0
    # Calculate the addition constant (do at start)
    def const_calc(self):
        self.addconst = 1 - (self.Ov + self.Or + self.Om)
    # Step the system forward by X years.
    def forstep(self, X):
        # Step R
        #print(self.R)
        self.R = self.R + X*(((self.R**2) * (self.H**2))**0.5)
        # Recalculate H
        Hsquared = (self.H0**2)*((self.Or/(self.R**4) + self.Om/(self.R**3) + self.Ov) + (self.addconst/(self.R**2)))
        self.H = np.sqrt(Hsquared)
        # Alter age GOOD
        self.age += X
    # Step the system forward by X years.
    def backstep(self, X):
        # Step R
        self.R = self.R - X*(((self.R**2) * (self.H**2))**0.5)
        # Recalculate H
        Hsquared = (self.H0**2)*((self.Or/(self.R**4) + self.Om/(self.R**3) + self.Ov) + (self.addconst/(self.R**2)))
        self.H = np.sqrt(Hsquared)

        # Alter age GOOD
        self.age -= X
    # Dump age/R
    def dump(self):
        return [self.age, self.R, self.H]

# Universe class for the w-model. MATTER - VACUUM - RADATION = order.
class Universe_w():
    def __init__(self,w, rho, R0, RPrime0, h0, T0, X):
        self.w = w
        self.rho = rho
        self.R = R0
        self.Rprime = RPrime0
        self.h = h0 # In the form of the integrator.
        self.t = T0
        self.step = X

    # Calculate R'' given parameters.
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



# main class that holds all stuff like integrators/etc.
# Group is save group identity. Step in years for integrators.
# This is for the finite step model that Andy said to use
class primary():
    def __init__(self, group, set, step, age_min, age_max):
        self.group, self.set = group, set
        self.X = step
        self.agemin = age_min
        self.agemax = age_max


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
        # Loop for FORWARD integrator
        while True:
            forwd.forstep(self.X)
            Fordump.append(forwd.dump())
            if forwd.age > self.agemax:
                break
        # Loop for BACKWARD integrator
        while True:
            backwd.backstep(self.X)
            Backdump.append(backwd.dump())
            if backwd.age <= self.agemin:
                break
        # Format + save data. Both have same initial value (same age/etc) so clip the first element of backdump.
        Backdump = Backdump[1:]
        Backdump.reverse()
        Backdump += Fordump
        table = Table()
        table['t'], table['R'], table['H'] = np.array(Backdump).T
        filer.write_table(self.group, self.set, table)


    # Lambda CDM
    def lambda_Cdm(self, t_list, density, oriH):
        lamfrac = (density[1] / density[0]) ** (1 / 3)
        hnought = oriH * 1e3 / mpc
        t_list = t_list * (365 * 86400)
        lamt = (2 / (3 * hnought)) / (density[1] ** 0.5)
        r_vals = [lamfrac * ((np.sinh(t / lamt)) ** (2 / 3)) for t in t_list]
        return r_vals

    # Visualize/Graph
    def R_grapher(self):
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
        x_norm = year*x*hub*1000/(mpc)
        logxnorm = np.log(x_norm)
        logy = np.log(y)

        # Generate analytic lambda CDM solution
        analytic_y = self.lambda_Cdm(x, den_array, hub)

        # Rescale by about 1.74
        analytic_y = [d/1.74 for d in analytic_y]

        fig, axes = plt.subplots(2, figsize=(10, 10), dpi=72)
        fig.suptitle("Evolution of R(t)")
        axes[0].plot(x,y, lw=1, color="black")
        axes[0].plot(x, analytic_y, lw=1, color="red")
        axes[0].set(xlabel="t / Yr",
                    ylabel="R(t)")

        legend_elements = [Patch(facecolor='black', edgecolor='black', label="Model"),
                           Patch(facecolor='red', edgecolor='black', label="Analytic")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[0].legend(handles=legend_elements, loc='upper right')

        # Calculate fractional
        analytic_fracs = []
        for num, val in enumerate(analytic_y):
            fractional = (analytic_y[num] - y[num]) / y[num]
            analytic_fracs.append(fractional)

        axes[1].plot(x, analytic_fracs, lw=1, color="green")
        axes[1].set(xlabel="t / Yr",
                    ylabel="Dif. Fract. R(t)")
        legend_elements = [Patch(facecolor='green', edgecolor='black', label="Analytic - Model")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[1].legend(handles=legend_elements, loc='upper right')




        plt.show()

# This is for the w-model simulator.
# Takes various w models (for matter, vac, rad)
# Also takes current density parameters, current Hubble. Uses density parameters to calculate the actual densities involved.
# This one works with densities directly (will be adapted in future to deal with density parameters once adequate equation is given.
# w and den in order of "matter, vacuum, radiation"
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

    # Returns value of R for flat radiationless lambda CDM, fed entire list of t
    # Den is mat vac rad
    def lambda_Cdm(self, t_list):
        lamfrac = (self.den[1]/self.den[0])**(1/3)
        hnought = self.OriH * 1e3 / mpc
        t_list = t_list * (365*86400)
        lamt = (2/(3 * hnought)) / (self.den[1]**0.5)
        r_vals = [lamfrac*((np.sinh(t/lamt))**(2/3)) for t in t_list]
        return r_vals


    # Visualize/Graph
    # Graphing/visualization
    def R_grapher(self):
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
                    ylabel="R(t)")
        legend_elements = [Patch(facecolor='black', edgecolor='black', label="Model"),
                           Patch(facecolor='red', edgecolor='black', label="Analytic")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[0].legend(handles=legend_elements, loc='upper right')

        axes[1].plot(x, self.lambda_Cdm(x) - y, lw=1, color="green")
        axes[1].set(xlabel="t / Yr",
                    ylabel="Dif. R(t)")
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

    # Run for initialization, integrating + graphing.
    def run(self):
        self.init_density()
        self.init_Rprime()
        self.integrator()
        self.R_grapher()

#w_universe = primary_w(w_array, den_array, hub, 1, current_age, 1e7, [0, 3*current_age], "testing", "test_run")
#w_universe.run()

main = primary("Testing", "Testset", 1e7, 0, 2*current_age)
main.integrator()
main.R_grapher()










