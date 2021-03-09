# Many imports
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
from astropy.io import ascii
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
# Own imports
import astrowrapper
import primary
from mpmath import mpf, nstr
from primary import font
import pickle
mpl.rc('font', **font)

# This will hold all the final R(t)/etc integrators.

# Canonical scalar field.
# init_psi = [psi0, psi'0, psi''0]
# Potentials are attached, so far massive_gamma, in units of mp^2.
# TODO: Units
class field(object):
    def __init__(self,  field_params): #massnorm, alpha, init_psi, X):
        self.massnorm, self.alpha = field_params[0], field_params[1]
        self.psi, self.psiprime = field_params[2]
        self.X = field_params[-1]

    # Evaluate the starting values of the constant (for if you aren't manually specifying it/etc, or just want to.)
    def start_eval(self, densmat_nought, R_start, H_start):
        lamden_nought = 1 - densmat_nought*(R_start**-3)
        first_bracket = lamden_nought*(3*(H_start**2)/(8*np.pi))
        next_bracket = first_bracket - (0.5*(self.psiprime**2))
        next_bracket = next_bracket*(self.psi**self.alpha)
        final_bracket = next_bracket**(1/(4+self.alpha))
        self.massnorm = final_bracket

    # Return value of the massive potential (GAMMA!!!), in unites of tp^-2/mp^2, alongside its derivative.
    def massive_gamma(self):
        return ((self.massnorm)**(4 + self.alpha))*(self.psi**(-1*self.alpha))
    def massive_gamma_prime(self):
        return (-self.alpha)*((self.massnorm) ** (4 + self.alpha)) * (self.psi ** ((-1 * self.alpha) - 1))

    # Return the dynamic quantity for the Friedmann (1/2prime^2 + gamma) in units of mp^2/tp^-2
    def total_e(self):
        return (0.5*(self.psiprime**2) + self.massive_gamma())

    # Step psi_prime in gigayears: this does both field primeprime, prime, and no prime
    def psiprime_forwd(self, H):
        psiprimechange = ((-self.massive_gamma_prime()) - 3*H*self.psiprime)*self.X
        self.psiprime += psiprimechange
        self.psiprimedot = psiprimechange/self.X
    def psiprime_backwd(self, H):
        psiprimechange = ((-self.massive_gamma_prime()) - 3*H*self.psiprime)*self.X
        self.psiprime -= psiprimechange
        self.psiprimedot = -psiprimechange / self.X

    # Step psi
    def psi_forwd(self):
        self.psi += self.psiprime*self.X
    def psi_backwd(self):
        self.psi -= self.psiprime*self.X

    # Dump
    def dump(self):
        return [self.psi, self.psiprime, self.massive_gamma()]

# Scalar-field Universe Class
class field_Universe(object):
    def __init__(self, massden, age, R0, H0, fieldparams, X):
        self.H0 = H0 # save H0
        self.massden, self.age, self.R, self.H, self.step = massden, age, R0, H0, X
        self.Rprime = R0*H0
        self.scalar = field(fieldparams)

    # DUMP IT ALL BABY!
    def dump(self):
        return [[self.age, self.R, self.H], self.scalar.dump()]

    # Step forward + dump return
    def forstep(self):
        # Step forward R, Rprime, etc
        self.R += (self.H*self.R*self.step)
        HSQRD = ((self.massden*(self.R**-3)*(self.H0**2)) + (((8*np.pi/3)*self.scalar.total_e())*((1/1)**2)))
        self.H = np.sqrt(HSQRD)

        # Step forward the scalar field + derivative
        self.scalar.psi_forwd()
        self.scalar.psiprime_forwd(self.H)

        self.age += self.step

        return [[self.age, self.R, self.H], self.scalar.dump()]

    # Step backward + dump return
    def backstep(self):
        # Step forward R, Rprime, etc

        self.R -= (self.H*self.R*self.step)
        HSQRD = ((self.massden * (self.R ** -3) * (self.H0 ** 2)) + (((8 * np.pi / 3) * self.scalar.total_e()) * ((1 / 1) ** 2)))
        self.H = np.sqrt(HSQRD)
        # Step forward the scalar field + derivative
        self.scalar.psi_backwd()
        self.scalar.psiprime_backwd(self.H)

        self.age -= self.step

        return [[self.age, self.R, self.Rprime], self.scalar.dump()]

# This is the finite step simulator but with the scalar field.
class primary_scalar():
    def __init__(self, group, set, age_min, age_max, UNIVERSEPARAMS):
        self.group, self.set = group, set
        self.agemin = age_min
        self.agemax = age_max
        self.UNIVERSEPARAMS = UNIVERSEPARAMS
        self.table = True
        self.massnorm = True
    # Integrator. Same as primary.
    def integrator(self):
        # Create forward/backward-integrating universes.
        forwd, backwd = field_Universe(*self.UNIVERSEPARAMS), field_Universe(*self.UNIVERSEPARAMS)
        # Quickly set the normalization of the potential/etc for flat universe.
        forwd.scalar.start_eval(self.UNIVERSEPARAMS[0], self.UNIVERSEPARAMS[2], forwd.H), backwd.scalar.start_eval(self.UNIVERSEPARAMS[0], self.UNIVERSEPARAMS[2], forwd.H)

        Backdump, Fordump = [forwd.dump()], [backwd.dump()]
        # Forward Loop
        while True:
            forwd.forstep()
            Fordump.append(forwd.dump())
            if forwd.age > self.agemax:
                break
        # Backward Loop
        while True:
            backwd.backstep()
            Backdump.append(backwd.dump())
            if backwd.age <= self.agemin:
                break

        # DUMP FORMAT
        """
        [[self.age, self.R, self.Rprime], self.scalar.dump()]
        [self.psi, self.psiprime, self.massive_gamma()]
        """
        # Format + save.
        Backdump = Backdump[1:]
        Backdump.reverse()
        Backdump += Fordump
        agedumps, psidumps = [],[]
        for dump in Backdump:
            agedumps.append(dump[0]), psidumps.append(dump[1])
        agedumps, psidumps = np.array(agedumps).T, np.array(psidumps).T
        ages, scales, hubbles = agedumps
        psi, psiprimes, massive_gammas = psidumps
        table = Table()
        table['t'], table['R'], table['H'], table['psi'], table['psiprime'], table['massive_gamma'] = ages, scales, hubbles, psi, psiprimes, massive_gammas

        # PICKLE THE BASTARD
        # Write to file
        f_myfile = open(self.set + '.pickle', 'wb')
        pickle.dump(table, f_myfile)
        f_myfile.close()

        # Save the mass normalization constant.
        massnorm = forwd.scalar.massnorm
        self.massnorm = massnorm




    # R_grapher, same as primary.
    # Visualize/Graph
    def grapher(self, tlims, rlims, t_shift, w_lims, field_lims, density_lims, gamma_lims):
        # Get data from the pickle!
        # Read from file
        f_myfile = open(self.set + '.pickle', 'rb')
        data = pickle.load(f_myfile)  # variables come out in the order you put them in
        f_myfile.close()

        x, y = data['t'], data['R']
        x, y = np.array(x), np.array(y)
        x = x*primary.planck_time/primary.gigayear
        x = np.array([float(d) for d in x])

        # Also generate analytic solution for this situation.
        lam_class = primary.lambda_cdm(primary.den_array, primary.hub)
        analytic_y = lam_class.lambda_timecalc(x, t_shift)
        # Generate rescale factor necessary for analytic_y
        scale_nought = lam_class.lambda_timecalc([primary.current_age], t_shift)[0]
        analytic_y = np.array(analytic_y) / scale_nought  # assumes R0 = 1
        data['R_ana'] = analytic_y

        # Create figure
        fig, axes = plt.subplots(6, figsize=(12,12),dpi=72, sharex="col",  constrained_layout=True)

        # Set up all the title parameters
        alpha = primary.field_params[1]
        psi, psidot = primary.field_params[2]
        timestep = primary.time_step_in_gigayears
        density_param = primary.dens_mat
        massnorm = self.massnorm

        psipsidotmass = r'$\psi, \dot{\psi}, M, \alpha$' + "=" + ("{0},{1},{2},{3:.1f}").format(nstr(psi, 2, show_zero_exponent=True),nstr(psidot, 2, show_zero_exponent=True),nstr(massnorm, 2, show_zero_exponent=True), alpha)
        timestep_normal = r'$\Delta{t}$' + "=" + ("{0:.2f} in Gyr").format(timestep)
        densparemnormal = r'$\Omega_{m,0}$' + "=" + ("{0:.2f}").format(density_param)
        full_title = psipsidotmass + "," + densparemnormal + "," + timestep_normal
        fig.suptitle(full_title)

        #fig.suptitle(("Evolution of R(t) for (t + shift) of {0:.3f}").format(t_shift))
        axes[0].plot(x,y, lw=1, color="red")
        axes[0].plot(x, analytic_y, lw=1, color="green")
        axes[0].set(ylabel="R(t) Scale")

        legend_elements = [Patch(facecolor='red', edgecolor='black', label="Model"),
                           Patch(facecolor='green', edgecolor='black', label="Analytic")]
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
                analytic_satisfied.append(frac), t_satisfied.append(x[num])


        # Error plotter
        axes[1].plot(t_satisfied, analytic_satisfied, lw=2, color="green")
        axes[1].plot(x, analytic_fracs, lw=1, color="red")
        axes[1].set(ylabel="Dif. Fract. R(t)")
        legend_elements = [Patch(facecolor='red', edgecolor='black', label="Analytic - Model"),
                           Patch(facecolor='green', edgecolor='black', label="Err <= 1%")]  # ,
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[1].legend(handles=legend_elements, loc='upper right')

        if tlims != 0:
            axes[0].set(xlim=tlims)
            axes[1].set(xlim=tlims)
        if rlims != 0:
            axes[0].set(ylim=rlims)
            axes[1].set(ylim=[-0.1,0.1])

        # Next up: grab the equation-of-state-parameter and also psi.
        # EOS = pressure divided by density. We have psi and gamma, which are just phi/V divided by mp^2
        EOS = (0.5 * data['psiprime'] ** 2 - data['massive_gamma']) / (
                    0.5 * data['psiprime'] ** 2 + data['massive_gamma'])

        axes[2].plot(x, EOS, lw=1, color="red")
        axes[2].set(ylabel="w(t) EOS")
        if tlims != 0:
            axes[2].set(xlim=tlims)
        axes[2].set(ylim=w_lims)

        # Also Psi.
        axes[3].plot(x, data['psi'], lw=1, color="red")
        axes[3].set(ylabel=r'$\psi$' + " Field (t)")
        if tlims != 0:
            axes[3].set(xlim=tlims,
                        ylim=field_lims)


        # Also want to plot the density parameter over time.
        # denmat = matter_0



        # We're also going for density parameter.
        density_scalar = (8*np.pi/(3*(data['H']**2)))*(0.5*data['psiprime']**2 + data['massive_gamma'])
        axes[4].plot(x, density_scalar, lw=1, color="red")
        axes[4].set(ylabel=r'$\Omega_\psi$' + " Field (t)")
        if tlims != 0:
            axes[4].set(xlim=tlims)
        if density_lims != False:
            axes[4].set(ylim=density_lims)

        # Also want to plot the density parameter over time.
        # denmat = matter_0

        # Massive potential, too.
        halfpsiprimesq, massive_gamma = 0.5*data['psiprime']**2, data['massive_gamma']
        axes[5].plot(x, halfpsiprimesq, lw=1, color="red")
        axes[5].plot(x, massive_gamma, lw=1, color="blue")
        axes[5].set(ylabel='Energies' + " (t)")
        if tlims != 0:
            axes[5].set(xlim=tlims)
        if gamma_lims != False:
            axes[5].set(ylim=gamma_lims)
        legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\psi$'+"'"),
                           Patch(facecolor='blue', edgecolor='black', label=r'$\Gamma_\psi$')]
        # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
        axes[5].legend(handles=legend_elements, loc='upper right')

        # Also want to plot the density parameter over time.
        # denmat = matter_0

        axes[5].set(xlabel="t / Gyr")

        plt.savefig(astrowrapper.rootdir + "\\Figures\\" + self.set + ".png", dpi=300)
        plt.show()

scalar_universe = primary_scalar(*primary.primary_scalar_params)
#scalar_universe.integrator()
scalar_universe.grapher([0, 1400], [0,1], 0, [-2,2], [-0.1,0.1], [0,1.1], False)  # t_lims, r_lims, t_shift_for_analytic, w_lims, field_lims