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
    def __init__(self,  field_params): #massnorm, alpha, init_psi):
        self.massnorm, self.alpha = field_params[0], field_params[1]
        self.psi, self.psiprime = field_params[2]

    # Evaluate the starting values of the constant (for if you aren't manually specifying it/etc, or just want to.)
    def start_eval(self, densmat_nought, H_start):
        lamden_nought = 1 - densmat_nought # assumes   has been accounted for.
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
    def psiprime_forwd(self, H, X):
        psiprimechange = ((-self.massive_gamma_prime()) - 3*H*self.psiprime)*X
        self.psiprime += psiprimechange
        self.psiprimedot = psiprimechange/X
    def psiprime_backwd(self, H, X):
        psiprimechange = ((-self.massive_gamma_prime()) - 3*H*self.psiprime)*X
        self.psiprime -= psiprimechange
        self.psiprimedot = -psiprimechange / X

    # Step psi
    def psi_forwd(self, X):
        self.psi += self.psiprime*X
    def psi_backwd(self, X):
        self.psi -= self.psiprime*X

    # Dump
    def dump(self):
        return [self.psi, self.psiprime, self.massive_gamma()]

# Scalar-field Universe Class
class field_Universe(object):
    def __init__(self, massden, age, R0, H0, fieldparams):
        self.H0 = H0 # save H0
        self.R0 = R0
        self.massden, self.age, self.R, self.H = massden, age, R0, H0
        self.Rprime = R0*H0
        self.scalar = field(fieldparams)

    # DUMP IT ALL BABY!
    def dump(self):
        return [[self.age, self.R, self.H], self.scalar.dump()]

    # Step forward + dump return
    def forstep(self, X):
        # Step forward R, Rprime, etc
        self.R += (self.H*self.R*X)
        HSQRD = ((self.massden*((self.R/self.R0)**-3)*(self.H0**2)) + (((8*np.pi/3)*self.scalar.total_e())*((1/1)**2)))
        self.H = np.sqrt(HSQRD)

        # Step forward the scalar field + derivative
        self.scalar.psi_forwd(X)
        self.scalar.psiprime_forwd(self.H,X)

        self.age += X

        return [[self.age, self.R, self.H], self.scalar.dump()]

    # Step backward + dump return
    def backstep(self, X):
        # Step forward R, Rprime, etc
        self.R -= (self.H*self.R*X)
        HSQRD = ((self.massden*((self.R/self.R0)**-3)*(self.H0**2)) + (((8*np.pi/3)*self.scalar.total_e())*((1/1)**2)))
        self.H = np.sqrt(HSQRD)
        # Step forward the scalar field + derivative
        self.scalar.psi_backwd(X)
        self.scalar.psiprime_backwd(self.H,X)

        self.age -= X

        return [[self.age, self.R, self.Rprime], self.scalar.dump()]

# This is the finite step simulator but with the scalar field.
class primary_scalar():
    def __init__(self, set, age_min, age_max, UNIVERSEPARAMS, minmaxR, multi_timestep, transition_R):
        self.set = set
        self.agemin = age_min
        self.agemax = age_max
        self.UNIVERSEPARAMS = UNIVERSEPARAMS
        self.table = True
        self.massnorm = True
        self.minmaxR = minmaxR
        self.multistep = multi_timestep
        self.transition_R = transition_R
    # Joe Kington https://stackoverflow.com/questions/36455083/how-to-plot-and-work-with-nan-values-in-matplotlib
    def interpolate_gaps(values, limit=None):
        """
        Fill gaps using linear interpolation, optionally only fill gaps up to a
        size of `limit`.
        """
        values = np.asarray(values)
        i = np.arange(values.size)
        valid = np.isfinite(values)
        filled = np.interp(i, i[valid], values[valid])

        if limit is not None:
            invalid = ~valid
            for n in range(1, limit + 1):
                invalid[:-n] &= invalid[n:]
            filled[invalid] = np.nan

        return filled

    # Integrator. Same as primary.
    def integrator(self):
        # Create forward/backward-integrating universes.
        forwd, backwd = field_Universe(*self.UNIVERSEPARAMS), field_Universe(*self.UNIVERSEPARAMS)
        # Quickly set the normalization of the potential/etc for flat universe.
        forwd.scalar.start_eval(self.UNIVERSEPARAMS[0], self.UNIVERSEPARAMS[3]), backwd.scalar.start_eval(self.UNIVERSEPARAMS[0], self.UNIVERSEPARAMS[3])

        Backdump, Fordump = [forwd.dump()], [backwd.dump()]

        # We consider two cases. Setting a minmaxR or not. In this case, age is considered also as is transition R, but minmaxR is false.
        if self.minmaxR == False:
            print("MinmaxR is False. Rolling by time.")
            # Forward Loop
            while True:
                if forwd.R >= self.transition_R:
                    forwd.forstep(self.multistep[1])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
                else:
                    forwd.forstep(self.multistep[0])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
            # Backward Loop
            while True:
                if backwd.R >= 0:
                    if backwd.R < self.transition_R:
                        backwd.backstep(self.multistep[0])
                        if backwd.R >= 0:
                            Backdump.append(backwd.dump())
                        else:
                            break
                    else:
                        backwd.backstep(self.multistep[1])
                        if backwd.R >= 0:
                            Backdump.append(backwd.dump())
                        else:
                            break

        # Here minmaxR isn't false. Integral will satisfy age conditions, R conditions, and transition R.
        else:
            # Forward Loop
            while True:
                # R conditions & Age Coonditions
                if forwd.R > self.minmaxR[1]:
                    break
                if forwd.age >= self.agemax:
                    break
                else:
                    # R exceeds transition. Use upper step.
                    if forwd.R >= self.transition_R:
                        forwd.forstep(self.multistep[1])
                        Fordump.append(forwd.dump())
                    # R is lower than transition, use lower step.
                    else:
                        forwd.forstep(self.multistep[0])
                        Fordump.append(forwd.dump())
            # Backward Loop. LEQ transition R.
            while True:
                # Satisfy the age requirement
                if backwd.age < self.agemin:
                    break
                # For consistency also apply the same as above. If lower than transition, use lower step.
                if backwd.R < self.transition_R:
                    backwd.backstep(self.multistep[0])
                # If higher than transition, use larger step.
                else:
                    backwd.backstep(self.multistep[1])
                # R boundary condition: must be greater than minimum
                if backwd.R >= self.minmaxR[0]:
                    Backdump.append(backwd.dump())
                else:
                    print("Breaking: R less than minimum target. Age is ", backwd.age)
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

        # Also create H_normal, as previously. H here is in units of planck time ^-1, convert to seconds.
        table['H_normal'] = (table['H']/primary.planck_time)*primary.mpc/1e3

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

        # Sort out the density before anything.
        rolling_density = primary.current_matter_density_absolute/(data['R']**3)
        rolling_critical_density = 3*((data['H_normal']*1e3/primary.mpc)**2)/(8 * np.pi * primary.G)
        rolling_parameter = rolling_density/rolling_critical_density

        x, y = data['t'], data['R']
        x, y = np.array(x), np.array(y)
        x = x*primary.planck_time/primary.gigayear
        x = np.array([float(d) for d in x])

        # Also generate analytic solution for this situation. X is in planck times. Need to sort that out.
        lam_class = primary.lambda_cdm(primary.den_array, primary.hub_R_ONE)
        analytic_y = lam_class.lambda_timecalc(x, t_shift)
        # Generate rescale factor necessary for analytic_y
        scale_nought = lam_class.lambda_timecalc([primary.age_R_ONE], t_shift)[0]
        analytic_y = np.array(analytic_y) / scale_nought  # assumes R0 = 1
        data['R_ana'] = analytic_y

        # Create figure
        fig, axes = plt.subplots(6, figsize=(12,12),dpi=72, sharex="col",  constrained_layout=True)

        # Set up all the title parameters
        UNIPARAMS = self.UNIVERSEPARAMS
        alpha = UNIPARAMS[4][1]
        psi, psidot = UNIPARAMS[4][2]
        density_param = primary.dens_mat
        massnorm = self.massnorm

        psipsidotmass = r'$\psi, \dot{\psi}, M, \alpha$' + "=" + ("{0},{1},{2},{3:.1f}").format(nstr(psi, 4, show_zero_exponent=True),nstr(psidot, 2, show_zero_exponent=True),nstr(massnorm, 2, show_zero_exponent=True), alpha)
        densparemnormal = r'$\Omega_{m,0}$' + "=" + ("{0:.2f}").format(density_param)
        full_title = psipsidotmass + "," + densparemnormal
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
            axes[3].set(xlim=tlims)
        if field_lims != False:
            axes[3].set(ylim=field_lims)


        # Also want to plot the density parameter over time.
        # denmat = matter_0



        # We're also going for density parameter.


        density_scalar = (8*np.pi/(3*(data['H']**2)))*(0.5*data['psiprime']**2 + data['massive_gamma'])
        axes[4].plot(x, density_scalar, lw=1, color="red")
        axes[4].plot(x, rolling_parameter, lw=1, color="blue")
        axes[4].set(ylabel=r'$\Omega$' + " (t)")
        if tlims != 0:
            axes[4].set(xlim=tlims)
        if density_lims != False:
            axes[4].set(ylim=density_lims)
        legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\Omega_\Lambda$'),
                           Patch(facecolor='blue', edgecolor='black', label=r'$\Omega_M$')]
        axes[4].legend(handles=legend_elements, loc='upper right')

        # Also want to plot the density parameter over time.
        # Massive potential, too.
        # We also need the normal density of the matter content, too.
        halfpsiprimesq, massive_gamma = 0.5*data['psiprime']**2, data['massive_gamma']
        total_lambda_density_norm = halfpsiprimesq + massive_gamma
        rolling_critical_density_natural = 3 * data['H']**2 / 8 * np.pi
        rolling_matter_density_natural = rolling_critical_density_natural*rolling_parameter
        quicktable = Table()
        quicktable['R'] = data['R']
        quicktable['lam'] = total_lambda_density_norm
        quicktable['mat'] = rolling_matter_density_natural
        quicktable['denmat'] = rolling_parameter
        quicktable['denvac'] = density_scalar
        filer = astrowrapper.hdf5_writer(astrowrapper.rootdir, "shp.hdf5")
        filer.write_table("testgroup", "testset", quicktable)

        axes[5].plot(x, halfpsiprimesq, lw=1, color="red")
        axes[5].plot(x, massive_gamma, lw=1, color="blue")
        axes[5].plot(x, rolling_matter_density_natural, lw=1, color="black")
        axes[5].set(ylabel='Densities' + " (t)")
        if tlims != 0:
            axes[5].set(xlim=tlims)
        if gamma_lims != False:
            axes[5].set(ylim=gamma_lims)
        legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\frac{1}{2}(\dot{\psi})^2$'),
                           Patch(facecolor='blue', edgecolor='black', label=r'$\Gamma_\psi$'),
                           Patch(facecolor='black', edgecolor='black', label=r'$\frac{\rho_m}{m_p^2}$')]
        axes[5].legend(handles=legend_elements, loc='upper right')

        # Also want to plot the density parameter over time.
        # denmat = matter_0

        axes[5].set(xlabel="t / Gyr")

        plt.savefig(astrowrapper.rootdir + "\\Figures\\" + self.set + ".png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()

    # Main runtime for multiprocessing
    def main_multi(self):
        self.integrator()
        self.grapher([-0.1727, 100], [0,20], 0.1726673736404837, [-2,2], False, False, [0, 10**-122]) # hard coded sorry mate


scalar_universe = primary_scalar(*primary.primary_scalar_params)
scalar_universe.integrator()
scalar_universe.grapher([-0.1727, 100], [0,20], 0.1726673736404837, [-2,2], False, False, [0, 10**-123])  # t_lims, r_lims, t_shift_for_analytic, w_lims, field_lims, gamma lims


# Misc multiprocessing stuff to make things run fast.
def masking(scalarparams):
    primary_scalar(*scalarparams).main_multi()
def miscellaneous():
    field_params = [1, 2, [0.7,0]] # Massnorm, Alpha, Psi and Psidot
    possible_masses = np.linspace(0.1,12.5,10)
    possible_alpha = [4, 3, 2, 1]

    parameter_arrays = []


    for mass in possible_masses:
        for alpha in possible_alpha:
            try:
                filename = ("{0:.3f}.{1:.2f}.png").format(mass, alpha)
                newfieldparams = [1, alpha, [mass,0]]
                newfieldparams.append(primary.time_step_in_planck_time)
                array = copy.deepcopy(primary.primary_scalar_params)
                array[0] = filename
                array[3][4] = newfieldparams
                parameter_arrays.append(array)
            except:
                print("fucked")
                pass



    if __name__ == '__main__':
        pool = multiprocessing.Pool(4)
        pool.map(masking, parameter_arrays)
#miscellaneous()