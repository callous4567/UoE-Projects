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
# Cba this just makes trials directory.
try:
    os.mkdir(astrowrapper.rootdir + "\\Trials")
except:
    pass
try:
    os.mkdir(astrowrapper.rootdir + "\\Figures")
except:
    pass

# This will hold all the final R(t)/etc integrators.

# Canonical scalar field.
# init_psi = [psi0, psi'0, psi''0]
# Potentials are attached, so far massive_gamma, in units of mp^2.
# TODO: Units
class field(object):
    def __init__(self,  field_params): #massnorm, alpha, init_psi):
        self.scalmatfrac, self.alpha = field_params[0], field_params[1]
        self.psi, self.psiprime = field_params[2]
        self.massnorm = 0 # Placeholder!!!

    # Go and solve for the initial field mass.
    def start_eval(self, matdeninit):
        scaldeninit = self.scalmatfrac*matdeninit
        scalden_natural = scaldeninit * primary.planck_time**2 * primary.G
        # Now solve for the field mass.
        potinit = scalden_natural - 0.5*(self.psiprime**2)
        init_mass = (potinit/(self.psi**-self.alpha))**(1/(4 + self.alpha))
        self.massnorm = init_mass

    # Return value of the massive potential (GAMMA!!!), in unites of tp^-2/mp^2, alongside its derivative.
    def massive_gamma(self):
        primer = ((self.massnorm)**(4 + self.alpha))*(self.psi**(-1*self.alpha))
        return primer
    def massive_gamma_prime(self):
        primer = (-self.alpha)*((self.massnorm) ** (4 + self.alpha)) * (self.psi ** ((-1 * self.alpha) - 1))
        return primer

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
    def __init__(self, model_start_age, model_start_R, fieldparams):
        # All the initializations from the PRESENT DAY! R=1/etc
        self.H0 = primary.scalar_hubble # For current day.
        self.R0 = primary.scalar_init_R_ONE  # For which we're defining the massden AND H0 AND initial densities.
        self.massden = primary.dens_mat  # For current day.

        # All the parameters for the actual Universe as it rolls along!
        self.H = "Not Defined Yet"
        self.age, self.R = model_start_age, model_start_R

        # The field parameters at initialization!!!
        self.scalar = field(fieldparams)

    # Runs scalar field start eval alongside getting an initial estimate for H0.
    def start_eval(self):
        # Initialize the scalar field
        self.scalar.start_eval(primary.integrator_matter_density_absolute)
        # Grab an initial hubble constant.
        HSQRD = (self.massden * ((self.R / self.R0) ** -3) * (self.H0 ** 2)) + ((8 * np.pi / 3) * self.scalar.total_e())
        self.H = np.sqrt(HSQRD)

    # DUMP IT ALL BABY!
    def dump(self):
        return [[self.age, self.R, self.H], self.scalar.dump()]

    # Step forward + dump return
    def forstep(self, X):
        # Step forward R, Rprime, etc
        self.R += (self.H*self.R*X)
        HSQRD = (self.massden*((self.R/self.R0)**-3)*(self.H0**2)) + ((8*np.pi/3)*self.scalar.total_e())
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

        return [[self.age, self.R], self.scalar.dump()]

# This is the finite step simulator but with the scalar field.
class primary_scalar():
    def __init__(self, set, age_min, age_max, UNIVERSEPARAMS, minmaxR, multi_timestep, transition_R):
        self.set = set
        self.agemin = age_min
        self.agemax = age_max
        self.UNIVERSEPARAMS = UNIVERSEPARAMS
        self.massnorm = True
        self.minmaxR = minmaxR
        self.multistep = multi_timestep
        self.transition_R = transition_R
        self.satisfied = "False" # Holder for a "True" or "False" for if the model is satisfied or not. Should be a 1 or a 0.
    # Joe Kington https://stackoverflow.com/questions/36455083/how-to-plot-and-work-with-nan-values-in-matplotlib
    def interpolate_gaps(self, values, limit=None):
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
        forwd.start_eval(), backwd.start_eval()
        Backdump, Fordump = [forwd.dump()], [backwd.dump()]
        # We consider two cases. Setting a minmaxR or not. In this case, age is considered also as is transition R, but minmaxR is false.
        if self.minmaxR == False:
            print("MinmaxR is False. Rolling by time.")
            # Forward Loop
            while True:
                if forwd.R >= self.transition_R:
                    forwd.forstep(self.multistep[0])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
                else:
                    forwd.forstep(self.multistep[1])
                    Fordump.append(forwd.dump())
                    if forwd.age > self.agemax:
                        break
            # Backward Loop
            while True:
                if backwd.R >= 0:
                    if backwd.R < self.transition_R:
                        backwd.backstep(self.multistep[1])
                        if backwd.R >= 0:
                            Backdump.append(backwd.dump())
                        else:
                            break
                    else:
                        backwd.backstep(self.multistep[0])
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
                if forwd.age > self.agemax:
                    break
                else:
                    # R exceeds transition. Use upper step.
                    if forwd.R >= self.transition_R:
                        forwd.forstep(self.multistep[0])
                        Fordump.append(forwd.dump())
                    # R is lower than transition, use lower step.
                    else:
                        forwd.forstep(self.multistep[1])
                        Fordump.append(forwd.dump())

            # Backward Loop. LEQ transition R.
            while True:
                # Satisfy the age requirement
                if backwd.age < self.agemin:
                    break
                # For consistency also apply the same as above. If lower than transition, use lower step.
                if backwd.R < self.transition_R:
                    backwd.backstep(self.multistep[1])
                # If higher than transition, use larger step.
                else:
                    backwd.backstep(self.multistep[0])
                # R boundary condition: must be greater than minimum
                if backwd.R >= self.minmaxR[0]:
                    Backdump.append(backwd.dump())
                else:
                    break

        # DUMP FORMAT. We're removing the RPrime.
        """
        [[self.age, self.R, self.H], self.scalar.dump()]
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
        #f_myfile = open(self.set + '.pickle', 'wb')
        #pickle.dump(table, f_myfile)
        #f_myfile.close()

        # Save the mass normalization constant.
        massnorm = forwd.scalar.massnorm
        self.massnorm = massnorm

        # We're going with returning the table, instead. But the option still remains to save it.
        return table

    # Visualize/graph the current primary_scalar object.
    # Also ticks off "self.satisfied" or not.
    def grapher(self, tlims, rlims, t_shift, w_lims, field_lims, density_lims, gamma_lims, table):
        data = table
        x, y = data['t'], data['R']
        x, y = np.array(x), np.array(y)
        x = x*primary.planck_time/primary.gigayear
        x = np.array([float(d) for d in x])

        # Sort out the matter density parameter before anything.
        rolling_density = primary.current_matter_density_absolute/(data['R']**3) # kgm^-3
        rolling_critical_density = 3*((data['H_normal']*1e3/primary.mpc)**2)/(8 * np.pi * primary.G) # kgm^-3
        rolling_parameter = rolling_density/rolling_critical_density

        # SETTING THE LIMITS FOR CONTINUING:
        # we want final omega M to be 0.3 around 13.7 Gyr. That's the only condition as of yet.
        # We want a density switchover somewhere between 8 and 10 Gyr, ish.
        # First verify correct density at our era.
        satisfactory_first = "Holder"
        for num, age in enumerate(x):
            if np.abs(age - primary.age_R_ONE) <= 0.5:
                if np.abs(rolling_parameter[num] - primary.dens_mat) <= 0.05:
                    satisfactory_first = "True"
                    break
        # Then verify the switch-over. 8->10 Gyr.
        if satisfactory_first == "True":
            for num, age in enumerate(x):
                if np.abs(age - 9) <= 1:
                    if np.abs(rolling_parameter[num] - 0.5) <= 0.01:
                        self.satisfied = "True"
                        break
        # Both density switchover between 8 and 10 WITH the actual density parameter agreeing in modern era.

        if self.satisfied == "True":
            print("Found a satisfactory model.")
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
            frac = UNIPARAMS[2][0]
            alpha = UNIPARAMS[2][1]
            psi, psidot = UNIPARAMS[2][2]
            massnorm = self.massnorm

            alpha_format = ("{0:.0f}").format(alpha)
            psiformat = ("{0:.1e}").format(psi)
            psiprimeformat = ("{0:.1e}").format(psidot)
            massnormformat = ("{0:.1e}").format(float(massnorm))
            fracformat = ("{0:.1e}").format(frac)

            full_title = r'$\alpha,\psi,\dot\psi,M,F$' + " = " + ("{0},{1},{2},{3},{4}").format(alpha_format,psiformat,psiprimeformat,massnormformat,fracformat)
            fig.suptitle(full_title)

            axes[0].plot(x,y, lw=1, color="red")
            axes[0].plot(x, analytic_y, lw=1, color="green")
            axes[0].set(ylabel="R(t) Scale")
            if rlims != False:
                axes[0].set(ylim=rlims)

            legend_elements = [Patch(facecolor='red', edgecolor='black', label="Model"),
                               Patch(facecolor='green', edgecolor='black', label="Analytic")]
            # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
            axes[0].legend(handles=legend_elements, loc='upper right')

            # THIS IS BEING HASHED OUT FOR NOW FOR THE SAKE OF PROCESSING SPEED!!!
            """
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
            axes[1].set(ylabel="Dif. Fractional. R(t)")
            legend_elements = [Patch(facecolor='red', edgecolor='black', label="Analytic - Model"),
                               Patch(facecolor='green', edgecolor='black', label="Err <= 1%")]  # ,
            # Patch(facecolor='green', edgecolor='black', label="Av Fit")]
            axes[1].legend(handles=legend_elements, loc='upper right')
    
            if tlims != 0:
                axes[0].set(xlim=tlims)
                axes[1].set(xlim=tlims)
            if rlims != 0:
                axes[0].set(ylim=rlims)
                axes[1].set(ylim=[-0.1,0.1]) """

            # Next up: grab the equation-of-state-parameter and also psi.
            # EOS = pressure divided by density. We have psi and gamma, which are just phi/V divided by mp^2
            EOS = (0.5 * data['psiprime'] ** 2 - data['massive_gamma']) / (
                        0.5 * data['psiprime'] ** 2 + data['massive_gamma'])


            axes[2].plot(x, EOS, lw=1, color="red")
            axes[2].set(ylabel=r'$w(t)$')
            if tlims != 0:
                axes[2].set(xlim=tlims)
            axes[2].set(ylim=w_lims)

            # Also Psi.
            axes[3].plot(x, data['psi'], lw=1, color="red")
            axes[3].set(ylabel=r'$\psi(t)$')
            if tlims != 0:
                axes[3].set(xlim=tlims)
            if field_lims != False:
                axes[3].set(ylim=field_lims)



            # We're also going for density-related stuff.
            halfpsiprimesq, massive_gamma = 0.5*data['psiprime']**2, data['massive_gamma']
            halfpsiprimesq_act, massive_gamma_act = halfpsiprimesq/(primary.planck_time**2), massive_gamma/(primary.planck_time**2)
            halfpsiprimesq_act, massive_gamma_act = halfpsiprimesq_act/primary.G, massive_gamma_act/primary.G
            total_act = halfpsiprimesq_act + massive_gamma_act
            scalar_parameter = total_act/(3*((data['H_normal']*1e3/primary.mpc)**2)/(8*np.pi*primary.G))
            axes[4].plot(x, scalar_parameter, lw=1, color="red")
            axes[4].plot(x, rolling_parameter, lw=1, color="blue")

            # Also plot down the current epoch results for density parameter (13.7/etc)
            axes[4].scatter(primary.age_R_ONE, primary.dens_vac, marker="x", s=2, color="red")
            axes[4].scatter(primary.age_R_ONE, primary.dens_mat, marker="x", s=2, color="blue")

            axes[4].set(ylabel=r'$\Omega$' + " (t)")
            if tlims != 0:
                axes[4].set(xlim=tlims)
            if density_lims != False:
                axes[4].set(ylim=density_lims)
            legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\Omega_\Lambda$'),
                               Patch(facecolor='blue', edgecolor='black', label=r'$\Omega_M$')]
            axes[4].legend(handles=legend_elements, loc='upper right')



            # Plot the actual densities, too.
            axes[5].plot(x, halfpsiprimesq_act, lw=1, color="red")
            axes[5].plot(x, massive_gamma_act, lw=1, color="blue")
            axes[5].plot(x, total_act, lw=1, color="purple")
            axes[5].plot(x, rolling_density, lw=1, color="black")
            axes[5].set(ylabel='Densities' + " (t)")
            if tlims != 0:
                axes[5].set(xlim=tlims)
            if gamma_lims != False:
                axes[5].set(ylim=gamma_lims)
            legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\frac{1}{2}(\dot{\psi})^2$'),
                               Patch(facecolor='blue', edgecolor='black', label=r'$\Gamma_\psi$'),
                               Patch(facecolor='purple', edgecolor='black', label=r'$\rho_\phi$'),
                               Patch(facecolor='black', edgecolor='black', label=r'$\frac{\rho_m}{m_p^2}$')]
            axes[5].legend(handles=legend_elements, loc='upper right')

            # Also want to plot the density parameter over time.
            # denmat = matter_0

            axes[5].set(xlabel="t / Gyr")

            try:
                os.mkdir(astrowrapper.rootdir + "\\Figures\\" + alpha_format)
            except:
                pass
            try:
                os.mkdir(astrowrapper.rootdir + "\\Figures\\" + alpha_format + "\\" + psiprimeformat)
            except:
                pass

            plt.savefig(astrowrapper.rootdir + "\\Figures\\" + alpha_format + "\\" + psiprimeformat + "\\" + self.set + ".png", dpi=300)
            #plt.show()
            plt.clf()
            plt.close()

            fig, axis = plt.subplots(1)
            fig.suptitle(full_title)
            redshiftplusone = (1/data['R'])
            axis.plot(redshiftplusone, total_act, lw=1, color="red")
            axis.plot(redshiftplusone, rolling_density, lw=1, color="black")
            axis.set(ylabel=r'$\log{\rho}$',
                     xlabel="z + 1")
            plt.gca().invert_xaxis()
            legend_elements = [Patch(facecolor='red', edgecolor='black', label=r'$\frac{1}{2}(\dot{\psi})^2 + \Gamma_\psi$'),
                               Patch(facecolor='black', edgecolor='black', label=r'$\frac{\rho_m}{m_p^2}$')]
            axis.legend(handles=legend_elements, loc='upper right')

            axis.set_xscale('log')
            axis.set_yscale('log')

            try:
                os.mkdir(astrowrapper.rootdir + "\\Figures\\" + alpha_format + "\\" + psiprimeformat + "\\Density")
            except:
                pass

            plt.savefig(astrowrapper.rootdir + "\\Figures\\" + alpha_format + "\\" + psiprimeformat + "\\Density\\" + self.set + "_DENSITY.png", dpi=300)
            #plt.show()
            plt.clf()
            plt.close()
        else:
            print("Model not satisfied.")


    # Main runtime for multiprocessing
    def main_multi(self):
        try:
            table = self.integrator()
            self.grapher([-0.3, 56], False, 0.1726673736404837, [-1.1,1.1], False, False, [0, 2e-26], table) # hard coded sorry mate
        except:
            pass


#scalar_universe = primary_scalar(*primary.primary_scalar_params)
#owo = scalar_universe.integrator()
#scalar_universe.grapher([-0.2, 14], [0,1], primary.t_shift, [-2,2], False, False, [0, 1e-26],owo)  # t_lims, r_lims, t_shift_for_analytic, w_lims, omega_lims, field_lims, gamma lims




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
miscellaneous_simrunner()

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



