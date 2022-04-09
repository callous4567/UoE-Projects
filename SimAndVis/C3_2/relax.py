import copy
import itertools
import pickle
import sys
import time
import timeit
import warnings
from astropy.table import Table
from matplotlib import rc, cm
import os
from matplotlib.patches import Patch
import hdfutils
import numpy as np
import matplotlib.pyplot as plt
from SimAndVis.C3_2.fast_relax import fast_gs, converged, potential


class charge_generator(object):
    def __init__(self, shape):
        self.shape = shape

    # Generate a point charge at the centre of charge q (we assume discretization is unity thus this equals density.)
    def point_charge(self, q):
        rho = np.zeros(self.shape)
        rho_mid = [int(d/2) for d in self.shape]
        rho[rho_mid[0], rho_mid[1], rho_mid[2]] = q
        return rho

    # Generate one instead at some point elsewhere
    def point_charge_point(self, q, coordinate):
        rho = np.zeros(self.shape)
        rho[coordinate[0], coordinate[1], coordinate[2]] = q
        return rho

    # Random spread of point charges
    def charge_spread(self, max_q, num):
        rho = np.zeros(self.shape)
        coords_x, coords_y, coords_z = np.random.randint(0, self.shape[0], num), \
                                       np.random.randint(0, self.shape[1], num), \
                                       np.random.randint(0, self.shape[2], num)
        coords = np.array([coords_x, coords_y, coords_z]).T
        chargs = np.random.normal(0, max_q, num)
        for coord, charg in zip(coords, chargs):
            rho[coord[0],coord[1],coord[2]] = charg

        return rho

class relaxed_poisson(object):
    # Generic Poisson solver for static charge distributions (inside cuboids- not necessarily cubes.)
    def __init__(self, boundary, rho, converge_frac, converge_rate, max_n):
        self.boundary, self.converge_frac = boundary, converge_frac # dirichlet condition & converge frac
        self.converge_rate = converge_rate # every X steps to check for convergence
        self.rho = rho # initial charge distribution, held static
        self.shape = np.shape(rho) # grid size
        self.pot_dist = np.zeros_like(rho) # potential grid.
        self.has_converged = False
        self.n = 1
        self.max_n = max_n

        """
        I have set potgrid to zero...
        Evaluating potgrid from rho would solve the problem in this case (the point is to solve this numerically)
        The way I understand it, we're trying to get the steady state potgrid- no analytic solution to it. 
        """

    # Plot x,y slice over time (ideal for point potential.)
    def run_sim_plot_point(self, binrate, sqr_dist_relax):

        # Interactive on
        plt.ion()

        # Evaluate the analytical solution
        analytical = potential(self.rho, sqr_dist_relax)

        # Analytically estimate the ideal colour limits

        # Set up figure, axes, etc
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
        zmid = int(self.shape[2]/2)
        im = axs[0].imshow(self.pot_dist[:,:,zmid], animated=True, cmap="bwr", aspect='equal', vmin=-0.1, vmax=0.1)
        axs[0].set(xlim=[0, self.shape[0] - 1],
                   ylim=[0, self.shape[1] - 1])

        # Set text.
        t1 = axs[0].text(1, 1, str(self.n), color="black", fontsize=20)

        # Also add in for the analytical solution
        im_analytic = axs[1].imshow(analytical[:,:,zmid], animated=True, cmap="bwr", aspect='equal', vmin=-0.1, vmax=0.1)
        axs[1].set(xlim=[0, self.shape[0] - 1],
                   ylim=[0, self.shape[1] - 1])
        axs[1].text(1,1,"Analytic",color='black',fontsize=20)

        # Until it's converged. Keep track of the "last converge set" to check
        last_converge_pot = copy.copy(self.pot_dist)

        # Iterate over. Terminate when new_converge_set near match to last, by fraction converge_frac
        while self.n < self.max_n:

            # If in a "converge set" run and  check for convergence, too
            if self.n % self.converge_rate == 0:

                # Get the current converge pot
                current_converge_pot, n = fast_gs(self.rho, self.pot_dist, self.shape, self.boundary, self.n)
                self.has_converged = converged(last_converge_pot, current_converge_pot, self.converge_frac)
                if self.has_converged == True:
                    print("Converged on ", self.n)
                    plt.close()
                    break
                else:
                    self.pot_dist = current_converge_pot
                    last_converge_pot = copy.copy(current_converge_pot)
                    self.n = n

            # Not in a "converge set" so run normally
            else:
                self.pot_dist, self.n = fast_gs(self.rho, self.pot_dist, self.shape, self.boundary, self.n)

            # Binrate check, too
            if self.n % binrate == 0:
                im.set_array(self.pot_dist[:,:,zmid])
                t1.set_text(str(self.n))
                fig.canvas.draw()
                fig.canvas.flush_events()

        if self.n > (self.max_n - 3):
            print("Failed to converge, sorry!")
            plt.close()


# Class for handling user input (i.e. checkpoint.)
class checkpoint(object):
    def __init__(self):
        self.delay_max = 2e-3
        self.rng = np.random.default_rng()

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        if __name__ == "__main__":
            print()
            for c in text:
                sys.stdout.write(c)
                sys.stdout.flush()
                resttime = self.rng.uniform(0.0001, self.delay_max)
                time.sleep(resttime)
            print()
        else:
            pass

    # User Input Specification
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to my Poisson Solving Thingy! \n" 
                        "This particular model will let you do a bunch of point charges \n"
                        "Please give parameters as outlined in the documentation!")
        print()
        self.time_delay("grid size x")
        nx = int(input())
        print()
        self.time_delay("grid size y")
        ny = int(input())
        print()
        self.time_delay("grid size z")
        nz = int(input())
        print()
        self.time_delay("number of point charges")
        ncharges = int(input())
        print()
        self.time_delay("magnitude of charges (recommended near unity.)")
        magcharge = float(input())
        print()
        self.time_delay("boundary value")
        boundary = float(input())
        print()
        self.time_delay("convergence fraction")
        confrac = float(input())
        print()
        self.time_delay("convergence rate")
        conrate = int(input())
        print()
        self.time_delay("max iterations")
        maxiter = int(input())
        print()
        self.time_delay("binrate for imaging (i.e. every X frames)")
        binrate = int(input())
        print()
        self.time_delay("relaxation distance for analytic solution (recommended 1e-3)")
        relaxation = float(input())
        print()
        shape= [nx,ny,nz]
        charges = [magcharge,ncharges]
        relpars = [boundary, confrac, conrate, maxiter]
        simpars = [binrate, relaxation]
        return shape, charges, relpars, simpars

    # Run
    def run(self):
        try:
            print("Do you just want to use some default settings? y/n")
            he_says = str(input())
            if he_says == "y":
                rho = charge_generator([20,20,5]).charge_spread(2, 20)
                poisson = relaxed_poisson(0, rho, 1e-12, 1000000, int(100e6))
                poisson.run_sim_plot_point(10000,0.00000001)
            else:
                pars = self.user_input()
                rho = charge_generator(pars[0]).charge_spread(*pars[1])
                poisson = relaxed_poisson(pars[2][0], rho, *pars[2][1:])
                poisson.run_sim_plot_point(*pars[3])
        except Exception as e:
            self.time_delay("An error occurred...")
            print(e)


"""
# Class for handling user input (i.e. checkpoint.)
class checkpoint(object):
    def __init__(self):
        self.delay_max = 2e-3
        self.rng = np.random.default_rng()

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        if __name__ == "__main__":
            print()
            for c in text:
                sys.stdout.write(c)
                sys.stdout.flush()
                resttime = self.rng.uniform(0.0001, self.delay_max)
                time.sleep(resttime)
            print()
        else:
            pass

    # User Input Specification
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to this Poisson Simulator!. \n"
                        "You will now be asked for a few parameters. Please give them. \n"
                        "Please note that the code is a bit optimized for multiple sequential runs (i.e. parallel) \n"
                        "Due to this, a one-time-run will incur a @jit compile cost, compared to regular python. \n"
                        "Now, onto the parameters!!!!")
        print()
        self.time_delay("Boundary Value")
        boundary = float(input())
        print()
        self.time_delay("Convergence Fraction (see documentation)")
        convergence_frac = float(input())
        print()
        self.time_delay("Convergence Rate (see documentation)")
        convergence_rate = int(input())
        print()
        self.time_delay("Max iterations before abandoning convergence")
        max_n = int(input())
        print()
        self.time_delay("")
        M = float(input())
        print()
        self.time_delay("iniphi")
        iniphi = float(input())
        print()
        self.time_delay("grid size")
        nx = int(input())
        print()
        self.time_delay("maximum number of sweeps")
        maxsweeps = int(input())
        print()
        self.time_delay("binrate for imaging (i.e. every X frames)")
        binrate = int(input())
        print()
        self.time_delay("Gaussian noise loc of 0 with dispersion of 0.01 will be used by default.")
        return a,k,dx,dt,M,iniphi,nx,0.1,maxsweeps,binrate

    # Run
    def run(self):
        try:
            pars = self.user_input()
            model = cahn(*pars)
            model.run()
        except Exception as e:
            self.time_delay("An error occurred... \n"
                            "You probably don't have the correct dependencies. \n"
                            "If the error regards the existence of LaTex: delete lines 21,22 \n"
                            "If the error regards missing packages, please install them.\n")
            print(e) """



checkpoint().run()
