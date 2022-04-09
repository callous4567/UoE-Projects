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
from fast_cahn import fast_cahn, phitegrate, free_energy


# Cahn-Hilliard Equation Simulator. a=b in this case.
class cahn(object):
    def __init__(self, a, k, dx, dt, M, iniphi, nx, noise, maxsweeps, binrate):
        """
        Cahn-Hilliard Process in 2D abiding the equation
        dphi/dt = M(grad**2)mu
        Where
        mu = -a(phi - phi**3) - k(grad**2)phi
        As detailed for Modelling/Visualization.

        :param a: see: ref
        :param k: see: ref
        :param dx: spatial discretization
        :param dt: temporal discretization
        :param M: see: ref
        :param iniphi: initial value (flat) for composite order matrix
        :param nx: grid size
        :param noise: gaussian noise scale (loc=0 default)
        :param maxsweeps: max sweep to run until
        :param binrate: binrate for viewing with run()

        """
        # Sim pars
        self.a, self.k, self.dx, self.dt, self.M = a, k, dx, dt, M
        # Matpars
        self.nx = nx
        # Initial condition for phi
        self.phimat = np.ones((nx,nx))*iniphi + np.random.normal(loc=0, scale=noise, size=(nx,nx))
        # Free Energy
        self.fe = None
        # Sweep
        self.sweep = 0
        self.max_sweeps = maxsweeps
        self.binrate = binrate

    # Iterate forward by one sweep
    def fast_cahn(self):
        self.phimat = fast_cahn(self.phimat, self.a, self.k, self.dx, self.dt, self.M)

    # Run the simulator graphically (plots the free energy/etc.)
    def run(self):
        # Interactive On
        plt.ion()

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8,8))
        im = ax.imshow(self.phimat, animated=True, cmap="bwr", aspect='equal', vmin=-1, vmax=1)
        ax.set(xlim=[0, self.nx - 1],
               ylim=[0, self.nx - 1])

        # Set text.
        t1 = ax.text(1, 1, str(self.sweep), color="black", fontsize=20)
        ax.set_title(("f = {0:.1f}").format(free_energy(self.phimat, self.a, self.k, self.dx)))

        while self.sweep < self.max_sweeps:
            self.fast_cahn()
            self.sweep += 1
            if self.sweep % self.binrate == 0:
                im.set_array(self.phimat)
                ax.set_title(("f = {0:.1f}").format(free_energy(self.phimat, self.a, self.k, self.dx)))
                t1.set_text(str(self.sweep))
                fig.canvas.draw()
                fig.canvas.flush_events()

        # All done.
        plt.close()

    # Get the free energy over time
    def freeplot(self):
        sweep = [self.sweep]
        free = [free_energy(self.phimat, self.a, self.k, self.dx)]

        while self.sweep < self.max_sweeps:
            self.fast_cahn()
            self.sweep += 1
            sweep.append(self.sweep)
            free.append(free_energy(self.phimat, self.a, self.k, self.dx))

        # plot
        sweep = [d*self.dt for d in sweep]
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.scatter(sweep, free, marker='x', color='red')
        ax.set(xlabel=r'$t$',
               ylabel='f')
        ax.grid(which='major', color='pink')
        plt.savefig("free-energy-test.png", dpi=300)


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
        self.time_delay("Yōkoso!!! Welcome to this Cahn-Hilliard Simulator!. \n"
                        "You will now be asked for a few parameters. Please give them. \n"
                        "Please note that the code is a bit optimized for multiple sequential runs (i.e. parallel) \n"
                        "Due to this, a one-time-run will incur a @jit compile cost, compared to regular python. \n"
                        "Now, onto the parameters!!!!")
        print()
        self.time_delay("a")
        a = float(input())
        print()
        self.time_delay("k")
        k = float(input())
        print()
        self.time_delay("dx")
        dx = float(input())
        print()
        self.time_delay("dt")
        dt = float(input())
        print()
        self.time_delay("M")
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
            print(e)

checkpoint().run()

