import copy
import multiprocessing
import os
import sys
import time

from astropy.table import Table
from matplotlib.patches import Patch
from scipy.optimize import curve_fit
from scipy.signal import convolve
import numpy as np
import matplotlib.pyplot as plt

import hdfutils
from SimAndVis.C3_2 import multiprocessing_functions
from SimAndVis.C3_2.fast_electrorelax import electrofield_nonnormal
from SimAndVis.C3_2.fast_relax import potential_stagger, potential, meshgrid, fast_gs, converged, fast_sor, \
    fast_sor_inplace
from SimAndVis.C3_2.fast_electrorelax import electrofield


class charge_generator(object):
    def __init__(self, shape):
        self.shape = shape
        self.charges = None
        self.charges_coords = None

    # Generate a point charge at the centre of charge q (we assume discretization is unity thus this equals density.)
    def point_charge(self, q):
        rho = np.zeros(self.shape)
        rho_mid = [int(d/2) for d in self.shape]
        rho[rho_mid[0], rho_mid[1], rho_mid[2]] = q
        self.charges = np.array([q])
        self.charges_coords = np.array([rho_mid])
        return rho

    # Generate one instead at some point elsewhere
    def point_charge_point(self, q, coordinate):
        rho = np.zeros(self.shape)
        rho[coordinate[0], coordinate[1], coordinate[2]] = q
        self.charges = np.array([q])
        self.charges_coords = np.array(coordinate)
        return rho

    # Random spread of point charges
    def charge_spread(self, max_q, num, pad):
        rho = np.zeros(self.shape)
        coords_x, coords_y, coords_z = np.random.randint(pad, self.shape[0] - pad, num), \
                                       np.random.randint(pad, self.shape[1] - pad, num), \
                                       np.random.randint(pad, self.shape[2] - pad, num)
        coords = np.array([coords_x, coords_y, coords_z]).T
        chargs = np.random.normal(0, max_q, num)
        self.charges = chargs
        self.charges_coords = coords
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

    # Analytical potential generator. Kernel width should be at least double the padding (to conserve charge.)
    def gen_analypot(self, gauss_kernel_width, zmid, sqr_dist_relax, stagger, nonnorm=False):

        if gauss_kernel_width == 0:

            analy_rho = self.rho

        else:

            # Specify empty gaussian. Make sure it's bigger than the largest binsize.
            kernel_width = gauss_kernel_width  # in number of cells
            gauss_array = np.zeros_like(self.rho)
            gaussi, gaussj, gaussk = np.shape(gauss_array)
            midi, midj, midk = int((gaussi + 1) / 2), int((gaussj + 1) / 2), int((gaussk + 1) / 2)
            # Craft the gaussian! Unnormalized.
            gauss = lambda ii, jj, kk: np.exp(
                (-1 / (2 * (kernel_width ** 2))) * ((ii - midi) ** 2 + (jj - midj) ** 2 + (kk - midk) ** 2))
            for i in range(gaussi):
                for j in range(gaussj):
                    for k in range(gaussk):
                        gauss_array[i, j, k] = gauss(i, j, k)
            gauss_array /= np.sum(gauss_array)
            gaussian = gauss_array

            # Convolve our density matrix with the Gaussian (with the gauss kernel width appropriate.)
            analy_rho = convolve(self.rho, gaussian, mode='same')

        if stagger != None:
            analytical = potential_stagger(analy_rho, stagger)
            if nonnorm == False:
                ana_field_vecs = electrofield(analytical, zmid, self.boundary)
            else:
                ana_field_vecs = electrofield_nonnormal(analytical, zmid, self.boundary)

        else:
            if sqr_dist_relax != None:
                analytical = potential(analy_rho, sqr_dist_relax)
                if nonnorm == False:
                    ana_field_vecs = electrofield(analytical, zmid, self.boundary)
                else:
                    ana_field_vecs = electrofield_nonnormal(analytical, zmid, self.boundary)

        # Return
        return analytical, ana_field_vecs

    # Plot x,y slice over time (ideal for point potential.)
    def run_sim_plot_point(self, binrate, gauss_kernel_width, sqr_dist_relax, stagger):


        # Interactive on
        plt.ion()

        # Midpoint
        zmid = int(self.shape[2]/2)

        # Get meshgrid for quiver
        shape = np.shape(self.rho)
        coord_slice = meshgrid(shape[0], shape[1])

        # Evaluate the analytical solution.
        analytical, ana_field_vecs = self.gen_analypot(gauss_kernel_width, zmid, sqr_dist_relax, stagger)

        # Set up figure, axes, etc
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(8,4))

        # Plot the electric field, also
        field_vecs = np.ones((2, shape[0], shape[1]))/np.sqrt(2)
        # Plot the potential start
        quiver = axs[0].quiver(coord_slice[0], coord_slice[1], field_vecs[0], field_vecs[1], color='black')
        im = axs[0].imshow(self.pot_dist[:,:,zmid], animated=True, cmap="bwr", aspect='equal', vmin=-0.1, vmax=0.1,
                           origin='lower')
        t1 = axs[0].text(0, 0, str(self.n), color="black", fontsize=20)

        plt.savefig(str(self.n) + ".png", dpi=150)

        # Also add in for the analytical solution
        # noinspection PyUnboundLocalVariable
        quiver_analytic = axs[1].quiver(coord_slice[0], coord_slice[1], ana_field_vecs[0], ana_field_vecs[1],
                                        color='black')
        im_analytic = axs[1].imshow(analytical[:,:,zmid], animated=True, cmap="bwr", aspect='equal', vmin=-0.1,
                                    vmax=0.1, origin='lower')
        axs[1].text(0,0,"Analytic",color='black',fontsize=20)

        # Also plot some residuals
        im_residual = axs[2].imshow(self.pot_dist[:,:,zmid]-analytical[:,:,zmid], animated=True,
                                    cmap="bwr", origin='lower',
                                    aspect='equal', vmin=-0.1, vmax=0.1)
        axs[2].text(0,0,"Residual",color='black',fontsize=20)

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
                field_vecs = electrofield(self.pot_dist, zmid, self.boundary)
                quiver.set_UVC(field_vecs[0], field_vecs[1])
                im.set_array(self.pot_dist[:,:,zmid])
                t1.set_text(str(self.n))
                im_residual.set_array(self.pot_dist[:,:,zmid]-analytical[:,:,zmid])
                fig.canvas.draw()
                fig.canvas.flush_events()
                plt.savefig(str(self.n) + ".png", dpi=150)

        if self.n > (self.max_n - 3):
            print("Failed to converge, sorry!")
            plt.close()

    # Run the sim without plotting (with/without SOR) just to produce a final slice
    def runsim(self, gauss_kernel_width, sqr_dist_relax, stagger, SOR_value, nonnorm=True):

        """
        The analytical potential is the same shape as the density field. The field vectors are a slice through the
        midplane of the potential grid, calculated appropriately. Set the SOR_value to unity for standard Gauss-Seidel.
        SOR update rule is in-place: too high will cause divergence/instability in simulation. Not-inplace exists-
        see fast_relax.

        :return: analytical, ana_field_vecs, self.pot_dist, field_vecs, self.n
        """

        # Midpoint
        zmid = int(self.shape[2] / 2)

        # Evaluate the analytical solution.
        analytical, ana_field_vecs = self.gen_analypot(gauss_kernel_width, zmid, sqr_dist_relax, stagger, nonnorm)

        # Until it's converged. Keep track of the "last converge set" to check
        last_converge_pot = copy.copy(self.pot_dist)

        # Iterate over. Terminate when new_converge_set near match to last, by fraction converge_frac
        while self.n < self.max_n:

            # If in a "converge set" run and  check for convergence, too
            if self.n % self.converge_rate == 0:

                # Get the current converge pot
                current_converge_pot, n = fast_sor_inplace(self.rho, self.pot_dist, self.shape, self.boundary, self.n, SOR_value)
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
                self.pot_dist, self.n = fast_sor_inplace(self.rho, self.pot_dist, self.shape, self.boundary, self.n, SOR_value)

        if self.n > (self.max_n - 3):
            print("Failed to converge, sorry!")
            plt.close()

        # Get the final field vectors for the final potential
        field_vecs = electrofield_nonnormal(self.pot_dist, zmid, self.boundary)

        # Returns
        return analytical, ana_field_vecs, self.pot_dist, field_vecs, self.n

    # Create plots radially from a point charge (assumed single charge) using charge_generator
    def radiplot(self, gauss_kernel_width, sqr_dist_relax, stagger, SOR_value, charge_object):

        # Generate
        analytical, ana_field_vecs, pot_dist, field_vecs, n = self.runsim(gauss_kernel_width,
                                                                          sqr_dist_relax,
                                                                          stagger, SOR_value)

        # Grab coord/charge
        charge, coord = charge_object.charges[0], charge_object.charges_coords[0]
        centre = coord[0]

        # Shape
        shape = np.shape(analytical)
        yy = np.arange(0, shape[0], 1)

        # Get the magnitude of the field vectors (this is an XY plane slice through the midplane, i.e. no z component.)
        ana_field_vecs, field_vecs = np.sqrt(ana_field_vecs[0]**2 + ana_field_vecs[1]**2), \
                                     np.sqrt(field_vecs[0]**2 + field_vecs[1]**2)


        # Get a radial slice for the charge along the i-axis for the potential
        anay, poty = analytical[:,coord[1],coord[2]], pot_dist[:,coord[1],coord[2]]

        # Get a radial slice for the charge along the i-axis for the electric field strength (direction skipped.)
        anafy, fieldy = ana_field_vecs[:,coord[1]], field_vecs[:,coord[1]]

        # Set up figure
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(16,8))

        fig.subplots_adjust(hspace=0.3)

        axs[0,0].scatter(yy, anay, color='green', marker='x')
        axs[0,0].scatter(yy, poty, color='red', marker='x')
        axs[0,1].scatter(yy, anafy, color='green', marker='x')
        axs[0,1].scatter(yy, fieldy, color='red', marker='x')
        axs[0,0].grid(which='major', color='pink')
        axs[0,1].grid(which='major', color='pink')
        axs[0,0].set(ylim=[0,0.6])
        axs[0,1].set(ylim=[0,0.6])
        axs[0,0].set(xlim=[2, shape[0]-2])
        axs[0,1].set(xlim=[2, shape[0]-2])

        # Also evaluate potential/etc analytically from a functional form
        def inverse_analytic(x, Q, power):
            return (Q/4)*(1/np.pi)*(1/np.abs((x - centre))**power)

        i_range = np.linspace(0, shape[0]-1, 1000)
        axs[0,0].plot(i_range, inverse_analytic(i_range, charge, 1), color='lime', lw=0.5)
        axs[0,1].plot(i_range, inverse_analytic(i_range, charge, 2), color='lime', lw=0.5)

        # Set up legend elements
        legend_elements = [Patch(edgecolor='black', facecolor='red', label='Sim'),
                           Patch(edgecolor='black', facecolor='green', label='Analytic (Numerical)'),
                           Patch(edgecolor='black', facecolor='lime', label='Analytic (Functional)')]
        axs[0,0].legend(handles=legend_elements, loc='upper right')
        axs[0,1].legend(handles=legend_elements, loc='upper right')
        axs[0,0].set(xlabel=r'$i$', ylabel=r'$\phi$')
        axs[0,1].set(xlabel=r'$i$', ylabel=r'$|\vec{E}|$')

        # Also set up a plot to try and see if the data follows square laws/etc as hoped for...
        yy_distances = yy[centre + 2:] - centre
        distance_ratios = yy_distances/yy_distances[0]
        distance_logs = np.log(distance_ratios)
        pot_edge = poty[centre + 2:]
        pot_ratios = pot_edge/pot_edge[0]
        potlogs = np.log(pot_ratios)
        field_edge = fieldy[centre + 2:]
        field_ratios = field_edge/field_edge[0]
        fieldlogs = np.log(field_ratios)

        # Standard polynomial fit
        def poly(r, mag, power):
            return mag*(r**power)

        # Test fit
        potfit = curve_fit(poly, distance_logs, potlogs)[0]
        fieldfit = curve_fit(poly, distance_logs, fieldlogs)[0]

        # Produce a plot of the logs/etc
        axs[1,0].scatter(distance_logs, potlogs, color='blue')
        axs[1,0].plot(distance_logs, -distance_logs, color='lime')
        axs[1,0].plot(distance_logs, poly(distance_logs, *potfit), color='red')
        axs[1,1].scatter(distance_logs, fieldlogs, color='blue')
        axs[1,1].plot(distance_logs, -2*distance_logs, color='lime')
        axs[1,1].plot(distance_logs, poly(distance_logs, *fieldfit))
        # Set up legend elements
        legend_elements = [Patch(edgecolor='black', facecolor='blue', label='Sim'),
                           Patch(edgecolor='black', facecolor='lime', label='Analytic (Functional)'),
                           Patch(edgecolor='black', facecolor='red', label='Polynomial Fit')]
        axs[1,0].legend(handles=legend_elements, loc='upper right')
        axs[1,1].legend(handles=legend_elements, loc='upper right')
        axs[1,0].set(xlabel=r'$\log(\frac{\Delta{x}}{\Delta{x}_0})$', ylabel=r'$\log(\frac{\phi}{\phi_0})$')
        axs[1,1].set(xlabel=r'$\log(\frac{\Delta{x}}{\Delta{x}_0})$', ylabel=r'$\log(\frac{|\vec{E}|}{|\vec{E_0}|})$')
        axs[1,0].grid(which='major', color='pink')
        axs[1,1].grid(which='major', color='pink')
        axs[1,0].set_title("Polynomial Fit of " + r'$m{r}^{a}$' + " where m,a " + r'$=$' + str(potfit))
        axs[1,1].set_title("Polynomial Fit of " + r'$m{r}^{a}$' + " where m,a " + r'$=$' + str(fieldfit))


        plt.savefig("example_fieldplot.png", dpi=300)
        plt.show()

        # Dump all the data
        writer = hdfutils.hdf5_writer(os.getcwd(), "electric_monopole.hdf5")
        columns = ["distance", "potential", "field"]
        data = np.array([yy_distances, pot_edge, field_edge]).T
        table = Table(data=data, names=columns)
        writer.write_table("CHARGE_EQUALS_TWO", "data", table)

    # Runsim but without any trimmings. Only returns the value of $n$
    def converge_runsim(self, SOR_value):

        # Midpoint
        zmid = int(self.shape[2] / 2)

        # Until it's converged. Keep track of the "last converge set" to check
        last_converge_pot = copy.copy(self.pot_dist)

        # Iterate over. Terminate when new_converge_set near match to last, by fraction converge_frac
        while self.n < self.max_n:

            # If in a "converge set" run and  check for convergence, too
            if self.n % self.converge_rate == 0:

                # Get the current converge pot
                current_converge_pot, n = fast_sor_inplace(self.rho, self.pot_dist, self.shape, self.boundary,
                                                           self.n, SOR_value)
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
                self.pot_dist, self.n = fast_sor_inplace(self.rho, self.pot_dist, self.shape, self.boundary, self.n,
                                                         SOR_value)

        if self.n > (self.max_n - 3):
            print("Failed to converge, sorry!")
            plt.close()

        # Returns
        return self.n

# Class for handling user input (i.e. checkpoint.)
class checkpoint(object):
    def __init__(self):
        self.delay_max = 2e-3
        self.rng = np.random.default_rng()

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        print()
        for c in text:
            sys.stdout.write(c)
            sys.stdout.flush()
            resttime = self.rng.uniform(0.0001, self.delay_max)
            time.sleep(resttime)
        print()


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
        self.time_delay("You can get analytical solution either by a staggered grid, or by relaxation distance. s or r")
        which = str(input())
        if which == "s":
            self.time_delay("grid stagger (recommended 0.01)")
            stagger = float(input())
            relaxation = None
            print()
        elif which == "r":
            self.time_delay("relaxation distance for analytic solution (recommended 1e-12)")
            relaxation = float(input())
            stagger = None
            print()
        else:
            print("You didn't select s or r ... will use a default relaxation of 1e-12")
            relaxation = 1e-12
            stagger = None


        shape= [nx,ny,nz+2]
        charges = [magcharge,ncharges,0]
        relpars = [boundary, confrac, conrate, maxiter]
        simpars = [binrate, 1, relaxation, stagger]
        return shape, charges, relpars, simpars

    # Run
    def run(self):
        print("Do you just want to use some default settings? y/n")
        he_says = str(input())
        if he_says == "y":
            rho = charge_generator([32,32,32]).charge_spread(1, 30, 0) # 2, 1)
            poisson = relaxed_poisson(0, rho, 1e-10, 10000000, int(10000e6))
            # noinspection PyTypeChecker
            poisson.run_sim_plot_point(1,0,None, 0.001)
        if he_says == "n":
            pars = self.user_input()
            rho = charge_generator(pars[0]).charge_spread(*pars[1])
            poisson = relaxed_poisson(pars[2][0], rho, *pars[2][1:])
            poisson.run_sim_plot_point(*pars[3])
        else:
            print("You said something wrong.")

# To replicate the radial plot in the PDF I gave
def do_radiplot():
    charge_dist = charge_generator([100, 100, 100])
    rho = charge_dist.point_charge(2)  # 2, 1)
    poisson = relaxed_poisson(0, rho, 1e-3, 2000, int(1e6))
    poisson.radiplot(0, None, 0.001, 1, charge_dist)

# Generate the convergence plot. We're doing it for a convergence percentage of 0.1%.
def do_sorplot():

    if __name__ == "__main__":

        # Charge distribution
        charge_dist = charge_generator([40, 40, 40])
        rho = charge_dist.point_charge(2)  # 2, 1)

        # Set up sorrange
        sor_range = np.linspace(1, 3, 50)
        sorsets = []
        for i in sor_range:
            sorsets.append([i, rho])

        pool = multiprocessing.Pool(8)
        n_vals = pool.map(multiprocessing_functions.do_sor, sorsets)
        pool.close()

        # Make plot
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.set(ylim=[0,1000])
        ax.grid(which='major', color='pink')
        ax.plot(sor_range, n_vals, color='red')
        ax.set(xlabel=r'$\omega$',
               ylabel="convergence / n")
        plt.savefig("sorvergence.png", dpi=300)
        plt.show()


# Run checkpoint
checkpoint().run()
#do_radiplot()
#do_sorplot()
