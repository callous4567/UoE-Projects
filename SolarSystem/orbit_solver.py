import sympy
from scipy.optimize import least_squares
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
# Can also be done using symfit and sympy. Will apply them to 3D fit if 2D fit isn't accurate enough.
from symfit import parameters, variables, Fit
from sympy import *

import numpy as np
import random

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D # Redundant here but will be used in future and thus is kept for my use.
import matplotlib.ticker as pltick

import sys
import time
import os

# Ellipse class.
# Ellipse object takes a, b, x_0, y_0, z_0 gamma, beta, alpha, noise_fraction
# gamma, beta, alpha are the extrinsic euler angles defined on the wikipedia page.
# This is just to generate noisy data with noise_fraction being a ratio of a to use as "noise" on all final x,y,z values.

"""
This is a note to future me. Some stuff to include on future rendition for report:
- The simulator takes a timestep. Modify the printing software to print timestep in, then read it.
- Modify for partial orbit simulation parameters: do for report. 
"""

# Modify for new import parameters at a later date.
# Returns the [x_list, y_list, z_list] and also [timestep, stepsplit] array of simulation parameters.
class file_import_tool(object):
    # Reads file in, skips first line, gives XYZ data AND the timestep, which is currently hard coded for some reason.

    def orbit_porter(self, filename_and_extension):
        file = open(filename_and_extension, 'r')
        filelines = file.readlines()
        # LINE ZERO ANALYSIS. Takes the timestep and "splitstep," i.e. how many steps that are skipped.
        x, y, z = [],[],[]

        # OPTIONAL AND FOR FUTURE USE. IGNORE FOR NOW (will modify later)
        #timeline = filelines[0]
        #timesplit = timeline.split("###") # Indicate splitting between constants/etc with ###.
        #timestep, stepsplit, orbitname = timesplit[0], timesplit[1], timesplit[2]

        # TEMPORARY HARD CODING
        orbitname = "Example"
        timestep, stepsplit = 1, 10 # Days and steps. :/
        filelines = filelines[1:]
        for line in filelines:
            linesplit = line.split(",")
            linesplit = [float(d) for d in linesplit]
            x.append(linesplit[0]), y.append(linesplit[1]), z.append(linesplit[2])
        xyz = [x,y,z]

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')

        graph = ax1.plot(x, y, z, 'o', lw=0.1, markersize=0.5)

        print(len(x))

        plt.show()

        file.close()

        return xyz, [timestep, stepsplit], orbitname


class ellipse(object):
    def __init__(self, semimajor, semiminor, x0, y0, z0, gamma, beta, alpha, noise_fraction):
        self.a = semimajor
        self.b = semiminor
        self.translation = [x0, y0, z0]
        self.gamma = gamma
        self.beta = beta
        self.alpha = alpha
        self.noise = noise_fraction
        # Note on angles.
        # first z by gamma
        # then x by beta
        # finally z by alpha

    # Generates 2D ellipse centred on 0,0,0
    def twod_gen(self):
        dp = 10000
        theta = np.linspace(0, 2*np.pi, dp)
        x , y ,z = [self.a*np.cos(d) for d in theta], [self.b*np.sin(d) for d in theta], [0 for d in theta]
        return [x, y, z]

    # 3D rotation about z by counterclockwise angle
    def threez_rot(self, angle):
        xprime, yprime, zprime = (np.cos(angle), -1 * np.sin(angle), 0), (np.sin(angle), np.cos(angle), 0), (0, 0, 1)
        z_trix = np.array((xprime, yprime, zprime), dtype=float)
        return z_trix

    # 3D rotation about x by counterclockwise angle
    def threex_rot(self, angle):
        x_prime, y_prime, z_prime = (1, 0, 0), (0, np.cos(angle), -1*np.sin(angle)), (0, np.sin(angle), np.cos(angle))
        x_trix = np.array((x_prime, y_prime, z_prime), dtype=float)
        return x_trix

    # Returns 3D ellipse centred on x0, y0, z0) with rotation and noise.
    def threed_gen(self):
        x, y, z = self.twod_gen()
        # Method for 3D rotation: standard matrix method.
        # First a ztrix by gamma, xtrix by beta, then another ztrix by alpha.
        # Need in vector form, though.
        xyz = np.column_stack((x, y, z))

        gammatrix = self.threez_rot(self.gamma)
        betatrix = self.threex_rot(self.beta)
        alphatrix = self.threez_rot(self.alpha)

        xyz_1 = [np.dot(gammatrix,d) for d in xyz]
        xyz_2 = [np.dot(betatrix,d) for d in xyz_1]
        xyz_3 = np.array([np.dot(alphatrix,d) for d in xyz_2])

        x, y, z = xyz_3.T

        # Add some noise. Add a hash if you want to omit this step. Noise is gaussian.
        a_fraction = self.a * self.noise
        x, y, z = [d + random.gauss(-a_fraction, a_fraction) for d in x], [d + random.gauss(-a_fraction, a_fraction) for d in y], [d + random.gauss(-a_fraction, a_fraction) for d in z]
        # Apply translation
        x, y, z = [d + self.translation[0] for d in x], [d + self.translation[1] for d in y], [d + self.translation[2] for d in z]
        return [x, y, z]

    def ellipse_plot(self):
        xyz = self.threed_gen()

        # Define figure and ax1/ax2 3D/2D representation
        fig = plt.figure(figsize=(10, 10))

        # ax1 for 3D.
        ax1 = fig.add_subplot(111, projection='3d')

        graph = ax1.plot(xyz[0], xyz[1], xyz[2], 'o', lw=0.1, markersize=0.5)

        # Sets limits and enables grid for axes. ax1.set_xlim(min, max) etc if you want lims.
        ax1.grid()

        plt.show()

# Okay. General method.
# First estimate x0, y0, z0 for the ellipse by vector subtraction over the entire list of data. This part lets you also estimate period.
# Translate to 0,0,0.
# Gram-schmidt inverse the ellipse by using it as an XY plane and force it into the 2D.
# Once in 2D, carry out a least-squares fit on the ellipse to estimate the semimajor and semiminor axis.
# The gram-schmidt inverse can be used to estimate the extrinsic euler angles.
# Input the (a,b), initial coordinates x0,y0,z0, and euler angles as first guess for a 3D least squares fit.
class ellipse_solver(object):
    def __init__(self, data, time_data, name):
        self.xyz = data # x y z transposed data, i.e. not vectorized
        self.period = 0 # orbit period in earth days
        self.xyzrow = 0 # vectorized xyz
        self.step = time_data[0]
        self.skip = time_data[1]
        self.transestimate = 0  # estimate for translation
        self.closure_num = 0  # estimate of closure number in xyzrow
        self.twodflat = 0
        self.ab = [0,0]
        self.name = name # For future name implementation when graphing/writing files. May not be used.


    # Calculates/estimates closure period. Sets that to ellipse period. Returns gram vectors.
    # Sets self.period to the ellipse orbit period, taken from the input time list.
    # Returns x',y',z' w/ gram1,gram2,gram3 in array [gram1,gram2,gram3]
    # gram1 x gram2 = gram3.
    def vecsep_calc(self):
        xyz = np.column_stack((self.xyz[0], self.xyz[1], self.xyz[2]))

        vecseps = []
        # Vector separations (indexed) vs. the original
        for d in range(len(xyz)):
            vecseps.append(np.linalg.norm(xyz[d] - xyz[0]))


        hey = np.arange(0, len(vecseps), 1)
        print(hey)

        plt.plot(hey, vecseps)
        plt.show()

        # Set a bound on the minimum vector separation beneath which we have closure.
        # Python searches list horizontally left-right. First point which this is satisfied will provide the estimate.
        # You can set it manually or to a lower index if you're confident in your integrator period +/+ want more accuracy.
        minsearch = vecseps[2]

        # Set range of values to infinity to prevent overlap when estimating first closure
        for i in range(0,150):
            vecseps[i] = np.inf

        closure_num = 0

        # Closure Index
        for num, item in enumerate(vecseps, 0):
            if item <= minsearch:
                self.period = self.step * num # Set the ellipse period :)
                closure_num = num
                break

        # Estimation of the translation for the ellipse.

        xyzsplit = xyz.T
        translation = np.array([np.mean(d) for d in xyzsplit])

        # In case we need it, closed xyz vector list, in row-vector format.
        closed_xyz = xyz[:closure_num]

        closure_num -= int(1)
        print(translation)
        print(closure_num)

        self.xyzrow = np.array(xyz)
        self.transestimate = translation
        self.closure_num = closure_num

    def gram_flattener(self):

        onemax, twomax, twomin = 0.2,0.4,0.3

        closed_xyz, closure_num, translation = self.xyzrow, self.closure_num, self.transestimate

        translated_xyz = [d - translation for d in closed_xyz]

        vector1, vector2 = translated_xyz[int(onemax * closure_num)] - translated_xyz[0], translated_xyz[int(twomax * closure_num)] - \
                           translated_xyz[int(twomin * closure_num)]

        print(vector1)
        print(vector2)

        gram2 = vector2 / np.linalg.norm(vector2)
        # v12_angle = np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
        mid1 = vector1 - (np.dot(vector1, gram2)*gram2)
        gram1 = mid1 / np.linalg.norm(mid1)
        gram3 = np.cross(gram1, gram2)


        print(gram1)
        print(gram2)
        print(gram3)
        # 1x2 = 3. 1,2 are xy, 3 iz z.

        grammatrix = np.array((gram1, gram2, gram3))
        invertrix = np.linalg.inv(grammatrix)
        flattened_xyz = np.array([np.dot(d, invertrix) for d in translated_xyz]).T # x,y,z form / transposed

        #plt.plot(flattened_xyz[0], flattened_xyz[1], 'o', lw=0.2, markersize=0.3)
        #plt.show()
        # We now have a (semi) flat 2D ellipse, rotated obviously.
        # Now then... to fit the ellipse...

        self.twodflat = flattened_xyz



    def twod_optimizer(self):
        x, y, z = self.twodflat

        # use method of ax^2 + by^2 + cxy + dx + ey = 1 and numpy linalg least square fit.
        # https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.lstsq.html
        # ax = b
        # x is solution matrix hence coefficients

        rowvects = []
        for i in range(len(x)):
            rowvects.append([x[i]**2, y[i]**2, x[i]*y[i], x[i], y[i]])
        rowvects = np.array(rowvects)
        univects = [1 for d in rowvects]
        lsq_solve = np.linalg.lstsq(rowvects, univects)[0]


        # solve for ellipse parameters + rotation angle gamma
        # tan(2*gamma) = c/(a-b)
        gamma = 0.5*np.arctan(lsq_solve[2]/(lsq_solve[0] - lsq_solve[1]))
        gammangle = gamma # Remnant for future use. Sorry not sorry for dirty code.

        # rotate ellipse to zero
        gammatrix = lambda gamma: np.array([(cos(gamma), -1*sin(gamma)),
                                            (sin(gamma), cos(gamma))])


        # RED IS UNROTATED
        # GREEN IS ROTATED
        x, y, z = self.twodflat
        xraw, yraw, zraw = self.xyz

        fig = plt.figure()
        ax1 = fig.add_subplot(111, projection='3d')


        graph = ax1.plot(xraw, yraw, zraw, 'o', lw=0.1, markersize=0.5)

        # Matrixify up a gammatrix
        gamma = np.array(gammatrix(gamma), dtype=float)

        # Rotate back the ellipse!
        xy = np.column_stack((x,y))
        xy_new = np.array([np.dot(d, gamma) for d in xy])
        xx, yy = xy_new.T
        graph2 = ax1.plot(xx, yy, [0 for d in xx], 'o', markersize=0.5, c='g')
        plt.show()


        gamma = gammangle
        # Solve for various geometric parameters and then for the ellipse parameters.
        # Format: x^2, y^2, xy, x, y ... all equal to unity (in matrix solve)
        aprime = lsq_solve[0]*(np.cos(gamma))**2 + lsq_solve[2]*np.cos(gamma)*np.sin(gamma) + lsq_solve[1]*(np.sin(gamma))**2
        cprime = lsq_solve[0]*(np.sin(gamma))**2 - lsq_solve[2]*np.cos(gamma)*np.sin(gamma) + lsq_solve[1]*(np.cos(gamma))**2
        dprime = lsq_solve[3]*np.cos(gamma) + lsq_solve[4]*np.sin(gamma)
        eprime = -1*lsq_solve[3]*np.sin(gamma) + lsq_solve[4]*np.cos(gamma)
        fprime = -1
        asquared = (-4*fprime*aprime*cprime + cprime*dprime**2 + aprime*eprime**2)/(4*cprime*aprime**2)
        bsquared = (-4*fprime*aprime*cprime + cprime*dprime**2 + aprime*eprime**2)/(4*aprime*cprime**2)
        a,b = np.sqrt(asquared), np.sqrt(bsquared)

        self.ab = [a,b]

    def main(self):
        self.vecsep_calc()
        self.gram_flattener()
        self.twod_optimizer()
        file = open((self.name + ".txt"), 'w+')
        file.write("Omitting txt, body is " + self.name + " Orbital a,b " + str(self.ab) + " with period " + str(self.period) + " days.")
        file.close()




# semimajor, semiminor, x0, y0, z0, gamma, beta, alpha, noise_fraction, Example ellipse.
# ellipsey = ellipse(333,111,5463,533,555.1,5.4, 2.8,13,0.001) # 4 4 4 should be x0 y0 z0
#plotter = ellipsey.ellipse_plot()
# xyz = ellipsey.threed_gen()
# uwu = ellipse_solver(xyz, [4,4], "daddy.txt")
# uwu.main() # JUST FOR BUG TESTING

# Import file.
filer = file_import_tool()
file_use = filer.orbit_porter("Earth.dat") # returns xyz, [timestep, stepsplit]
orbit_object = ellipse_solver(file_use[0], file_use[1], file_use[2])
orbit_object.main()