from sympy import *
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
# Proudly coded by Callicious <3
from mpl_toolkits.mplot3d import Axes3D

"""
This is a note to future me. Some stuff to include on future rendition for report:
- The simulator takes a timestep. Modify the printing software to print timestep in, then read it. CHECK! 
- Modify for partial orbit simulation parameters: do for report. Works in theory. Can't do period though. So half a check. 
"""

# Class with all import tools for the files in this project.
class file_import_tool(object):
    # Gives transposed XYZ data, i.e. [xlist, ylist, zlist] from read in file.
    def orbit_porter(self, filename_and_extension):
        xlist, ylist, zlist = [],[],[]

        with open(filename_and_extension) as f:
            for line in f:
                linesplit = line.split(",")
                linesplit = [float(d) for d in linesplit]
                x,y,z = linesplit[0], linesplit[1], linesplit[2]
                xlist.append(x), ylist.append(y), zlist.append(z)

        #fig = plt.figure()
        #ax1 = fig.add_subplot(111, projection='3d')

        #graph = ax1.plot(xlist, ylist, zlist, 'o', lw=0.1, markersize=0.5)

        #plt.show()

        xyz_list = [xlist, ylist, zlist]

        return xyz_list
    # Returns timestep (days) for each step in file, modified to handle specific parameters file that we use.
    # FORMAT:
    # timestep, number of steps, some_random_constant_I_don't_know
    # This takes the timestep.
    def paramgetter(self, paramname_and_extension):
        owo = open(paramname_and_extension)
        # Timestep, Total Steps, Graviconstant
        uwu = owo.readlines()
        step = uwu[0].split(",")[0]
        owo.close()
        return step
    # Gets list of names + extensions for bodies involved.
    # Rips these from the "particles.txt" file.
    def namesgetter(self, names_and_extension):
        hairimasuyo = open(names_and_extension)
        ryokai = hairimasuyo.readlines()
        names = []
        for line in ryokai:
            linesplit = line.split(",")
            names.append(linesplit[0] + ".dat")
        hairimasuyo.close()
        return names


# Ellipse class.
# Ellipse object takes a, b, x_0, y_0, z_0 gamma, beta, alpha, noise_fraction
# gamma, beta, alpha are the extrinsic euler angles defined on the wikipedia page.
# This is just to generate noisy data with noise_fraction being a ratio of a to use as "noise" on all final x,y,z values.
# I needed an ellipse to practice fitting with so I wrote this to generate it.
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

# ELLIPSE SOLVE INFORMATION, ROUGH GUIDE, AND INTUITION. Recommend reading this. If you don't, you might be at a loss as to why I'm doing some things.
#################################################################################################################################################################
# Note.
# Data is [xlist, ylist, zlist]
# Name is example: "Earth.dat".... the filename and extension
#################################################################################################################################################################
# vecsepmin selects the value of vector separation magnitude relative to the initial position.
# This is used in calculating orbit period by picking off the next point that vector separation magnitude relative to initial position dips beneath this value.
#################################################################################################################################################################
# infinitycllosurelimit sets all elements from 0 up to this value to infinity for closure calculation
# ensures that no closure is done just for values near the origin.
# this should be set large enough such that only the first closure is done, but not such that it passes first closure.
#################################################################################################################################################################
# The general method utilised here for the PERIOD is as follows:
# A list of vector separation magnitudes is calculated for all elements in the provided orbit, relative to initial position.
# This varies periodically in time, and thus we select for the first minimum (where it should tend to zero) as our period point.
# This is done by setting "vecsepmin" as the upper limit for vector separation at first closure and setting the first "infinityclosurelimit" values to infinity.
# Once we pass vecsepmin, we select this point as the point at which we have orbit closure.
# We then remove all elements of the list before this closure point, and then instead of searching for values less than vecsepmin, we search for values above.
# This takes us across the trough in vector separation as a function of time, passing vector separation magnitude equal to zero until we hit vecsepmin again.
# The length across this point is halved and added to the previous closure_num, giving us a better estimate than otherwise would have been had.
#################################################################################################################################################################
# The general method utilised here FOR GEOMETRIC PARAMETERS is as follows:
# First we translate ellipse to origin, using a translation vector from the average X Y Z coordinate over the closed loop.
# This isn't entirely necessary (we can get a,b using conic mathematics without even translating) but it makes it easier to debug and makes graphing better.
# Flatten the ellipse via gram-schmidt orthogonalization of the coordinate system, using the ellipse as a plane. We cut XYZ data down to the "closed orbit".
# This should (I imagine...?) save on processing, though it may reduce how accurate the ellipse is solved for... worth writing about in the report?
# 2D least-squares solve using numpy.linalg.leastsq and matrix solve for equation ax^2 + by^2 + cxy + dx + ey (or whatever symbols I used)
# Ellipse parameters are determined using conic section mathematics as outlined in this convenient Stack Exchange post:
# https://math.stackexchange.com/questions/280937/finding-the-angle-of-rotation-of-an-ellipse-from-its-general-equation-and-the-ot
# A more formal math derivation can be had in the "Mathematical Methods for Physics and Engineering" book, you know the one, the big one, but I dislike it...
#################################################################################################################################################################
# The general method for LUNAR SOLVING is as follows:
# Exact same functional method as above.
# The data for the moon is instead the relative vector from the moon to the earth, i.e. moonvectors - earthvectors.
# Then just solve for the lunar ellipse instead.
#################################################################################################################################################################
# Further note on recommended vecsepmin, infinityclosurelimit, and timestep.
# The timestep used initially was 1 day. This was not satisfactory. 0.1 gives reasonable results (~0.5% accurate on average... which is okay)
# vecsepmin can be pushed down to 30 timesteps, equivalent to an error in 3 days, which is reduced further by the method used in averaging closure num for the
# period, down to less than 1/2 a day for smaller orbits.
# Halley requires special treatment: if you examine the graphing plot for it, you'll notice that vector separation magnitude never actually reaches 0 (or even
# near 0) at all later on, thus the custom loop.
# If you've read this far, thanks for reading this far! You might ask "Why the stupid long method... why not just do a 3D least squares fit on the ellipse?"...
# My answer to that is that while I'm sure it's doable, this seemed more fun and less... painful... even if it's hit-or-miss and less accurate.
################################################################################################################################################################
# Callicious = Sebastian = Me. That's my nickname online! c:
class ellipse_solver(object):
    def __init__(self, data, step, name, vecsepmin, infinityclosurelimit):
        self.xyz = data # x y z transposed data, i.e. not vectorized
        self.period = 0 # orbit period in earth days
        self.xyzrow = 0 # vectorized xyz
        self.step = float(step) # self explanatory... timestep
        self.transestimate = 0  # estimate for translation
        self.closure_num = 0  # estimate of closure number in xyzrow
        self.twodflat = 0 # placeholder
        self.ab = [0,0] # major-minor placeholder
        self.name = name # For future name implementation when graphing/writing files. May not be used.
        self.closed_xyz = 0 # Closed orbit placeholder
        self.vecsepmin = vecsepmin # Sets minimum tolerance for closure
        self.infinityclosurelimit = infinityclosurelimit # Sets elements of list to set to infinity (to avoid closure_num being near start/before full orbit)

    # Calculates period, and returns closed orbit data.
    def vecsep_calc(self):
        xyz = np.column_stack((self.xyz[0], self.xyz[1], self.xyz[2]))

        vecseps = []
        # Vector separations (indexed) vs. the original
        for d in range(len(xyz)):
            vecseps.append(np.linalg.norm(xyz[d] - xyz[0]))



        hey = np.arange(0, len(vecseps), 1)

        #plt.plot(hey, vecseps)
        #plt.show()

        # Set a bound on the minimum vector separation beneath which we have closure.
        # Python searches list horizontally left-right. First point which this is satisfied will provide the estimate.
        # You can set it manually or to a lower index if you're confident in your integrator period +/+ want more accuracy.
        minsearch = vecseps[self.vecsepmin]

        # Set range of values to infinity to prevent overlap when estimating first closure. Make sure to set to less than orbit of smallest body.
        for i in range(0, self.infinityclosurelimit):
            vecseps[i] = np.inf

        closure_num = 0

        # New closure method. Assuming symmetry in the graphs, which seems to be the case.
        # Closure Index
        for num, item in enumerate(vecseps, 0):
            if item <= minsearch:
                closure_num = num
                break

        # Clip vecseps down.
        vecseps = vecseps[int(closure_num):]

        # Find the point where the search exceeds minsearch and average for the middle, the bottom, giving better closure estimate. 10/10 works. Fantastic.
        # Makes sense if you examine graphically... sorry for not being so intuitive :(
        for num, item in enumerate(vecseps, 0):
            if item >= minsearch:
                closure_num += 0.5*num
                closure_num = int(closure_num)
                self.period = closure_num*self.step # Set period.
                break

        # Estimation of the translation for the ellipse.
        xyzsplit = xyz.T
        translation = np.array([np.mean(d) for d in xyzsplit])

        # In case we need it, closed xyz vector list, in row-vector format.
        closed_xyz = np.array(xyz[:closure_num])

        closure_num -= int(1)

        # Set various placeholders for onward.
        self.closed_xyz = closed_xyz
        self.xyzrow = np.array(xyz)
        self.transestimate = translation
        self.closure_num = closure_num

    # Gram-schmidt orbit flattener.
    def gram_flattener(self):

        onemax, twomax, twomin = 0.2,0.4,0.3

        closed_xyz, closure_num, translation = self.closed_xyz, self.closure_num, self.transestimate

        translated_xyz = [d - translation for d in closed_xyz]

        print(self.closure_num)
        print(len(translated_xyz))

        vector1, vector2 = translated_xyz[int(onemax * closure_num)] - translated_xyz[0], translated_xyz[int(twomax * closure_num)] - \
                           translated_xyz[int(twomin * closure_num)]



        gram2 = vector2 / np.linalg.norm(vector2)
        # v12_angle = np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
        mid1 = vector1 - (np.dot(vector1, gram2)*gram2)
        gram1 = mid1 / np.linalg.norm(mid1)
        gram3 = np.cross(gram1, gram2)

        # 1x2 = 3. 1,2 are xy, 3 iz z.

        grammatrix = np.array((gram1, gram2, gram3))
        invertrix = np.linalg.inv(grammatrix)
        flattened_xyz = np.array([np.dot(d, invertrix) for d in translated_xyz]).T # x,y,z form / transposed

        #plt.plot(flattened_xyz[0], flattened_xyz[1], 'o', lw=0.2, markersize=0.3)
        #plt.show()
        # We now have a (semi) flat 2D ellipse, rotated obviously.
        # Now then... to fit the ellipse...

        self.twodflat = flattened_xyz

    # Least-squares fitting and orbit solving.
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


        x, y, z = self.twodflat
        xraw, yraw, zraw = self.xyz


        # Matrixify up a gammatrix
        gamma = np.array(gammatrix(gamma), dtype=float)

        # Rotate back the ellipse!
        xy = np.column_stack((x,y))
        xy_new = np.array([np.dot(d, gamma) for d in xy])
        xx, yy = xy_new.T

        try:
            # GRAPH EVERYTHING FOR THE FINAL PRODUCT
            fig3 = plt.figure()
            ax3 = fig3.add_subplot(111, projection='3d')
            ax3.plot(self.xyz[0], self.xyz[1], self.xyz[2], color='red')
            ax3.plot(xx, yy, [0 for d in xx], color='blue')
            namesplit = self.name.split('.')
            ax3.set(title=('Orbit for ' + namesplit[0]), xlabel="AU", ylabel="AU", zlabel="AU")
            legend_elements = [Patch(facecolor='red', edgecolor='black', label='Original Orbit'),
                               Patch(facecolor='blue', edgecolor='black', label="Processed Orbit")]
            ax3.legend(handles=legend_elements, loc='upper right', title="Legend")
            plt.show()
            fig3.savefig(namesplit[0] + ".png", dpi=600)
        except:
            pass


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
        ab = [np.sqrt(asquared), np.sqrt(bsquared)]

        self.ab = [max(ab), min(ab)]



    # Main. Execute everything for orbit solve.
    def main(self):
        self.vecsep_calc()
        self.gram_flattener()
        self.twod_optimizer()

        # Okay. That gives the first estimate using a large minimum vecsep, hardcoded for the biggest radius we have (Halley has problems with it)

        # Write parameters to files.
        # This is a custom choice and I prefer it, but I realise it might not be preferable (and having one file for the whole system might have been better)
        # Anyway, I like it this way... don't know why.
        # OH. Apoapsis and periapsis... Okidoke...

        eccentricity = np.sqrt(1 - ((self.ab[1]**2)/(self.ab[0]**2)))
        apo, peri = self.ab[0]*(1 + eccentricity), self.ab[0]*(1 - eccentricity)
        file = open((self.name + "-Solved-Parameters-" + ".txt"), 'w+')
        file.write("Omitting txt and dat... body is " + self.name + "\n")
        file.write("Orbital period is " + str(self.period) + " days. \n")
        file.write("Semimajor, semiminor, eccentricty, apoapsis, and periapsis are respectively " + str(self.ab) + " " + str(eccentricity) + " " + str(apo) + " " + str(peri))




        file.close()

# Solves for the moon specifically (the 27.5 day orbit, not 365.whatever)
class moon_solver(object):
    # Note that earth_input and moon_input are unvectorized, i.e. the transposed form of the XYZ data.
    def __init__(self, timestep, earthdata, lunadata):
        self.relvectata = [] # Relative vector placeholder
        self.earthdata = earthdata
        self.lunadata = lunadata
        self.timestep = timestep

    # Calculate relative vectors, moon - earth.
    def relvects(self):
        earthrow, lunarow = np.column_stack((self.earthdata[0], self.earthdata[1], self.earthdata[2])), np.column_stack((self.lunadata[0], self.lunadata[1], self.lunadata[2]))
        relvects = lunarow - earthrow
        self.relvectata = relvects.T

    # Orbit solve as for other orbits.
    def mainethestate(self):
        self.relvects()
        lunasolve = ellipse_solver(self.relvectata, self.timestep, "Luna.txt", 200, 1500)
        lunasolve.main()



# BUG TESTING FOR ELLIPSES. You can ignore this.
# semimajor, semiminor, x0, y0, z0, gamma, beta, alpha, noise_fraction, Example ellipse.
# ellipsey = ellipse(333,111,5463,533,555.1,5.4, 2.8,13,0.001) # 4 4 4 should be x0 y0 z0
#plotter = ellipsey.ellipse_plot()
# xyz = ellipsey.threed_gen()
# uwu = ellipse_solver(xyz, [4,4], "daddy.txt")
# uwu.main() # JUST FOR BUG TESTING. Needed some ellipse data to fit.



# NOTE URGENT NOTE READ THIS NOTE.
# THIS SECTION IS "OPTIMISED" SO TO SAY FOR A TIMESTEP OF 0.1 WITH 1,000,000 DATAPOINTS.
# Data is not included (duh... 50 MB limit? I can't fit like 2 GB of files in that!)

# Set up timestep, names of orbits/etc, and also the import tool.
import_tool = file_import_tool()
timestep = import_tool.paramgetter("parameters.dat")
filedatalist = import_tool.namesgetter("particles.dat")

# Precursors for lunar solve (Earth and Moon data. I could have coded this into the loop to save double-processing, but the time is negligible anyway and I don't think we're graded on the efficiency of the file)
earthdata = import_tool.orbit_porter("Earth.dat")
lunadata = import_tool.orbit_porter("Moon.dat")

# Lunar solve
moonsolver = moon_solver(timestep, earthdata, lunadata)
moonsolver.mainethestate()

# Run solve for all orbits.
for d in range(len(filedatalist)):
    if filedatalist[d] == "Sun.dat":
        pass
    elif filedatalist[d] == "Halley.dat":
        print(filedatalist[d])
        solver_object = ellipse_solver(import_tool.orbit_porter(filedatalist[d]), timestep, filedatalist[d], 60, 400)
        # 1, 109500 are the step and number of steps for 300 years (Halley Closure, supposedly.)
        # 1, 100,000 with 1500, 10000
        # 0.1, 1,000,000 with 15000, 100000
        # 0.01, 10,000,000 with 150000, 1000000

        solver_object.main()
    else:
        print(filedatalist[d])
        solver_object = ellipse_solver(import_tool.orbit_porter(filedatalist[d]), timestep, filedatalist[d], 1, 3)
        # 1, 109500 are the step and number of steps for 300 years (Halley Closure, supposedly.)
        # 1, 100,000 with 30, 25
        # 0.1, 1,000,000 with 30, 250
        # 0.01, 10,000,000 with 300, 2500

        solver_object.main()



