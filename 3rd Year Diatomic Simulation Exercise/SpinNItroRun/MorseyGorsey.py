from __init__ import Particle3D
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from textwrap import wrap
# from scipy import fftpack
# from pydsm import ft
# from mpmath import nsum, exp, inf

# Quick documentation notes.
# There are lot of random pieces of code in here that may seem irrelevant. That's because they are.
# I originally started (and did a lot of this) in line for working on an array of loads of particles, which is why the particles are even in an array.
# There are also loads of amendments and changes I did along the way.
# Cleaning them up is a headache, so there are loads of quick fixes instead, i.e. the "base" quantities and whatnot, alongside the inconsistency in graphs.
# Anyway, it seems to work and gives the correct-ish frequencies, so that's good enough for me.
# One last note. SciPy curve fitting was used for getting frequency from the curves. IT ISN'T THAT GREAT IN THE LONG RUN. Only use it for single oscillation.
# The consequence of this is that you can't fit long curves, but instead have to rely on a single oscillation with a very low timestep. It works "well enough."
# Originally I attempted to do a manual DTFT, to no avail. The summation broke.
# Then I tried to use pydsm for a DFT: no dice.
# SciPy FFT's didn't really seem to help either.
# Least squares fitting was nice but was a pain in the *** and you'd need to give it initial parameters to fit.
# In the end, I just used the curve_fit.



# Timemax is the maximum time to do this for. The units of time are given by 1.018050571x10^-14
# i.e. roughly 10.18 femtoseconds, but don't hold me to that.
timemax = 2.3   # Set maximum time for integration. The one here is the one I fitted.
timestep = 0.00001   # Set timestep of any integrators. 0.00001 seems to be a *very* good spot. 0.000001 isn't doable. I used 0.0001 as recommended.
basestep = timestep   # Placeholder to reset the timestep.
whichatompointstart = 35   # Defined by the conditions data file. See "morseconditions.dat" for example, or particle3d.
which_quantity = 0 # 0->3. Separation, Potential, Total KE, Total E.
base_quantity = which_quantity # Placeholder to reset the quantity to default.
filename = "morseconditions.dat"   # Give the name of the file that has the data of the particles. Form is given below.

# This file reader assumes the same form as previously,
# X Y Z VX VY VZ MASS INDEX
# The constants should be given in the order D_E, R_E, ALPHA directly beneath the particle information.
# Returns an array of the two particles and also a trio of constants. (Start Conditions)
def dp_file_reader(namextension):
    pointstart = whichatompointstart
    file = open(namextension, "r")
    liner = file.readlines()
    parray = np.ndarray(shape=(1,0),dtype=Particle3D.particle)
    for d in range(pointstart, pointstart + int(2)):
        linered = liner[d].split(" ")
        pos = np.array([float(linered[0]), float(linered[1]), float(linered[2])])
        vel = np.array([float(linered[3]), float(linered[4]), float(linered[5])])
        mass = float(linered[6])
        index = str(linered[7])
        particle_d = Particle3D.particle(pos, vel, mass, index, timestep)
        parray = np.append(parray, particle_d)
    constant = liner[pointstart + int(2)].split(" ")
    constants = [float(d) for d in constant]
    file.close()
    return parray, constants

# Some quick lore. I spent about an hour wondering "Why is the force zero?" before I realised it was probably the equilibrium or something similar.
# Takes two arguments: The particle array, containing particles [0, 1] and also the constants array, [D_E, R_E, alpha].
# Returns two lists.
# List #1 provides the forces on [0,1] respectively.
# List #2 provides the [separation magnitude, potential, ke, total energy] respectively.
def morse_calc(p_0, p_1, constants):
    separation = p_0 - p_1 # Array. Zeroth element is vector separation. 1st is magnitude. Gives the vector from 1 to 0.
    current_potential = constants[0]*((1 - np.exp(-constants[2]*(separation[1] - constants[1])))**2 - 1) # These are self explanatory
    current_ke_tot = p_0.kinetic() + p_1.kinetic()
    current_energy_tot = current_potential + current_ke_tot

    # Quick note. The vector deduced above points from 1 to 0. The form of the force that is to be used is:
    # F ON 0 = 2*alpha*D_e*(1 - exp(-a(r_10 - r_e)))*exp(-a(r_10 - r_e)) * -r10 (unit) , where r10 is the vector FROM 1 pointing TOWARD 0.
    # Goes against the convention in the notes, but that's how I was taught so that's why I altered the vector direction ;-;
    # FORCE ON ZERO (particle vector above points toward)
    current_force_0 = 2*constants[2]*constants[0]*(1 - np.exp(-constants[2]*(separation[1] - constants[1])))*np.exp(-constants[2]*(separation[1] - constants[1])) * (-separation[0]/separation[1])
    # FORCE ON ONE
    current_force_1 = -current_force_0
    return [current_force_0, current_force_1], [separation[1], current_potential, current_ke_tot, current_energy_tot]

# Symplectic iterator. Returns various graphing stuff (separation, potential, KE, total E, and time in random units from notes)
# Also writes a file with all the data.
# Returns array of SEPARATION, POTENTIAL, KINETIC ENERGY, TOTAL ENERGY, TIME, over the entire integration.
def euler_symplectic(file_read_concat):
    file = open("symplectic_data.dat", "w") # Open a file to write.
    array = file_read_concat[0] # Start array, straight from the file.
    constants = file_read_concat[1] # Constants, just incase.
    file.write("Symplectic Euler Script. Separation, Potential, KE, Total Energy and also the number of steps.\n" + time.ctime())
    file.write("The step interval is the default, 0.0001 units of time, the units given in the code.")
    step = 0 # Step counter.
    separation = [] # Various things for graphing. I could make an animated one but... well... I left this last second, quite literally, so I won't ;-;
    stepper = []
    total_e = []
    pot = []
    ke = []

    while True:
        # Splits array.
        p_0, p_1 = array[0], array[1]

        # Gathers energies/potentials/etc (before any new step, i.e. at t)
        morse_result = morse_calc(p_0, p_1, constants)
        value = morse_result[1]
        value_string = str(("{} + {} + {} + {} + {} \n").format(value[0], value[1], value[2], value[3], step)) # String consists of current separation, potential, KE, and energy total
        file.write(value_string)
        separation.append(value[0]), stepper.append(step*p_0.step), total_e.append(value[3]), pot.append(value[1]), ke.append(value[2])

        # Position step
        p_0.euler_r_step(), p_1.euler_r_step()

        # Gathers forces (after position step, before velocity step)
        morse_resultnew = morse_calc(p_0, p_1, constants)
        forcenew = morse_resultnew[0]

        # Velocity step and step index. Also remake array.
        p_0.euler_v_step(forcenew[0]), p_1.euler_v_step(forcenew[1])
        array[0], array[1] = p_0, p_1
        step += int(1)

        # Lets you set the maximum time (in the units that the sim runs in, which are in the handbook)
        if step*p_0.step >= timemax:
            break

    return [separation, pot, ke, total_e, stepper]

# Euler iterator. Made because Symplectic and Verlet are too similar to debug. Didn't help, In.. the... SLIGHTEST.
# Same returns as above.
def euler(file_read_concat):
    file = open("euler_data.dat", "w") # Open a file to write.
    array = file_read_concat[0] # Start array, straight from the file.
    constants = file_read_concat[1] # Constants, just incase.
    file.write("Symplectic Euler Script. Separation, Potential, KE, Total Energy and also the number of steps.\n" + time.ctime())
    file.write("The step interval is the default, 0.0001 units of time, the units given in the code.")
    step = 0 # Step counter.
    separation = [] # Various things for graphing. I could make an animated one but... well... I left this last second, quite literally, so I won't ;-;
    stepper = []
    total_e = []
    pot = []
    ke = []

    while True:
        # Splits array.
        p_0, p_1 = array[0], array[1]

        # Gathers energies/potentials/etc (before any new step, i.e. at t)
        morse_result = morse_calc(p_0, p_1, constants)
        value = morse_result[1]
        value_string = str(("{} + {} + {} + {} + {} \n").format(value[0], value[1], value[2], value[3], step)) # String consists of current separation, potential, KE, and energy total
        file.write(value_string)
        separation.append(value[0]), stepper.append(step*p_0.step), total_e.append(value[3]), pot.append(value[1]), ke.append(value[2])
        force = morse_result[0]

        # Position step
        p_0.euler_r_step(), p_1.euler_r_step()

        # velocity step
        p_0.euler_v_step(force[0]), p_1.euler_v_step(force[1])
        array[0], array[1] = p_0, p_1
        step += int(1)

        # Les you set the maximum time (in the units that the sim runs in, which are in the handbook)
        if step*p_0.step >= timemax:
            break

    return [separation, pot, ke, total_e, stepper]

# Verlet iterator.
# Same returns as above.
def verlererator(file_read_concat):
    file = open("verlet.dat", "w") # Open a file to write.
    array = file_read_concat[0] # Start array, straight from the file.
    constants = file_read_concat[1] # Constants, just incase.
    file.write("Velocity Verlet Script. Separation, Potential, KE, Total Energy and also the number of steps.\n" + time.ctime())
    file.write("The step interval is the default, 0.0001 units of time, the units given in the code.")
    step = 0 # Step counter.
    separation = [] # Various things for graphing. I could make an animated one but... well... I left this last second, quite literally, so I won't ;-;
    stepper = []
    total_e = []
    pot = []
    ke = []

    while True:
        # Split array
        p_0, p_1 = array[0], array[1]

        # Morse Force/Values
        morse_result = morse_calc(p_0, p_1, constants)
        value = morse_result[1]
        value_string = str(("{} + {} + {} + {} + {} \n").format(value[0], value[1], value[2], value[3], step)) # String consists of current separation, potential, KE, and energy total
        file.write(value_string)
        separation.append(value[0]), stepper.append(np.float64(step*p_0.step)), total_e.append(value[3]), pot.append(value[1]), ke.append(value[2])

        # Run position step
        force = morse_result[0]
        p_0.verlet_r_step(force[0]), p_1.verlet_r_step(force[1])

        # Morse Force/Values ROUND TWO FIGHT
        morse_resultnew = morse_calc(p_0, p_1, constants)

        # Run velocity step
        forcenew = morse_resultnew[0]
        p_0.verlet_v_step(force[0], forcenew[0]), p_1.verlet_v_step(force[1], forcenew[1])


        step += int(1)

        # Les you set the maximum time (in the units that the sim runs in, which are in the handbook)
        if step*p_0.step >= timemax:
            break
    return [separation, pot, ke, total_e, stepper]

# Data format is the same as output from the two iterators. Just plug the iterator straight in! <3
"""  **** doesn't even work 
# Discrete time fourier transform
def dtft(xdata, k):
    datan = lambda n: xdata[n]
    summation = nsum(lambda n: datan(int(n))*exp(-np.complex(0,1)*n*k*timestep), [-inf,inf])
    return summation
"""

# CURVE FITTING DETAILS
# URGENT NOTE <--------------------------------------------------------------------------------
# Only run for ~1 oscillation tops, with an insanely low timestep, otherwise the curve fit totally ****'s up. Seems to be an issue with SciPy of some sort, according to the SE Forums.
# Anyways, yeah. Trial and error with the "max_time" argument at the top of the page until you get a reasonable fit.
# Functional form for the sinusoid we expect from the oscillating atoms.
def funct(t,a,b,c,d):
    return (a*np.sin(b*t - c) + d)
# Parameter getter for funct. Returns array [a,b,c,d]. Uses SciPy numerical curve fit.
def param(xdata, ydata):
    param,misc = curve_fit(funct, xdata, ydata)
    return param
# Just a misc item for graphs. Identifies the quantity for the titles/etc.
def identifier(what_quantity):
    if what_quantity == 0:
        return "Radial Separation"
    if what_quantity == 1:
        return "Potential Energy"
    if what_quantity == 2:
        return "Total Kinetic Energy"
    if what_quantity == 3:
        return "Total Energy"
# Graph labels. Again, miscellaneous.
def identilabel(what_quantity):
    if what_quantity == 0:
        return ("Radial Separation in " + r'$\AA$')
    if what_quantity == 1:
        return "Potential Energy in eV"
    if what_quantity == 2:
        return "Total Kinetic Energy in eV"
    if what_quantity == 3:
        return "Total Energy in eV"

# Graphs the integrator output. Takes entire integrator as argument. Index is what you want the file to be called. Overlays the fitted curve.
# As is before, only ever run this for several oscillations, max of 3-5. The curve fit cocks up past that.
# Late addition, "what_quantity" varies from 0 to 3, same indices as the integrator output. Ex: 0 would be separation. 3 total energy.
# I recommend only using this for short integrations, ~1 oscillation tops. Otherwise use the 69 edition.
def grapherator_9000(dataset, index, what_quantity):
    figure = plt.figure(111, figsize=(8,8))

    x, y = dataset[4], dataset[what_quantity]
    ax1 = plt.subplot()  # Subplot on figure
    ax1.set(title=("\n".join(wrap("Plot of {} vs. time of simulation for atoms in diatomic covalent configuration, {}."))).format(identifier(what_quantity), index))
    ax1.plot(x, y, color='black', linewidth=0.5, label=("{} Data").format(index))  # Plot graph
    ax1.grid(True, which='major', color="blue", alpha=0.1)  # Enable grids on subplot
    ax1.grid(True, which='minor', color="pink", alpha=0.4)
    ax1.set(xlabel="Time elapsed, " + r'$\approx 10.18$' + " fs per unit time", ylabel=identilabel(what_quantity))


    # Curve optimization instead...
    # asin(bx + c) + d
    # b is equal to k, which is omega, which is 2*pi*f. The oscillatory frequency can be had by dividing omega by 2pi.
    # There's a try-except just incase of a nasty error with the curve fitter.
    xdata, ydata = x,y
    try:
        parameters = param(xdata, ydata)
        freq = np.abs(parameters[1]/(2*np.pi)) # Take abs. It's negative.
        # Spectroscopic frequency calculation. The units of time (calculated from CODATA 2018 are 1.018050571x10^-14
        true_freq = freq * 1/(10**-10 * ((1.66053906660*(10**(-27)))/(1.602176634*(10**(-19))))**(1/2))
        lightwavelengthincentimetres = 100*(2.99792458*(10**8))/true_freq
        inversewavelength = lightwavelengthincentimetres**-1 # IN CM^-1


        recalc_y = [funct(d, *parameters) for d in xdata]
        ax3 = plt.subplot()
        ax3.plot(xdata, recalc_y, 'r-', label="Sinusoid fit, asin(bt + c) + d")
        ax3.annotate(s=(("Frequency is {0:.5e}, in units of cm").format(inversewavelength) + r'$^{-1}$'), xy=(0,min(recalc_y)))
        plt.legend(loc="upper right")

        figure.savefig(("{}{}.png").format(index, identilabel(what_quantity))) # Saves figure. Always do it (idk, the old code that I used in SciProg always had this since I needed to save it and it's habit.)

        plt.show()
    except:
        print("Wasn't able to fit curve.")
        plt.show()
        plt.savefig(("{}{}nofit.png").format(index, identilabel(what_quantity)))

# Doesn't provide curve fit. Same output, though.
def grapherator_69(dataset, index, what_quantity):
    # file_read_concat = dp_file_reader("morseconditions.dat", 13) # Or change to 19, your choice. Just check the file for the two start points c: 19 is for Nitrogen, 13 oxygen.

    figure = plt.figure(111, figsize=(8,8))

    x, y = dataset[4], dataset[what_quantity]
    ax1 = plt.subplot()  # Subplot on figure
    ax1.set(title=("\n".join(wrap("Plot of {} vs. time of simulation for atoms in diatomic covalent configuration, {}."))).format(identifier(what_quantity),index))
    ax1.plot(x, y, color='black', linewidth=0.5, label=("{} Data").format(index))  # Plot graph
    ax1.grid(True, which='major', color="blue", alpha=0.1)  # Enable grids on subplot
    ax1.grid(True, which='minor', color="pink", alpha=0.4)
    ax1.set(xlabel="Time elapsed, " + r'$\approx 10.18$' + " fs per unit time", ylabel=identilabel(what_quantity))

    xdata = dataset[4]
    ydata = dataset[0]

    figure.savefig(("{}{}.png").format(index, identifier(what_quantity))) # Saves figure. Always do it (idk, the old code that I used in SciProg always had this since I needed to save it and it's habit.)
    plt.show()

# Frequency getter (from raw data): takes entire integrator as argument, once more. Calculates frequency for the chosen quantity.
# Once more, only really suitable for the radial separation. It gets messy otherwise.
def frequetter(integrator):
    xdata, ydata = integrator[4], integrator[which_quantity]
    popt = param(xdata, ydata)
    freq = np.abs(popt[1] / (2 * np.pi))  # Take abs. It's negative.
    # Spectroscopic frequency calculation. The units of time (calculated from CODATA 2018 are 1.018050571x10^-14
    true_freq = freq * 1 / (10 ** -10 * ((1.66053906660 * (10 ** (-27))) / (1.602176634 * (10 ** (-19)))) ** (1 / 2))
    lightwavelengthincentimetres = 100 * (2.99792458 * (10 ** 8)) / true_freq
    inversewavelength = lightwavelengthincentimetres ** -1  # IN CM^-1
    return inversewavelength

# Relative error in timestep and energy of the system.
# Returns 3 lists.
# Each list is for euler, symplectic euler, or verlet respectvely.
# Form of list: [timestep_for_integration_frequency, relative_frequency_error, timestep_for_integration_energy, relative_energy_error]
# Relative error for frequency is given by the absolute value of the difference from the "best" frequency divided by that best frequency. "Best" is for the timestep chosen at the top of the page.
# Energy error is given by half the range of the values of total energy, normalised by best energy. Best in both cases is chosen for Verlet integrator.
# Urgent point: the error in energy can be had perpetually, up to the 0.5 requested for the timestep. Error in frequency is limited by the point that SciPy can no longer use curve_fit().
# It appears to not favor a fit on long curves, to be expected for a numerical sim I suppose, but yeah.

def convergerator(separation):
    global timestep
    global which_quantity
    if separation == "base":
        separation = basestep

    eulersteplist, eulerfreq, eulerEsteplist, eulerenergy = [],[],[],[]
    sympeulersteplist, sympeulerfreq, sympeulerEsteplist, sympeulerenergy = [],[],[],[]
    verletsteplist, verlerfreq, verletEsteplist, verletenergy = [],[],[],[]

    timestep = basestep
    while (timestep <= 0.5):
        print("euler")
        print(timestep)
        eulerdata = euler(dp_file_reader("morseconditions.dat"))
        try:
            which_quantity = 0
            eulerfreqt = frequetter(eulerdata)
            eulerfreq.append(eulerfreqt)
            eulersteplist.append(timestep)
            which_quantity = base_quantity
        except:
            which_quantity = base_quantity
        energyerror = (max(eulerdata[3]) - min(eulerdata[3]))/eulerdata[3][0]
        eulerenergy.append(energyerror)
        eulerEsteplist.append(timestep)
        timestep += separation
    timestep = basestep
    while (timestep <= 0.5):
        print("sympeuler")
        print(timestep)
        sympeulerdata = euler_symplectic(dp_file_reader("morseconditions.dat"))
        try:
            which_quantity = 0
            sympeulerfreqt = frequetter(sympeulerdata)
            sympeulerfreq.append(sympeulerfreqt)
            sympeulersteplist.append(timestep)
            which_quantity = base_quantity
        except:
            which_quantity = base_quantity
        energyerror = (max(sympeulerdata[3]) - min(sympeulerdata[3]))/sympeulerdata[3][0]
        sympeulerenergy.append(energyerror)
        sympeulerEsteplist.append(timestep)
        timestep += separation
    timestep = basestep
    while (timestep <= 0.5):
        print("verlet")
        print(timestep)
        verleta = verlererator(dp_file_reader("morseconditions.dat"))
        try:
            which_quantity = 0
            verletfreqt = frequetter(verleta)
            verlerfreq.append(verletfreqt)
            verletsteplist.append(timestep)
            which_quantity = base_quantity
        except:
            which_quantity = base_quantity
        energyerror = (max(verleta[3]) - min(verleta[3]))/verleta[3][0]
        verletenergy.append(energyerror)
        verletEsteplist.append(timestep)
        timestep += separation

    eulerenergy, sympeulerenergy, verletenergy = [np.abs(d) for d in eulerenergy], [np.abs(d) for d in sympeulerenergy], [np.abs(d) for d in verletenergy]
    timestep = basestep
    normeuler = [np.abs(d/verlerfreq[0]) for d in eulerfreq]
    normeuler = [np.abs(1-d) for d in normeuler]
    normsympeuler = [np.abs(d/verlerfreq[0]) for d in sympeulerfreq]
    normsympeuler = [np.abs(1-d) for d in normsympeuler]
    normverlet = [np.abs(d/verlerfreq[0]) for d in verlerfreq]
    normverlet = [np.abs(1-d) for d in normverlet]
    print(len(normverlet))
    print(len(verletenergy))
    print(len(verletsteplist))
    print(len(verletEsteplist))

    return [eulersteplist,normeuler, eulerEsteplist, eulerenergy], [sympeulersteplist,normsympeuler, sympeulerEsteplist, sympeulerenergy] , [verletsteplist,normverlet, verletEsteplist, verletenergy]

# Graphs the data straight from convergerator. Which is 1 or 3. 1 is frequency error, 3 is energy error.
# Strange sharps occur in the timesteps... not sure why. In any case, I'll eyeball an interpolation for when I decide the max timestep for error.
def convergerapherator(sim_step_norm_energy, which):
    type = "dummy"
    steplist = ""
    if which == 1:
        type = "Frequency"
        steplist = 0
    else:
        type = "Energy"
        steplist = 2

    figure = plt.figure(111, figsize=(8,8))
    ax1 = plt.subplot()  # Subplot on figure
    ax1.set(title=("Plot of {} Fractional Error vs. Timestep").format(type))
    ax1.set(ylim=(0, 0.1), xlim=(0, 0.5))
    ax1.plot(sim_step_norm_energy[0][steplist], sim_step_norm_energy[0][which], color='red', linewidth=1, label="Euler Error")
    ax1.grid(True, which='major', color="blue", alpha=0.1)  # Enable grids on subplot
    ax1.grid(True, which='minor', color="pink", alpha=0.4)
    ax1.set(xlabel=("Timestep, " + r'$\approx$' + " 10.18 fs per unit time"), ylabel=("{} Fractional Error").format(type))
    ax2 = plt.subplot()
    ax2.plot(sim_step_norm_energy[1][steplist], sim_step_norm_energy[1][which], color='blue', linewidth=1, label="Symplectic Euler Error")
    ax2.grid(True, which='major', color="blue", alpha=0.1)  # Enable grids on subplot
    ax2.grid(True, which='minor', color="pink", alpha=0.4)
    ax3 = plt.subplot()
    ax3.plot(sim_step_norm_energy[2][steplist], sim_step_norm_energy[2][which], color='green', linewidth=1, label="Verlet Error")
    ax3.grid(True, which='major', color="blue", alpha=0.1)  # Enable grids on subplot
    ax3.grid(True, which='minor', color="pink", alpha=0.4)

    # Plot to illuminate 1%/0.01 error. Doesn't extend the entire graph, but it's useful enough to indicate which points are the limits.
    ax4 = plt.subplot()
    x_range = np.linspace(0, 0.5, 5000)
    y = [0.01 for d in x_range]
    ax4.plot(x_range, y, color='black', linewidth=0.8)
    plt.legend()
    plt.show()
    figure.savefig(("Error_plot_{}.png").format(type))



# If you purely want the data for a single run, just run the iterator of choice.
# Other than that that, this gives various graphs/data/etc.
if which_quantity == 0:
    grapherator_9000(verlererator(dp_file_reader(filename)), "Verlet", which_quantity)
    grapherator_9000(euler_symplectic(dp_file_reader(filename)), "Symplectic Euler", which_quantity)
    grapherator_9000(euler(dp_file_reader(filename)), "Normal Euler", which_quantity)
    convergedata = convergerator("base")
    convergerapherator(convergedata, 1)
    convergerapherator(convergedata, 3)
else:
    grapherator_69(verlererator(dp_file_reader(filename)), "Verlet", which_quantity)
    grapherator_69(euler_symplectic(dp_file_reader(filename)), "Symplectic Euler", which_quantity)
    grapherator_69(euler(dp_file_reader(filename)), "Normal Euler", which_quantity)










