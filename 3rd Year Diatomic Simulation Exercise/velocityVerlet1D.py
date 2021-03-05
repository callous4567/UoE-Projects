"""
CMod Ex2: velocity Verlet time integration of
a particle moving in a double well potential.

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

The potential is V(x) = a*x^4 - b*x^2, where
a and b are hard-coded in the main() method
and passed to the functions that
calculate force and potential energy.
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Particle1D import Particle1D

def force_dw(particle, a, b):
    """
    Method to return the force on a particle
    in a double well potential.
    Force is given by
    F(x) = -dV/dx = -4*a*x^3 + 2*b*x

    :param particle: Particle1D instance
    :param a: parameter a from potential
    :param b: parameter b from potential
    :return: force acting on particle as Numpy array
    """
    force = -4*a*particle.position**3 + 2*b*particle.position
    return force


def pot_energy_dw(particle, a, b):
    """
    Method to return potential energy 
    of particle in double-well potential
    V(x) = a*x^4 - b*x^2

    :param particle: Particle1D instance
    :param a: parameter a from potential
    :param b: parameter b from potential
    :return: potential energy of particle as float
    """
    potential = a*particle.position**4 - b*particle.position**2
    return potential


# Begin main code
def main():
    # Read name of output file from command line
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[1]

    # Open output file
    outfile = open(outfile_name, "w")

    # Set up simulation parameters
    dt = 0.01
    numstep = 2000
    time = 0.0
    a = 0.1
    b = 1.0

    # Set up particle initial conditions:
    #  position x0 = 0.0
    #  velocity v0 = 1.0
    #  mass      m = 1.0
    p1 = Particle1D(0.0, 1.0, 1.0)

    # Write out initial conditions
    energy = p1.kinetic_energy() + pot_energy_dw(p1, a, b)
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,p1.position,energy))

    # Get initial force
    force = force_dw(p1, a, b)

    # Initialise data lists for plotting later
    time_list = [time]
    pos_list = [p1.position]
    energy_list = [energy]

    # Start the time integration loop
    for i in range(numstep):
        # Update particle position
        p1.leap_pos2nd(dt, force)
        
        # Update force
        force_new = force_dw(p1, a, b)
        # Update particle velocity by averaging
        # current and new forces
        p1.leap_velocity(dt, 0.5*(force+force_new))
        
        # Re-define force value
        force = force_new

        # Increase time
        time += dt
        
        # Output particle information
        energy = p1.kinetic_energy() + pot_energy_dw(p1, a, b)
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,p1.position,energy))

        # Append information to data lists
        time_list.append(time)
        pos_list.append(p1.position)
        energy_list.append(energy)
    

    # Post-simulation:
    # Close output file
    outfile.close()

    # Plot particle trajectory to screen
    pyplot.title('Velocity Verlet: position vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Position')
    pyplot.plot(time_list, pos_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Velocity Verlet: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()

