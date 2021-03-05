"""
 CMod Ex2: Particle3D, a class to describe 3D particles
"""
import math

import numpy as np

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    Label(string) - The particle's name
    position(numpy array) - Position in three dimensions (x, y, z)
    velocity(numpy array) - Velocity in three dimensions (v(x), v(y), v(z))
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first order velocity update
    * first and second order position updates
    * define a new particle from a file
    * vector separation
    """

    def __init__(self, par):
        """
        Initialise a Particle3D instance
        
        :param label: label as string
        :param position: position as numpy array
        :param velocity: velocity as numpy array
        :param mass: mass as float
        """
        self.label = str(par[0])
        self.position = np.array([float(par[1]), float(par[2]), float(par[3])])
        self.velocity = np.array([float(par[4]), float(par[5]), float(par[6])])
        self.mass = float(par[7])
    

    def stringy(self):
        """
        Define output format.
        Prints the particle name, followed by position and mass
        "Particle: x = 1.0, y = 1.0, z = 1.0, v(x) =1.0 v(y) = 1.0 v(z) = 1.0, m = 1.0"
        """
        return str(self.position[0]) + ", " + str(self.position[1]) + ", " + str(self.position[2]) + ", " + str(self.velocity[0]) + ", " + str(self.velocity[1]) + ", " + str(self.velocity[2])

    
    def kinetic_energy(self):
        """
        velocity total then
        Return kinetic energy as
        1/2*mass*vel^2
        
        using vel^2 as the velocity vector squared
        
        """
        vsqr = self.velocity[0]**2 +self.velocity[1]**2 +self.velocity[2]**2
        
        return 0.5*self.mass*vsqr
        

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as numpy array
        """
        self.velocity += dt*force/self.mass


    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position += dt*self.velocity


    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as numpy array
        """
        self.position += dt*self.velocity + 0.5*(dt**2)*(force/self.mass)
        
    @staticmethod
    def new_particle (file_handle):
        """
        New particle definition,
        Takes a file parameter and splits it into a list of numpy arrays
        Each element of the list is a different particle
        does not return particle if it is not defined to specifications

        :param file_handle: a file with 8 numbers separated by commas in each row
        """
        Part3D = []
        
        for line in file_handle.readlines():
            partline = line.split(", ")
            if len(partline) != 8:
                break
            particle = np.array([partline[0], partline[1], partline[2], partline[3], partline[4], partline[5], partline[6], partline[7] ])
            Part3D.append(particle)
            
        return Part3D
    @staticmethod
    def vect_sep (p1, p2):
        """
        The Vector separation of two particles,
        Finds the separation in three dimentions of two particles
        Position of particle 2 (in 3D) - position of particle 1 (in 3D)

        :param p1: particle number 1
        :param p2: particle number 2
        """

        vectsep = p2.position-p1.position
        
        
        return vectsep