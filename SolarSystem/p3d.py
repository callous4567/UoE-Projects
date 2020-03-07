import numpy as np
import random

# 3D Particle Object.
# Contains euler, symplectic euler, and verlet integrators.
# Contains single-line file reading functionality.
# pos, vel, mass, index, timestep must be defined as arguments, the first two as ARRAYS and NOT lists.
# Separation calculator via subtraction magic method.
class particle(object):
    # Define various properties. Pos, vel, accel_past (for integrators), mass, index and timestep.
    def __init__(self, pos, vel, mass, index, timestep):
        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.index = index
        self.step = timestep

    # Subscriptability. 0,1,2,3 is xyz, vxvyvz, mass and name.
    def __getitem__(self, item):
        if item == int(0):
            return self.pos
        if item == int(1):
            return self.vel
        if item == int(2):
            return self.mass
        if item == int(3):
            return self.index

    # String magic method
    def __str__(self):
        return "Particle: " + str(self.index) + " Position: " + str(self.pos) + " Velocity: " + str(self.vel) + " Mass: " + str(self.mass)

    # Generates a difference, returning relative position vector and magnitude.
    # Example: parti_a - parti_b = v_a - v_b. Returns [vector, distance].
    # The return is a numpy array.
    def __sub__(self, other):
        vector_difference = self.pos - other.pos
        magnitude = np.linalg.norm(vector_difference)
        return np.array([vector_difference, magnitude])

    # Functionality for momentum summation. Useful for CoM correction in making the code cleaner to deal with.
    def __add__(self, other):
        momentum_sum = self.mass*self.vel + other.mass*other.vel
        return momentum_sum


    # Returns kinetic energy 1/2*m*v**2
    def kinetic(self):
        ke = (1/2)*self.mass*(np.linalg.norm(self.vel)**2)
        return ke


    # EULER INTEGRATORS. For standard euler, run both simultaneously. For symplectic:
    # Run the r_step process. x(t+dt) = x(t) + v(t)*dt
    # This runs on current position and velocity.
    # Run velocity update after updating all positions forward to x(t+dt): calculate forces for the new position.
    # v(t+dt) = v(t) + a(t+dt)*dt (uses the new acceleration, not the older one from the past step)
    def euler_r_step(self):
        self.pos += self.vel*self.step
    def euler_v_step(self, f):
        self.vel += self.step*f/self.mass
    # VERLET INTEGRATOR
    # First do second order approximation for position
    # Velocity approximation uses the time average of the acceleration one step ahead and the current acceleration
    # 0.5*(a(t) + a(t+dt))*timestep + v(t) = v(t+dt)
    def verlet_r_step(self, f):
        self.pos += self.vel*self.step + (1/2)*(f/self.mass)*(self.step**2)
        self.accel_past = (f/self.mass) # Saves acceleration from (t) to use in the next step where new force is calculated.
    def verlet_v_step(self, f, f_n):
        self.vel += (1/2)*((f/self.mass) + (f_n/self.mass))*self.step




# Contains the reader_former method.
# Reads file containing particle data.
# Produces ndarray of particles, unidimensional.
# File format:
# X Y Z VX VY VZ MASS INDEX
# Takes arguments:
# - NAMETENSION... i.e. "hey.txt". Must be in string form.
# - POINTSTART... i.e. the pythonic line that it starts on. If your data begins on the 2nd line, pointstart equals 1.
# - TIMESTEP... the iterator timestep you desire, which is defined for all the particels involved.
class filereader(object):
    def reader_former(self, nametension, pointstart, timestep):
        file = open(nametension, "r")
        liner = file.readlines()
        parray = np.ndarray(shape=(1,0),dtype=particle)
        for d in range(pointstart, len(liner)):
            linered = liner[d].split(" ")
            pos = np.array([float(linered[0]), float(linered[1]), float(linered[2])])
            vel = np.array([float(linered[3]), float(linered[4]), float(linered[5])])
            mass = float(linered[6])
            index = str(linered[7])
            particle_d = particle(pos, vel, mass, index, timestep)
            parray = np.append(parray, particle_d)
        file.close()
        return parray

# Class which contains a random particle generator and also an array generator of said random particles.
# Bit of a preemptive for when I inevitably have to test the scripts/etc to see if any nasty stuff crops up.
class randparticlearray(object):
    def __init__(self, number):
        self.number = number
        self.dumbstep = 1
    def generator(self, d):
        randmax = 50
        pos = np.array([random.randint(0, randmax), random.randint(0,randmax), random.randint(0, randmax)])
        vel = np.array([random.randint(0, randmax), random.randint(0,randmax), random.randint(0, randmax)])
        mass = random.randint(0, randmax)
        index = str("{0:d}").format(d)
        return particle(pos, vel, mass, index, self.dumbstep)
    def array_gen(self):
        parray = np.ndarray(shape=(1,0), dtype=particle)
        for d in range(0, self.number):
            particle_generated = self.generator(d)
            parray = np.append(parray, particle_generated)
        return parray


