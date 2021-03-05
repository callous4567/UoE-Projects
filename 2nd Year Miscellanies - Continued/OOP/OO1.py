import matplotlib.pyplot as plt
import numpy as np
import time
import random
import sys
import matplotlib.ticker as pltick

class particle:
    def __init__(self, mass, velocity):  # __init__ defines the variables with respect to self
        self.mass = mass # define each variable as self.mass = mass, w/ the mass given by the variable
        self.velocity = velocity

    def get_mass(self):
        print(self.mass)

    def get_velocity(self):
        print(self.velocity)

# When calling the particle:
"""
hello = particle(mass, velocity)
i.e.
hello = particle(2, [1,2,3])
this calls the init with the mass and velocity called, with the new values of the particle : self.mass = mass, and self.velocity = velocity... so...

self.velocity = velocity;

when you call it as hello = particle(mass, velocity), you specify the variables
then

self.velocity = velocity; calls velocity as self.velocity = velocity, with self as the name of the particle, i.e. hello. 
"""
# So yeah.

a = particle(2,[1,2,3])
b = particle(4,[2,3,-1])

