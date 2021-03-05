#Let's make a program for Gravity and scaling with a planet,
#Radius R initially, mass M and gravity = 6.67E-11
#First, define parameters

import math
import cmath

print("This program will let you define the gravitational field strength\
of a body at a radius from centre, smaller than HS equilibrium radius")

def gField(mass, radius, constantG):
    fieldstrength = mass*constantG*radius**-2
    return fieldstrength

mass = float(input("Mass:"))
radius = float(input("Radius:"))
constantG = float(6.67)*float(10)**-11

gField(mass, radius, constantG)

#Okay. That worked great.
#Now then, at a radius below. We have gField()...

print("Next step will be to calculate it at a radius below.")

def radiusCorrection(gField, depth):
    gField = gField(mass, radius, constantG)
    depthFactor = float(radius - depth)/radius
    newg = depthFactor*gField
    print(newg)
    return(newg)

depth = float(input("Depth below surface:"))

radiusCorrection(gField, depth)

# Done! <3 
          
