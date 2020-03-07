# Hello Folks! Here's my radius calculator!
# The goal is to get the volume in mm^-3 and output in metres & millimetres!
# First let's import math.

import math
import cmath
import random
import sys
import time

def timeDelay(t):
    print()
    for c in t:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.05)
        time.sleep(numbertime)
    print()

def quantityGetter(quantityText, quantity):
    timeDelay(quantityText)
    print()
    quantity = float(input(quantity + " Value : "))
    return quantity

def vconverterToM(v, e):
    # Converting from exponent, i.e. cm^2, is (1E-2)^3 = 10^-6 metres cubed.
    newe = float(3)*e
    newv = v*float(10)**newe
    return newv
    # A quick test for this... take v = 1 and e = -3, i.e. 1 mm^3...
    # We go forward and take e = -3 to newe = -9,
    # 1x10^-9 metres cubed.
    # This gives us the volume expressed in metres cubed.

# The next step is to take V in metres cubed to radius in metres.
# v = 4/3 * pi * r^3, hence 3v/4pi to the 1/3 is equal to r.

def volumeToRadius(vol):
    # Let's start up with solving for the radius^3. I could do it all in one entire sentence swing but then I would have to convert back...
    radiusCubed = float(3)*vol/(float(4)*math.pi)
    radius = radiusCubed**float(1/3)
    return radius

# So, we have the radius in metres! So... next step to do... We need to question and see what exponent the user wants it in, i.e. to metres, millimetres, etc...

def rconvertertoe(radiusM):
    print()
    whatexponent = "What exponent in metres would you like your quantity to be in? mm = -3, m = 0, fm = -15, you get the idea :)"
    timeDelay(whatexponent)
    print()
    exponentR = float(input("Printing Exponent : "))
    radiusExponentiated = radiusM*((float(10)**exponentR)**-1)
    print()
    printtext = "Your final radius is... *drumroll*... {0:4.2e}! The units are given by the exponent you chose, i.e. mm if you chose -3.".format(radiusExponentiated)
    timeDelay(printtext)

# So we have all the necessary functions to do this. First up then... we need the volume to some exponent, i.e. millimetres, -3. Let's do it!
# -----------------------------------------------------------------------------

volume = "Volume Significand"
volumeText = "Please give the volume significand, i.e. 2mm^3 → 2 <3"
volumeSignificand = quantityGetter(volumeText, volume)

Exponent = "Exponent"
exponentText = "Also, please give the superscript exponent for the unit of length that defines the volume, i.e. mm^3 → -3, as a millimetre is 1x10^-3 metres"
exponentValue = quantityGetter(exponentText, Exponent)
Volume = vconverterToM(volumeSignificand, exponentValue)
RadiusM = volumeToRadius(Volume)
rconvertertoe(RadiusM)


