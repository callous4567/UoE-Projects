# Volume to Radius calculator.

import time
import math
import cmath
import sys
import random

# This function is just for printing text with delayed release for characters
def timeDelay(t):
    print()
    for c in t:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.005)
        time.sleep(numbertime)
    print()

# Volume = 4/3 * pi * r^3 and hence r = 3v/4pi to the 1/3. 1 mm^3 = 1E-9 m^3.
convert_to_r = lambda e: ((float(3)*(e*float(10)**-9))/(float(4)*math.pi))**float(1/3)
# This takes the volume in mm^3 and takes it to a radius in m
area_calc = lambda e: (float(4)*math.pi*(e**2))
# This gets the area of the sphere.
# --------------------------------------------------------------------------
timeDelay("Please give the volume of the sphere in mm^3...: ")
try:
    volume_mm = float(input("Volume: "))
    radius = convert_to_r(volume_mm)
    area = area_calc(radius)
    timeDelay(("Your radius is {0:.2e} metres and your area is {1:.2e} metres^2").format(radius, area))
except:
    timeDelay("You gave a non-float value that can't be used for calculation... sorry!")
