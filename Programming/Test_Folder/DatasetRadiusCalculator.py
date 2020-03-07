# This script will calculate the radius for multiple spheres of a volume V.
# Should account for exponent of base unit, and allow custom base unit output.
import numpy
import sys
import time
import random
import math

def time_delay(text):
    print()
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.005)
        time.sleep(numbertime)
    
def volume_significand_getter():
    volume_list = [] 
    n_o_text = ("How many spheres are we calculating for? : ")
    time_delay(n_o_text)
    n_o = int(input())
    for i in range( 0, n_o):
        significand_i = float(input(("Sphere {0:d} volume significand: ").\
        format(i + int(1))))
        volume_list.append(significand_i)
    return volume_list

def v_s_l_rounder(v_s_list):
    v_s_l_rounded = numpy.round_(v_s_list, 2)
    return v_s_l_rounded
    
# Numpy signifies a list, rounds it all, and returns the list. It doesn't auto-
# correct the list you provide.
# To install Numpy, sudo install apt-get python3.6-numpy, or something like it!

# So we have a function to get the volume significands, next up is the exponent
# and exponent conversion. There is also the rounding for print-purposes c;

def v_and_e(v_s_list):
    time_delay("Please give your base unit for the unit of volume, i.e. for mm^3, mm = -3... m^3? then 0! : ")
    b_u_e = float(input())
    v_s_e_list = [d*float(10)**(float(3)*b_u_e) for d in v_s_list] # Gives the v    olume in m^3
    # We have the list for the volume and exponential in m...
    radius_list = [(float(3)*e/(float(4)*math.pi))**(float(1/3)) for e in v_s_e_list]
    return v_s_e_list

# Next step is volume to radius, then defining that function and using it with
# the list.

def v_to_r(v):
    # v = 4/3 * pi * r^3 and hence r = 3v/4pi to the 1/3
    r = ((float(3)*v)/(float(4)*math.pi))**(float(1/3))
    return r

# So... we have the function to generate a radius from V... next up will be to
# apply this to the entire list.
"""
Let's start up our program here!"""

volume_significands = volume_significand_getter()
volume_exponentiated = v_and_e(volume_significands)
radius_lister = lambda e: v_to_r(e)
# Converting volume to radius with a function isn't working... OH! We forgot to
# do the return r... so now it works <3
new_radius_list = [radius_lister(h) for h in volume_exponentiated]
# Note, a lambda function:
# function = lambda e: function(e), or anything else that transforms e, will
# give you a new function <3
# Google "Lambda functions" Python

e_c_f_text = "What exponent would you like your radius to? -3 = mm, 0 = m, etc...: "
time_delay(e_c_f_text)
expo = int(input())

def e_c_for_r(r):
    converter = lambda e: e*float(10)**expo
    e_r_list = [converter(h) for h in r]
    return e_r_list

final_list = e_c_for_r(new_radius_list)

for x in final_list:
    haeitch = final_list.index(x) + int(1)
    devi = ("{0} X 10^{1:d} metres").format(x, expo)
    print(("Your final radius for sphere {0:d} is {1}").format(haeitch, devi))

