import numpy as np
import random 

# We can define infinite loops.
def loop_ex():
    while True:
        x = random.uniform(-radius, radius)
        y = random.uniform(-radius, radius)
        z = random,uniform(-radius, radius)
        if x**2 + y**2 + z**2 <= radius**2:
            break
    return x, y, z
# This loop gets you x, y, z for random points in a sphere...
# while True defines an infinite loop, i.e. while True, True is always True, hence an infinite loop.
# This loop gets a random x, random y, and random z
# It then checks if x**2 + y**2 + z**2 is less than or equal to radius**2
# If it is, then it breaks the loop
# If the random x y and z squared isn't, the loop continuous and cycles back to the start of the loop, and re-attains a new x y and z, and returns that random x,y,z, provided it's within the condition we specified.
"""
So...
while True: infinite loop
x, y, z = random.uniform(etc) : gets random x y and z
if x**2 + y**2 + z**2 <= radius**2 : makes sure that x, y and z satisfy the condition to terminate the loop
break : ends the infinite loop
return x, y, z : returns the values! <3 
"""
