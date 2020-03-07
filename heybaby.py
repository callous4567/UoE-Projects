import numpy as np

print("Hey, what's the volume in mm^3?")
print()
try:
    volume = float(input())
    print()
    print("Arigatou Gozaimashta")
except:
    print("Sorry, an error has occurred. Give a float you dolt.")
print()
thanks = ("Your volume is {} mm^3").format(volume)
print(thanks)
# Okay. That's that done.

newvolume = volume*10**-9
radiuscalc = lambda v: ((3*v)/(4*np.pi))**(1/3)
radius = radiuscalc(newvolume)
print(radius)
print("Enjoy")