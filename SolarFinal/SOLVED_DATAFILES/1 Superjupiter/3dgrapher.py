import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D # Redundant here but will be used in future and thus is kept for my use.
import matplotlib.ticker as pltick
from matplotlib.patches import Patch
import matplotlib as mpl

import sys
import time
import random
import pandas as pd
import io
import threading as thread
# Proudly callicioused by coded.


# EXAMPLE FORMAT
"""

12 # NUMBER OF OBJECTS. USEFUL. 
Point = 0 # Step. PYTHONIC step! 
Sun -0.004138911767258502 0.00735004172952756 3.241351872221244e-05
Mercury 0.1028291409161863 0.2958105653291369 0.01379108978698249
Venus 0.2452359016821033 0.6838213962971565 -0.005076132664660649
Earth -0.7756973325652013 0.6224693028784704 3.743540802014699e-06
Mars -0.9070804697534424 -1.240660408192739 -0.003964594547707404
Jupiter 0.8253100420406844 -5.138826959489569 0.002848999155655939
Saturn 3.99154270746388 -9.19249314048187 0.0009362653163950606
Uranus 16.1273500761803 11.5105824229213 -0.1661812043097535
Neptune 29.26524292268337 -6.233162987600283 -0.5460868878676082
Pluto 13.09469898830053 -31.34070948610209 -0.4341095185617982
Moon -0.7780163096409389 0.6230929317983712 0.0002034048255254121
Halley -20.30214358826866 26.61812318922205 -9.972768864700402

"""


# Class to take a VMD XYZ file format with known number of objects, decode it and plot it using 3D MPL toolkit.
# Not very efficient. Do *NOT* use when your filesize exceeds however much RAM you are willing to give.
class dataporter(object):
    def __init__(self, filetension):
        self.filename = filetension # Dummy object. Only way I know >.<
        self.bodies = 12 # Number of bodies involved, non-pythonic.

    # Return lists for each of the bodies, including their names.
    def vmd_porter(self):
        names = []
        # Generate a list of names from the VMD file. Needed for later.
        token = 0
        with open(self.filename) as f:
            for line in f:
                if token >= (self.bodies + 2):
                    break
                token += 1
                liner = line.split()
                if len(liner) == 4:
                    names.append(liner[0])
        datalists = np.array([[] for d in names])
        # Generate data from VMD file.
        token = 0
        baselist = np.array([])
        with open(self.filename) as f:
            for line in f:
                linesplit = line.split()
                if len(linesplit) == 4:
                    linesplit = linesplit[1:]
                    linesplit = np.array([float(d) for d in linesplit])
                    np.append(baselist, linesplit)









        print(names)




hey = dataporter("vmdoutput.xyz")
hey.vmd_porter()


