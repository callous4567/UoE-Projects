import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
from matplotlib.patches import Patch
import matplotlib as mpl

from scipy.optimize import curve_fit
from scipy import constants
from lmfit import minimize, Parameter, Parameters
import lmfit as lmf
from scipy.signal import savgol_filter

import time
import pandas as pd
import io
import sys
from OOP1 import utilities as utils
# Proudly callicioused by coded.

class parameterfile_handler(object):
    def __init__(self):
        self.dummy = 0

    # EXAMPLE DATA
    """
    Omitting txt and dat... body is Moon.dat
    Orbital period is 366.0 days. 
    Semimajor, semiminor, eccentricty, apoapsis, and periapsis are respectively [0.999938756485758, 0.9990307579276912] 0.04260614719279441 1.04254229432837 0.957335218643146
    """
    def pararipper(self, file):
        fileopened = open(file)
        filelines = fileopened.readlines()

        periodoline = filelines[1]
        heyo = periodoline.split("Orbital period is")
        heyo = heyo[1]
        heyo = heyo.split("days")
        period = float(heyo[0])

        importoline = filelines[2]
        importosplit = importoline.split("respectively")
        data = importosplit[1]
        hey = data.replace("[", " ")
        hey = hey.replace("]", " ")
        hey = hey.replace(",", " ")
        hey = hey.split()
        hey.append(period)
        hey = [float(d) for d in hey]
        fileopened.close()
        return hey

    def listgetter(self, bodyfile):
        fileopen = open(bodyfile)
        filelines = fileopen.readlines()
        namelist = []
        for line in filelines:
            linesplit = line.split(",")
            if linesplit[0] == "Sun":
                pass
            else:
                namelist.append(linesplit[0])
        filenamelist = [d + ".dat-Solved-Parameters-.txt" for d in namelist]
        filenamelist.append("Luna.txt-Solved-Parameters-.txt")
        namelist.append("Luna")
        return filenamelist, namelist

    def bulkhandler(self, namelist):
        datalist = []
        for d in namelist[0]:
            data = self.pararipper(d)
            data = [str(d) for d in data]
            data = data
            datalist.append(data)
        datatranspose = np.array(datalist).T.tolist()
        datatranspose.append(namelist[1])
        file = open("finaldata.csv", "w+")
        for d in datatranspose:
            file.write(str(d) + "\n")
        file.close()


    def main(self):
        filelist = self.listgetter("particles.dat")
        self.bulkhandler(filelist)

instance = parameterfile_handler()
instance.main()



