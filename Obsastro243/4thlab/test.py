import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants

h, c, k, T = constants.h, constants.c, constants.Boltzmann, 5000
x = np.linspace(0, 1*10**-6, 1*10**-9)

def plancker(wl):
    return ((h*c**2)/(wl**5))/(np.exp(h*c/(wl*k*T)) -1)

