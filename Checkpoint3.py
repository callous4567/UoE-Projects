import numpy as np
import matplotlib.pyplot as plt
import time
import random
import math
import sys
import scipy

def time_delay(text):
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.001, 0.005)
        time.sleep(numbertime)
    print()
def value_getter(name):
    time_text = ("Please provide the value for {}").format(name)
    time_delay(time_text)
    value = float(input())
    print(("The value for {} is {}").format(name, value))
    return(value)
def ode_mode(): # Returns m, b, k, gamma, omega.
    time_delay("Provide coefficients involved... mx** + bx* + kx = 0...")
    m = value_getter("m")
    b = value_getter("b")
    k = value_getter("k")
    if k == 0:
        print("You are going to get an error with k = 0.")
    else:
        print("k =/= 0 :)")
    gamma = b/m
    omegasq = k/m
    omega = np.sqrt(omegasq)
    a = 1
    p = np.sqrt(((gamma**2)/4) - omegasq)
    omega_true = omegasq - (gamma**2)/4
    return(gamma, omega, a, p, omega_true)
def COEF_mode(): # Returns [Gamma, Omega, a, p]
    time_delay("Provide gamma and omega...")
    gamma = value_getter("Gamma")
    omega = value_getter("Omega")
    a = 1
    p = np.sqrt(((gamma**2)/4) - omega**2)
    omega_true = omega**2 - (gamma**2)/4
    return(gamma, omega, a, p, omega_true)
# Both modes return [GAMMA, OMEGA, A, P, OMEGA_TRUE]

# Functions 1,2,3 for situations 1,2,3.
eq1 = lambda gamma, t, a, p, b: np.exp((-1)*(gamma*t)/2)*(a*np.cosh(p*t) + b*np.sinh(p*t))

"""
def eq1(gamma, t, a, p, b):
    p1 = np.exp(-(1/2)*gamma*t)
    p2 = a*np.cosh(p*t) + b*np.sinh(p*t)
    p3 = p1*p2
    return p3
"""

eq2 = lambda gamma, t, a, b: np.exp((-1)*(gamma*t)/2)*(a + b*t)
eq3 = lambda gamma, t, a, omega_true, b: np.exp((-1)*(gamma*t)/2)*(a*np.cos(omega_true*t) + b*np.sin(omega_true*t))





def gomegachecker(lists):
    if lists[0] > 2*lists[1]:
        print("1")
        b = lists[0]/(2*lists[3])
        timespan = np.linspace(0, 5*np.pi/lists[1], 201)
        itemspan = [eq1(lists[0], d, 1, lists[3], b) for d in timespan]
        return(timespan, itemspan)

    elif lists[0] == 2*lists[1]:
        print("2")
        b = lists[0]/2
        print(b)
        timespan = np.linspace(0, 5*np.pi/lists[1], 201)
        itemspan = [eq2(lists[0], d, 1, b) for d in timespan]
        return (timespan, itemspan)


    else:
        print("3")
        b = lists[0]/(2*lists[4])
        print(b)
        timespan = np.linspace(0, 5*np.pi/lists[1], 201)
        itemspan = [eq3(lists[0], d, 1, lists[4], b) for d in timespan]
        return (timespan, itemspan)
# Returns the timespan and itemspan respectively in list form. Time is X and itemspan is Y



# Decides the comparisons. 1,2,3 decide gamma>2omega, gamma=2omega, or gamma<2omega. Debug purposes.
def comparitor():
    time_delay("Are you going to be directly inputting the ODE, or the omega and gamma coefficients?")
    time_delay("ODE or COEF")
    type = str(input())
    if type == ("ODE"):
        print("ODE Mode.")
        lists = ode_mode()
        ranges = gomegachecker(lists)
        plt.plot(ranges[0], ranges[1])
        plt.show()

    elif type == ("COEF"):
        print("COEF Mode.")
        lists = COEF_mode()
        ranges = gomegachecker(lists)
        plt.plot(ranges[0], ranges[1])
        plt.show()

    else:
        print("Sorry. That fucked up.")

comparitor()

