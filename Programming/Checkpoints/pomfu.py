import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys
import time
import random
import matplotlib.ticker as pltick

# This is the one that uses odeint

# Old handy time-delay function...
def time_delay(text):
    t = text
    for c in t:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.002)
        time.sleep(numbertime)



# We have the equation mx'' + bx' + k = 0
# Hence x'' + b/m x' + k/m = 0
# Let gamma equal b/m and omega_0^2 equal k/m
# Hence x'' + gamma x' + o^2 = 0
# We need to dissolve this into first order linear ODEs and use scipy odeint
# So...
# First I'll test out using odeint on a normal equation... i.e.
# " y'(t) = (-4 + k)y(t)

"""
def funct(y,t):
    time_delay("Input your k value")
    print()
    k = float(input())
    print()
    dy_dt = (float(-4) + k)*y
    return dy_dt
# That gives a function func(y,t) and returns dy/dt for each t in the range
# So... Initial conditions

y_0 = float(10)
t = np.linspace(0, 100, 500)

# The solver...
y = odeint(funct, y_0, t) # You don't specify the y or t as funct(t)... funct is just the actual equation you use. Odeint puts in the y_0 as the value at t=0 and pushes through all the t values, takes all the resultant equations it gets for y: it takes differential of y'(t) you supply (-ky) and the initial condition for y, and integrates across the time range to give you y(t)
# Then you can plot it.

# Another example...
# y'(t) + y**2 + y + 4 = 0

def new_funct(y, t): # Y as a function of Y
    new_function = -y**2 -y + float(-4)
    return new_function # You want to define the differential of y, i.e. new_funct equals dy/dt for y, t...

y_1 = float(50) # Initial Y condition
time = np.arange(0, 100, 0.2) # All the range
new_y = odeint(new_funct, y_1, time) # Initial Y condition, and time period.

def funco(y, t):
    y_deriv = y**2 + y 
    return y_deriv
y0 = 10
t = np.linspace(0, 50, 400)
all_y = odeint(funco, y0, t)
"""



# Okay. So... we have (x'') + gamma(x') + o^2(x) = 0 and need to solve that for some initial gamma and o. We need to split this into linear first order ODE's or go find another function to use other than odeint.

# Letting u = x', we have:
# x '' + gamma x' + o^2 x = 0 and u = x'
# u' + gamma u + o^2 x = 0 so u' = - gamma u - o^2 x
# Okidoke. Let's go.
# x'' + gamma x' + km x
def condition_getter(condition):
    time_delay(("Please give {} value...: ").format(condition))
    print()
    try:
        condition_Value = float(input())
        return condition_Value
    except:
        print("Try giving a numerical value.")
        quit()
        
def x_funct(X, t): # X is an array of [X, X']
    return [X[1], float(-1)*gamma*X[1] + float(-1)*(omega**2)*X[0]]

# So derivatives of elements of X[0] equals X[1], i.e. x't and for X[1] equal to x''t. 2nd order linear equation solved by substitution to two first orders...
# So [y, y'] derivatives to [y', y''] and the function solves through y'' to y' and y' to y, vicariously getting y. odeint knows how to handle the dependencies between its provided lists to work its way back to the original variable you want.

# Define the initial conditions!
gamma = condition_getter("Gamma")
omega = condition_getter("Omega")
initial_x,initial_v = 1,0
maxi =  float(5)*np.pi*float(1/(omega)) # Add condition_getter(time) to get.
time_range = np.linspace(0, maxi, 500)
y_range = odeint(x_funct, [initial_x, initial_v], time_range) # Tuple of lists
# The tuple has [y values over range, y' values over range]
# Next up is to isolate the y.
y_true = y_range[:,0] #(gets you the left side of the list, i.e. the y values.)

fig = plt.figure(1, figsize=(80, 20))
ax = plt.subplot(111)
ax.set(xlim=[0, maxi], ylim=[float(-1), float(1)], xlabel=r'$\it\bf{t}$', ylabel=r'$\it\bf{x}$')
ax.plot(time_range, y_true, linestyle="--", color="lightblue", linewidth=1.5, label=("Graph of particle displacement {0} versus time {1}").format(r'$\it\bf{x}$', r'$\it\bf{t}$'))
plt.legend()
ax.grid(True, which='major')
ax.grid(True, which='minor')
ax.minorticks_on()
int_round_maxi = int(round(maxi))
majorticks_t = pltick.MultipleLocator(int_round_maxi/20)
majorticks_x = pltick.MultipleLocator(0.2)
minorticks_t = pltick.MultipleLocator(int_round_maxi/40)
minorticks_x = pltick.MultipleLocator(0.1)
ax.margins(x=0, y=0)
ax.xaxis.set_major_locator(majorticks_t)
ax.yaxis.set_major_locator(majorticks_x)
ax.xaxis.set_minor_locator(minorticks_t)
ax.yaxis.set_minor_locator(minorticks_x)
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig("latest.png")
plt.show()

