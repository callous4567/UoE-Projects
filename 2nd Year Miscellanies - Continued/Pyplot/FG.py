# Okay. This will be our graphing document.
import matplotlib.pyplot as plt
import numpy as np
import sys
import random
import time
import tkinter
import math
import cmath
import matplotlib.ticker as pltick

def time_delay(text):
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform( 0.001, 0.05)
        time.sleep(numbertime)
# Just our bog-standard time-delay text format.
def list_getter(variable_name):
    print()
    time_delay("How many variables in your list?")
    try:
        variable_lem = int(input())
    except:
        print()
        time_delay("Integers only, thanks.")
        quit()
    true_length = variable_lem - int(1)
    time_delay("Please give the value for each variable in the list...")
    list = []
    for d in range(int(1), variable_lem + int(1)):
        print(("Please give the value of element {0:d}: ").format(d))
        input_value = float(input())
        list.append(input_value)
    return list, variable_lem
# This returns us a list of the variables we needed alongside the list length. Just a key note to have and to show how multiple variables work. You get a tupled result, i.e. [[list],[variable_lem]], and have to call each one individually.
"""
def func(t):
    func = float(10)*np.exp(-t/4)*np.cos(2*np.pi*t)
    return func
# Just a random function for us to use...
# Note, to get a list of values, np.arange(start, end, print_every)... i.e. np.arange(0, 11, 1) and graphed will give 0 to 10, with a graph point every 1... If you make it 0, then a smooth curve.
# So...
# Let's define tdata and ydata first.
tdata = np.arange(0, 11, 0.2) # An array, 0 to 10, every 0.2
ydata = [func(t) for t in tdata]
plt.figure(1, figsize=(9,3)) # Defines the size ratio, x to y... 9 long, 3 high
plt.subplot(1, 2, 1) # Defines subplot location, i.e. (1, 2, 1): occupies entire row height, but there are 2 columns, it's in the 1st location.
plt.scatter(tdata, ydata, label="Graph") # Plots tdata on x against ydata on y. We need to define the tdata and ydata.
plt.xlabel("Time")
plt.ylabel("Y Position")
plt.legend()
plt.title("Graph of Y versus T")
plt.grid(True)
plt.subplot(1, 2, 2) # Second position, i.e. index = 2
plt.plot(tdata, ydata)
plt.show()
"""


# Yep. That all works fine. The next steps of the tutorial is text.
"""
func = lambda e: np.pi*np.exp(-0.5*np.cos(2*np.pi*e)) # A fucked up func
x = np.arange(0, 20, 0.2) # Spacing of 0.2 for points.
plt.plot(x, [func(x) for x in x], 'b', label="Weird")
plt.xlabel("X")
plt.ylabel("Fucked up func(x)")
plt.title("Graph")
plt.legend()
plt.text(5, 200, "Hello World") # This is the first example. Next up, symbols.
plt.text(10, 400, r'$\alpha > \beta$')
# You use r to denote the start. $ is the syntax to encapsulate the texworks variables, and \alpha and \beta just signify the actual codes for the variables in latex.
plt.grid(True)
plt.axis([0, 20, 0, 500])
plt.show()
"""

# The next up will be more elaborate symbol use. So..
"""
plt.figure(1, figsize=(10,10))
plt.subplot(2, 1, 1)
func = lambda e: np.sin(e**2 + np.pi*e)
sigma = np.arange(0, 1000, 1)
func_sigma = [func(x) for x in sigma]
plt.plot(sigma, func_sigma, 'g', label="Strange Graph")
plt.title(("Graph of {}").format(r'$\sigma_\mu$ versus $funct(\sigma_\mu)$'))
plt.xlabel(r'$\sigma_\mu$')
plt.ylabel(r'$func(\sigma_\mu)$')
plt.grid(True)
plt.legend()
plt.show()
"""

# Now, using the ax commands. This is more complicated.
"""
fig = plt.figure(1) # Define a figure as fig.
ax = fig.add_subplot(1, 1, 1) # Defines our axis within this figure.
ax.set(xlim=[0, 10], ylim=[0, 50], title=("Example Graph {}").format(r'$\mu$'), ylabel="Y Axis", xlabel="X Axis") # Defines variabels of our axis.
# We used plt.xlabel and plt.ylabel, or plt.title, and plt.axis([xstart, xfinish]) and etc.... this just puts all that into ax.set, setting the variables of ax, hence our subplot.
# We can also plot on this axis. So, example. ax.plot([x], [y], label, etc) and ax.scatter...
x = np.arange(0, 11, 0.2)
ax.plot(x, [d**2 for d in x], linewidth=1, color="red", label="Normal Plot")
ax.scatter(x, [d**1/2 for d in x], marker='^', label="Scatter")
# We can also set xticks and yticks. ax.set(xticks=[], yticks=[]) sets ticks at each point defined, i.e. on the sides, but doesn't set the grid. So, a if xticks=[1] you get one singular tick on the entire x axis, only at one, LABELLED as 1.
# So... to get custom minor ticks on our axis...
ax.minorticks_on() # Enables the minor ticks.
ax.grid(True, which='minor') # Enables minor ticks on the grid
ax.grid(True, which='major') # Enables major ticks on the grid
# To set custom minor and major locators... First define the locator(s) using pltick, matplotlib.ticker!
x_major_locator = pltick.MultipleLocator(base=5)
x_minor_locator = pltick.MultipleLocator(base=0.2)
y_major_locator = pltick.MultipleLocator(base=5)
y_minor_locator = pltick.MultipleLocator(base=1)
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
plt.show()
"""
"""
# There you have it. Custom ticks, a grid, fancy labels, etc... all sorted
# So. To recount the method!
fig = plt.figure(1)
ax = plt.subplot(1,1,1) # One row, one column, all just one big graph.
ax.set(xlim=[0, 50], ylim=[20, 100], title="X v. Y", xlabel="X", ylabel="Y")
ax.grid(True, which='minor')
ax.grid(True, which='major')
ax.minorticks_on()
x = np.arange(0, 100, 1)
y = [d**2*np.cos(2*d*np.pi) for d in x]
ax.plot(x, x + 10, label="X versus Y", color='red', linewidth='1')
ax.legend()
locator_major = pltick.MultipleLocator(5)
locator_minor = pltick.MultipleLocator(1)
ax.xaxis.set_major_locator(locator_major)
ax.yaxis.set_major_locator(locator_major)
ax.xaxis.set_minor_locator(locator_minor)
ax.yaxis.set_minor_locator(locator_minor)
ax.margins(x=0, y=0)
plt.show()
"""
