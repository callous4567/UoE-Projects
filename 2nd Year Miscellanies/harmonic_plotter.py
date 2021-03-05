"""
Simple Python code that plots the cosine function
"""

# Import relevant python modules
import math
import matplotlib.pyplot as pyplot
from OOP1 import utilities
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import matplotlib as mpl
import sys

# Define the cosine function
def my_function(x):
    return math.cos(x)

# Summerator. The summation iterator. Give it an x and an N and it'll give you the sum to the Nth term.
def summerator(x, N):
    funct = lambda x, i: (1 / (2*i - 1))*np.sin((2*i - 1)*x) # Lambda funct for a given x and i based on the fourier series for the square wave provided
    sum = 0
    for d in range(1, N + int(1)): # Iterates using lambda function from summation bounds 1 -> N
        sum+= funct(x, d) # Adds the calculated sub-sum to the total sum
    return sum

# Grapher that makes it look nicer.
def grapher(x, y):
    figure = plt.figure(111) # Define figure
    ax1 = plt.subplot(111) # Subplot on figure
    ax1_label = "f(x) =" + " " + r'$\sum_{1}^N {\frac{sin[(2\it{i} - 1)x]}{2\it{i} - 1}}$' + " " + "from -$\pi$ to +$\pi$" # LaTeX label showing function
    ax1.plot(x, y, color='black', linewidth=0.5, label=ax1_label) # Plot graph

    mpl.rc('font', size=8) # Redefine matplotlib rcParams to change font size

    major_ticks_x = pltick.MultipleLocator(np.pi/4) # Decide/calculate locators for the major and minor ticks on subplot
    major_ticks_y = pltick.MultipleLocator(max(y) / 4)
    minor_ticks_x = pltick.MultipleLocator(np.pi/8)
    minor_ticks_y = pltick.MultipleLocator(max(y)/8)

    ax1.xaxis.set_major_locator(major_ticks_x) # Apply locators to ax1 subplot
 #   ax1.yaxis.set_major_locator(major_ticks_y)
    ax1.xaxis.set_minor_locator(minor_ticks_x)
 #   ax1.yaxis.set_minor_locator(minor_ticks_y)

    ax1.grid(True, which='major', color="blue", alpha=0.1) # Enable grids on subplot
    ax1.grid(True, which='minor', color="pink", alpha=0.4)


    plt.legend() # Shows the legend/label
    plt.show()
    figure.savefig("heyo.png") # Saves figure. Always do it (idk, the old code that I used in SciProg always had this since I needed to save it and it's habit.)

# Main method
def main():
    # number of data points and the number of iterations respectively
    n_loop = int(utilities.value_getter("Number of datapoints")) # Known bug in pycharm.
    N = int(utilities.value_getter("Number of iterable loops"))
    # Define filename
    file_name = utilities.string_getter("output file name")

    # open output file with desired filename
    out_file = open(file_name + ".dat","w")

    # prepare data lists
    x_values = []
    y_values = []

    # obtain function values and write them to file
    for i in range(n_loop):
        x = 2*math.pi*i/n_loop - math.pi
        f = summerator(x, N)
    
        # append data to lists and output file
        x_values.append(x)
        y_values.append(f)

        out_file.write(str(x) + " " + str(f) + "\n")

    # close output file
    out_file.close()

    # plot result
    grapher(x_values, y_values)
    print(len(sys.argv))
    """
    pyplot.plot(x_values,y_values)
    pyplot.suptitle('Plotting the cosine function')
    pyplot.xlabel('X')
    pyplot.ylabel('Cos(X)')
    pyplot.show()
    """

# Execute main method
if __name__ == "__main__": main()