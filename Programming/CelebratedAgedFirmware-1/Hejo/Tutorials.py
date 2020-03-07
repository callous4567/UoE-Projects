import math
import cmath
import time
import random
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import numpy as np

def time_delay(text):
  for c in text:
    sys.stdout.write(c)
    sys.stdout.flush()
    timer = random.uniform(0.002, 0.005)
    time.sleep(timer)

# That's our time delay script. Nice and important to have.
# Okay. Let's test out the elif, else and if statements.

def float_variable_getter(variable_name):
  variable_text = ("Please give your {} value...: ").format(variable_name)
  print()
  time_delay(variable_text)
  print()
  try:
    variable_value = float(input())
    return variable_value
  except:
    time_delay("Sorry. Invalid value. Try again <3")
    quit()

# That's our float variable obtainer. Checks to see if the input is a float. If not, try again.
# Next up is going to be the "List" obtainer. It obtains a list.

def list_obtainer(list_name):
  print()
  time_delay("Please give the number of elements in your list... Integer Value...: ")
  print()
  try:
    elem_number = int(input())
  except:
    time_delay("Sorry. Not an integer value :/... try again.")
    quit()
  time_delay("Now please give the value for each element of the list.")
  list_list = []
  print()
  for d in range(0, elem_number):
    value_d = float(input(("Value of the element {0:d}..:").format(d + int(1))))
    list_list.append(value_d)
  print(("Your final list is {}").format(list_list))
  return list_list

# That gets a list of floats. Problem solved there. So... yeah.
# Let's try to map each element of that list with a function.

def list_function(list_list):
  function = lambda d: d**2 + float(4)
  new_list = [function(elem) for elem in list_list]
  return new_list

# That list_function is just to put the list we obtained to a certain function and then return the list with
# the function applied. So... yeah.
# Next step... plotting graphs...
"""
The way the graphing works is by taking two lists, i.e. xdata and ydata, and then plotting the two lists by doing it via this .... doing.... plt.plot( xdata, ydata)
Once it plots the list data you end up with a nice looking graph. So, let's make a function, take an original range and list, and solve it!.

Note: Here is how to put matplotlib as a png...

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from random import randint

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

points = []
last = 0
bound = 100
for i in range(0, 100):
  last += randint(-bound, bound)
  points.append(last)
  
ax.plot(points)
fig.savefig('graph.png')
"""
###################################################################################################################
def t_function():
  time_delay("Please give the start of the interval, time-wise, and the end of the interval.")
  print()
  time_start = int(input("The start of the interval, Integer-wise...: "))
  time_end = int(input("The end of the interval, Integer-wise...: "))
  time_function = lambda t: t**2
  # This is s = ut + 0.5at^2 with initial velocity of 5, up, and 9.81 = g downwards.
  y_value_list = []
  t_value_list = []
  for d in range(time_start, time_end + int(1)):
    y_value = time_function(d)
    y_value_list.append(y_value)
    t_value = d
    t_value_list.append(t_value)
  y_value_list.sort()
  t_value_list.sort()
  return y_value_list, t_value_list
  
# Graphing time. # First let's get the Y and T lists.
list_tuple = t_function()
y_list = list_tuple[0]
t_list = list_tuple[1]
# We have the lists for y and t. Next, the easy way to do it would be this...
# plt.plot(t_list, y_list)
# plt.savefig("plot.png")
# But, we need to actually put on axes, a grid, etc... and add minor spacings, etc. So...
# Let's do it the "Technical" way.
# First, the plot.
plt.plot(t_list, y_list) # This plots X against Y, T against Y Position
plt.xlabel("Time")
plt.ylabel("Y Position")
plt.title("Graph of Y versus T")
# Next up is to add a grid.
plt.grid(b=True, which="major", linestyle="-") # Defines the major grid, w/ a const. styl
plt.grid(b=True, which="minor", linestyle="--") # Minor grid. Dashed/dotted style.
# Easy enough. Now then. Tutorial time.
