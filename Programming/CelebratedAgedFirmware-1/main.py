import numpy as np
import matplotlib.pyplot as plt

"""
defined_x = np.arange(0., 51., 1) # 0 to 50. Spacing of one, i.e. 0, 1, 2...
x = defined_x # For any single x...
plt.plot(x, x**2, 'bs', label="Quadratic") # Blue square notation. 
plt.grid()
plt.savefig("Quadratic.png")
plt.legend()
# To plot multiple things on this axis, we can also do this.
plt.plot(x, x**2, "bs", x,  x**3, "ro") 
plt.savefig("Quadratic.png")
# That puts in two sets of data.
# Note that all of the above is put back int othis one graph. With plt, if it isn't inside a function, you'll end up getting everything in one graph. So... yeah.
"""

# The next step is to define figures. Example:
"""
names = ['group_a', 'group_b', 'group_c'] 
values = [1, 10, 100]

plt.figure(1, figsize=(9, 3))

plt.subplot(131)
plt.bar(names, values)
plt.subplot(132)
plt.scatter(names, values)
plt.subplot(133)
plt.plot(names, values)
plt.suptitle('Categorical Plotting')
plt.show()
"""
# We define the names and values, just some random set.
# subplot(nrows, ncols, index, **kwargs)
# For each subplot, we specify the number of rows columns and the index position for the columns. 
# So, for 1, 3, 3,... the graph has one row, three columns, and this subplot is for the third column.
"""
x = []
for d in range(0, 50):
  x.append(d)
plt.figure(1, figsize=(10, 4))

plt.subplot(1, 10, 1)
plt.scatter(x, [x**2 for x in x], label="Quadratic")
plt.subplot(2, 10, 2)
plt.bar(x, [x**2 for x in x], label="Bar Quadratic")
plt.subplot(3, 10, 3)
plt.plot(x, [x**2 for x in x])

plt.savefig("frick.png")
"""
# We can see that the first graph occupies the full height, 10, and is on the first row.
# The second is on the 2nd row, with full height, but the row limits it to half that.
# Etc.

# Okay. Working with multiple figures and axes.
# Let's define some functions.
"""
def func(t):
  return np.cos(t)*np.exp(-t) # Function of time. Exponentially decaying cosine.

t_one = np.arange(0, 50, 0.05)
t_two = np.arange(0, 25, 0.05)

plt.figure(1, figsize=(9,4))
plt.subplot(2, 2, 1)
plt.plot(t_two, [func(t) for t in t_two])
plt.subplot(1, 2, 2)
plt.plot(t_one, [func(t) for t in t_one])


plt.savefig("hello.png")
"""
# 2, 2, 1 specifies it as being in the 2nd row up, with there being two columns, with it being in the first position per the index.
# 1, 2, 2 specifies one row, hence why it occupies entire height, 2 columns, with 2nd position per the 2nd index.