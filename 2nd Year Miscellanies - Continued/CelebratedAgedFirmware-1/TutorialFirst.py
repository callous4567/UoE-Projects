import numpy as np
import matplotlib.pyplot as plt

# plt.plot(x_set, y_set) defines the plot of x against y
# plt.axis(x min, x max, y min, y max) defines axis limits. 
x = np.arange(0, 20, 50)
plt.plot(x, x**2)
plt.axis(0, 10, 0, 50)
plt.xlabel("x")
plt.ylabel("y")
plt.title("The graph")
plt.legend()
plt.savefig("plot.png")