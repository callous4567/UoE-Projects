import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

x_range = [0,4]
y_range = [-4,4]
z_range = [-4,4]
t_range = np.linspace(0, 10000, 10000)
xx = np.linspace(x_range[0], x_range[1], 500)
yy = np.linspace(y_range[0], y_range[1], 500)
zz = np.linspace(z_range[0], z_range[1], 500)
k1, w1, d1 = 1, 0.1, np.pi/2
k2, w2, d2 = 1, 0.2, 0
dt = t_range[1] - t_range[0]



# Run the simulator graphically (plots the free energy/etc.)
def run():
    # Interactive On
    plt.ion()

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection='3d')
    ax.set(xlabel='x', ylabel='y', zlabel='z')

    for t in t_range:
        ax.set(xlabel='x', ylabel='y', zlabel='z')
        ax.plot([x_range[0] for d in t_range], np.sin(k1*x_range[0] - w1*t_range + d1), np.sin(k2*x_range[0] - w2*t_range + d2), color='green')
        ax.plot(xx, np.sin(k1*xx - w1*t + d1), np.zeros_like(zz), color='red')
        ax.plot(xx, np.zeros_like(yy), np.sin(k2*xx - w2*t + d2), color='blue')
        ax.plot(xx, np.sin(k1*xx - w1*t + d1), np.sin(k2*xx - w2*t + d2), color='green')
        fig.canvas.draw()
        fig.canvas.flush_events()
        ax.clear()
    # All done.
    plt.close()

run()

