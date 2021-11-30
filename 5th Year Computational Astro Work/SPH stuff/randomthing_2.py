import numpy as np
import matplotlib.pyplot as plt

# Define function. Four variables: one dependent (t)
x = lambda amp, omega, t, phi: amp*np.sin(omega*t - phi)

# Set spring constant and mass for our "generator"
k = 10
m = 2
omega = np.sqrt(k/m)
time_for_two = 4*np.pi/omega
phi_deg = 60
amplitude = 5

# Set time domain/number of measurement points.
t = np.linspace(0, time_for_two, 10000)

# Generate noisy data (fractional noise.)
x_of_t = x(5, omega, t, np.deg2rad(phi_deg))
for num, xval in enumerate(x_of_t):
    noise = np.random.default_rng().normal(loc=xval, scale=1*amplitude, size=1)
    x_of_t[num] += noise

# Plot data example
plt.scatter(t, x_of_t, s=1)
plt.xlabel(r'$t/s$')
plt.ylabel(r'$x(t)/m$')
plt.savefig("example_1.png")
plt.show()