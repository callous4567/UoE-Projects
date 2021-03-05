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

def twod_iterator(B, v_0, t_0, timestep):
    velocity = [[v_0*np.cos(t_0), v_0*np.sin(t_0)]]
    v_mag = lambda v_x, v_y: np.sqrt(v_x**2 + v_y**2)
    g = 9.8106 # in metres per second^2
    accel_x = lambda v_x, B, v: -B*v_x*v
    accel_y = lambda v_y, B, v: -B*v_y*v - g
    position = [[0,0]] # List of lists... like the velocities.
    times = [0]
    token = []
    v_magnitudes = [v_0]
    print(velocity)
    while True:
        currenttoken = len(token)
        currentposition = position[currenttoken]
        currentvelocity = velocity[currenttoken]
        currentmagnitude = v_mag(currentvelocity[0], currentvelocity[1])
        currentacceleration = [accel_x(currentvelocity[0], B, currentmagnitude), accel_y(currentvelocity[1], B, currentmagnitude)]
        token.append(1) # Increase token count by 1 for the sake of the times. Length of token list is number of steps.

        new_velocity = [currentvelocity[0] + currentacceleration[0]*timestep, currentvelocity[1] + currentacceleration[1]*timestep]
        new_position = [currentposition[0] + new_velocity[0]*timestep, currentposition[1] + new_velocity[1]*timestep]

        velocity.append(new_velocity)
        position.append(new_position)
        v_magnitudes.append(v_mag(new_velocity[0], new_velocity[1]))


        newtoken = len(token)

        times.append(newtoken*timestep)
        if new_position[1] < 0:
            print(len(token))
            break
    print(len(position))
    print(len(times))
    print(len(velocity))
    return position, times, velocity, v_magnitudes
# returns position time and velocity and magnitudes

def main():
    time_delay("Welcome to the 2D projectile drag checkpoint.")
    v_0 = value_getter("Initial Velocity v_0")
    theta_0_deg = value_getter("Initial Trajectory theta_0 in degrees")
    theta_0 = np.deg2rad(theta_0_deg)
    beta = value_getter("Normalised Drag  Coefficient B")
    step = value_getter("Step interval in seconds")
    values = twod_iterator(beta, v_0, theta_0, step)
    x_vals = []
    y_vals = []
    x_vels = []
    y_vels = []
    v_mags = []
    combinedcoords = values[0]
    combinedvelocities = values[2]
    for d in combinedcoords:
        x = d[0]
        y = d[1]
        x_vals.append(x)
        y_vals.append(y)
    for d in combinedvelocities:
        vx = d[0]
        vy = d[1]
        x_vels.append(vx)
        y_vels.append(vy)

    plt.plot(x_vals, y_vals)
    plt.show()
    print(len(values[3]))
    print(len(values[1]))

    # So apparently we need to iterate the iterator. Fucking long day.
    angle_list = np.linspace(0, np.pi/2, 300)
    velocity_list_angles = []
    for d in angle_list:
        anglesvalues = twod_iterator(beta, v_0, d, step)
        velociterinos = anglesvalues[2]
        final_velocity = velociterinos[len(velociterinos) - 1]
        magnitudefinal = np.sqrt(final_velocity[0]**2 + final_velocity[1]**2)
        velocity_list_angles.append(magnitudefinal)
    v_squared_list = [d**2 for d in velocity_list_angles]
    ratio_list = [d/(v_0**2) for d in v_squared_list]
    plt.plot(angle_list, ratio_list)
    plt.show()

main()