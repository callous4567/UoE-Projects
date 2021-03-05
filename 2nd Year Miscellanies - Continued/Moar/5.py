import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import time
import random
import sys
import sympy

def time_delay(text, delay_max):
    print()
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        resttime = random.uniform(0.001, delay_max)
        time.sleep(resttime)

def quantity_get_print(quant):
    time_delay(("Value for {0}:..").format(quant), 0.005)
    try:
        quant_val = float(input())
        return quant_val
    except:
        time_delay("You probably have a value error in there somewhere. I need a float!")
        quit()
    
def value_getter():
    time_delay(("Please provide values for initial velocity magnitude, angular deviation from the positive x axis, normalised drag coefficient, and desired step interval."),0.005)
    v_m = quantity_get_print("Velocity Magnitude")
    a_m = quantity_get_print("Angular Deviation from Horizontal +X")
    n_d_c = quantity_get_print("Normalised Drag Coefficient")
    s_i = quantity_get_print("Step Interval")
    return v_m, a_m, n_d_c, s_i

# So we have to generate an iterative method.
############################################################################
"""
We have a = -bv^2 * unit_v 
ax = -bv^2 cos(theta) 
ay = -bv^2 sin(theta) - g
vx0 = v cos theta0 vy0 = v sin theta0

so
x0 = 0 y0 = 0
x1 = x0 + vx0(dt) y1 = y0 + vy0(dt) vx1 = vx0 + ax(dt) vy1 = vy0 + ay(dt)
xi = x(i-1) + vx(i-1)(dt), etc
and the theta(i) equals arctan{vy(i)/vx(i)}
"""
############################################################################

def iterator(v0, a0, B, step):
    x_time, x_velocity, y_time, y_velocity, time, theta = [],[],[],[],[],[]
    x_time.append(0)
    theta.append(a0)
    y_time.append(0)
    x_velocity.append(v0*np.cos(a0))
    y_velocity.append(v0*np.sin(a0))
    time.append(0)
    time_token = 0
    gravity = float(9.81) # Configurable if you'd like. Opted not to.
    
    while True:
        vel_mag = np.sqrt(x_velocity[time_token]**2 + y_velocity[time_token]**2)
        x_accel = lambda L,v,angle: float(-1)*L*(v**2)*np.cos(angle)
        y_accel = lambda L,v,angle,grav: float(-1)*(L*(v**2)*np.sin(angle) + grav)
        new_vx = x_velocity[time_token] + (x_accel(B, vel_mag, theta[time_token]))*step
        new_vy = y_velocity[time_token] + (y_accel(B, vel_mag, theta[time_token], gravity))*step
        new_x = x_time[time_token] + x_velocity[time_token]*step
        new_y = y_time[time_token] + y_velocity[time_token]*step
        new_theta = np.arctan(new_vy/new_vx)
        theta.append(new_theta)
        time.append(step)
        x_time.append(new_x)
        y_time.append(new_y)
        x_velocity.append(new_vx)
        y_velocity.append(new_vy)
        time_token = time_token + int(1)
        if new_y == 0:
            break

    return x_time, x_velocity, y_time, y_velocity, time


new_example = iterator(4, np.pi/6, 4, 0.1)
print(new_example)
