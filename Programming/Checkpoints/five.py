import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import time
import random
import sys

def time_delay(text, maxi): # Reliable time_delay function to delay text writing, It looks super fancy.
    print()
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        time_rest = random.uniform(0.001, maxi)
        time.sleep(time_rest)
    print()

def quantity_getter(name, typesan): # Gets the quantity you want. Modified to adapt the angle from degrees to radians, if that's your tonic.
    if typesan == int(0):
        time_delay(("Please state the value for your {0}...  :").format(name), 0.005)
        try:
            name_value = float(input())
            return name_value
        except:
            time_delay("Sorry. An error occurred. Probably a value error (give a real value for any quantity.", 0.005)
            quit()

    elif typesan == int(1):
        time_delay("Is your angle in degrees or radians... ?", 0.005)
        typeo = str(input())
        list_san = list(typeo)
        try:
            if list_san[2] == "g":
                time_delay("Please give the value of the angle... :", 0.005)
                try:
                    value = float(input())
                    radian_value = value*(1/180)*np.pi
                    return radian_value
                except:
                    time_delay("Real values, please and thanks.", 0.005)
                    quit()
            elif list_san[2] == "d":
                time_delay("Please give the value of the angle... :", 0.005)
                try:
                    value = float(input())
                    return value
                except:
                    time_delay("Real values, please and thanks.", 0.005)
                    quit()
                else:
                    time_delay("Sorry to say that the type of angle you specified was invalid... answer lowercase, degrees or radians.", 0.005)
                    quit()
        except:
            time_delay("Sorry. Probably a value error. Degrees and radians, not 0 and 2's.", 0.005)
            quit()

    else:
        time_delay("Sorry. A coding error has occurred. Please consult your nearest pythonic programmer.", 0.005)
        quit()
            

def constants_list(): # Returns all the specified constants.
    initial_v = quantity_getter("initial Velocity Magnitude", int(0))
    initial_theta = quantity_getter("initial angle above horizontal positive X axis", int(1))
    initial_B = quantity_getter("Normalised Drag Coefficient", int(0))
    initial_step = quantity_getter("desired step interval", int(0))
    return initial_v, initial_theta, initial_B, initial_step

def iterator(constants): # Iterates, Euler-style

    initial_v = constants[0] # Define initial values... makes it easier on the eyes.
    initial_theta = constants[1]
    initial_B = constants[2]
    initial_step = constants[3]
    
    x_list = [0] # Initial lists
    y_list = [0] #
    
    vel_y = [initial_v*np.sin(initial_theta)] #
    vel_x = [initial_v*np.cos(initial_theta)] #

    vel_mag_list = [initial_v] #
    
    step_list = [0] #
    time_token = int(0) # Time Token! Lets us keep track of which iteration we're in. Also useful for debugging to see if anything is actually happening.

    grav = float(9.81) # A fundamental constant.

    accel_y_inst = lambda B, vel, vely, g: (float(-1)*B*vel*vely - g)
    accel_x_inst = lambda B, vel, velx: (float(-1)*B*vel*velx)
    vel_square = lambda vely, velx: (np.sqrt(vely**2 + velx**2))

    while True: # Infinite loop using while True, the actual iterator.
        new_x = x_list[time_token] + initial_step*vel_x[time_token]
        new_y = y_list[time_token] + initial_step*vel_y[time_token]
        
        v_mag = vel_mag_list[time_token]
        new_velx = vel_x[time_token] + initial_step*accel_x_inst(initial_B, v_mag, vel_x[time_token])
        new_vely = vel_y[time_token] + initial_step*accel_y_inst(initial_B, v_mag, vel_y[time_token], grav)
        new_v_mag = vel_square(new_vely, new_velx)

        x_list.append(new_x)
        y_list.append(new_y)

        vel_x.append(new_velx)
        vel_y.append(new_vely)
        vel_mag_list.append(new_v_mag)

        time_token = time_token + int(1)
        
        step_list.append(time_token*initial_step)

        if new_y <= 0: # break is just to stop the loop. All the prints are just for debugging purposes, so, ignore them. 
         #   print(new_y)
         #   print(len(y_list))
         #   print(len(x_list))
         #   print(len(vel_x))
         #   print(len(vel_y))
         #   print(len(step_list))
         #   print(len(vel_mag_list))
            break

    kinetic_factor = [d**2/vel_mag_list[0]**2 for d in vel_mag_list]
    # Was originally going to plot K_factor(time), but I misread the actual question.
    # We're not concerned with processing power, but I will acknowledge that just calculating it outright using the last element versus first element of the vel_mag array would be better.
    
    return x_list, y_list, vel_y, vel_x, step_list, vel_mag_list, kinetic_factor, constants

def kinetic_iterator(constants): # Iterates the iterator to get kinetic factor as a function of theta.
    angles = [] # Define lists
    factors = [] #
    for d in range(0, 401): # Iterate. 400 points makes for a smooth curve.
        angular_equivalent = (d/400)*np.pi*1/2
        angles.append(angular_equivalent)
        new_constants = [constants[0], angular_equivalent, constants[2], constants[3]]
        angle_iterable = iterator(new_constants)[6]
        factors.append(angle_iterable[len(angle_iterable) - 1]) 
    return angles,factors # You are returned a list of angles and final kinetic energy factor.
        
def rounder(value, b): # Round value to a base b with accuracy defined by b
    return round(value/b)*b 

def log_rounder(value): # Round and get logarithmic base of value, i.e. 450 -> 100. 950 -> 1000. 
    log_factor = np.log10(value)
    return int(10)**rounder(log_factor, 1)

def graph_maker(iterable, kinetirable): # Need to graph y(x) and K(t) # V theta Beta step
    
    constants = iterable[7]
    
    fig1 = plt.figure(111, figsize=(24, 48))
    ax1 = plt.subplot(211) # y(x)
    ax2 = plt.subplot(212) # k(theta)

    # Plot the trajectory
    ax1.set(xlim=[0, max(iterable[0])], ylim=[0, max(iterable[1]) + 0.1*log_rounder(max(iterable[1])) ], title=("Trajectory of particle for specified initial conditions: {0} = {1:.1e}, {2} = {3:.1e}, {4} = {5:.1e}, {6} = {7:.1e}").format(r'$\it{\beta}$', constants[2], r'$\it{v_0}$', constants[0], r'$\it{\theta_0}$', constants[1], r'$\it{Step}$', constants[3]), xlabel="$\it{x}$ (m)", ylabel="$\it{y}$ (m)")
    ax1.plot(iterable[0], iterable[1], label="y(x)", color="lightblue")

    x_loc = [pltick.MultipleLocator(log_rounder(max(iterable[0]))/10), pltick.MultipleLocator(log_rounder(max(iterable[0]))/20)] # Set custom ticks to make a prettier graph.
    y_loc = [pltick.MultipleLocator(log_rounder(max(iterable[1]))/5), pltick.MultipleLocator(log_rounder(max(iterable[1]))/10)] #

    ax1.xaxis.set_major_locator(x_loc[0]) #
    ax1.xaxis.set_minor_locator(x_loc[1]) #
    ax1.yaxis.set_major_locator(y_loc[0]) #
    ax1.yaxis.set_minor_locator(y_loc[1]) #

    ax1.grid(True, which='major', color="red") # Enable the gridlines... for a prettier graph.
    ax1.grid(True, which='minor', color="pink")
    ax1.minorticks_on()

    # Plot the K_factor(angle)
    ax2.set(xlim=[0, np.pi/2], ylim=[0, 1.1], title=("KE Ratio vs. {4} for specified initial conditions: {0} = {1:.1e}, {2} = {3:.1e}, {6} = {7:.1e}").format(r'$\it{\beta}$', constants[2], r'$\it{v_0}$', constants[0], r'$\it{\theta_0}$', constants[1], r'$\it{Step}$', constants[3]), xlabel=(r'$\it{\theta_0}$' + ' ' + r'($\it{^o}$)'), ylabel="$\it{K}$ Ratio")
    ax2.plot(kinetirable[0], kinetirable[1], label="K($\it{\theta}$)", color="blue")

    plt.xticks(np.linspace(0, np.pi/2, 19), np.arange(0, 95, 5)) # Custom ticks. Angle was in radians. So, nineteen majors from 0 to pi/2 results in every five degrees. 

    ax2.grid(True, which='major', color="red") # Grids again.
    ax2.grid(True, which='minor', color="pink")
    ax2.minorticks_on()

    fig1.savefig("Hejoooo.png") # I always do this incase I need to open the file using an online editor, which may or may not support plt.show(). I know some don't.
    plt.show()

# Run the functions in order. 
hejo = constants_list()
iterable = iterator(hejo)
kinetirable = kinetic_iterator(hejo)
graph_maker(iterable, kinetirable)
