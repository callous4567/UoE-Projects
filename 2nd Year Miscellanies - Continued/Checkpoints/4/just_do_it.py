import sys
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import math

# Delays text release for more natural use.
def time_delay(text, max_delay):
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.001, max_delay)
        time.sleep(numbertime)
    print()

# Gathers the lists of voltage, current, and time.
def file_former(file_name, file_extension, linestart): # Gets the file
    file_file = open(("{0}.{1}").format(file_name, file_extension), "r")
    file_file_lines = file_file.readlines() # A tuple of strings. 0,1,2 = no.
    new_file_file = []
    for d in range(linestart - int(1), len(file_file_lines)):
        new_file_file.append(file_file_lines[d]) # remove first 3 with loop


    v_l = [] # Voltage
    c_l = [] # Current

    for line in new_file_file:
        new_line = line.split(",")
        try:
            v_l.append(float(new_line[0]))
            c_l.append(float(new_line[1]))
        except:
            time_delay("Integers, Numbers, and everything String... string?", 0.005)
            time.sleep(1)
            time_delay("... Hm ...", 0.5)
            time.sleep(1)
            time_delay("I think you may have mixed some words into your datastream. Try a different linestart.", 0.005)
            quit()

    file_file.close()

    t_l = np.arange(0, len(v_l)*float(4e-5), 4e-5) # Time w/ sampling rate
    # Sampling rate was 25khz and hence 1/25000 ~ 4e-5

    if len(v_l) == len(c_l):
        return v_l, c_l, t_l
    else:
        time_delay("Sorry. Your original file was a bit messed up.", 0.005)
        quit()

# Logarithmic function that was described (log of v*i)
def log_function(v_list, c_list):
    p = []
    for d in range(0, len(v_list)):
        p_d = np.log(v_list[d]*c_list[d])
        p.append(p_d)
    return p

def rounder(value):
    log10 = np.log10(value) # Gets the log of the value, to base 10
    log10_floored = math.floor(log10) # Rounds down to closest int
    return 10**log10_floored # Returns 10^the log base of the value.

# Graphs p(t)
def grapher(p_list, t_list):
    fig = plt.figure(111, figsize=(24,12))
    ax = plt.subplot(111)
    label = "Graph of $\it{p(t)}$  = $\ln(V(t)*I(t))$ versus time $\it{t}$ over specified interval. Sampling rate of 25 kHz"
    ax.set(ylim=[min(p_list), max(p_list)], xlim=[0, max(t_list)], title=label, xlabel=r'time / ($\it{s}$)', ylabel='p(t) / $\ln(V(t)*I(t))$')
    ax.plot(t_list, p_list, color='black', linewidth=0.5, label="$\it{p(t)}$ for provided Dataset")

    # Okay. We need the ticks to vary dependent on the values.

    t_base = rounder(max(t_list))
    p_base = rounder(max(p_list) - min(p_list))

    major_ticks_x = pltick.MultipleLocator(t_base/8)
    major_ticks_y = pltick.MultipleLocator(p_base/2)
    minor_ticks_x = pltick.MultipleLocator(t_base/10)
    minor_ticks_y = pltick.MultipleLocator(p_base/5)
    
    ax.xaxis.set_major_locator(major_ticks_x)
    ax.xaxis.set_minor_locator(minor_ticks_x)
    ax.yaxis.set_major_locator(major_ticks_y)
    ax.yaxis.set_minor_locator(minor_ticks_y)

    ax.grid(True, which='major', color="blue")
    ax.grid(True, which='minor', color="pink")
    ax.minorticks_on()

    plt.legend()
    plt.show()
    fig.savefig("hello.png")

# Used to get the name (and file extension) used for this
def name_getter():
    time_delay("Hey! Thanks for running me. Could you give me the filename of the file you want to use? You can optionally specify the extension... :)", 0.005)
    list_string = input()
    print()
    time_delay("Could you also say on which line the datastream begins? I.e. if the first three are text, say that you get data on line four. Etc...:", 0.005)
    try:
        line_start = int(input())
    except:
        time_delay("Needs to be an integer, buster brown.", 0.005)
        quit()
    typeo = []
    extension = []
    new_list = list_string.split(".")
    if len(new_list) == int(1):
        return new_list[0], line_start
    elif len(new_list) == int(2):
        return new_list[0], new_list[1], line_start
    else:
        time_delay("Something's fishy.", 0.005)
        quit()

# Execute to get names/extensions
hey = name_getter() ##################
# Try and except for the name, and if no extension exists, defaults to .txt as the task was specified. 
try:
    new_file = file_former(str(hey[0]), str(hey[1]), hey[2])
    logged_p = log_function(new_file[0], new_file[1])
    grapher(logged_p, new_file[2])
except:
    new_file = file_former(hey[0], "txt", hey[1])
    logged_p = log_function(new_file[0], new_file[1])
    grapher(logged_p, new_file[2])
