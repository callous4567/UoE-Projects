import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.axes3d import Axes3D # Redundant here but will be used in future and thus is kept for my use.
import matplotlib.ticker as pltick
import sys
import time
import random

"""
QUICK README
------------
I was unsure as to what was requested... at first I wondered if you wanted a 2D array simulation, so I did that. 
Then I realised "3D...?" so I did that...
Then I realised that none of that was even wanted and that algorithm was just meant to be written down... so...

ALGORITHM:
- Take length scale of cube
- Take position
- Calculate scale of position against length of cube and floor it for a scale
- Subtract (scale * length scale)  vector from position vector
- This gives an image inside the original space 

Cube = (1,1,1)
Position: (2,2.5,3.5) Scale: (2, 2, 3) -> (0,0.5,0.5)
Position: (2,2.5,-3.5) Scale: (2, 2, -4) -> (0,0.5,0.5)

The 2D/3D simulation is there to show the periodic boundary condition of a particle flying out/coming back in...
Hopefully it counts for some sort of grade on this or something like that since it took so long to learn how to do... >.<

METHODS REQUESTED FOR IMAGE AND "CLOSEST" IMAGE:
------------------------------------------------
They're wrapped in with the particle class... the image_closest and image functions both run via the console and print.
That should count as a test?
"""
# print() except it delays release of characters individually and makes the script look smoother on display.
# Change arguments for numbertime to 0,0 if you intend to do this vector stuff for decimal vectors instead of integer vectors... printing can take a while otherwise.
def time_delay(text):
    print()
    text = str(text)
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.005, 0.01)
        time.sleep(numbertime)
    print()

# Hard-code box parameters. Pos_max is size/lengthscale, vel_max is particle velocity, step is timestep and number is number of particles.
pos_max = 10000
vel_max = 500
step = 0.1
number = 1

# To change size and such for the sake of "custom l" condition
time_delay("Please provide box size, integers only.")
try:
    pos_max = int(input())
except:
    while True:
        print("Should have given an integer.")

# Defines a particle via: hey = particle([x,y,z], [vx, vy, vz])
# particle.iterate() iterates by a hard-defined timestep, assuming no acceleration (for simulation)
# particle.image() provides an image within the cube (pos_max, pos_max, pos_max)
# particle.distance() provides distance to origin 0,0,0
# particle.image_closest() returns the closest image to the origin (by distance)
# URGENT NOTE WORTH READING --->>> all objects returned are new particle objects.
# Old object editing is not done due to issues in the graphing iterator (object NoneType error after exactly three iterations.)
class particle():
    def __init__(self, xyz, vel):
        self.position = np.array(xyz)
        self.velocity = np.array(vel)
        self.distance = np.sqrt(np.inner(self.position, self.position))
    # Images particle to produce new one. Functionality is that: new_particle = particle.image()
    # For each vector component, takes ratio of component to pos_max (length of cube) and floors this value...
    # Once floored, it subtracts the product of this "scale factor" and the pos_max, forcing an image in the cube.
    def image(self):
        xyz = self.position
        newxyz = []
        for d in xyz:
            scale = d/pos_max
            image_coord = d - np.floor(scale)*pos_max
            newxyz.append(image_coord)

        return particle(newxyz, self.velocity)
    # Finds the image closest to the origin and provides a new particle object.
    def image_closest(self):
        image_array = np.ndarray(shape=(0,8), dtype=object)
        image_base = self.image()
        xyz_base = image_base.position
        array_of_arrays = [xyz_base]
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    new_xyz = xyz_base - np.array([i*pos_max, j*pos_max, k*pos_max])
                    print(new_xyz) # For debug purposes/you guys want to see that the final image choice is correct in console?
                    array_of_arrays.append(new_xyz)
        length_array = [np.sqrt(np.inner(d,d)) for d in array_of_arrays]
        mini = min(length_array)
        minidex = length_array.index(mini)

        return particle(array_of_arrays[minidex], self.velocity)
    # Quick iterate step based on velocity and timestep. New pos = pos_0 + v*timestep.
    # Note, depreciated use of "if/for/etc" loops instead of using the added "image" method.
    # Works fine for simulation, so it's been kept. I'm lazy :-:
    def iterate(self):
        index = int(0)
        new_particle = np.array([0,0,0])
        while True:
            new_pos_comp = self.position[index] + self.velocity[index]*step
            if new_pos_comp > pos_max:
                new_pos_comp -= pos_max
            if new_pos_comp < 0:
                new_pos_comp += pos_max
            new_particle[index] = (new_pos_comp)
            index += int(1)
            if index == int(3):
                break
        return particle(new_particle, self.velocity)
# Generates random particle, x between 0 and pos_max, and velocity 0 and +-vel_max. Integer wise.
def randgen():
    test_particle = "False"
    if test_particle == "True":
        global number
        number = int(1)
        norpos = pos_max/2
        norvel = vel_max/8
        return particle([norpos, norpos, norpos],[norvel, 0, 0])
    else:
        z = lambda r: np.random.randint(0, r)
        v = lambda r: np.random.randint(-r, r)
        xyz = pos_max
        vel = vel_max
        return particle([z(xyz), z(xyz), z(xyz)], [v(vel), v(vel), v(vel)])
# Generates "number" random particles in an ndarray, with element 0 being the particles.
def parti_gen():
    pa = np.ndarray(shape=[1,number], dtype=object) # Object array
    for i in range(0, number):
            pa[0,i] = randgen()
    return(pa)
# Iterates particle matrix according to iterate function in particle class
def parti_iterator(matrix):
    hey = np.ndarray(shape=(1,number), dtype=object)
    for i in range(0, number):
        u = matrix[0, i]
        u_new = u.iterate()
        hey[0,i] = u_new
    return hey
# Forms list of XY (and z) values for the particle matrix.
def xyz_getter(pamat):
    x_list = [0]*number
    y_list = [0]*number
    z_list = [0]*number
    for i in range(0, number):
        u = pamat[0, i]
        xyz = u.position
        x_list[i] = (xyz[0])
        y_list[i] = (xyz[1])
        z_list[i] = (xyz[2])
    return x_list, y_list, z_list


# Take x and return image in box and closest to origin. Quick and dirty (since I just added it at the end and really just wanna sleep! >.<
time_delay("Please provide vector coordinates for position to be 'imaged'... provide in x,y,z form without braacketing.")
coords = str(input())
coordsplit = coords.split(',')
try:
    coord = [float(d) for d in coordsplit]
except:
    while True:
        print("Give me the right format!")
coord = np.array(coord)
particle_test = particle(coord, [0,0,0])
p_img = particle_test.image()
p_origin = particle_test.image_closest()
time_delay("Your image in the positive box, and also your closest image to the origin, is:")
print(p_img.position)
print(p_origin.position)




# Here's the script for the particle cell boundary simulation that was in the exercise... or something like it?
# 3D simulation
"""
# Define figure and ax1/ax2 3D/2D representation
fig = plt.figure(figsize=(10,10))
# ax1 for 3D.
ax1 = fig.add_subplot(111, projection='3d')

# Defines particle matrix and also quickly gets the initial state from it using xyz_getter()
pamat = parti_gen()
heyo = xyz_getter(pamat)
graph, = ax1.plot(heyo[0], heyo[1], heyo[2], 'o')

# Sets limits and enables grid for axes
ax1.set_xlim(0, pos_max)
ax1.set_ylim(0, pos_max)
ax1.set_zlim(0, pos_max)
ax1.grid()

# Defines ticks on axes
x_majors = pltick.MultipleLocator(1000)
ax1.xaxis.set_major_locator(x_majors)
ax1.yaxis.set_major_locator(x_majors)
ax1.zaxis.set_major_locator(x_majors)

# Animator 3D updater for FuncAnimation (updates plot data)
def updater(i):
    global pamat
    data = xyz_getter(pamat)
    graph.set_data(data[0], data[1])
    graph.set_3d_properties(data[2])
    pamat = parti_iterator(pamat)
    return graph,
# Animation function from matplotlib. Figure, update function, number of frames in animation, timestep.
# Blit defines visual optimisation. Disabling causes issues on my setup and thus it is enabled by default.
# If plot gives list index error, attempt to disable blit. If index error persists, your setup is to blame...
ani = animation.FuncAnimation(fig, updater, frames=6000, interval=step, blit=True)

plt.show()
"""

# 2D simulation
"""
# 2D Plot on ax2. Aspect scaling is set for square, autoscale fixed to false to prevent nasty scaling glitch.
ax2 = fig.add_subplot(211, aspect='equal', autoscale_on=False, xlim=(0, pos_max), ylim=(0, pos_max))

# Scatter plot in 2D: size of points ~1**2 and red.
ax2.scatter([], [], s=1, c='r')
ax2.grid()

# Create matrix of random particles
pamat = parti_gen()

# Define 2D line for ax2 animations
line, = ax2.plot([], [], 'bo', ms=6)

# Set majors and axis ticks
x_majors = pltick.MultipleLocator(1000)
ax2.xaxis.set_major_locator(x_majors)
ax2.yaxis.set_major_locator(x_majors)


# Initialisation for anim 2D
def init():
    line.set_data([], [])
    return line,

# Animator 2D on xy axis (by default. If you want xz, zy, etc, adjust indices for x_list and y_list.
# Index 0,1,2 are x,y,z from xyz_getter() function.
def animate(i):
    global pamat
    x_list = xyz_getter(pamat)[0]
    y_list = xyz_getter(pamat)[1]
    concat_list = [x_list, y_list]
    line.set_data(concat_list)
    pamat = parti_iterator(pamat)
    return line,
    
# Initiate animation.
ani = animation.FuncAnimation(fig, animate, frames=6000, interval=step, blit=True, init_func=init)

plt.show()
"""

# If you want to save, remove the hash. Added for the sake of my laptop, on which Blit (or no blit) doesn't help and it doesn't work.
# ani.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])



