import copy
import math
import pickle
import re
import sys
import time
import pandas
from astropy.table import Table
from matplotlib import rc, cm
from matplotlib import animation as anim
import os
from matplotlib.patches import Patch
from numba import njit, jit, types, int32
from numba.experimental import jitclass
import fast_gol
import hdfutils
import numpy as np
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Same conventions as twod_ising, for the most part.
"""
Note that here "not_up" allows for various patterns, either from file or from pre-defined set
There are a few old bits from twod_ising, like measurements/equilibration: rest assured that sims end when absorbing.
"""
class twod_gol(object):
    # TODO: make sure the current init is fine. Remember we added input instead of two separate inits for multirun.
    def __init__(self, lx=None, binrate=None,
                 equilibration=None, measurements=None, not_up=True, identifier=None, checkpoint=False):
        # This just lets you avoid dealing with __init__ if you want to run a method that doesn't need all this stuff.
        if lx!=None:
            # Initialization parameters.
            self.not_up = not_up
            self.lx = lx
            # When running this myself, I don't want to specify this, so this is a dirty fix for the checkpoint.
            if binrate == None:
                self.binrate = 1e6  # 1e6 # binrate for animations. Every i'th sweep is selected as frames.
            else:
                self.binrate = binrate
            self.equilibration = equilibration # In number of flips.
            self.measurements = measurements # In number of flips.
            self.extras = 100 # some extra sweeps to run after reaching equilibrium state (to prove equilibrium.)
            self.autocorrelation_length = 25000 # In number of flips.
            # Derived self parameters.
            self.rng = np.random.default_rng()
            self.A = None
            self.all_A = None
            self.absorbing = False
            self.checkpoint = checkpoint
            # If you want totally random, set not_up to True
            if not_up == True:
                self.mat = self.rng.choice([False, True], size=(self.lx, self.lx))

            # Else, create an empty array, then put in the pattern (or load it) that you selected in not_up
            else:
                # Empty slate
                self.mat = np.zeros(shape=(self.lx, self.lx))
                try:
                    # Load pattern from file
                    if not_up.split(".")[1] == "txt":
                        self.plaintext(not_up)
                # Or load from the predefined presets I've slapped in.
                except:
                    self.pattern(not_up)
            self.sweep = 0 # in NUMBER OF FLIPS
            self.max_sweeps = self.equilibration + self.measurements + 1
            self.delay_max = 5e-3
            self.directory = os.getcwd()
            self.filename = "datafile.hdf5"
            self.set = "data"
            self.numset = "results"
            self.distinct = ("{0}_{1}_{2}").format(self.lx, self.max_sweeps, identifier)
            if checkpoint==False:
                if __name__ == "__main__":
                    self.writer = hdfutils.hdf5_writer(self.directory, self.filename)
                imgdir = self.directory + "\\" + "img"
                try:
                    os.mkdir(imgdir)
                except:
                    pass
                self.imgdir = imgdir + "\\" + self.distinct
                try:
                    os.mkdir(self.imgdir)
                except:
                    pass
                self.animation_name = self.imgdir + "\\animation.mp4" # saving.
            if __name__ == "__main__":
                self.draw = True # set to True to produce a live animation.
                self.framerate = 15  # for animated mp4 saved to file, if applicable. Set to 0 for no animation.
            else:
                self.draw = False # set to True to produce a live animation.
                self.framerate = 0  # for animated mp4 saved to file, if applicable. Set to 0 for no animation.
            self.saveattempts = 5
            self.fast = fast_gol.fast_gol(lx)
            self.all_I = None
        else:
            pass
    def run_high(self):
        # Calculate I for initial step
        start = np.sum(self.mat)
        self.A = start

        # Custom colormap (blue and red.)
        cmap = cm.get_cmap('bwr', 2)

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(self.mat, animated=True, cmap=cmap, aspect='equal', vmin=0, vmax=1)
        ax.set(xlim=[0, self.lx - 1],
               ylim=[0, self.lx - 1])
        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='Alive', edgecolor='black'),
                           Patch(facecolor='blue', label='Dead', edgecolor='black')]
        ax.legend(handles=legend_elements, loc='upper right')

        # Set title.
        ax.set_title(("A = {0:.1f}").format(self.A))
        t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)

        # Run sim
        if self.not_up == True:
            while self.sweep < self.max_sweeps:
                self.fast_gol()
                # Plot every binrate (normally.)
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    ax.set_title(("A = {0:.1f}").format(self.A))
                    t1.set_text(str(self.sweep))
                    # Preliminary Save
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)

    # Load in a Plaintext following https://conwaylife.com/wiki/Plaintext guidelines. Filename WITH .txt.
    # Supports trailing zeros and empty lines.
    def plaintext(self, filename):
        # Load the pattern raw
        lines = []
        with open(os.getcwd() + "\\" + "plaintext" + "\\" + filename, 'rb') as f:
            for line in f:
                line = line.decode("utf-8")
                newline = line.replace("\n","")
                newline = newline.replace("O","1")
                newline = newline.replace(".","0")
                newline = newline.strip()
                newline = list(newline)
                newline = [int(float(d)) for d in newline]
                lines.append(newline)
        # Check the max length and fill in empty lines/trailing zeros.
        maxlength = max([len(d) for d in lines])
        for num, line in enumerate(lines):
            length = len(line)
            zeros_to_add = maxlength - length
            zeros_to_add = [0 for d in range(zeros_to_add)]
            lines[num] += zeros_to_add

        lines = np.array(lines)

        # Get the shape of the pattern
        patshape = np.shape(lines)

        # Verify it's smaller than the matrix
        if np.max(patshape) > self.lx:
            self.time_delay("Your pattern is larger than your Life array. Please increase lx, or more errors will soon come.")
            time.sleep(30)
            self.time_delay("Your loss. Here come the errors.")

        # Get the midpoint of the main array
        mid = int(self.lx/2)

        # Use Numpy Slicing to embed the pattern
        self.mat[mid - int(patshape[0]/2): mid - int(patshape[0]/2) +patshape[0],
        mid - int(patshape[1]/2): mid+patshape[1]-int(patshape[1]/2)] = lines

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        if __name__ == "__main__":
            print()
            for c in text:
                sys.stdout.write(c)
                sys.stdout.flush()
                resttime = self.rng.uniform(0.001, self.delay_max)
                time.sleep(resttime)
        else:
            pass

    # Run the simulator (this is for the actual checkpoint rather than run_high for fun.)
    """
    Given that this simulation is very small, saving everything to file is perfectly reasonable.
    I've removed save functionality from run_high (made for 200, 400x matrices, with 10^6 sweeps- too big.) 
    """
    def run(self):
        # Calculate I for initial step
        start = np.sum(self.mat)
        self.A = start
        all_A = [start]

        # Custom colormap (blue and red.)
        cmap = cm.get_cmap('bwr', 2)

        # Interactive On
        plt.ion()

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8,8))
        im = ax.imshow(self.mat, animated=True, cmap=cmap, aspect='equal', vmin=0, vmax=1)
        ax.set(xlim = [0, self.lx - 1],
               ylim = [0, self.lx - 1])
        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='Alive', edgecolor='black'),
                           Patch(facecolor='blue', label='Dead', edgecolor='black')]
        ax.legend(handles=legend_elements, loc='upper right')

        # Set title.
        ax.set_title(("A = {0:.1f}").format(self.A))
        t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

        # Run Simulation.
        start = time.time()
        mpl = 0
        self.time_delay("Starting Simulation..." + self.distinct)

        # One loop specific for dealing with when we're after checking for "absorbing states"
        if self.not_up == True:
            # Run it until we hit an absorbing state.
            while self.sweep < self.max_sweeps:
                self.fast_gol()
                all_A.append(self.A)
                # Do a check on the previous few states to see if A has been constant (i.e. stdev is zero.)
                # If A has been constant for the last few states, then it's absorbing.
                if int(np.std(all_A[-11:-1])) == 0 and self.sweep > 10:
                    equilibrium_sweep = self.sweep - 10
                    mins = time.time()
                    im.set_array(self.mat)
                    ax.set_title(("Equilibrium state reached, sweep {0:.0f} with A = {1:.1f}.").format(equilibrium_sweep,
                                                                                                       self.A))
                    t1.set_text(self.sweep)
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    ens = time.time()
                    mpl += (ens-mins)
                    break
                # Plot every binrate (normally.)
                if self.sweep % self.binrate == 0:
                    mins = time.time()
                    im.set_array(self.mat)
                    ax.set_title(("A = {0:.1f}").format(self.A))
                    t1.set_text(str(self.sweep))
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    ens = time.time()
                    mpl += (ens-mins)

            # Run it a bit more, to show it truly is "equilibrated." If not, throw a warning.
            self.time_delay(("Running {} extra sweeps to prove we're equilibrated.").format(self.extras))
            for i in range(self.extras):
                self.fast_gol()
                all_A.append(self.A)
                if self.absorbing != True:
                    print("Error occurred and we stopped at a non-absorbing state previously. Balls.")
                # Check to see if an absorbing state has been reached: terminate and do a special plot.
                mins = time.time()
                im.set_array(self.mat)
                ax.set_title(("Equilibrium state reached, sweep {0:.0f}.").format(equilibrium_sweep))
                t1.set_text(self.sweep)
                fig.canvas.draw()
                fig.canvas.flush_events()
                ens = time.time()
                mpl += (ens - mins)
        # One loop specific for dealing with generic patterns (and not worrying about absorbing states and the like.
        else:
            while self.sweep < self.max_sweeps:
                self.fast_gol()
                all_A.append(self.A)
                # Plot every binrate (normally.)
                if self.sweep % self.binrate == 0:
                    mins = time.time()
                    im.set_array(self.mat)
                    ax.set_title(("A = {0:.1f}").format(self.A))
                    t1.set_text(str(self.sweep))
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    ens = time.time()
                    mpl += (ens - mins)
        # All done.
        end = time.time()
        self.time_delay(("Simulation finished. Took {0:.1f} seconds. Plotting stole approximately {1:.1f} seconds.").format(end - start, mpl))
        plt.close()

    # This is for the glider (specifically, tracking the glider, returning metrics, etc.)
    def run_glider(self):
        # just like run, but plots the gliders path along the way.
        if self.checkpoint == True:
            # Calculate I for initial step
            start = np.sum(self.mat)
            self.A = start
            all_A = [start]

            # Get centre of mass for the glider at start
            com = self.fast.glide_com(self.mat)
            all_com = [com]

            # Custom colormap (blue and red.)
            cmap = cm.get_cmap('bwr', 2)

            # Interactive On
            plt.ion()

            # Set up figure, axes, etc
            fig, ax = plt.subplots(figsize=(8, 8))
            im = ax.imshow(self.mat, animated=True, cmap=cmap, aspect='equal', vmin=0, vmax=1)
            ax.set(xlim=[0, self.lx - 1],
                   ylim=[0, self.lx - 1])
            # Set up image plot
            legend_elements = [Patch(facecolor='red', label='Alive', edgecolor='black'),
                               Patch(facecolor='blue', label='Dead', edgecolor='black')]
            ax.legend(handles=legend_elements, loc='upper right')

            # Set title.
            ax.set_title(("A = {0:.1f}").format(self.A))
            t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

            # Also draw on the path of the glider as a scatter plot
            scat = plt.scatter(*zip(*all_com), color='white', marker='x')

            # Run Simulation.
            start = time.time()
            mpl = 0
            self.time_delay("Starting Simulation..." + self.distinct)

            # Run it until we hit an absorbing state.
            while self.sweep < self.max_sweeps:
                self.fast_gol()
                all_A.append(self.A)
                # Get centre of mass for the glider
                com = self.fast.glide_com(self.mat)
                all_com.append(com)
                # Plot every binrate (normally.)
                if self.sweep % self.binrate == 0:
                    mins = time.time()
                    im.set_array(self.mat)
                    scat.set_offsets(all_com)
                    ax.set_title(("A = {0:.1f}").format(self.A))
                    t1.set_text(str(self.sweep))
                    fig.canvas.draw()
                    fig.canvas.flush_events()
                    ens = time.time()
                    mpl += (ens - mins)

            # All done.
            end = time.time()
            self.time_delay(
                ("Simulation finished. Took {0:.1f} seconds. Plotting stole approximately {1:.1f} seconds.").format(
                    end - start, mpl))
            plt.close()
        # designed to produce a graph, specifically, to show glider position and estimate the speed.
        else:
            # Calculate I for initial step
            start = np.sum(self.mat)
            self.A = start
            all_A = [start]

            # Get centre of mass for the glider at start
            com = self.fast.glide_com(self.mat)
            all_com = [com]

            # Run Simulation.
            while self.sweep < self.max_sweeps:
                self.fast_gol()
                all_A.append(self.A)
                # Get centre of mass for the glider
                com = self.fast.glide_com(self.mat)
                all_com.append(com)

            # Create a scatter plot
            fig = plt.figure(figsize=(8,8))
            plt.scatter(*zip(*all_com), color='red', marker='x')
            plt.grid(which='major', color='pink')

            # Get x/y displacements for speedsteps. Use a grid size of 200 for a good result with speedsteps of 100.
            speedsteps = 100
            x_disp, y_disp = all_com[speedsteps] - all_com[0]
            x_mag, y_mag = np.abs(x_disp), np.abs(y_disp)
            print(x_mag, y_mag)

            # Also ascertain the speed according to https://conwaylife.com/wiki/Speed (i.e. in x or y terms.)
            speed = np.max([x_mag,y_mag])/speedsteps

            # Set the title up
            plt.title(("The speed of this glider was {0:.2f}c (see https://conwaylife.com/wiki/Speed definition.)").format(speed))
            plt.savefig("gliderspeed.png", dpi=300)
            plt.show()

    # For multiprocessed version: no thrills attached.
    def run_multi(self):
        # Calculate I for initial step
        self.A = np.sum(self.mat)
        all_A = [self.A]

        # Run Simulation.
        self.time_delay("Starting Simulation..." + self.distinct)

        # Run it.
        while self.sweep < self.max_sweeps:
            self.fast_gol()
            all_A.append(self.A)
            # Check to see if an equilibrium state has been reached: terminate and do a special plot.
            if int(np.std(all_A[-11:-1])) == 0 and self.sweep > 10:
                break
        # Create datasave format (ignoring "non-equilibrium states.")
        self.all_A = np.array(all_A[:-11])

        # Save it
        self.save()

    # Iterate Sequential (i.e. monte-carlo the points one at a time.)
    def fast_gol(self):
        self.mat, self.A, self.sweep = self.fast.fast_gol(self.mat, self.sweep) # , self.absorbing

    # Some loops to generate "Examples of patterns."
    # See Wikipedia https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life
    def pattern(self, which):
        mid = int(self.lx/2)
        if which == "block":
            self.mat[mid, mid], \
            self.mat[mid + 1, mid], \
            self.mat[mid, mid + 1], \
            self.mat[mid + 1, mid + 1] = True, \
                                         True, \
                                         True, \
                                         True
        if which == "beehive":
            self.mat[mid,mid], self.mat[mid,mid+1], self.mat[mid+1,mid-1], \
            self.mat[mid+1,mid+2], self.mat[mid+2, mid], self.mat[mid+2, mid+1] = True, True, True, \
                                                                                  True, True, True
        if which == "blinker":
            self.mat[mid, mid], self.mat[mid+1, mid], self.mat[mid+2, mid] = True, True, True
        if which == "pentalon":
            column = np.arange(mid, mid+3, 1)
            row = np.arange(mid, mid+8, 1)
            for i in row:
                for j in column:
                    self.mat[i,j] = True
            self.mat[mid+1, mid+1], self.mat[mid+6, mid+1] = False, False
        if which == "glider":
            self.mat[mid,mid], self.mat[mid, mid+1], self.mat[mid, mid+2], \
            self.mat[mid-1, mid+2], self.mat[mid-2, mid+1] = True, True, \
                                                            True, True, True
        if which == "rentomino":
            self.mat[mid,mid], self.mat[mid-1,mid], self.mat[mid+1,mid], \
            self.mat[mid,mid-1], self.mat[mid-1,mid+1] = True, True, \
                                                         True, True, True
        if which == 'queenbee':
            self.mat[mid,mid], self.mat[mid+1, mid], self.mat[mid+5,mid], \
            self.mat[mid+6, mid], self.mat[mid+2, mid+1], self.mat[mid+3, mid+1], \
            self.mat[mid+4, mid+1], self.mat[mid+1, mid+2], self.mat[mid+5, mid+2], \
            self.mat[mid+2, mid+3], self.mat[mid+4, mid+3], self.mat[mid+3, mid+4] = True, True, True, \
                                                                                     True, True, True, \
                                                                                     True, True, True, \
                                                                                     True, True, True

    # Save data.
    """
    Note on the form of all_A
    The zeroth element corresponds to the initial condition, i.e. SWEEP ZERO
    Nice and pythonic, right? 
    """
    def save(self, table=None):
        if table == None:
            i = 0
            while i <= self.saveattempts:
                try:
                    # Dump it
                    with open(self.imgdir + "\\" + "data.txt", 'wb') as f:
                        pickle.dump(obj=np.array(self.all_A), file=f)
                    break
                except:
                    i += 1
                    time.sleep(self.rng.integers(0, 30, 1)[0])
                    pass
        else:
            i = 0
            while i <= self.saveattempts:
                try:
                    self.writer.write_table(self.distinct, self.numset, table)
                    break
                except:
                    i += 1
                    time.sleep(self.rng.integers(0, 30, 1)[0])
                    pass

    # Load data.
    def load(self, table=None):
        try:
            if table == None:
                # Load it
                with open(self.imgdir + "\\" + "data.txt", 'rb') as f:
                    self.all_A = pickle.load(f)
            else:
                return self.writer.read_table(self.distinct, self.numset)
            self.time_delay("Loaded successfully.")
        except:
            self.time_delay("Load failed.")

    # Graphing tool for part (2) to generate the distribution for reaching absorbing states
    def multigraph_absorbing(self):
        # Set "identifiers" and get "distincts"
        lx = 50
        nrun=10000
        equilibration = int(0)
        measurementss = int(20e3)
        max_sweepssss = equilibration + measurementss + 1
        identifiers = np.arange(0, nrun, 1)
        distincts = [("{0}_{1}_{2}").format(lx, max_sweepssss, identifier) for identifier in identifiers]
        imgdir = os.getcwd() + "\\" + "img"
        distinct_imgdirs = [imgdir + "\\" + distinct for distinct in distincts]
        all_A_list = []
        for distinct_imgdir in distinct_imgdirs:
            with open(distinct_imgdir + "\\data.txt", 'rb') as f:
                all_A_list.append(pickle.load(f))
        # Get lengths of all_A (and hence the maximum sweeps before absorbing states are reached.)
        max_sweeps = [len(all_A) for all_A in all_A_list]
        # Create Histogram Data
        nbins = 30
        range = (0, 5000)
        binned = np.histogram(max_sweeps, bins=nbins, range=range, density=True)
        binned *= 100 # in percent and not in pure probability
        # Produce plot
        fig = plt.figure()
        plt.plot(np.linspace(range[0], range[1], nbins), 100*binned[0], color='red')
        plt.grid(which='major', color='pink')
        plt.xlabel("Final Sweep")
        plt.ylabel("Probability \%")
        plt.savefig("gol_probability.png", dpi=300)
        plt.show()

# Class for handling user input (i.e. checkpoint.)
class checkpoint(object):
    def __init__(self):
        self.delay_max = 2e-3
        self.rng = np.random.default_rng()

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text, nospacing=False):
        if __name__ == "__main__":
            if nospacing==False:
                print()
                for c in text:
                    sys.stdout.write(c)
                    sys.stdout.flush()
                    resttime = self.rng.uniform(0.0001, self.delay_max)
                    time.sleep(resttime)
                print()
            else:
                for c in text:
                    sys.stdout.write(c)
                    sys.stdout.flush()
                    resttime = self.rng.uniform(0.0001, self.delay_max)
                    time.sleep(resttime)
        else:
            pass

    # User Input Specification
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to the Game of Life. \n"
                        "You will now be asked for a few parameters. Please give them. \n"
                        "Please note that the code is a bit optimized for multiple sequential runs (i.e. parallel) \n"
                        "Due to this, a one-time-run will incur a @jit compile cost, compared to regular python. \n"
                        "Now, onto the parameters!!!!")
        print()
        self.time_delay("Grid size. This is a square simulator, so just one integer will suffice.")
        lx = int(input())
        self.time_delay("Do you want a random state? y/n.")
        randyesno = str(input())
        if randyesno == "y":
            not_up = True
        elif randyesno == 'n':
            # Get the list of plaintexts
            files = os.listdir(os.getcwd() + "\\" + "plaintext")
            stringformat = ["- " + file + "\n" for file in files]
            self.time_delay("You don't want it random. Your options are:\n"
                            "- block \n"
                            "- beehive \n"
                            "- blinker \n"
                            "- pentalon \n"
                            "- rentomino \n"
                            "- queenbee \n"
                            "\n"
                            "You also have the option (drawn from plaintext) of: \n")
            for string in stringformat:
                self.time_delay(string, nospacing=True)
            self.time_delay("If selecting a plaintext, please include the file extension. \n"
                            "Note that some patterns may larger grids than the status quo of 50x50.")
            not_up = str(input())
        if not_up == True:
            self.time_delay("Since you selected a random state, we will run until we hit an absorbing state.")
            measure = int(1e6)
        else:
            self.time_delay("You selected a pattern- how long do you want to run the simulation for? Give an integer.")
            measure = int(input())
        equi = 0
        return lx, equi, measure, not_up

    # Run
    def run(self):
        try:
            lx, equi, measure, not_up = self.user_input()
            model = twod_gol(lx=lx, binrate=1, equilibration=equi, measurements=measure, not_up=not_up, identifier=0,
                             checkpoint=True)
            model.run()
        except Exception as e:
            self.time_delay("An error occurred... \n"
                            "You probably don't have the correct dependencies. \n"
                            "If the error regards the existence of LaTex: delete lines 21,22 \n" 
                            "If the error regards missing packages, please install them.\n")
            print(e)

# Run the checkpoint/etc.
check = checkpoint()
check.run()