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
import fast_sirs
import hdfutils
import numpy as np
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Same conventions as twod_ising, for the most part.
class twod_sirs(object):
    # TODO: make sure the current init is fine. Remember we added input instead of two separate inits for multirun.
    def __init__(self, lx=None, p1=None, p2=None, p3=None, binrate=None,
                 equilibration=None, measurements=None, not_up=True):
        # This just lets you avoid dealing with __init__ if you want to run a method that doesn't need all this stuff.
        if lx!=None:
            # Initialization parameters.
            self.lx = lx
            self.p = np.array([p1, p2, p3])
            # When running this myself, I don't want to specify this, so this is a dirty fix for the checkpoint.
            if binrate == None:
                self.binrate = 1e6  # 1e6 # binrate for animations. Every i'th sweep is selected as frames.
            else:
                self.binrate = binrate
            self.equilibration = equilibration # In number of flips.
            self.measurements = measurements # In number of flips.
            self.autocorrelation_length = 25000 # In number of flips.
            # Derived self parameters.
            self.rng = np.random.default_rng()
            self.I = None
            if not_up == True:
                self.mat = self.rng.choice([False, True, 2], size=(self.lx, self.lx))  # Random S I or R matrix.
            else:
                self.mat = np.ones(shape=(self.lx, self.lx)) # all-up spin matrix (all infected.)
            self.sweep = 0 # in NUMBER OF FLIPS
            self.max_sweeps = self.equilibration + self.measurements + 1
            self.delay_max = 5e-3
            self.directory = os.getcwd()
            self.filename = "datafile.hdf5"
            self.set = "data"
            self.numset = "results"
            if __name__ == "__main__":
                self.writer = hdfutils.hdf5_writer(self.directory, self.filename)
            imgdir = self.directory + "\\" + "img"
            try:
                os.mkdir(imgdir)
            except:
                pass
            self.distinct = ("{0}_{1:.4f}_{2:.4f}_{3:.4f}_{4}").format(self.lx, p1, p2, p3, self.max_sweeps)
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
            self.fast = fast_sirs.fast_sirs(lx, p1, p2, p3)
            self.all_I = None
        else:
            pass

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
    def run(self, checkpoint=False):
        # Calculate I for initial step
        self.I = self.fast.fast_infsum(self.mat)
        all_I = [self.I]

        # Custom colormap (blue and red.)
        cmap = cm.get_cmap('bwr', 3)

        # Interactive On
        plt.ion()

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8,8))
        im = ax.imshow(self.mat, animated=True, cmap=cmap, aspect='equal', vmin=0, vmax=2)
        ax.set(xlim = [0, self.lx - 1],
               ylim = [0, self.lx - 1])
        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='R', edgecolor='black'),
                           Patch(facecolor='white', label='I', edgecolor='black'),
                           Patch(facecolor='blue', label='S', edgecolor='black')]
        ax.legend(handles=legend_elements, loc='upper right')

        # Set title.
        ax.set_title(("I = {0:.1f}").format(self.I))
        t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

        # Run Simulation.
        start = time.time()
        mpl = 0
        self.time_delay("Starting Simulation..." + self.distinct)
        while self.sweep < self.max_sweeps:
            self.fast_sequential()
            all_I.append(self.I)
            if self.sweep % self.binrate == 0:
                mins = time.time()
                im.set_array(self.mat)
                ax.set_title(("I = {0:.1f}").format(self.I))
                t1.set_text(str(self.sweep))
                fig.canvas.draw()
                fig.canvas.flush_events()
                ens = time.time()
                mpl += (ens-mins)
        end = time.time()

        # All done.
        self.time_delay(("Simulation finished. Took {0:.1f} seconds. Plotting stole approximately {1:.1f} seconds.").format(end - start, mpl))
        plt.close()

        # Create datasave format
        self.all_I = all_I

        if checkpoint==False:
            # Save it
            self.save()

    # Iterate Sequential (i.e. monte-carlo the points one at a time.)
    def fast_sequential(self):
        self.mat, self.I, self.sweep = self.fast.fast_sequential(self.mat, self.I, self.sweep)

    # Iterate parallel (i.e. sweep entire lattice at once.)
    def fast_parallel(self):
        self.mat, self.I, self.sweep = self.fast.fast_parallel(self.mat, self.I, self.sweep)

    # Saving files.
    def save(self, table=None):
        if table == None:
            i = 0
            while i <= self.saveattempts:
                try:
                    # Dump it
                    with open(self.imgdir + "\\" + "data.txt", 'wb') as f:
                        pickle.dump(obj=self.all_I, file=f)
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

    # Loading files.
    def load(self, table=None):
        try:
            if table == None:
                # Load it
                with open(self.imgdir + "\\" + "data.txt", 'rb') as f:
                    self.all_I = pickle.load(f)
            else:
                return self.writer.read_table(self.distinct, self.numset)
            self.time_delay("Loaded successfully.")
        except:
            self.time_delay("Load failed.")

    # For multiprocessed version: no thrills attached.
    def run_multi(self):
        # Calculate I for initial step
        self.I = self.fast.fast_infsum(self.mat)
        all_I = [self.I]

        # Run Simulation.
        while self.sweep < self.max_sweeps:
            self.fast_sequential()
            all_I.append(self.I)

        # Create datasave format
        self.all_I = np.array(all_I)
        self.save()

    # To generate runs for graphs/etc. Set run=False to re-generate averages but not do anything else.
    def main_multi(self, run=True):
        # To run the simulation
        if run == True:
            # Check if data already exists
            if os.path.isfile(self.imgdir + "\\" + "data.txt") == True:
                pass
            else:
                start = time.time()
                self.run_multi()
                mid = time.time()
                self.fast_averages_errors()
                end = time.time()
                print(("Simulation took (including saving) {0:.1f} and averages/errors took {1:.1f}").format(mid-start, end-mid))
        # Re-generate averages with new parameters you've added to multigraph
        else:
            self.load()
            self.fast_averages_errors()

    # Averages/errors
    def fast_averages_errors(self):
        # Calculate averages and errors
        avg_I, \
        avg_I_err, \
        chi_true, \
        chi_error = self.fast.averages_errors(self.all_I,
                                              self.equilibration,
                                              self.measurements,
                                              self.autocorrelation_length)
        # Only save if not multiprocessing
        if __name__ == "__main__":
            # Produce a table and save it.
            table = Table()
            table["I"], table["chi"] = np.array([avg_I, avg_I_err]), np.array([chi_true, chi_error])
            self.save(table)
        else:
            array_dump = np.array([avg_I, avg_I_err, chi_true, chi_error, self.p[0], self.p[1], self.p[2]])
            with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
                pickle.dump(obj=array_dump, file=f)

    # Produce graphs for multiprocessed runs selected.
    def multigraph(self):
        # Sim Params.
        lx = 50
        equilibration = int(2.5e6)
        measurements = int(25e6)
        max_sweeps = equilibration + measurements + 1

        # Generate the parameter space for the runs.
        one_range = np.linspace(0, 1, 50)
        two_range = np.linspace(0, 1, 50)
        tre_range = np.linspace(0, 1, 50)
        zipped = list(zip(one_range, two_range, tre_range))
        distincts = [("{0}_{1:.4f}_{2:.4f}_{3:.4f}_{4}").format(self.lx,
                                                                p1,
                                                                p2,
                                                                p3,
                                                                self.max_sweeps) for p1, p2, p3 in zipped]
        imgdirs = [os.getcwd() + "\\" + "img" + "\\" + distinct for distinct in distincts]
        files = [imgdir + "\\" + "results.txt" for imgdir in imgdirs]
        arrays = [] # avg_I, avg_I_err, chi, chi_err, p1, p2, p3 (chi is susceptibility of I. I'm a lazy sod!)
        # Grab dumps and make table.
        for file in files:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass
        arrays=np.array(arrays)
        labels = ['avg_I', 'avg_I_err', 'chi_true', 'chi_error', 'p1', 'p2', 'p3']
        table = Table(arrays, names=labels)

        # Select which planes to use as "x" or "y" and which to split along the "z"
        x_plane, y_plane, z_plane = 'p1', 'p2', 'p3'

# Class for handling user input (i.e. checkpoint.)
class checkpoint(object):
    def __init__(self):
        self.delay_max = 2e-3
        self.rng = np.random.default_rng()

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        if __name__ == "__main__":
            print()
            for c in text:
                sys.stdout.write(c)
                sys.stdout.flush()
                resttime = self.rng.uniform(0.0001, self.delay_max)
                time.sleep(resttime)
            print()
        else:
            pass

    # User Input Specification
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to this 2D SIRS Simulation Suite. \n"
                        "You will now be asked for a few parameters. Please give them. \n"
                        "Please note that the code is a bit optimized for multiple sequential runs (i.e. parallel) \n"
                        "Due to this, a one-time-run will incur a @jit compile cost, compared to regular python. \n"
                        "Now, onto the parameters!!!!")
        print()
        self.time_delay("Grid size. This is a square simulator, so just one integer will suffice.")
        lx = int(input())
        print()
        self.time_delay("Probability of S to I (p1).")
        p1 = int(input())
        print()
        self.time_delay("Probability of I to R (p2).")
        p2 = int(input())
        print()
        self.time_delay("Probability of R to S (p3).")
        p3 = int(input())
        print()
        not_up = True
        self.time_delay("What binrate do you want for the sim? (spacing between animation prints. Code bottleneck.)")
        binrate = int(input())
        print()
        self.time_delay("We will use the default equilibration and measurement periods for this run. Enjoy!")
        equi, measure = int(0.25e6), int(25e6)
        return lx, equi, measure, not_up, binrate, p1, p2, p3

    # Run
    def run(self):
        try:
            lx, equi, measure, not_up, binrate, p1, p2, p3 = self.user_input()
            model = twod_sirs(lx, p1, p2, p3, binrate, equi, measure, not_up)
            model.run(checkpoint=True)
        except Exception as e:
            self.time_delay("An error occurred... \n"
                            "You probably don't have the correct dependencies. \n"
                            "If the error regards the existence of LaTex: delete lines 21,22 \n" 
                            "If the error regards missing packages, please install them.\n")
            print(e)


uwu = twod_sirs(lx=50, p1=0.3, p2=0.4, p3=0.5, binrate=5000, equilibration=1000000, measurements=10000000, not_up=True)
uwu.main_multi(run=True)