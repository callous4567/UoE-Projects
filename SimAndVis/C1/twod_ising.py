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
from fast_ising import fast_ising
import hdfutils
import numpy as np
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Important Note regarding comments on Sweeps/Flips.
""" 
As mentioned in the email, I have accidentally misunderstood a sweep as being equal to a flip.
Consequently, this code uses them interchangeably. 
I have, however, altered the user-input definition to represent the conventional definition of sweep
I.e., one sweep = N^2 flips (as mentioned in an email I sent to Davide.)
If that's you (who is marking this) then Hi! :D 
"""
# If you want to see completed animations of this code in action (rather than deal with live-run.)
"""
https://www.youtube.com/watch?v=X4qfKsjoVds 800x800 Kawasaki and Glauber for 160,000,000 Flips T = 0.5->3.5
https://www.youtube.com/watch?v=MBItgHSFCic 200x200 Glauber T = 0.5->3.5
https://www.youtube.com/watch?v=qrdOGegpwr4 200x200 Glauber T = 1
https://www.youtube.com/watch?v=E9T4MoEZstU 2001x2001 Glauber (uses fast_ising instead.) 
"""
class twod_ising(object):
    # TODO: make sure the current init is fine. Remember we added input instead of two separate inits for multirun.
    def __init__(self, lx=None, T=None, dyn=None, equilibration=None, measurements=None, not_up=False, binrate=None):
        # This just lets you avoid dealing with __init__ if you want to run a method that doesn't need all this stuff.
        if lx!=None:
            # Initialization parameters.
            self.lx = lx
            self.T = T
            self.dyn = dyn
            # When running this myself, I don't want to specify this, so this is a dirty fix for the checkpoint.
            if binrate == None:
                self.binrate = 1e6  # 1e6 # binrate for animations. Every i'th sweep is selected as frames.
            else:
                self.binrate = binrate
            self.equilibration = equilibration # int(0.25e6) # int(4000e6) # in NUMBER OF FLIPS
            self.measurements = measurements # int(25e6) # in NUMBER OF FLIPS
            self.autocorrelation_length = 25000
            # Derived self parameters.
            self.rng = np.random.default_rng()
            self.M, self.E = None, None
            self.sweep = 0 # in NUMBER OF FLIPS
            self.max_sweeps = self.equilibration + self.measurements + 1
            self.delay_max = 5e-3
            self.directory = os.getcwd()
            self.filename = "datafile.hdf5"
            self.set = "data"
            self.numset = "results"
            if __name__ == "__main__":
                self.writer = hdfutils.hdf5_writer(self.directory, self.filename)
            if self.dyn == "g": # glauber dynamic
                # When running this myself, I don't want to specify this, so this is a dirty fix for the checkpoint.
                if not_up == True:
                    self.mat = self.rng.choice([-1, 1], size=(self.lx, self.lx))  # random spin matrix
                else:
                    self.mat = np.ones(shape=(self.lx, self.lx)) # all-up spin matrix
            if self.dyn == "k": # kawasaki dynamic
                self.mat = self.rng.choice([-1, 1], size=(self.lx, self.lx)) # random spin matrix
            imgdir = self.directory + "\\" + "img"
            try:
                os.mkdir(imgdir)
            except:
                pass
            self.distinct = ("{0}_{1:.4f}_{2}_{3}").format(self.lx, self.T, self.dyn, self.max_sweeps)
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
            self.temprate = 1000e6
            self.temprise = 0.6
            self.all_M, self.all_E = None, None
            self.saveattempts = 5
            self.fast = fast_ising(self.lx)
        else:
            pass

    # Run the simulator (this is for the actual checkpoint rather than run_high for fun.)
    """
    Given that this simulation is very small, saving everything to file is perfectly reasonable.
    I've removed save functionality from run_high (made for 200, 400x matrices, with 10^6 sweeps- too big.) 
    """
    def run(self, checkpoint=False):
        # Calculate M and E for initial step
        self.M, self.E = self.fast.fast_magenergy(self.mat)
        all_M, all_E = [self.M], [self.E]

        # Custom colormap (blue and red.)
        cmap = cm.get_cmap('bwr', 2)

        # Interactive On
        plt.ion()

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8,8))
        im = ax.imshow(self.mat, animated=True, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
        ax.set(xlim = [0, self.lx - 1],
               ylim = [0, self.lx - 1])
        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                           Patch(facecolor='blue', label='-1', edgecolor='black')]
        ax.legend(handles=legend_elements, loc='upper right')

        # Set title.
        ax.set_title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,
                                                                    self.M,
                                                                    self.E))
        t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

        # Run Simulation.
        self.time_delay("Starting Simulation..." + self.distinct)
        if self.dyn == 'g':
            while self.sweep < self.max_sweeps:
                self.fast_glauber()
                all_M.append(self.M), all_E.append(self.E)
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    ax.set_title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,
                                                                                self.M,
                                                                                self.E))
                    t1.set_text(str(self.sweep))
                    fig.canvas.draw()
                    fig.canvas.flush_events()
        if self.dyn == 'k':
            while self.sweep < self.max_sweeps:
                self.fast_kawasaki()
                all_M.append(self.M), all_E.append(self.E)
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    ax.set_title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,
                                                                                self.M,
                                                                                self.E))
                    t1.set_text(str(self.sweep))
                    fig.canvas.draw()
                    fig.canvas.flush_events()

        # All done.
        plt.close()

        # Create datasave format
        self.all_M, self.all_E = all_M, all_E

        if checkpoint==False:
            # Save it
            self.save()

    # For multiprocessed version: no thrills attached.
    def run_multi(self):
        # Calculate M and E for initial step
        self.M, self.E = self.fast.fast_magenergy(self.mat)
        all_M, all_E = [self.M], [self.E]

        # Run Sim
        if self.dyn == 'g':
            while self.sweep < self.max_sweeps:
                self.fast_glauber()
                all_M.append(self.M), all_E.append(self.E)
        # For Kawasaki scenario
        if self.dyn == 'k':
            while self.sweep < self.max_sweeps:
                self.fast_kawasaki()
                all_M.append(self.M), all_E.append(self.E)

        # Create datasave format
        self.all_M, self.all_E = np.array(all_M), np.array(all_E)

        # Save it
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

    # Produce graphs for multiprocessed runs selected.
    def multigraph(self):
        # Same "settings" as used by Multirun (dynamics are "k" or "g" obviously.)
        equilibration = int(0.25e6) # int(2e3)
        measurements = int(25e6) # int(20e3)
        max_sweeps = equilibration + measurements + 1
        lx = 50
        dyn = 'k'
        # Check dynamics type. Note: make sure all temps attached here are same as multirun.py.
        if dyn == 'g':
            temps_1 = np.linspace(1, 2, 10)
            temps_2 = np.linspace(2, 2.2, 20)
            temps_4 = np.linspace(2.2, 2.4, 40)
            temps_5 = np.linspace(2.4, 2.5, 10)
            temps_6 = np.linspace(2.5, 3.5, 10)
            temperatures = np.concatenate([temps_1, temps_2, temps_4, temps_5, temps_6])
        if dyn == 'k':
            temps_1 = np.linspace(0.1, 3, 30)
            temps_2 = np.linspace(-10, 0, 100)
            temps_2 = np.array([10 ** x for x in temps_2])
            temperatures = np.concatenate([temps_1, temps_2])  # , temps_2, temps_4, temps_5, temps_6])
        distincts = [("{0}_{1:.4f}_{2}_{3}").format(lx, T, dyn, max_sweeps) for T in temperatures]
        imgdirs = [os.getcwd() + "\\" + "img" + "\\" + distinct for distinct in distincts]
        files = [imgdir + "\\" + "results.txt" for imgdir in imgdirs]
        arrays = [] # avg_M, avg_E, avg_M_err, avg_E_err, chi_true, chi_error, c_true, c_error, self.T
        # Grab dumps and make table.
        for file in files:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass
        arrays=np.array(arrays)
        labels = ['avg_M', 'avg_E', 'avg_M_err', 'avg_E_err', 'chi_true', 'chi_error', 'c_true', 'c_error', 'T']
        table = Table(arrays, names=labels)
        table.sort('T')
        # Do figure for glauber
        if dyn == 'g' or "k":
            fig, axs = plt.subplots(nrows=2, ncols=2, sharex='col', figsize=(15,7))
            plt.subplots_adjust(wspace=0.15, hspace=0.1)
            x = table['T']
            colors = ['blue', 'red', 'blue', 'red']
            y_index = ['avg_M', 'avg_E', 'chi_true', 'c_true']
            y_error = ['avg_M_err', 'avg_E_err', 'chi_error', 'c_error']
            x_label = "Temperature"
            y_labels = ["Magnetisation", "Energy", "Magnetic Susceptibility", "Specific Heat p. Spin"]
            for i in range(2):
                for j in range(2):
                    ax = axs[i,j]
                    table_index = j
                    if i == 1:
                        table_index += 2
                    y = table[y_index[table_index]]
                    y_err = table[y_error[table_index]]
                    ax.errorbar(x=x, y=y, yerr=y_err, color=colors[table_index], ecolor='pink', lw=0.5)
                    ax.grid(True, which='major', alpha=0, linestyle='dotted')  # Enable grids on subplot
                    ax.grid(True, which='minor', alpha=0, linestyle='dotted')
                    ax.tick_params(axis="x", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                                   top=True)
                    ax.tick_params(axis="y", which="both", direction="in", length=4, bottom=True, left=True, right=True,
                                   top=True)
                    ax.set(xlabel=x_label,
                           ylabel=y_labels[table_index])

                    if table_index == 3:
                        if dyn == 'k':
                            ax.set(ylim=[0, 50])


        # Dump figure
        plt.savefig(dyn + "_multi.png", dpi=600)
        plt.show()
        # Save the table so that people can have a gander at it.
        writer = hdfutils.hdf5_writer(os.getcwd(), "multirun_data.hdf5")
        writer.write_table(dyn, "results", table)

    # High-scale run (no saving) that will output images to imgdir rather than keeping in memory/writing.
    """
    Estimated run time: 
    3.6 GHz Single Thread
    160x10^6 flips on an 800x800 
    7 hours (Glauber) 8 hours (Kawasaki) 
    """
    def run_high(self):
        # Define initial condition
        self.M, self.E = self.fast.fast_magenergy(self.mat)
        # Various fire prelim stuff
        fig = plt.figure(figsize=(8,8))
        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                           Patch(facecolor='blue', label='-1', edgecolor='black')]
        plt.legend(handles=legend_elements, loc='upper right')
        cmap = cm.get_cmap('bwr', 2)
        plt.xlim([0, self.lx-1])
        plt.ylim([0, self.lx - 1])
        im = plt.imshow(self.mat, animated=True, cmap=cmap, aspect='equal')
        plt.clim(-1, 1)
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        plt.title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,self.M,self.E))

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        # Run sim
        if self.dyn == "g":
            while self.sweep < self.max_sweeps:
                self.fast_glauber()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    plt.title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,self.M,self.E))
                    if self.sweep % self.temprate == 0:
                        self.T += self.temprise
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        if self.dyn == "k":
            while self.sweep < self.max_sweeps:
                self.fast_kawasaki()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    plt.title(("T = {0:.1f} M = {1:.1f} E = {2:.1f}").format(self.T,self.M,self.E))
                    if self.sweep % self.temprate == 0:
                        self.T += self.temprise
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)

    # Save data.
    """
    all_M and all_E are saved to data.txt as a pickled numpy array np.array([all_M, all_E])
    Deprecated format is pandas DataFrame (much too slow) so some code may reflect this. 
    """
    def save(self, table=None):
        if table == None:
            i = 0
            while i <= self.saveattempts:
                try:
                    # Dump it
                    with open(self.imgdir + "\\" + "data.txt", 'wb') as f:
                        pickle.dump(obj=np.array([self.all_M, self.all_E]), file=f)
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

    # Load data from file
    def load(self, table=None):
        try:
            if table == None:
                # Load it
                with open(self.imgdir + "\\" + "data.txt", 'rb') as f:
                    self.all_M, self.all_E = pickle.load(f)
            else:
                return self.writer.read_table(self.distinct, self.numset)
            self.time_delay("Loaded successfully.")
        except:
            self.time_delay("Load failed.")

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

    # Glauber, but using fast_ising.
    def fast_glauber(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_glauber(self.mat,
                                                                                       self.M,
                                                                                       self.E,
                                                                                       self.T,
                                                                                       self.sweep)

    # Kawasaki, but using fast_ising
    def fast_kawasaki(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_kawasaki(self.mat,
                                                                               self.M,
                                                                               self.E,
                                                                               self.T,
                                                                               self.sweep)

    # Averages and Errors, but fast_ising.
    def fast_averages_errors(self):
        # Calculate averages and errors
        avg_M, \
        avg_E, \
        avg_M_err, \
        avg_E_err, \
        chi_true, \
        chi_error, \
        c_true, \
        c_error = self.fast.averages_errors(self.all_M,
                                        self.all_E,
                                        self.T,
                                        self.equilibration,
                                        self.measurements,
                                        self.autocorrelation_length)


        # Only save if not multiprocessing
        if __name__ == "__main__":
            # Produce a table and save it.
            table = Table()
            table["M"], table["E"], table["chi"], table["c"] = np.array([avg_M, avg_M_err]), \
                                                               np.array([avg_E, avg_E_err]), \
                                                               np.array([chi_true, chi_error]), \
                                                               np.array([c_true, c_error])
            self.save(table)
        else:
            array_dump = np.array([avg_M, avg_E, avg_M_err, avg_E_err, chi_true, chi_error, c_true, c_error, self.T])
            with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
                pickle.dump(obj=array_dump, file=f)

    # TODO: DEPRECATED! See fast_ising.
    """
    ***This should reduce computation time when having to scour for periodic BCs.*** 
    Go through array of all points [i,j] 
    Lots and lots of if/else loops
    Save nearest-neighbour identifications for later use.
    """
    def build_nns(self):
        for i in range(self.lx):
            for j in range(self.lx):
                if 0 < i < self.lx - 1 and 0 < j < self.lx - 1:
                    """
                    nns format is [[top bot left right i_indices],[j_indices]]
                    """
                    self.nns[i,j] = np.array([[i-1,j],[i+1,j],[i,j-1],[i,j+1]]).T
                else:
                    # Bottom
                    if i == self.lx-1:
                        top, bot = [i-1,j], [0,j]
                    # Top
                    elif i == 0:
                        top, bot = [self.lx-1,j], [i+1,j]
                    # All other i
                    else:
                        top, bot = [i - 1, j], [i + 1, j]
                    # Right
                    if j == self.lx-1:
                        left, right = [i,j-1], [i,0]
                    # Left
                    elif j == 0:
                        left, right = [i,self.lx-1], [i,j+1]
                    # All other j
                    else:
                        left, right = [i, j - 1], [i, j + 1]
                    self.nns[i, j] = np.array([top, bot, left, right]).T # top bot left right

    # TODO: DEPRECATED! See fast_ising.
    def glauber_nonns(self):
        i, j = self.rng.integers(0, self.lx, 2)
        nnsum = np.sum(self.mat[[self.pbc(i-1), self.pbc(i+1), i, j], [j, j, self.pbc(j-1), self.pbc(j+1)]])
        # Calculate change in magnetisation (mag new - mag old)
        mag_cost = -2 * self.mat[i, j]
        # Generate an energy cost for doing this flip (energy new - energy old)
        energy_cost = -1 * mag_cost * nnsum
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[i, j] *= -1
            # Add change in mag/energy
            self.M += mag_cost
            self.E += energy_cost
        else:
            probability = math.exp(-1 * energy_cost / self.T)
            if self.rng.random(1)[0] <= probability:
                self.mat[i, j] *= -1
                # Add change in mag/energy
                self.M += mag_cost
                self.E += energy_cost
        # Sweep complete.
        self.sweep += 1

    # TODO: Report Issue Numba. See fast_ising information. Deprecated.
    def fast_kawaglauber(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_kawaglauber(self.mat,
                                                                                  self.M,
                                                                                  self.E,
                                                                                  self.T,
                                                                                  self.sweep)

    # TODO: DEPRECATED! See fast_ising.
    def glauber(self):
        # Generate a candidate location to flip and gather near neighbour sum
        candidate = self.rng.integers(0, self.lx, 2)
        nns = self.nns[candidate[0], candidate[1]]
        nnsum = np.sum(self.mat[nns[0], nns[1]])
        # Calculate change in magnetisation (mag new - mag old)
        mag_cost = -2 * self.mat[candidate[0], candidate[1]]
        # Generate an energy cost for doing this flip (energy new - energy old)
        energy_cost = -1 * mag_cost * nnsum
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidate[0], candidate[1]] *= -1
            # Add change in mag/energy
            self.M += mag_cost
            self.E += energy_cost
        else:
            probability = math.exp(-1 * energy_cost / self.T)
            if self.rng.random(1)[0] <= probability:
                self.mat[candidate[0], candidate[1]] *= -1
                # Add change in mag/energy
                self.M += mag_cost
                self.E += energy_cost
        # Sweep complete.
        self.sweep += 1

    # TODO: DEPRECATED! See fast_ising.
    def kawasaki(self):
        # Generate two candidate locations AB (and their spins) to flip, and gather near-neighbour sums.
        candidates = self.rng.integers(0, self.lx, (2, 2))
        canditrans = candidates.T
        spindidates = self.mat[canditrans[0], canditrans[1]]
        B_min_A = spindidates[1] - spindidates[0]
        # Gather neighbours and their spin sum for site B
        nns = self.nns[candidates[1][0], candidates[1][1]]
        nnsum = np.sum(self.mat[nns[0], nns[1]])
        # Gather change in magnetisation from setting B to A
        mag_cost_ba = -1 * B_min_A
        # Get energy cost
        energy_cost_ba = -1 * mag_cost_ba * nnsum
        # Create new array (N^2 operation I think- idk for sure?...) and set B to A.
        mat_baab = copy.copy(self.mat)  # copy to avoid any risk of altering original matrix (unsure.)
        mat_baab[candidates[1][0], candidates[1][1]] = spindidates[0]
        # Gather neighbours and their spin sum for site A
        nns = self.nns[candidates[0][0], candidates[0][1]]
        nnsum = np.sum(mat_baab[nns[0], nns[1]])
        # Gather change in magnetisation from setting A to B
        mag_cost_ab = B_min_A
        energy_cost_ab = -1 * mag_cost_ab * nnsum
        # Total magnetism cost (zero for kawasaki.)
        # mag_cost = 0
        # Total energy cost
        energy_cost = energy_cost_ba + energy_cost_ab
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidates[1][0], candidates[1][1]] = spindidates[0]
            self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
            self.E += energy_cost
        else:
            probability = math.exp(-1 * energy_cost / self.T)
            if self.rng.random(1)[0] <= probability:
                self.mat[candidates[1][0], candidates[1][1]] = spindidates[0]
                self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
                self.E += energy_cost
        # Sweep complete.
        self.sweep += 1

    # TODO: DEPRECATED! See fast_ising.
    """
    List functionality isn't needed anymore (this was from when I was only dealing with 20,000 or so "sweeps." 
    See information regarding how I misconstrued sweep=flip.
    """
    def magenergy(self, sweep_list):
        all_M = [np.sum(i) for i in sweep_list]

        def energy(matrix):
            E = 0
            for i in range(self.lx):
                for j in range(self.lx):
                    nns = self.nns[i, j]
                    nnsum = np.sum(matrix[nns[0], nns[1]])
                    E += -1 * matrix[i, j] * nnsum
            E = (1 / 2) * E  # doublecounting
            return E

        all_E = [energy(i) for i in sweep_list]
        return all_M, all_E

    # TODO: DEPRECATED! See fast_ising.
    """
    The following is done for the range equilibration -> equilibration + measurements 
    Calculate the average values of magnetisation, energy, and their errors
    Calculate susceptibility, scaled heat capacity, and their errors (specify error type.)
    ONE SWEEP IS 2,500 ATTEMPTED FLIPPEROOS! DON'T YOU DARE FORGET!
    """
    def averages_errors(self):
        # Trim data to measurement range
        raw_dataframe = self.data[self.equilibration:self.equilibration + self.measurements]
        # Sample every n'th measurement, autolength from autocorrelation
        autolength = self.autocorrelation_length
        samplelength = self.autocorrelation_length
        df = raw_dataframe[raw_dataframe.index % samplelength == 0]
        num_samples = len(df)
        # Grab Magnetism/Energy samples
        all_M, all_E = df['all_M'].to_numpy(), df['all_E'].to_numpy()
        all_MM, all_EE = all_M ** 2, all_E ** 2
        # Get averages
        avg_M, avg_E = np.average(all_M), np.average(all_E)
        avg_MM, avg_EE = np.average(all_MM), np.average(all_EE)
        # Get Errors
        avg_M_err, avg_E_err = np.sqrt(((avg_MM - avg_M ** 2) * (2 * autolength)) / (num_samples * samplelength)), \
                               np.sqrt(((avg_EE - avg_E ** 2) * (2 * autolength)) / (num_samples * samplelength))
        # Estimate susceptibility and specific heat/spin
        chi_true, c_true = (1 / num_samples) * (1 / self.T) * (avg_MM - avg_M ** 2), \
                           (1 / num_samples) * (1 / self.T ** 2) * (avg_EE - avg_E ** 2)
        # Error estimation for chi and c via the Bootstrap method
        # Bootstrap error as defined by https://thestatsgeek.com/2013/07/02/the-miracle-of-the-bootstrap/
        # The error in the notes didn't match up to the error methods online so I just went with these...
        chi_list, c_list = [], []
        number_of_resamples = 1000
        for i in range(number_of_resamples):
            # Select num_samples random from num_samples
            resample = self.rng.integers(0, num_samples, num_samples)
            # Grab Magnetism/Energy samples
            Qall_M, Qall_E = all_M[resample], all_E[resample]
            Qall_MM, Qall_EE = all_M ** 2, all_E ** 2
            # Get averages
            Qavg_M, Qavg_E = np.average(Qall_M), np.average(Qall_E)
            Qavg_MM, Qavg_EE = np.average(Qall_MM), np.average(Qall_EE)
            # Estimate susceptibility and specific heat/spin
            Qchi = (1 / num_samples) * (1 / self.T) * (Qavg_MM - Qavg_M ** 2)
            Qc = (1 / num_samples) * (1 / self.T ** 2) * (Qavg_EE - Qavg_E ** 2)
            # Append
            chi_list.append(Qchi), c_list.append(Qc)
        chi_list, c_list = np.array(chi_list), np.array(c_list)
        chi_average, c_average = np.average(chi_list), np.average(c_list)
        boot_chi, boot_c = np.sqrt((1 / number_of_resamples) * np.sum((chi_list - chi_average) ** 2)), \
                           np.sqrt((1 / number_of_resamples) * np.sum((c_list - c_average) ** 2))
        chi_error, c_error = boot_chi, boot_c
        # Only save if not multiprocessing
        if __name__ == "__main__":
            # Produce a table and save it.
            table = Table()
            table["M"], table["E"], table["chi"], table["c"] = np.array([avg_M, avg_M_err]), \
                                                               np.array([avg_E, avg_E_err]), \
                                                               np.array([chi_true, chi_error]), \
                                                               np.array([c_true, c_error])
            self.save(table)
        else:
            array_dump = np.array([avg_M, avg_E, avg_M_err, avg_E_err, chi_true, chi_error, c_true, c_error, self.T])
            with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
                pickle.dump(obj=array_dump, file=f)

    # TODO: Deprecated. See "sweeps=flips" info. See "self.save" info.
    def autocorrelation(self):
        # Load the energy/magnetisation
        all_M, all_E = list(self.data['all_M']), list(self.data['all_E'])
        # Clip to measurements
        all_M, all_E = all_M[self.equilibration:self.equilibration + self.measurements],\
                       all_E[self.equilibration:self.equilibration + self.measurements]
        # Generate a correlation given
        pandas.plotting.autocorrelation_plot(all_M)
        plt.savefig("autocorrelation_" + self.distinct + ".png")

    # TODO: Deprecated
    def anigen_high(self):
        image_files = [os.path.join(self.imgdir, img)
                       for img in os.listdir(self.imgdir)
                       if img.endswith(".png")]
        image_files.sort(key=lambda f: int(re.sub('\D', '', f)))
        clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=self.framerate)
        clip.write_videofile(self.animation_name)

    # TODO: Deprecated
    def anigen(self):
        if self.framerate == 0:
            pass
        else:
            # Grab data
            all_sweeps = list(self.data['all_sweeps'])

            # Trim down by binrate
            selected_sweeps = np.arange(0, len(all_sweeps), self.binrate)
            all_sweeps = [all_sweeps[d] for d in selected_sweeps]

            # Set up figure
            fig, ax = plt.subplots()

            # Set up image plot
            legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                               Patch(facecolor='blue', label='-1', edgecolor='black')]
            plt.legend(handles=legend_elements, loc='upper right')
            cmap = cm.get_cmap('bwr', 2)
            plt.xlim([0, self.lx - 1])
            plt.ylim([0, self.lx - 1])
            im = ax.imshow(all_sweeps[0], animated=True, cmap=cmap, aspect='equal', vmin=-1, vmax=1)
            t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)
            ax.text(1, self.lx - 5, str("T = ") + str(self.T), va='top', fontsize=20, color='white')
            self.time_delay("Animation Progress...")
            def updatefig(i):
                im.set_array(all_sweeps[i])
                t1.set_text(str(selected_sweeps[i]))
                # Print is included inside here to suppress spam in the console if multiprocessing.
                if self.draw == True:
                    print(("{0:.2f}%").format(100 * (i + 1) / len(selected_sweeps)))
                    fig.canvas.draw()  # draw
                    fig.canvas.flush_events()  # deal with resize
                return im,
            ani = anim.FuncAnimation(fig, updatefig, interval=1e3/self.framerate, frames=len(all_sweeps), blit=True, repeat=True)
            ani.save(self.animation_name, 'ffmpeg_file')
            plt.close()

    # TODO: Deprecated
    def delete(self):
        try:
            os.remove(self.filename)
        except:
            print("Failed to delete file. Do it manually, mate.")

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
        self.time_delay("Yōkoso!!! Welcome to this 2D Ising Simulation Suite. \n"
                        "You will now be asked for a few parameters. Please give them. \n"
                        "Please note that the code is a bit optimized for multiple sequential runs (i.e. parallel) \n"
                        "Due to this, a one-time-run will incur a @jit compile cost, compared to regular python. \n"
                        "Now, onto the parameters!!!!")
        print()
        self.time_delay("Grid size. This is a square simulator, so just one integer will suffice.")
        lx = int(input())
        self.time_delay("Temperature (note that J=KB=unity in this model.")
        T = float(input())
        self.time_delay("Dynamical mode: g for glauber, k for kawasaki.")
        dyn = str(input())
        if dyn == "g":
            self.time_delay("You said Glauber. Do you want initial spins entirely random? y/n")
            not_up = str(input())
            if not_up == 'y':
                not_up = True
        else:
            not_up = False
        self.time_delay("What binrate do you want for the sim? (spacing between animation prints. Code bottleneck.)")
        binrate = int(input())
        print()
        self.time_delay("We will use the default equilibration and measurement periods for this run. Enjoy!")
        equi, measure = int(0.25e6), int(25e6)
        return lx, T, dyn, equi, measure, not_up, binrate

    # Run
    def run(self):
        try:
            model = twod_ising(*self.user_input())
            model.run(checkpoint=True)
        except Exception as e:
            self.time_delay("An error occurred... \n"
                            "You probably don't have the correct dependencies. \n"
                            "If the error regards the existence of LaTex: delete lines 21,22 \n" 
                            "If the error regards missing packages, please install them.\n")
            print(e)

