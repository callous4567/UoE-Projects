import copy
import math
import re
import sys
import time
import pandas
from matplotlib import rc, cm
from matplotlib import animation as anim
import os
from matplotlib.patches import Patch
import hdfutils
import numpy as np
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\Callicious\\Documents\\Prog\\pycharm\\venv\\ffmpeg\\bin\\ffmpeg.exe'
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Square-lattice twod Ising following MVP specification (Checkpoint 1).
class twod_ising(object):
    def __init__(self):
        # Specify the lattice dimension, system size, temperature, and dynamics. We take J=Kb=1.
        self.lx = 50
        self.T = 1
        self.dyn = "g"
        # Placeholder to recognize the current "sweep" that has just been carried out
        self.sweep = 0
        self.equilibration = int(1e3)
        self.measurements = int(20e3)
        self.max_sweeps = self.equilibration + self.measurements + 1
        # This holds a graph [0,1,2,3...] for rolling variable plotting, etc. Cheaper than generating it.
        self.graph_sweeps = []
        # Use hard-coded parameters or ask the user?
        self.hard_coded = True
        self.delay_max = 5e-3 # for user input typography, delay rate of text
        self.framerate = 15 # for animated mp4 produced, if applicable
        self.binrate = 2000 # only produce an image from every "binrate" in the sweeps.
        # For writing the storage file for this particular instance. Overwrite if desired. Group is 'distinct.'
        self.directory = os.getcwd()
        self.filename = "datafile.hdf5"
        self.set = "data"
        self.writer = hdfutils.hdf5_writer(self.directory, self.filename)
        # Placeholder for the indices for nearest neighbours
        self.nns = np.ndarray(shape=(self.lx, self.lx), dtype=list)
        # Generate a randomized initial spin matrix
        self.mat = np.random.default_rng().choice([-1, 1], size=(self.lx, self.lx))
        # Build near neighbours
        self.build_nns()
        # For storing/writing images
        imgdir = self.directory + "\\" + "img"
        try:
            os.mkdir(imgdir)
        except:
            pass
        # Unique string stating important parameters ***used as group for storage inside file.)
        self.distinct = str(self.lx) + "_" + str(self.T) + "_" + str(self.dyn) + "_" + str(self.max_sweeps)
        self.imgdir = imgdir + "\\" + self.distinct
        try:
            os.mkdir(self.imgdir)
        except:
            pass
        # Placeholder for dataframe
        self.data = None

    # Save data to file
    def save(self):
        self.writer.write_df(self.distinct, self.set, self.data)

    # Load data from file
    def load(self):
        self.data = self.writer.read_df(self.distinct, self.set)

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        print()
        for c in text:
            sys.stdout.write(c)
            sys.stdout.flush()
            resttime = np.random.default_rng().uniform(0.001, self.delay_max)
            time.sleep(resttime)

    # Build a catalogue of nearest-neighbours for each and every [i,j].
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
                    self.nns[i,j] = np.array([[i-1,j],[i+1,j],[i,j-1],[i,j+1]]) # top bot left right
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
                    self.nns[i, j] = np.array([top, bot, left, right]) # top bot left right

    # Visualize the matrix, and optionally plot the neighbours for [i,j] (to see if nns works... it seems to.)
    def visualize(self, nns=None):
        if nns == None:
            cmap = cm.get_cmap('bwr',2)
            plt.imshow(self.mat, cmap=cmap, aspect='equal')
            legend_elements = [Patch(facecolor='blue', label='-1', edgecolor='black'),
                               Patch(facecolor='red', label='+1', edgecolor='black')]
            plt.legend(handles=legend_elements, loc='upper right')
            plt.show()
        # plot a 5-point star for the nns selected (debug provision.)
        # Also verify that all nns lists are of cardinality 4
        else:
            for i in range(self.lx):
                for j in range(self.lx):
                    if len(self.nns[i,j]) != 4:
                        print("nns build was bad.")
            cmap = cm.get_cmap('bwr', 3)
            # Generate points that are the centre and near neighbour and mark as white
            NNS = self.nns[nns[0],nns[1]]
            NNS = np.append(NNS, [np.array(nns)], axis=0).T
            mat_copy = copy.deepcopy(self.mat)
            mat_copy[NNS[0],NNS[1]] = 0
            plt.imshow(mat_copy, cmap=cmap, aspect='equal')
            legend_elements = [Patch(facecolor='blue', label='-1', edgecolor='black'),
                               Patch(facecolor='red', label='+1', edgecolor='black')]
            plt.legend(handles=legend_elements, loc='upper right')
            plt.show()

    # Glauber-dynamic change of the matrix.
    """
    Randomly select point in matrix
    Flip the point (multi by -1.) 
    Calculate energy cost of doing so 
    Return energy cost
    """
    def glauber(self):
        # Generate a candidate location to flip and gather near neighbour sum
        candidate = np.random.default_rng().integers(0, self.lx, 2)
        nns = self.nns[candidate[0],candidate[1]].T
        nnsum = np.sum(self.mat[nns[0],nns[1]])
        # Generate an energy cost for doing this flip
        energy_cost = 2*self.mat[candidate[0],candidate[1]]*nnsum
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidate[0],candidate[1]] *= -1
        else:
            probability = math.exp(-1*energy_cost/self.T)
            if np.random.default_rng().random(1)[0] <= probability:
                self.mat[candidate[0], candidate[1]] *= -1
        # Sweep complete.
        self.sweep += 1

    # Kawasaki-dynamic change of the matrix.
    """
    If they're not next to each-other, I'm doing it as a "simultaneous change" (ii) in the MVP specification.
    They're next to each-other if the distance between them is less than sqrt(2) which I will approx as 1.4.
    If they're next to each-other, I'm doing it via changing B to A in the entire file (which is slower.)
    """
    def kawasaki(self):
        # Generate two candidate locations AB (and their spins) to flip, and gather near-neighbour sums.
        candidates = np.random.default_rng().integers(0, self.lx, (2,2))
        distance = np.linalg.norm(candidates[1]-candidates[0])
        canditrans = candidates.T
        spindidates = self.mat[canditrans[0], canditrans[1]]
        B_min_A = spindidates[1]-spindidates[0]
        # For nearest-neighbours, deal with copying entire matrix.
        if distance < 1.4:
            # Gather neighbours and their spin sum for site B, get a cost for setting B to A.
            nns = self.nns[candidates[1][0],candidates[1][1]].T
            nnsum = np.sum(self.mat[nns[0],nns[1]])
            energy_cost_ba = B_min_A*nnsum
            # Create new array (N^2 operation I think- idk for sure?...) and set B to A.
            mat_baab = self.mat
            mat_baab[candidates[1][0],candidates[1][1]] = spindidates[0]
            # Gather neighbours and their spin sum for site A, get a cost for setting A to B in new array.
            nns = self.nns[candidates[0][0], candidates[0][1]].T
            nnsum = np.sum(mat_baab[nns[0], nns[1]])
            energy_cost_ab = -1*B_min_A*nnsum
            # Total energy cost
            energy_cost = energy_cost_ba + energy_cost_ab
        # For non-nearest neighbours
        else:
            # Literally just two glauber energy calculations added up: see report pdf.
            nnss = [self.nns[candidate[0], candidate[1]].T for candidate in candidates]
            nnsums = [np.sum(self.mat[nns[0], nns[1]]) for nns in nnss]
            energy_cost = B_min_A*nnsums[1] - B_min_A*nnsums[0]
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidates[1][0],candidates[1][1]] = spindidates[0]
            self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
        else:
            probability = math.exp(-1*energy_cost/self.T)
            if np.random.default_rng().random(1)[0] <= probability:
                self.mat[candidates[1][0], candidates[1][1]] = spindidates[0]
                self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
        # Sweep complete.
        self.sweep += 1

    # User Input
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to this 2D Ising Simulation Suite. Please provide various parameters.")
        self.time_delay("Grid size. This is a square simulator, so just one integer will suffice.")
        self.lx = int(input())
        self.time_delay("Temperature (note that J=KB=unity in this model.")
        self.T = float(input())
        self.time_delay("Maximum number of sweeps to run for")
        self.max_sweeps = int(input())
        self.time_delay("Dynamical mode: g for glauber, k for kawasaki.")
        self.dyn = str(input())
        self.time_delay("If you do not wish to produce an animation for this, respond with a 0."
                        "Framerate in FPS for output video of simulation. "
                        "Note that by default, this program will produce an animation."
                        "You should further note that it takes literally forever to produce.")
        self.framerate = int(input())
    # Run the simulator (this is for the actual checkpoint rather than run_high for fun.)
    """
    Given that this simulation is very small, saving everything to file is perfectly reasonable.
    I've removed save functionality from run_high (made for 200, 400x matrices, with 10^6 sweeps- too big.) 
    """
    def run(self):
        # User input specific that was desired.
        if self.hard_coded == True:
            pass
        else:
            self.user_input()

        # List that holds all the individual steps. Only suitable for small cases.
        all_sweeps = [copy.copy(self.mat)]

        # Run Sim
        self.time_delay("Starting Simulation..." + self.distinct)
        if self.dyn == 'g':
            while self.sweep < self.max_sweeps:
                self.glauber()
                all_sweeps.append(copy.copy(self.mat))
        # For Kawasaki scenario
        if self.dyn == 'k':
            while self.sweep < self.max_sweeps:
                self.kawasaki()
                all_sweeps.append(copy.copy(self.mat))

        # Complete Sim
        self.time_delay("Simulation is done. Grabbing Energy/Magnetisation. Takes a while.")

        # Now that simulation is done, grab magnetisation and energy for all sweeps.***takes a while!***
        all_M, all_E = self.magenergy(all_sweeps)

        # Create df and save it to file.
        self.time_delay("Creating dataframe...")
        df = pandas.DataFrame()
        df['all_sweeps'] = all_sweeps
        df['all_M'], df['all_E'] = all_M, all_E
        self.data = df
        self.save()
        self.time_delay("Data save successful. No need to re-run simulation: use load() to reload dataframe.")

    # High-scale run (no saving) that will output images to imgdir rather than keeping in memory/writing.
    """
    This is for the case of 200x200+ arrays with 1e6+ sweeps.
    This is just for producing images at the moment, but you could conceivably produce scientific results, too.
    I just wanted to produce large-array movies- no science data output. 
    """
    def run_high(self):
        # Various fire prelim stuff
        fig = plt.figure(figsize=(8,8))
        # Set up image plot
        legend_elements = [Patch(facecolor='blue', label='-1', edgecolor='black'),
                           Patch(facecolor='red', label='+1', edgecolor='black')]
        plt.legend(handles=legend_elements, loc='upper right')
        cmap = cm.get_cmap('bwr', 2)
        plt.xlim([0, self.lx-1])
        plt.ylim([0, self.lx - 1])
        im = plt.imshow(self.mat, animated=True, cmap=cmap, aspect='equal')
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        plt.text(1, self.lx-5, str("T = ") + str(self.T), va='top', fontsize=20, color='white')

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")
        # Run sim
        if self.dyn == "g":
            while self.sweep < self.max_sweeps:
                self.glauber()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")
        if self.dyn == "k":
            while self.sweep < self.max_sweeps:
                self.kawasaki()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")

    # Snip all images in img together. This is faster than anigen by miles.
    def anigen_high(self):
        fps = self.framerate
        image_files = [os.path.join(self.imgdir, img)
                       for img in os.listdir(self.imgdir)
                       if img.endswith(".png")]
        image_files.sort(key=lambda f: int(re.sub('\D', '', f)))
        clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
        clip.write_videofile('simulation.mp4')

    # Generate (and save) animation using all_sweeps.
    def anigen(self):
        if self.framerate == 0:
            pass
        else:
            # Grab data
            all_sweeps = self.all_sweeps
            # Trim down by binrate
            selected_sweeps = np.arange(0, len(all_sweeps), self.binrate)
            all_sweeps = [all_sweeps[d] for d in selected_sweeps]

            # Set up figure
            fig = plt.figure()
            legend_elements = [Patch(facecolor='blue', label='-1', edgecolor='black'),
                               Patch(facecolor='red', label='+1', edgecolor='black')]
            plt.xlim([0, self.lx - 1])
            plt.ylim([0, self.lx - 1])
            plt.legend(handles=legend_elements, loc='upper right')
            cmap = cm.get_cmap('bwr', 2)
            im = plt.imshow(all_sweeps[0], animated=True, cmap=cmap, aspect='equal')
            def updatefig(i):
                print(i)
                im.set_array(all_sweeps[i])
            ani = anim.FuncAnimation(fig, updatefig, interval=1e3/self.framerate, frames=len(all_sweeps))
            ani.save("test.mp4", 'ffmpeg_file')

    # Calculate Magnetisation/Energy for each and every sweep, and save it to "all_M" "all_E"
    def magenergy(self, sweep_list):
        all_M = [np.sum(i) for i in sweep_list]
        def energy(matrix):
            E = 0
            for i in range(self.lx):
                for j in range(self.lx):
                    nns = self.nns[i,j].T
                    nnsum = np.sum(matrix[nns[0], nns[1]])
                    E += -1*matrix[i,j]*nnsum
            E = (1/2)*E # doublecounting
            return E
        all_E = [energy(i) for i in sweep_list]
        return all_M, all_E

    # Create an autocorrelation plot (to estimate autocorrelation time.)
    def autocorrelation(self):
        # Load the energy/magnetisation
        all_M, all_E = self.data['all_M'], self.data['all_E']
        # Clip to measurements
        all_M, all_E = all_M[self.equilibration:self.equilibration + self.measurements],\
                       all_E[self.equilibration:self.equilibration + self.measurements]
        #


    # Calculate various things in measurement range
    """
    The following is done for the range equilibration -> equilibration + measurements 
    Calculate the average values of magnetisation, energy, and their errors
    Calculate susceptibility, scaled heat capacity, and their errors (specify error type.)
    """
    def averages_errors(self):
        # Get measurement sweeps (i.e. ignore equilibration)
        #sweep_measurements = self.all_sweeps[self.equilibration:self.equilibration + self.measurements]
        hey = None


uwu = twod_ising()
uwu.load()
#uwu.visualize(nns=[100, 0])
#uwu.run_high()
#uwu.anigen_high()
#uwu.run()
#uwu.autocorrelation()