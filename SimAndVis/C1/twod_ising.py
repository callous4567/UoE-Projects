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
from scipy.signal import savgol_filter

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
        self.T = 2.2
        self.dyn = "g"
        # Placeholder to recognize the current "sweep" that has just been carried out
        self.sweep = 0
        self.M, self.E = None, None # current E and M: only applicable for run()
        # Specify equilibration sweeps, measurement sweeps, and simulation maximum
        self.equilibration = int(0.25e6) # int(2e3)
        self.measurements = int(25e6) # int(20e3)
        self.max_sweeps = self.equilibration + self.measurements + 1
        # Use hard-coded parameters or ask the user?
        self.hard_coded = True
        self.delay_max = 5e-3 # for user input typography, delay rate of text
        # For writing the storage file for this particular instance. Overwrite if desired. Group is 'distinct.'
        self.directory = os.getcwd()
        self.filename = "datafile.hdf5"
        self.set = "data" # for dataframe
        self.numset = "results" # for numerical results in an astropy table: M, E, CHI, C, and errors.
        if __name__ == "__main__":
            self.writer = hdfutils.hdf5_writer(self.directory, self.filename)
        # Placeholder for the nearest-neighbour array
        self.nns = np.ndarray(shape=(self.lx, self.lx), dtype=list)
        # Specify if you want all up, all down, or totally random (set to None for random.) +-1 for Up/Down.
        all_up = 1
        if all_up != None:
            self.mat = all_up*np.ones(shape=(self.lx, self.lx))
        else:
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
        self.distinct = ("{0}_{1:.2f}_{2}_{3}").format(self.lx, self.T, self.dyn, self.max_sweeps)
        self.imgdir = imgdir + "\\" + self.distinct
        # Only create this directory here if main. If not main, create directory in init_multiprocess.
        if __name__ == "__main__":
            try:
                os.mkdir(self.imgdir)
            except:
                pass
        # Various animation parameters
        self.animation_name = self.imgdir + "\\animation.mp4"
        self.draw = False # set to True to produce a live animation.
        self.framerate = 0  # 15 # for animated mp4 produced, if applicable. Set to 0 for no animation.
        self.binrate = 200 # binrate for animations. Every i'th is selected as frames.
        self.temprate = 4e6 # raise temperature every temprate: only applicable for run_high (for 80 million, 20 rises.)
        self.temprise = 0.5 # by this amount: only applicable for run_high
        # Placeholder for dataframe
        self.data = None
        # Save attempts for saving the file. After this, saving will be abandoned and stuff goes bad.
        self.saveattempts = 5

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
        # Calculate change in magnetisation (mag new - mag old)
        mag_cost = -2*self.mat[candidate[0], candidate[1]]
        # Generate an energy cost for doing this flip (energy new - energy old)
        energy_cost = -1*mag_cost*nnsum
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidate[0],candidate[1]] *= -1
            # Add change in mag/energy
            self.M += mag_cost
            self.E += energy_cost
        else:
            probability = math.exp(-1*energy_cost/self.T)
            if np.random.default_rng().random(1)[0] <= probability:
                self.mat[candidate[0], candidate[1]] *= -1
                # Add change in mag/energy
                self.M += mag_cost
                self.E += energy_cost
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
            # Gather neighbours and their spin sum for site B
            nns = self.nns[candidates[1][0],candidates[1][1]].T
            nnsum = np.sum(self.mat[nns[0],nns[1]])
            # Gather change in magnetisation from setting B to A
            mag_cost_ba = -1*B_min_A
            # Get energy cost
            energy_cost_ba = -1*mag_cost_ba*nnsum
            # Create new array (N^2 operation I think- idk for sure?...) and set B to A.
            mat_baab = self.mat
            mat_baab[candidates[1][0],candidates[1][1]] = spindidates[0]
            # Gather neighbours and their spin sum for site A
            nns = self.nns[candidates[0][0], candidates[0][1]].T
            nnsum = np.sum(mat_baab[nns[0], nns[1]])
            # Gather change in magnetisation from setting A to B
            mag_cost_ab = B_min_A
            energy_cost_ab = -1*mag_cost_ab*nnsum
            # Total magnetism cost (zero for kawasaki.)
            # mag_cost = 0
            # Total energy cost
            energy_cost = energy_cost_ba + energy_cost_ab
        # For non-nearest neighbours
        else:
            # Literally just two glauber energy calculations added up: see report pdf.
            # Gather neighbours/etc.
            nnss = [self.nns[candidate[0], candidate[1]].T for candidate in candidates]
            nnsums = [np.sum(self.mat[nns[0], nns[1]]) for nns in nnss]
            # Get costs.
            # mag_cost = 0
            energy_cost = B_min_A*nnsums[1] - B_min_A*nnsums[0]
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            self.mat[candidates[1][0],candidates[1][1]] = spindidates[0]
            self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
            self.M += 0
            self.E += energy_cost
        else:
            probability = math.exp(-1*energy_cost/self.T)
            if np.random.default_rng().random(1)[0] <= probability:
                self.mat[candidates[1][0], candidates[1][1]] = spindidates[0]
                self.mat[candidates[0][0], candidates[0][1]] = spindidates[1]
                self.M += 0
                self.E += energy_cost
        # Sweep complete.
        self.sweep += 1

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

        # List that holds all the individual steps. Only suitable for small cases. Used to animate and save sim.
        all_sweeps = [copy.copy(self.mat)]
        # Calculate M and E for initial step
        all_M, all_E = self.magenergy(all_sweeps)
        self.M, self.E = all_M[0], all_E[0]

        # Run Sim
        self.time_delay("Starting Simulation..." + self.distinct)
        if self.dyn == 'g':
            while self.sweep < self.max_sweeps:
                self.glauber()
                all_sweeps.append(copy.copy(self.mat))
                all_M.append(self.M), all_E.append(self.E)
        # For Kawasaki scenario
        if self.dyn == 'k':
            while self.sweep < self.max_sweeps:
                self.kawasaki()
                all_sweeps.append(copy.copy(self.mat))
                all_M.append(self.M), all_E.append(self.E)

        # Complete Sim
        self.time_delay("Simulation is done. Running quick final checks...")

        # Generate final magnetisation/energy (if you want to double check that rolling calcs match doing it manually.)
        #fin_M, fin_E = self.magenergy([all_sweeps[-1]])
        #print(all_M[-1], fin_M[0], all_E[-1], fin_E[0])

        # Now that simulation is done, grab magnetisation and energy for all sweeps.***takes a while!***
        # Update: We've made the sim calculate and sum changes on the fly. Much cheaper computationally.
        # all_M, all_E = self.magenergy(all_sweeps)

        # Create df
        self.time_delay("Creating dataframe...")
        df = pandas.DataFrame()
        df['all_sweeps'] = all_sweeps
        df['all_M'], df['all_E'] = all_M, all_E
        self.data = df
        # Only save if not multiprocessing.
        if __name__ == "__main__":
            self.save()
            self.time_delay("Data save successful. No need to re-run simulation: use load() to reload dataframe.")

    # For multiprocessed version: no sweep saving.
    def run_multi(self):
        # Calculate M and E for initial step
        all_M, all_E = self.magenergy([self.mat])
        self.M, self.E = all_M[0], all_E[0]

        # Run Sim
        self.time_delay("Starting Simulation..." + self.distinct)
        if self.dyn == 'g':
            while self.sweep < self.max_sweeps:
                self.glauber()
                all_M.append(self.M), all_E.append(self.E)
        # For Kawasaki scenario
        if self.dyn == 'k':
            while self.sweep < self.max_sweeps:
                self.kawasaki()
                all_M.append(self.M), all_E.append(self.E)

        # Create df
        self.time_delay("Creating dataframe...")
        df = pandas.DataFrame()
        df['all_M'], df['all_E'] = all_M, all_E
        self.data = df

        # Dump it
        with open(self.imgdir + "\\" + "data.txt", 'wb') as f:
            pickle.dump(obj=self.data, file=f)

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

    # Calculate various things in measurement range
    """
    The following is done for the range equilibration -> equilibration + measurements 
    Calculate the average values of magnetisation, energy, and their errors
    Calculate susceptibility, scaled heat capacity, and their errors (specify error type.)
    ONE SWEEP IS 2,500 ATTEMPTED FLIPPEROOS!
    """
    def averages_errors(self):
        # Trim data to measurement range
        raw_dataframe = self.data[self.equilibration:self.equilibration + self.measurements]
        # Sample every n'th measurement, autolength from autocorrelation
        autolength = 250000
        samplelength = autolength
        df = raw_dataframe[raw_dataframe.index % samplelength == 0]
        num_samples = len(df)
        # Grab Magnetism/Energy samples
        all_M, all_E = df['all_M'].to_numpy(), df['all_E'].to_numpy()
        all_MM, all_EE = all_M**2, all_E**2
        # Get averages
        avg_M, avg_E = np.average(all_M), np.average(all_E)
        avg_MM, avg_EE = np.average(all_MM), np.average(all_EE)
        # Get Errors
        avg_M_err, avg_E_err = np.sqrt(((avg_MM - avg_M**2)*(2*autolength))/(num_samples*samplelength)), \
                               np.sqrt(((avg_EE - avg_E**2)*(2*autolength))/(num_samples*samplelength))
        # Estimate susceptibility and specific heat/spin
        chi_true, c_true = (1/num_samples)*(1/self.T)*(avg_MM - avg_M**2),\
                           (1/num_samples)*(1/self.T**2)*(avg_EE - avg_E**2)
        # Error estimation for chi and c via the Bootstrap method
        chi_list, c_list = [],[]
        number_of_resamples = 1000
        for i in range(number_of_resamples):
            # Select num_samples random from num_samples
            resample = np.random.default_rng().integers(0, num_samples, num_samples)
            # Grab Magnetism/Energy samples
            Qall_M, Qall_E = all_M[resample], all_E[resample]
            Qall_MM, Qall_EE = all_M ** 2, all_E ** 2
            # Get averages
            Qavg_M, Qavg_E = np.average(Qall_M), np.average(Qall_E)
            Qavg_MM, Qavg_EE = np.average(Qall_MM), np.average(Qall_EE)
            # Estimate susceptibility and specific heat/spin
            Qchi = (1/num_samples)*(1/self.T)*(Qavg_MM - Qavg_M**2)
            Qc = (1/num_samples)*(1/self.T**2)*(Qavg_EE - Qavg_E**2)
            # Append
            chi_list.append(Qchi), c_list.append(Qc)
        chi_list, c_list = np.array(chi_list), np.array(c_list)
        chi_average, c_average = np.average(chi_list), np.average(c_list)
        # Bootstrap error as defined by https://thestatsgeek.com/2013/07/02/the-miracle-of-the-bootstrap/
        # The error in the notes didn't match up to the error methods online so I just went with these...
        boot_chi, boot_c = np.sqrt((1/number_of_resamples)*np.sum((chi_list - chi_average)**2)), \
                           np.sqrt((1/number_of_resamples)*np.sum((c_list - c_average)**2))
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
            # Confusing as hell but I forgot to relabel. Anyway, "g.txt" is the array that has the above.
            with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
                pickle.dump(obj=array_dump, file=f)


    # Full Run
    def main(self):
        self.run()
        self.anigen()
        self.averages_errors()

    # Multi
    def main_multi(self, run=True):
        # To run the simulation
        if run == True:
            self.run_multi()
            self.averages_errors()
        # Re-generate averages with new parameters you've added to multigraph
        else:
            with open(self.imgdir + "\\" + "data.txt", 'rb') as f:
                self.data = pickle.load(f)
            self.averages_errors()

    # Produce graphs for multiprocessed runs. Only applicable for multiprocessed runs.
    def multigraph(self):
        # Same "settings" as used by Multirun (dynamics are "k" or "g" obviously.)
        self.equilibration = int(1e5) # int(2e3)
        self.measurements = int(1e6) # int(20e3)
        max_sweeps = self.equilibration + self.measurements + 1
        lx = 50
        dyn = 'g'

        if dyn == 'g':
            temps_1 = np.linspace(1, 2, 10)
            temps_2 = np.linspace(2, 2.2, 20)
            temps_4 = np.linspace(2.2, 2.4, 40)
            temps_5 = np.linspace(2.4, 2.5, 10)
            temps_6 = np.linspace(2.5, 3.5, 10)
            temperatures = np.concatenate([temps_1, temps_2, temps_4, temps_5, temps_6])
        if dyn == 'k':
            temps_1 = np.linspace(0.0001, 0.1, 30)
            temps_2 = np.linspace(0.1, 2, 5)
            temps_4 = np.linspace(2.2, 2.4, 10)
            temps_5 = np.linspace(2.4, 2.5, 10)
            temps_6 = np.linspace(2.5, 3.5, 10)
            temperatures = np.concatenate([temps_1, temps_2, temps_4, temps_5, temps_6])


        distincts = [("{0}_{1:.2f}_{2}_{3}").format(lx, T, dyn, max_sweeps) for T in temperatures]
        imgdirs = [os.getcwd() + "\\" + "img" + "\\" + distinct for distinct in distincts]
        files = [imgdir + "\\" + "results.txt" for imgdir in imgdirs]
        arrays = [] # avg_M, avg_E, avg_M_err, avg_E_err, chi_true, chi_error, c_true, c_error, self.T
        for file in files:
            with open(file, 'rb') as f:
                loaded = pickle.load(f)
                arrays.append(loaded)
        arrays=np.array(arrays)
        labels = ['avg_M', 'avg_E', 'avg_M_err', 'avg_E_err', 'chi_true', 'chi_error', 'c_true', 'c_error', 'T']
        table = Table(arrays, names=labels)
        table.sort('T')

        fig, axs = plt.subplots(nrows=2, ncols=2, sharex='col', figsize=(15,7))
        plt.subplots_adjust(wspace=0.15, hspace=0)
        x = table['T']
        table['avg_M'] = np.abs(table['avg_M'])
        colors = ['blue', 'red', 'blue', 'red']
        y_index = ['avg_M', 'avg_E', 'chi_true', 'c_true']
        y_error = ['avg_M_err', 'avg_E_err', 'chi_error', 'c_error']
        x_label = "Temperature"
        y_labels = ["Magnetisation", "Energy", "Magnetic Susceptibility", "Specific Head p. Spin"]
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
                #if table_index == 0:
                #    ax.set(ylim=[-0.1,1.1])
        plt.savefig(dyn + "_multi.png", dpi=600)
        #plt.show()
        # Save table for people to look at.
        writer = hdfutils.hdf5_writer(os.getcwd(), "multirun_data.hdf5")
        writer.write_table(dyn, "results", table)

    # Snip all images in img together. This is faster than anigen by miles.
    def anigen_high(self):
        image_files = [os.path.join(self.imgdir, img)
                       for img in os.listdir(self.imgdir)
                       if img.endswith(".png")]
        image_files.sort(key=lambda f: int(re.sub('\D', '', f)))
        clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=self.framerate)
        clip.write_videofile(self.animation_name)

    # Generate (and save) animation using all_sweeps. Only applicable for small sims.
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

    # Delete storage file (not multiprocessed one.)
    def delete(self):
        try:
            os.remove(self.filename)
        except:
            print("Failed to delete file. Do it manually, mate.")

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
        legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                           Patch(facecolor='blue', label='-1', edgecolor='black')]
        plt.legend(handles=legend_elements, loc='upper right')
        cmap = cm.get_cmap('bwr', 2)
        plt.xlim([0, self.lx-1])
        plt.ylim([0, self.lx - 1])
        im = plt.imshow(self.mat, animated=True, cmap=cmap, aspect='equal')
        plt.clim(-1, 1)
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        t2 = plt.text(1, self.lx-5, str("T = ") + str(self.T), va='top', fontsize=20, color='white')

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")
        # Run sim
        if self.dyn == "g":
            while self.sweep < self.max_sweeps:
                self.glauber()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    if self.sweep % self.temprate == 0:
                        self.T += self.temprise
                        t2.set_text(("T = {0:.2f}").format(self.T))
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")
        if self.dyn == "k":
            while self.sweep < self.max_sweeps:
                self.kawasaki()
                if self.sweep % self.binrate == 0:
                    im.set_array(self.mat)
                    t1.set_text(str(self.sweep))
                    if self.sweep % self.temprate == 0:
                        self.T += self.temprise
                        t2.set_text(("T = {0:.2f}").format(self.T))
                    plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png")

    # Create an autocorrelation plot (to estimate autocorrelation time.) TODO: Figure out why time >> 10
    def autocorrelation(self):
        # Load the energy/magnetisation
        all_M, all_E = list(self.data['all_M']), list(self.data['all_E'])
        # Clip to measurements
        all_M, all_E = all_M[self.equilibration:self.equilibration + self.measurements],\
                       all_E[self.equilibration:self.equilibration + self.measurements]
        # Generate a correlation given
        pandas.plotting.autocorrelation_plot(all_M)
        plt.savefig("autocorrelation_" + self.distinct + ".png")

    # User Input
    def user_input(self):
        self.time_delay("Yōkoso!!! Welcome to this 2D Ising Simulation Suite. "
                        "Please be aware that there is a lot of customisation available."
                        "However, for brevity, here you should only provide a small few parameters."
                        "Then, the simulation will run, a live animation will be produced, and that's that."
                        "See the report document if you desire results. Without further ado... provide the parameters!")
        print()
        print()
        self.time_delay("Grid size. This is a square simulator, so just one integer will suffice.")
        self.lx = int(input())
        self.time_delay("Temperature (note that J=KB=unity in this model.")
        self.T = float(input())
        self.time_delay("Maximum number of sweeps to run for")
        self.max_sweeps = int(input())
        self.time_delay("Dynamical mode: g for glauber, k for kawasaki.")
        self.dyn = str(input())
        self.time_delay("Please enjoy!")

    # Regenerate the __init__ section with the new variables provided
    """
    This is for multiprocessing-related stuff where you want to change dynamics or temperature on the fly.
    Regenerates all the dependencies for these two parameters.
    """
    def init_multiprocess(self):
        # For storing/writing images
        imgdir = self.directory + "\\" + "img"
        # Unique string stating important parameters ***used as group for storage inside file.)
        self.distinct = ("{0}_{1:.2f}_{2}_{3}").format(self.lx, self.T, self.dyn, self.max_sweeps)
        self.imgdir = imgdir + "\\" + self.distinct
        try:
            os.mkdir(self.imgdir)
        except:
            pass
        self.animation_name = self.imgdir + "\\animation.mp4"

    # Save data to file.
    def save(self, table=None):
        if table == None:
            i = 0
            while i <= self.saveattempts:
                try:
                    self.writer.write_df(self.distinct, self.set, self.data)
                    break
                except:
                    i += 1
                    time.sleep(np.random.default_rng().integers(0, 30, 1)[0])
                    pass
        else:
            i = 0
            while i <= self.saveattempts:
                try:
                    self.writer.write_table(self.distinct, self.numset, table)
                    break
                except:
                    i += 1
                    time.sleep(np.random.default_rng().integers(0, 30, 1)[0])
                    pass

    # Load data from file
    def load(self, table=None):
        if table == None:
            self.data = self.writer.read_df(self.distinct, self.set)
        else:
            return self.writer.read_table(self.distinct, self.numset)

    # Old bit of code to make text in the console appear slower and crisper (2nd year???)
    def time_delay(self, text):
        if __name__ == "__main__":
            print()
            for c in text:
                sys.stdout.write(c)
                sys.stdout.flush()
                resttime = np.random.default_rng().uniform(0.001, self.delay_max)
                time.sleep(resttime)
        else:
            pass

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
            legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                               Patch(facecolor='blue', label='-1', edgecolor='black')]
            plt.legend(handles=legend_elements, loc='upper right')
            plt.clim(-1, 1)
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
            plt.clim(-1,1)
            legend_elements = [Patch(facecolor='red', label='+1', edgecolor='black'),
                               Patch(facecolor='blue', label='-1', edgecolor='black')]
            plt.legend(handles=legend_elements, loc='upper right')
            plt.show()
