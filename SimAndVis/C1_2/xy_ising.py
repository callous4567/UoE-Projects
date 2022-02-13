import pickle
import time

import astropy.convolution
import matplotlib
import numba.core.types
from astropy.table import Table
from matplotlib import rc, cm
import os
from matplotlib.patches import Patch
from fast_xy import fast_xy, fast_glauber, fast_clustglauber
import hdfutils
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.pyplot as plt
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
class xy_ising(object):
    # TODO: make sure the current init is fine. Remember we added input instead of two separate inits for multirun.
    def __init__(self, T=None):
        # Initialization parameters.
        self.lx = 50
        if T != None:
            self.T = T
        else:
            self.T=1
        self.binrate = 1  # 1e6 # binrate for animations. Every i'th sweep is selected as frames.
        self.rng = np.random.default_rng()
        self.M, self.E = None, None
        self.sweep = 0 # in NUMBER OF FLIPS
        self.max_sweeps = int(1e6)
        self.directory = os.getcwd()
        self.filename = "datafile.hdf5"
        self.set = "data"
        self.numset = "results"
        not_up = False # True if you want random False if else
        if not_up == True:
            mat = np.empty((self.lx, self.lx, 2))
            for i in range(self.lx):
                for j in range(self.lx):
                    random_vector = np.random.normal(size=2)
                    random_vector = random_vector/np.linalg.norm(random_vector)
                    mat[i, j] = random_vector
            self.mat = mat
        else:
            # Generate a uniform initial state for the XY model
            mat = np.empty((self.lx, self.lx, 2))
            for i in range(self.lx):
                for j in range(self.lx):
                    mat[i, j] = np.array([0,1])
            self.mat = mat
        imgdir = self.directory + "\\" + "img"
        try:
            os.mkdir(imgdir)
        except:
            pass
        self.distinct = ("{0}_{1:.4f}_{2}_{3}").format(self.lx, self.T, "g", self.max_sweeps)
        self.imgdir = imgdir + "\\" + self.distinct
        try:
            os.mkdir(self.imgdir)
        except:
            pass
        self.all_M, self.all_E = None, None
        self.saveattempts = 5
        self.fast = fast_xy(self.lx)
        self.convolve_radius = 2
        self.temprate = int(300000000)
        self.temprise = 0.1
        self.autocor = int(1e3)

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
        plt.xlim([0, self.lx-1])
        plt.ylim([0, self.lx - 1])
        # Set up the X,Y coordinates for each and every point in our array
        #ax = plt.quiver(self.mat[:,:,0], self.mat[:,:,1], headlength=16, scale=128, headwidth=8)
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        plt.title(("T = {} M = {} E = {}").format(self.T,self.M,self.E))
        # Compute the divergence field (first order) to start with
        #convolution_kernel = astropy.convolution.Gaussian2DKernel(self.convolve_radius,
        #                                                          self.convolve_radius,
        #                                                          theta=0)
        #div_field = self.fast.divergence_first(self.mat)
        #col_field = astropy.convolution.convolve_fft(div_field, kernel=convolution_kernel)
        col_field = self.fast.angle(self.mat)
        im = plt.imshow(col_field, cmap='plasma', interpolation='nearest')#'bwr')
        plt.colorbar(label="Phase")

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        # Run sim
        while self.sweep < self.max_sweeps:
            self.fast_clustglauber()
            if self.sweep % self.binrate == 0:
                #ax.set_UVC(self.mat[:,:,0], self.mat[:,:,1])
                t1.set_text(str(self.sweep))
                plt.title(("T = {} M = {} E = {}").format(self.T, self.M, self.E))
                #div_field = self.fast.divergence_first(self.mat)
                #col_field = astropy.convolution.convolve_fft(div_field, kernel=convolution_kernel)
                col_field = self.fast.angle(self.mat)
                im.set_array(col_field)
                if self.sweep % self.temprate == 0:
                    self.T += self.temprise
                plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)

    # Save data.
    """
    all_M and all_E are saved to data.txt as a pickled numpy array np.array([all_M, all_E])
    Deprecated format is pandas DataFrame (much too slow) so some code may reflect this. 
    """
    def save(self):
        i = 0
        while i <= self.saveattempts:
            try:
                # Dump it
                with open(self.imgdir + "\\" + "data.txt", 'wb') as f:
                    pickle.dump(obj=np.array([self.all_E]), file=f)
                break
            except:
                i += 1
                time.sleep(self.rng.integers(0, 30, 1)[0])
                pass

    # Produce graphs for multiprocessed runs selected.
    def multigraph(self):
        max_sweeps = int(1e6)
        lx = 50
        dyn = 'g'
        # Check dynamics type. Note: make sure all temps attached here are same as multirun.py.
        if dyn == 'g':
            temperatures = np.linspace(0.5, 2, 30)
        distincts = [("{0}_{1:.4f}_{2}_{3}").format(lx, T, dyn, max_sweeps) for T in temperatures]
        imgdirs = [os.getcwd() + "\\" + "img" + "\\" + distinct for distinct in distincts]
        files = [imgdir + "\\" + "results.txt" for imgdir in imgdirs]
        arrays = []  # avg_M, avg_E, avg_M_err, avg_E_err, chi_true, chi_error, c_true, c_error, self.T
        # Grab dumps and make table.
        for file in files:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass
        arrays = np.array(arrays)
        labels = ["avg_E", "avg_E_err", "avg_C", "avg_C_err", 'T']
        table = Table(arrays, names=labels)
        table.sort('T')

        # Plot
        fig = plt.figure()
        plt.plot(table['T'], table['avg_C'], color='red')
        plt.title("Specific Heat/Spin Fluctuation vs. Temperature")
        plt.xlabel("Temperature")
        plt.ylabel("C/Spin")

        # Dump figure
        plt.savefig(dyn + "_multi.png", dpi=600)
        plt.show()

        # Save the table so that people can have a gander at it.
        writer = hdfutils.hdf5_writer(os.getcwd(), "multirun_data.hdf5")
        writer.write_table(dyn, "results", table)

    # For multiprocessed version: no thrills attached.
    def run_multi(self):
        # Calculate M and E for initial step
        self.M, self.E = self.fast.fast_magenergy(self.mat)
        all_M, all_E = [self.M], [self.E]

        # Run Sim
        while self.sweep < self.max_sweeps:
            self.fast_clustglauber()
            all_M.append(self.M), all_E.append(self.E)

        # Create datasave format
        self.all_E = np.array(all_E)

        # Trim it (for equilibration/etc.)
        self.all_E = self.all_E[int(1e5):]

        # Save it
        self.save()
        self.fast_averages_errors()

    # Averages and Errors, but fast_ising.
    def fast_averages_errors(self):
        # Calculate averages and errors
        avg_E, avg_E_err, avg_C, avg_C_err = self.fast.averages_errors(self.all_E, self.T, self.autocor)

        array_dump = np.array([avg_E, avg_E_err, avg_C, avg_C_err, self.T])
        with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
            pickle.dump(obj=array_dump, file=f)

    """
    Estimated run time: 
    3.6 GHz Single Thread
    160x10^6 flips on an 800x800 
    7 hours (Glauber) 8 hours (Kawasaki) 
    """
    def run_high_streamplot(self):
        # Define initial condition
        self.M, self.E = self.fast.fast_magenergy(self.mat)
        # Various fire prelim stuff
        fig = plt.figure(figsize=(16,16))
        # Create x,y
        x, y = np.arange(0, self.lx, 1), np.arange(0, self.lx, 1)
        # Set up image plot
        plt.xlim([0, self.lx-1])
        plt.ylim([0, self.lx - 1])
        # Set up the X,Y coordinates for each and every point in our array
        ax = plt.streamplot(x, y, self.mat[:,:,0], self.mat[:,:,1], color='black', density=self.lx/30)
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        plt.title(("T = {}").format(self.T))

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        # Run sim
        while self.sweep < self.max_sweeps:
            self.fast_glauber()
            if self.sweep % self.binrate == 0:
                plt.cla()
                ax = plt.streamplot(x, y, self.mat[:, :, 0], self.mat[:, :, 1], color='black', density=self.lx / 30)
                t1.set_text(str(self.sweep))
                plt.title(("T = {}").format(self.T))
                if self.sweep % self.temprate == 0:
                    self.T += self.temprise
                plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)

    # Glauber, but using fast_ising.
    def fast_glauber(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_glauber(self.mat, self.M, self.E,
                                                                              self.T, self.sweep)

    # Glauber, but using fast_ising. Clustered.
    def fast_clustglauber(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_clustglauber(self.mat, self.M, self.E,
                                                                                   self.T,self.sweep)

    # Testing fastmath
    def fast_glauber_fastmath(self):
        self.mat, self.M, self.E, self.T, self.sweep = fast_glauber(self.lx, self.mat, self.M, self.E,
                                                                              self.T, self.sweep)

    # Glauber, but using fast_ising. Clustered.
    def fast_clustglauber_fastmath(self):
        self.mat, self.M, self.E, self.T, self.sweep = fast_clustglauber(self.lx, self.mat, self.T,self.sweep)

