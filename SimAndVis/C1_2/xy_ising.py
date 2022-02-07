import astropy.convolution
from matplotlib import rc, cm
import os
from matplotlib.patches import Patch
from fast_xy import fast_xy
import hdfutils
import numpy as np
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
    def __init__(self):
        # Initialization parameters.
        self.lx = 800
        self.T = 0.2
        self.binrate = int(0.5e6)  # 1e6 # binrate for animations. Every i'th sweep is selected as frames.
        self.rng = np.random.default_rng()
        self.M, self.E = None, None
        self.sweep = 0 # in NUMBER OF FLIPS
        self.max_sweeps = int(400e6)
        self.directory = os.getcwd()
        self.filename = "datafile.hdf5"
        self.set = "data"
        self.numset = "results"
        not_up = True # True if you want random False if else
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
        self.temprate = int(50e6)
        self.temprise = 0.5
        self.all_M, self.all_E = None, None
        self.saveattempts = 5
        self.fast = fast_xy(self.lx)
        self.convolve_radius = 2


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
        ax = plt.quiver(self.mat[:,:,0], self.mat[:,:,1], headlength=16, scale=32, headwidth=8)
        t1 = plt.text(1, 1, str(self.sweep), color="white", fontsize=20)
        plt.title(("T = {} M = {} E = {}").format(self.T,self.M,self.E))
        # Compute the divergence field (first order) to start with
        convolution_kernel = astropy.convolution.Gaussian2DKernel(self.convolve_radius,
                                                                  self.convolve_radius,
                                                                  theta=0)
        div_field = self.fast.divergence_first(self.mat)
        col_field = astropy.convolution.convolve_fft(div_field, kernel=convolution_kernel)
        im = plt.imshow(col_field, cmap='bwr')

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        # Run sim
        while self.sweep < self.max_sweeps:
            self.fast_glauber()
            if self.sweep % self.binrate == 0:
                ax.set_UVC(self.mat[:,:,0], self.mat[:,:,1])
                t1.set_text(str(self.sweep))
                plt.title(("T = {} M = {} E = {}").format(self.T, self.M, self.E))
                div_field = self.fast.divergence_first(self.mat)
                col_field = astropy.convolution.convolve_fft(div_field, kernel=convolution_kernel)
                im.set_array(col_field)
                if self.sweep % self.temprate == 0:
                    self.T += self.temprise
                plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
    # High-scale run (no saving) that will output images to imgdir rather than keeping in memory/writing.

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
        plt.title(("T = {} M = {} E = {}").format(self.T,self.M,self.E))

        # Preliminary Save
        plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
        # Run sim
        while self.sweep < self.max_sweeps:
            self.fast_glauber()
            if self.sweep % self.binrate == 0:
                plt.cla()
                ax = plt.streamplot(x, y, self.mat[:, :, 0], self.mat[:, :, 1], color='black', density=self.lx / 30)
                t1.set_text(str(self.sweep))
                plt.title(("T = {} M = {} E = {}").format(self.T, self.M, self.E))
                if self.sweep % self.temprate == 0:
                    self.T += self.temprise
                plt.savefig(self.imgdir + "\\" + str(self.sweep) + ".png", dpi=300)
    # Glauber, but using fast_ising.
    def fast_glauber(self):
        self.mat, self.M, self.E, self.T, self.sweep = self.fast.fast_glauber(self.mat,
                                                                              self.M,
                                                                              self.E,
                                                                              self.T,
                                                                              self.sweep)

uwu = xy_ising()
uwu.run_high()









