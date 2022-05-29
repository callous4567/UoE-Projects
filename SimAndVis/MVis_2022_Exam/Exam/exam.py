import pickle
import sys
import time

import numpy as np
from matplotlib import cm, pyplot as plt, colors
from matplotlib.patches import Patch

from SimAndVis.Exam.fast_exam import fast_model, typefield, get_typefraction


class model(object):
    def __init__(self, dx=None, dt=None, nx=None, D=None, q=None, p=None, max_sweeps=None, binrate=None, filedirname=None):

        # Holder for variables
        self.nx, self.D, self.q, self.p = nx, D, q, p
        self.dx, self.dt = dx,dt

        # RNG (faster.)
        self.rng = np.random.default_rng()

        # Matrices a,b,c
        self.a, self.b, self.c = self.initial_condition_random()

        # The extra array d (which is held in memory for faster plotting/less calculation + iterated forward.)
        self.d = 1 - (self.a + self.b + self.c)

        # Keep track of sweeps
        self.sweep = 0

        # Set max number before auto-terminate
        self.max_sweeps = max_sweeps

        # For plotting
        self.binrate = binrate

        # For saving/etc
        self.filedirname = filedirname
        self.savedata = None

    def initial_condition_random(self):
        a,b,c = self.rng.uniform(0, 1/3, (self.nx,self.nx)), \
                self.rng.uniform(0, 1/3, (self.nx,self.nx)),\
                self.rng.uniform(0, 1/3, (self.nx,self.nx))
        return a,b,c

    def fast_model(self):
        self.a, self.b, self.c, self.d, self.sweep = fast_model(self.a, self.b, self.c, self.d,
                                                                self.D, self.q, self.p,
                                                                self.dx, self.dt,
                                                                self.sweep)
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

    def run(self):

        # Custom colormap (gray/R/G/B - d/a/b/c)
        cmap = colors.ListedColormap(['gray','r','g','b'])
        bounds=[-0.5,0.5,1.5,2.5,3.5]
        norm = colors.BoundaryNorm(bounds, cmap.N)

        # Interactive On
        plt.ion()

        # Set up figure, axes, etc
        fig, ax = plt.subplots(figsize=(8,8))
        im = ax.imshow(typefield(self.a, self.b, self.c, self.d), animated=True, cmap=cmap, norm=norm, aspect='equal')
        ax.set(xlim = [0, self.nx - 1],
               ylim = [0, self.nx - 1])

        # Set up image plot
        legend_elements = [Patch(facecolor='red', label='a', edgecolor='black'),
                           Patch(facecolor='green', label='b', edgecolor='black'),
                           Patch(facecolor='blue', label='c', edgecolor='black'),
                           Patch(facecolor='grey', label='d', edgecolor='black')]
        ax.legend(handles=legend_elements, loc='upper right')

        # Set title.
        t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)
        # Run Simulation.
        start = time.time()
        mpl = 0
        self.time_delay("Starting Simulation...")
        while self.sweep < self.max_sweeps:
            self.fast_model()
            if self.sweep % self.binrate == 0:
                mins = time.time()
                im.set_array(typefield(self.a, self.b, self.c, self.d))
                t1.set_text(str(self.sweep))
                fig.canvas.draw()
                fig.canvas.flush_events()
                ens = time.time()
                mpl += (ens-mins)
        end = time.time()

        # All done.
        self.time_delay(("Simulation finished. Took {0:.1f} seconds. "
                         " Plotting stole approximately {1:.1f} seconds.").format(end - start, mpl))
        plt.close()

    def type_over_time(self):
        """
        Run the simulation until max_sweeps
        Get the type fractions over time of the simulation
        Returns array [[a1, a2, a3, ...], [b1, b2, b3, ...], ...] for a,b,c
        Runs until max_sweeps.
        """
        typefractions = []
        while self.sweep < self.max_sweeps:
            if self.sweep == int(self.max_sweeps/4):

                # Custom colormap (gray/R/G/B - d/a/b/c)
                cmap = colors.ListedColormap(['gray', 'r', 'g', 'b'])
                bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
                norm = colors.BoundaryNorm(bounds, cmap.N)

                # Set up figure, axes, etc
                fig, ax = plt.subplots(figsize=(8, 8))
                im = ax.imshow(typefield(self.a, self.b, self.c, self.d), animated=True, cmap=cmap, norm=norm,
                               aspect='equal')
                ax.set(xlim=[0, self.nx - 1],
                       ylim=[0, self.nx - 1])

                # Set up image plot
                legend_elements = [Patch(facecolor='red', label='a', edgecolor='black'),
                                   Patch(facecolor='green', label='b', edgecolor='black'),
                                   Patch(facecolor='blue', label='c', edgecolor='black'),
                                   Patch(facecolor='grey', label='d', edgecolor='black')]
                ax.legend(handles=legend_elements, loc='upper right')

                # Set title.
                t1 = ax.text(1, 1, str(self.sweep), color="white", fontsize=20)

                # Save
                plt.savefig("sweep_5000_example_Q2", dpi=300)
            typefraction = get_typefraction(typefield(self.a, self.b, self.c, self.d))
            typefractions.append(typefraction)
            self.fast_model()
        return np.array(typefractions).T

    def run_until_absorption(self):

        """
        Run the simulation until one of the typefractions np.isclose to unity
        Will return boolean False is failed to reach this point within time
        Else will return the time T at which we reached this "absorbing" state
        :return:
        """
        ones = np.ones(3)
        while self.sweep < self.max_sweeps:
            typefraction = get_typefraction(typefield(self.a, self.b, self.c, self.d))
            if True in np.isclose(ones, typefraction):
                break
            self.fast_model()
        if self.sweep >= self.max_sweeps - 10:
            return False
        else:
            return self.sweep*self.dt

    def a_over_time_twopoints(self, coord_1, coord_2):
        as_one, as_two = [],[]
        while self.sweep < self.max_sweeps:
            as_one.append(self.a[coord_1[0],coord_1[1]])
            as_two.append(self.a[coord_2[0],coord_2[1]])
            self.fast_model()
        return as_one, as_two

    def type_over_time_twopoints(self, coord_1, coord_2):
        types_1, types_2 = [],[]
        while self.sweep < self.max_sweeps:
            typematrix = typefield(self.a, self.b, self.c, self.d)
            types_1.append(typematrix[coord_1[0],coord_1[1]])
            types_2.append(typematrix[coord_2[0],coord_2[1]])
            self.fast_model()
        return types_1, types_2



    def save(self):
        with open(self.filedirname, "wb") as f:
            pickle.dump(obj=self.savedata, file=f)

    def load(self):
        with open(self.filedirname, "rb") as f:
            self.savedata = pickle.load(file=f)

    def multigraph_first(self):
        none = None

