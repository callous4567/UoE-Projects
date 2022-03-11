import itertools
import pickle
import sys
import time
from astropy.table import Table
from matplotlib import rc, cm
import os
from matplotlib.patches import Patch
import fast_sirs
import hdfutils
import numpy as np
import matplotlib.pyplot as plt


class twod_sirs(object):
    # TODO: make sure the current init is fine. Remember we added input instead of two separate inits for multirun.
    # Optional "identifier" exists to allow for multirun saving of many instances. Default to zero.
    def __init__(self, lx=None, p1=None, p2=None, p3=None, binrate=None,
                 equilibration=None, measurements=None, not_up=True, identifier=0, immune_frac=0):
        # This just lets you avoid dealing with __init__ if you want to run a method that doesn't need all this stuff.
        if lx!=None:
            # Initialization parameters.
            self.lx = lx
            self.p = np.array([p1, p2, p3])
            self.equilibration = equilibration # In number of flips.
            self.measurements = measurements # In number of flips.
            self.autocorrelation_length = 25000 # In number of flips.
            self.binrate = binrate
            # Derived self parameters.
            self.rng = np.random.default_rng()
            self.I = None
            if not_up == True:
                self.mat = self.rng.choice([False, True, 2], size=(self.lx, self.lx))  # Random S I or R matrix.
            else:
                self.mat = np.ones(shape=(self.lx, self.lx)) # all-up spin matrix (all infected.)
                # Go and set a fraction of them to recovered (for immunity)
                random_i = np.random.randint(0, lx, int(immune_frac*(lx**2)))
                random_j = np.random.randint(0, lx, int(immune_frac*(lx**2)))
                immu = np.zeros_like(self.mat)
                for i,j in zip(random_i,random_j):
                    immu[i,j] = True
                self.immu = immu

            self.sweep = 0 # in NUMBER OF FLIPS
            self.max_sweeps = self.equilibration + self.measurements
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
                self.distinct = ("{0}_{1:.4f}_{2:.4f}_{3:.4f}_{4}_{5}").format(self.lx, p1, p2, p3, self.max_sweeps, identifier)
                self.imgdir = imgdir + "\\" + self.distinct
                try:
                    os.mkdir(self.imgdir)
                except:
                    pass
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
        legend_elements = [Patch(facecolor='blue', label='S', edgecolor='black'),
                           Patch(facecolor='white', label='I', edgecolor='black'),
                           Patch(facecolor='red', label='R', edgecolor='black')]
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

        # Create datasave format and save if applicable.
        self.all_I = all_I
        if checkpoint==False:
            # Save it
            self.save()

    # Iterate Sequential (i.e. monte-carlo the points one at a time.)
    def fast_sequential(self):
        self.mat, self.I, self.sweep = self.fast.fast_sequential(self.mat, self.I, self.sweep)

    # Iterate Sequential (i.e. monte-carlo the points one at a time.) but with immu
    def fast_immu(self):
        self.mat, self.I, self.sweep = self.fast.fast_sequential_immu(self.mat, self.I, self.sweep, self.immu)

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
        breakpoint = self.max_sweeps
        # Run Simulation. Will run until self.max_sweeps = self.equilibration + self.measurements. start at 0.
        while self.sweep < self.max_sweeps:
            self.fast_immu()
            all_I.append(self.I)
            # If absorbing, then just put "I=0" for the rest (saves computational effort.)
            if self.I == 0:
                empties = list(np.zeros(self.max_sweeps - self.sweep))
                all_I += empties
                #breakpoint = self.sweep
                break

        # Create datasave format
        self.all_I = np.array(all_I)
        #self.save()

        # Save a graph, if graph is true
        #fig = plt.figure()
        #plt.plot(np.arange(0, 12000, 1), all_I[0:12000], color='red', lw=0.5)
        #plt.title(self.distinct)
        #plt.savefig(self.imgdir + "\\" + "I_of_t_plot.png", dpi=72)
        #plt.close()

    # To generate runs for graphs/etc. Set run=False to re-generate averages but not do anything else.
    def main_multi(self):
        self.run_multi()
        array_dump = self.fast_averages_errors()
        return array_dump

    # Averages/errors. Note that we're combining several runs (0, identifier) and consequently need to return.
    def fast_averages_errors(self):
        # Calculate averages and errors
        avg_I, \
        avg_I_err, \
        chi_true, \
        chi_error = self.fast.averages_errors(self.all_I,
                                              self.equilibration,
                                              self.measurements,
                                              self.autocorrelation_length)
        # Only save if not multiprocessing. Saving disabled for now.
        if __name__ == "__main__":
            # Produce a table and save it.
            table = Table()
            table["I"], table["chi"] = np.array([avg_I, avg_I_err]), np.array([chi_true, chi_error])
            self.save(table)
        else:
            array_dump = np.array([avg_I, avg_I_err, chi_true, chi_error, self.p[0], self.p[1], self.p[2]])
            return array_dump
            #with open(self.imgdir + "\\" + "results.txt", 'wb') as f:
            #    pickle.dump(obj=array_dump, file=f)


    # Produce graphs for multiprocessed runs selected.
    def multigraph(self):
        # Sim Params.
        lx = 50
        equilibration = 100*int(lx**2)
        measurements = 1000*int(lx**2)
        max_sweeps = equilibration + measurements

        # Generate the parameter space for the runs (for probability.)
        one_range = np.linspace(0, 1, 20)
        two_range = np.array([0.5])  # set p2 = 0.5
        tre_range = np.linspace(0, 1, 20)
        zipped = list(itertools.product(one_range, two_range, tre_range))

        # Generate file ids
        rootdir = os.getcwd()
        savedir = "part_3"
        unique_strings = [("{0:.3f}_{1:.3f}_{2:.3f}").format(p1, p2, p3) for p1,p2,p3 in zipped]
        fileids = [rootdir + "\\" + savedir + "\\" + unique_string + ".txt" for unique_string in unique_strings]
        arrays = [] # avg_I, error, chi, error (note that I and error are not normalized, but chi and error are per cell
        for file in fileids:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass

        # Arrays is as a vector at the moment.
        # Create an empty, with three extra columns
        array_new = np.empty((len(arrays), 7))
        array_new[:, 0:4] = np.array(arrays)
        # Put probability in, too
        p1, p2, p3 = np.array(zipped).T
        array_new[:, 4], array_new[:, 5], array_new[:, 6] = p1,p2,p3

        # Create table
        labels = ['avg_I', 'avg_I_err', 'chi_true', 'chi_error', 'p1', 'p2', 'p3']
        table = Table(data=array_new, names=labels) # table that has all the data.

        # Dump the table.
        writer = hdfutils.hdf5_writer(os.getcwd(), "datafile.hdf5")
        writer.write_table("part_3_table", "table", table)

        # Divide by the size of the array for I and I error
        table['avg_I'] /= lx**2
        table['avg_I_err'] /= lx**2

        # Get contour plot data
        XYZ = ['p1','p3','chi_true']
        X, Y, Z = [table[d] for d in XYZ]

        # Create the plot
        fig, ax1 = plt.subplots(nrows=1,ncols=1)
        ax1.tricontour(X, Y, Z, levels=15, linewidths=0.5, colors='k')
        cntr_fill = ax1.tricontourf(X, Y, Z, levels=15, cmap="RdBu_r")
        fig.colorbar(cntr_fill, ax=ax1, label=r'$<\chi>$')
        ax1.set(xlabel=XYZ[0], ylabel=XYZ[1])
        ax1.set(xlim=[0,1], ylim=[0,1])
        plt.savefig(savedir + "average_chi.png", dpi=300)
        plt.show()

    # Produce graphs for multiprocessed runs selected. For part 4,
    def multigraph_line(self):
        # Sim Params.
        lx = 50
        equilibration = 100*int(lx**2)
        measurements = 10000*int(lx**2)
        max_sweeps = equilibration + measurements

        # Generate the parameter space for the runs (for probability.)
        one_range = np.linspace(0.2, 0.5, 20)
        two_range = [0.5 for d in one_range]  # set p2 = 0.5
        tre_range = two_range  # set p3 = 0.5, too
        zipped = list(itertools.product(one_range, two_range, tre_range))

        # Generate file ids
        rootdir = os.getcwd()
        savedir = "part_4" # all of it for p2 = p3 = 0.5 and p1 variable.
        unique_strings = [("{0:.3f}_{1:.3f}_{2:.3f}").format(p1, p2, p3) for p1,p2,p3 in zipped]
        fileids = [rootdir + "\\" + savedir + "\\" + unique_string + ".txt" for unique_string in unique_strings]
        arrays = [] # avg_I, error, chi, error (note that I and error are not normalized, but chi and error are per cell
        for file in fileids:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass

        # Arrays is as a vector at the moment.
        # Create an empty, with three extra columns
        array_new = np.empty((len(arrays), 7))
        array_new[:, 0:4] = np.array(arrays)
        # Put probability in, too
        p1, p2, p3 = np.array(zipped).T
        array_new[:, 4], array_new[:, 5], array_new[:, 6] = p1,p2,p3

        # Create table
        labels = ['avg_I', 'avg_I_err', 'chi_true', 'chi_error', 'p1', 'p2', 'p3']
        table = Table(data=array_new, names=labels) # table that has all the data.

        # Divide by the size of the array for I and I error
        table['avg_I'] /= lx**2
        table['avg_I_err'] /= lx**2
        # Dump the table.
        writer = hdfutils.hdf5_writer(os.getcwd(), "datafile.hdf5")
        writer.write_table("part_4_table", "table", table)

        # We're not contouring this time. We're doing a line plot.
        fig = plt.figure()
        plt.errorbar(table['p1'], table['chi_true'], table['chi_error'], color='blue', ecolor='red')
        plt.xlim([0.2,0.5])
        plt.xlabel(r'$p_1$')
        plt.ylabel(r'$\chi$')
        plt.grid(which='major', color='pink')
        plt.savefig(rootdir + "\\" + "_" + savedir + "_chiplot.png", dpi=300)
        plt.show()

    # Produce graphs for multiprocessed runs selected. For part 4,
    def multigraph_immu(self):
        # Sim Params.
        lx = 50
        equilibration = 100*int(lx**2)
        measurements = 1000*int(lx**2)
        max_sweeps = equilibration + measurements

        # Generate the parameter space for the runs (for probability.)
        one_range = np.array([0.5])
        two_range = np.array([0.5])  # set p2 = 0.5
        tre_range = np.array([0.5])  # set p3 = 0.5, too
        immu_range = np.linspace(0, 1, 20)
        zipped = list(itertools.product(one_range, two_range, tre_range, immu_range))

        # Generate file ids
        rootdir = os.getcwd()
        savedir = "part_5" # all of it for p2 = p3 = 0.5 and p1 variable.
        unique_strings = [("{0:.3f}").format(immu) for p1,p2,p3,immu in zipped]
        fileids = [rootdir + "\\" + savedir + "\\" + unique_string + ".txt" for unique_string in unique_strings]
        arrays = [] # avg_I, error, chi, error (note that I and error are not normalized, but chi and error are per cell
        for file in fileids:
            try:
                with open(file, 'rb') as f:
                    loaded = pickle.load(f)
                    arrays.append(loaded)
            except:
                pass

        # Arrays is as a vector at the moment.
        arrays = np.array(arrays)

        # Get immus
        zipped = np.array(zipped).T
        immus = zipped[3]

        # Create table
        labels = ['avg_I', 'avg_I_err', 'chi_true', 'chi_error']
        table = Table(data=arrays, names=labels) # table that has all the data.
        table['immus'] = immus

        # Divide by the size of the array for I and I error
        table['avg_I'] /= lx**2
        table['avg_I_err'] /= lx**2
        # Dump the table.
        writer = hdfutils.hdf5_writer(os.getcwd(), "datafile.hdf5")
        writer.write_table("part_5_table", "table", table)

        # We're not contouring this time. We're doing a line plot.
        fig = plt.figure()
        plt.errorbar(table['immus'], table['avg_I'], table['avg_I_err'], color='blue', ecolor='red')
        plt.xlabel('R')
        plt.ylabel(r'$\langle\psi\rangle$')
        plt.grid(which='major', color='pink')
        plt.savefig(rootdir + "\\" + "_" + savedir + "immuplot.png", dpi=300)
        plt.show()

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
        p1 = float(input())
        print()
        self.time_delay("Probability of I to R (p2).")
        p2 = float(input())
        print()
        self.time_delay("Probability of R to S (p3).")
        p3 = float(input())
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

checkpoint().run()