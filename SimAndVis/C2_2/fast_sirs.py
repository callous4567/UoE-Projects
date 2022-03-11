import numpy as np
from numba.experimental import jitclass
from numba import types

# Some notes.
"""
Here we use the values:
0 SUSCEPTIBLE (False)
1 INFECTED (True)
2 RECOVERED (Else) 
Here we use the Ising Definition of NN's (only 4- TOP BOTTOM LEFT RIGHT.) 
In testing, using the "Parallel Update" method is much MUCH faster than sequential updates. 
By a factor of ~200.
"""
variable_specification = [
    ('good_luck', types.unicode_type),
    ('lx', types.int32),
    ('p1', types.float64),
    ('p2', types.float64),
    ('p3', types.float64),
    ('I', types.int32),
    ('sweep', types.int32),
    ('immu', types.int32[:,:]),
]
@jitclass(variable_specification)
class fast_sirs():
    def __init__(self, lx, p1, p2, p3):
        good_luck = "Godspeed, mate. It's the second round!"
        self.lx = lx
        self.p1, self.p2, self.p3 = p1, p2, p3
    # Periodic boundary conditions.
    def pbc(self, i):
        if i >= self.lx:
            return 0
        if i < 0:
            return self.lx - 1
        else:
            return i
    # "Glauber" version of SIRS (Random Sequential Updates.)
    def fast_sequential(self, mat, I, sweep):
        # Pick a spot
        i, j = np.random.randint(0, self.lx, 2)
        # If susceptible
        if mat[i,j] == False:
            possible_infections = np.array([mat[self.pbc(i - 1), j],
                                            mat[self.pbc(i + 1), j],
                                            mat[i, self.pbc(j - 1)],
                                            mat[i, self.pbc(j + 1)]])
            # If a neighbour is infected, consider infecting it.
            if True in possible_infections:
                if np.random.random(1)[0] <= self.p1:
                    mat[i,j] = True
                    I += int(1)
        # If infected
        elif mat[i,j] == True:
            # Consider recovery
            if np.random.random(1)[0] <= self.p2:
                mat[i, j] = int(2)
                I -= int(1)
        # If recovered
        else:
            # Consider transition to susceptible
            if np.random.random(1)[0] <= self.p3:
                mat[i, j] = False
        # Sweep complete
        sweep += int(1)
        return mat, I, sweep
    # "Glauber" version of SIRS (Random Sequential Updates) with Immunity
    def fast_sequential_immu(self, mat, I, sweep, immu):
        # Pick a spot
        i, j = np.random.randint(0, self.lx, 2)
        # If susceptible
        if mat[i,j] == False:
            possible_infections = np.array([mat[self.pbc(i - 1), j],
                                            mat[self.pbc(i + 1), j],
                                            mat[i, self.pbc(j - 1)],
                                            mat[i, self.pbc(j + 1)]])
            # If a neighbour is infected, consider infecting it.
            if True in possible_infections:
                if np.random.random(1)[0] <= self.p1:
                    mat[i,j] = True
                    I += int(1)
        # If infected
        elif mat[i,j] == True:
            # Consider recovery
            if np.random.random(1)[0] <= self.p2:
                mat[i, j] = int(2)
                I -= int(1)
        # If recovered
        else:
            # Verify it isn't immune:
            if immu[i,j] == True:
                pass
            # Not immune, hence
            else:
                # Consider transition to susceptible
                if np.random.random(1)[0] <= self.p3:
                    mat[i, j] = False
        # Sweep complete
        sweep += int(1)
        return mat, I, sweep
    # Do the entire lattice at once.
    # "Parallel" version of SIRS (Entire lattice at once.)
    def fast_parallel(self, mat, I, sweep):
        # Define a new matrix/I
        new_mat = np.copy(mat)
        I_change = 0
        # Go over the old matrix all-cells and put new values in new_mat
        for i in range(self.lx):
            for j in range(self.lx):
                # If susceptible
                if mat[i, j] == False:
                    possible_infections = np.array([mat[self.pbc(i - 1), j],
                                                    mat[self.pbc(i + 1), j],
                                                    mat[i, self.pbc(j - 1)],
                                                    mat[i, self.pbc(j + 1)]])
                    # If a neighbour is infected, consider infecting it.
                    if True in possible_infections:
                        if np.random.random(1)[0] <= self.p1:
                            new_mat[i, j] = True
                            I_change += int(1)
                # If infected
                elif mat[i, j] == True:
                    # Consider recovery
                    if np.random.random(1)[0] <= self.p2:
                        new_mat[i, j] = int(2)
                        I_change -= int(1)
                # If recovered
                else:
                    # Consider transition to susceptible
                    if np.random.random(1)[0] <= self.p3:
                        new_mat[i, j] = False
        # Sweep complete
        sweep += int(1)
        return new_mat, I + I_change, sweep
    # Calculate infection per cell
    def fast_infsum(self, mat):
        sum = 0
        for i in range(self.lx):
            for j in range(self.lx):
                if mat[i,j] == True:
                    sum += 1
        return sum
    # Averages and errors calculation, but faster!
    def averages_errors(self, all_I, equil, measure, autocor):
        # Trim data to measurement range
        all_I = all_I[equil:equil+measure]
        # Sample every autocor'th measurement
        all_I = all_I[0::autocor+1]
        num_samples = len(all_I)
        # Get squares
        all_II = all_I**2
        # Get averages of uncorrelated measurements of I
        avg_I = np.mean(all_I)
        avg_II = np.mean(all_II)
        # Get Errors. Note sample length = Autocorrelation length, so cancels out.
        avg_I_err = np.sqrt(((avg_II - avg_I**2)*(2))/(num_samples))
        # Estimate infection per site average variance (susceptibility.)
        chi_true = (1/self.lx**2)*(avg_II - avg_I**2)
        # Error estimation for chi and c via the Bootstrap method
        number_of_resamples = 2000
        chi_list = np.empty(number_of_resamples)
        for i in range(number_of_resamples):
            # Select num_samples random from num_samples
            resample = np.random.randint(0, num_samples, num_samples)
            # Grab Magnetism/Energy samples
            Qall_I = all_I[resample]
            Qall_II = Qall_I ** 2
            # Get averages
            Qavg_I = np.mean(Qall_I)
            Qavg_II = np.mean(Qall_II)
            # Estimate susceptibility and specific heat/spin
            Qchi = (1/self.lx**2)*(Qavg_II - Qavg_I**2)
            # Append
            chi_list[i] = Qchi
        chi_average = np.mean(chi_list)
        chichi_average = np.mean(chi_list**2)
        boot_chi = np.sqrt(chichi_average - chi_average**2)
        return avg_I, avg_I_err, chi_true, boot_chi # Average value of I for the array, error, Variance in (I/N), error
