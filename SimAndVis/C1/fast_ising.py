import numpy as np
from numba.experimental import jitclass
from numba import types

variable_specification = [
    ('good_luck', types.unicode_type),
    ('lx', types.int32)
]
@jitclass(variable_specification)
class fast_ising():
    def __init__(self, lx):
        good_luck = "Godspeed, mate. It actually works!"
        self.lx = lx
    def pbc(self, i):
        if i >= self.lx:
            return 0
        if i < 0:
            return self.lx - 1
        else:
            return i

    # Glauber Dynamics Stepping
    def fast_glauber(self, mat, M, E, T, sweep):
        i, j = np.random.randint(0, self.lx, 2)
        nnsum = mat[self.pbc(i - 1), j] + mat[self.pbc(i + 1), j] + mat[i, self.pbc(j - 1)] + mat[i, self.pbc(j + 1)]
        # Calculate change in magnetisation (mag new - mag old)
        mag_cost = -2 * mat[i, j]
        # Generate an energy cost for doing this flip (energy new - energy old)
        energy_cost = -1 * mag_cost * nnsum
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            mat[i, j] *= -1
            # Add change in mag/energy
            M += mag_cost
            E += energy_cost
        else:
            probability = np.exp(-1 * energy_cost / T)
            if np.random.random(1)[0] <= probability:
                mat[i, j] *= -1
                # Add change in mag/energy
                M += mag_cost
                E += energy_cost
        # Sweep complete.
        sweep += 1
        return mat, M, E, T, sweep

    # Kawasaki Dynamics Stepping
    def fast_kawasaki(self, mat, M, E, T, sweep):
        i, j, k, l = np.random.randint(0, self.lx, 4)
        a, b = mat[i, j], mat[k,l]
        # Gather neighbours and their spin sum for site B
        nnsum = mat[self.pbc(k - 1), l] \
                + mat[self.pbc(k + 1), l] \
                + mat[k, self.pbc(l - 1)] \
                + mat[k, self.pbc(l + 1)]
        # Get energy cost
        energy_cost_ba = -1*(a-b)*nnsum
        # Create new array (N^2 operation I think- idk for sure?...)
        mat_baab = mat.copy() # copy to avoid any risk of altering original matrix (unsure.)
        mat_baab[k,l] = a # set b to a
        # Gather neighbours and their spin sum for site A
        nnsum = mat_baab[self.pbc(i - 1), j] \
                + mat_baab[self.pbc(i + 1), j] \
                + mat_baab[i, self.pbc(j - 1)] \
                + mat_baab[i, self.pbc(j + 1)]
        # Gather change in magnetisation from setting A to B
        energy_cost_ab = -1*(b-a)*nnsum
        # Total energy cost
        energy_cost = energy_cost_ba + energy_cost_ab
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            mat[i,j] = b
            mat[k,l] = a
            E += energy_cost
        else:
            probability = np.exp(-1 * energy_cost / T)
            if np.random.random(1)[0] <= probability:
                mat[i, j] = b
                mat[k, l] = a
                E += energy_cost
        # Sweep complete.
        sweep += 1
        return mat, M, E, T, sweep

    # This is about 25% faster than the above: it has problems, though, for some reason.
    # TODO: DO NOT USE THIS
    # TODO: File report on https://github.com/numba/numba
    # TODO: Unexpected behaviour in handling arrays.
    # TODO: Indices flip unexpectedly: see https://github.com/numba/numba/issues/2747
    """
    Code produces unexpected behaviour.
    Set up two separate fast_ising objects
    Feed them the same matrix
    Ask them for the same a and b (replace ijkl by static values)
    50% of the the time you receive different a,b: the values switch: causes noise and makes issues.
    The issue is caused by the lower half of the code: doing two steps at once (setting mat[i,j]=mat[k,l] and etc.)
    This backfires, and switches the values for (a,b) for some stupid reason (probably how Numba compiles it?)
    Tested VIGOROUSLY ****ing christ I spent an entire day trying to fix this.
    In any cases, that's the summary. The lower half affects the upper.
    Defining copies of the arrays works instead of feeding the matrix in raw
    Running this code will incur you accuracy costs due to this noise error affecting near-neighbour swaps 50% time
    But it will give you a 25-33% performance saving.
    R E P O R T THE ERROR TO N U M B A ON THEIR G I T H U B !!!! 
    """
    def fast_kawaglauber(self, mat, M, E, T, sweep):
        # Generate (i,j)(k,l) for ab
        i, j, k, l = np.random.randint(0, self.lx, 4)

        # Grab spins of a and b
        a = mat[i,j]
        b = mat[k,l]

        # Calculate distance
        dist = np.sqrt((k-i)**2 + (l-j)**2)

        # they're neighbours (not diagonal but horizontal/vertical.)
        if dist < 1.4:

            # Gather neighbours and their spin sum for site B
            nnsum = mat[self.pbc(k - 1), l] \
                    + mat[self.pbc(k + 1), l] \
                    + mat[k, self.pbc(l - 1)] \
                    + mat[k, self.pbc(l + 1)]

            # Get energy cost of switching b to a
            energy_cost_ba = -1*(a-b)*nnsum

            # Create new array via copy: costly operation.
            mat_baab = mat.copy() # copy to avoid any risk of altering original matrix (unsure.)
            mat_baab[k,l] = a # set b to a

            # Gather neighbours and their spin sum for site A
            nnsum = mat_baab[self.pbc(i - 1), j] \
                    + mat_baab[self.pbc(i + 1), j] \
                    + mat_baab[i, self.pbc(j - 1)] \
                    + mat_baab[i, self.pbc(j + 1)]

            # Gather change in energy from setting A to B
            energy_cost_ab = -1*(b-a)*nnsum

        # they're not neighbours: we can avoid costly array copying.
        else:

            # Gather neighbours and their spin sum for site B
            nnsum = mat[self.pbc(k - 1), l] \
                    + mat[self.pbc(k + 1), l] \
                    + mat[k, self.pbc(l - 1)] \
                    + mat[k, self.pbc(l + 1)]

            # Get energy cost of switching b to a
            energy_cost_ba = -1 * (a - b) * nnsum

            # Gather neighbours and their spin sum for site A
            nnsum = mat[self.pbc(i - 1), j] \
                    + mat[self.pbc(i + 1), j] \
                    + mat[i, self.pbc(j - 1)] \
                    + mat[i, self.pbc(j + 1)]

            # Gather change in magnetisation from setting A to B
            energy_cost_ab = -1 * (b - a) * nnsum

        # Total energy cost
        energy_cost = energy_cost_ba + energy_cost_ab

        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            mat[i,j] = b
            mat[k,l] = a
            E += energy_cost
        else:
            probability = 1
            if np.random.random(1)[0] <= probability:
                mat[i, j] = b
                mat[k, l] = a
                E += energy_cost
        # Sweep complete.
        sweep += 1
        return a, b # mat, M, E, T, sweep

    # Get energy for matrix
    def energy(self, mat):
        E = 0
        for i in range(self.lx):
            for j in range(self.lx):
                nnsum = mat[self.pbc(i - 1), j] \
                        + mat[self.pbc(i + 1), j] \
                        + mat[i, self.pbc(j - 1)] \
                        + mat[i, self.pbc(j + 1)]
                E += -1 * mat[i, j] * nnsum
        E = (1 / 2) * E  # doublecounting
        return E

    # Both magnetisation and energy in one
    def fast_magenergy(self, mat):
        M = np.sum(mat)
        E = self.energy(mat)
        return M, E

    # Averages and errors calculation, but faster!
    def averages_errors(self, all_M, all_E, T, equil, measure, autocor):
        # Trim data to measurement range
        all_M, all_E = all_M[equil:equil+measure], all_E[equil:equil+measure]
        # Sample every autocor'th measurement
        all_M, all_E = all_M[0::autocor+1], all_E[0::autocor+1]
        num_samples = len(all_M)
        # Get squares
        all_MM, all_EE = all_M**2, all_E**2
        # Get averages
        avg_M, avg_E = np.mean(all_M), np.mean(all_E)
        avg_MM, avg_EE = np.mean(all_MM), np.mean(all_EE)
        # Get Errors. Note sample length = Autocorrelation length, so cancels out.
        avg_M_err, avg_E_err = np.sqrt(((avg_MM - avg_M**2)*(2))/(num_samples)), \
                               np.sqrt(((avg_EE - avg_E**2)*(2))/(num_samples))
        # Estimate susceptibility and specific heat/spin
        chi_true, c_true = (1/num_samples)*(1/T)*(avg_MM - avg_M**2),\
                           (1/num_samples)*(1/T**2)*(avg_EE - avg_E**2)
        # Error estimation for chi and c via the Bootstrap method
        number_of_resamples = 4000
        chi_list, c_list = np.empty(number_of_resamples), \
                           np.empty(number_of_resamples)
        for i in range(number_of_resamples):
            # Select num_samples random from num_samples
            resample = np.random.randint(0, num_samples, num_samples)
            # Grab Magnetism/Energy samples
            Qall_M, Qall_E = all_M[resample], all_E[resample]
            Qall_MM, Qall_EE = Qall_M ** 2, Qall_E ** 2
            # Get averages
            Qavg_M, Qavg_E = np.mean(Qall_M), np.mean(Qall_E)
            Qavg_MM, Qavg_EE = np.mean(Qall_MM), np.mean(Qall_EE)
            # Estimate susceptibility and specific heat/spin
            Qchi = (1/num_samples)*(1/T)*(Qavg_MM - Qavg_M**2)
            Qc = (1/num_samples)*(1/T**2)*(Qavg_EE - Qavg_E**2)
            # Append
            chi_list[i], c_list[i] = Qchi, Qc
        chi_average, c_average = np.mean(chi_list), np.mean(c_list)
        chichi_average, cc_average = np.mean(chi_list**2), np.mean(c_list**2)
        boot_chi, boot_c = np.sqrt(chichi_average - chi_average**2), np.sqrt(cc_average - c_average**2)
        return avg_M, avg_E, avg_M_err, avg_E_err, chi_true, boot_chi, c_true, boot_c


    # Cluster version of the Glauber. Note that this will just ignore changes in M and E (lazy.)
    def fast_clustglauber(self, mat, M, E, T, sweep):
        # Create a True/False array of spins (for the ones that we decide to flip.) False is unflipped.
        mark_array = np.zeros((self.lx, self.lx), dtype=types.boolean)

        # Pick a random point to start the cluster (the seed)
        start = np.random.randint(0, self.lx, 2)

        # Grab the spin state of this
        clust_spin = mat[start[0], start[1]]
        # Mark it "to flip" too
        mark_array[start[0], start[1]] = True

        # Set the flat probability of addition
        flatprob = 1 - np.exp(-2/T)

        # Set the tree
        current_centres = [start]

        # If len is zero, then it breaks.
        while len(current_centres) > 0:
            # Placeholder for the newly-flipped points we generate below
            new_centres = []

            # Work each recently-flipped centre in current_centres
            for centre in current_centres:
                # Set up centre index
                i, j = types.int32(centre[0]), types.int32(centre[1])
                # Grab indices for periodic edges
                leftrighttopbot_x = np.array([self.pbc(i - 1), self.pbc(i + 1), i, i])
                leftrighttopbot_y = np.array([j, j, self.pbc(j - 1), self.pbc(j + 1)])
                # Iterate over each of the edges
                for ii, jj in zip(leftrighttopbot_x, leftrighttopbot_y):
                    # Only work if not previously flipped
                    if mark_array[ii, jj] == False:
                        # Check whether the spin of the candidate is equal to the spin of the centre
                        if mat[ii,jj] == clust_spin:
                            # Generate random number [0,1]
                            u = np.random.random()
                            # Decide whether to flip
                            if flatprob > u:
                                ## Change the mark to True
                                mark_array[ii, jj] = True
                                ## Add this newly-flipped point as a new centre
                                new_centres.append(np.array([ii, jj]))
                        # If the spin of the candidate is not, then do not flip it.
                        else:
                            pass
                    else:
                        pass
            # Replace the old centres with the new ones
            current_centres = new_centres

        # Sweep complete.
        sweep += 1

        # Flip the entire cluster in one fell swoop
        mat = mat - 2*mark_array*mat

        # Calculate the magnetization and energy of the array
        M, E = self.fast_magenergy(mat)

        return mat, M, E, T, sweep
