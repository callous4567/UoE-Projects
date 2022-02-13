import numpy as np
from numba.experimental import jitclass
from numba import types, njit

# Somewhat-deprecated non-fastmath class that holds all the above njitted methods. Slower, but less noisy/no floats.
variable_specification = [
    ('good_luck', types.unicode_type),
    ('lx', types.int32),
    ('mat', types.float32[:,:,:]),
    ('leftrighttopbot_x', types.int32),
    ('leftrighttopbot_y', types.int32),
    ('mark_array', types.boolean),
    ('ii', types.int32),
    ('jj', types.int32)
]
@jitclass(variable_specification)
class fast_xy():
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
        # Pick random point
        i, j = np.random.randint(0, self.lx, 2)
        # Get near neighbour vectorial sum
        nnsum = mat[self.pbc(i - 1), j] + mat[self.pbc(i + 1), j] + mat[i, self.pbc(j - 1)] + mat[i, self.pbc(j + 1)]
        # Select a new angle for the spin i,j within 0, 2pi
        new = np.array([np.random.normal(), np.random.normal()])
        new_ij = new/np.linalg.norm(new)
        # Generate difference in magnetisation
        mag_cost = new_ij - mat[i,j]
        # Generate an energy cost for doing this flip (energy new - energy old)
        energy_cost = -1 * np.dot(mag_cost, nnsum)
        # If/else. If leq 0, do the flip. Else? Decide probabilistically.
        if energy_cost <= 0:
            mat[i, j] = new_ij
            # Add change in mag/energy
            M += mag_cost
            E += energy_cost
        else:
            probability = np.exp(-1 * energy_cost / T)
            if np.random.random(1)[0] <= probability:
                mat[i, j] = new_ij
                # Add change in mag/energy
                M += mag_cost
                E += energy_cost
        # Sweep complete.
        sweep += 1
        return mat, M, E, T, sweep

    # Get energy for matrix for the XY Model (uses dot products instead.)
    def energy(self, mat):
        E = 0
        for i in range(self.lx):
            for j in range(self.lx):
                # This is a vector sum
                nnsum = mat[self.pbc(i - 1), j] \
                        + mat[self.pbc(i + 1), j] \
                        + mat[i, self.pbc(j - 1)] \
                        + mat[i, self.pbc(j + 1)]
                # This is a dot product
                E -= np.dot(mat[i, j],nnsum)
        E /= 2
        return E

    # Both magnetisation and energy in one for the XY Model
    def fast_magenergy(self, mat):
        M = np.zeros(2)
        for i in range(self.lx):
            for j in range(self.lx):
                M += mat[i,j]
        E = self.energy(mat)
        return M, E

    # Compute the first-order divergence field
    def divergence_first(self, mat):
        dxy = 1/self.lx
        div_field = np.zeros((self.lx, self.lx))
        for i in range(self.lx):
            for j in range(self.lx):
                top, bottom, left, right = mat[self.pbc(i - 1), j], \
                                           mat[self.pbc(i + 1), j], \
                                           mat[i, self.pbc(j - 1)], \
                                           mat[i, self.pbc(j + 1)]
                div_x = (right[0] - left[0])/(2*dxy) # in increasing "x" to the right since j=x
                div_y = (bottom[1] - top[1])/(2*dxy) # in increasing "y" to the bottom since i=y
                div = div_x + div_y
                div_field[i,j] = div
        return div_field

    # Compute the angles for all elements of the array
    def angle(self, mat):
        angles = np.zeros((self.lx, self.lx))
        for i in range(self.lx):
            for j in range(self.lx):
                vec = mat[i,j]
                angles[i,j] = np.arctan2(vec[1],vec[0])
        # reverse the negatives
        for i in range(self.lx):
            for j in range(self.lx):
                if angles[i,j] > 0:
                    pass
                else:
                    angles[i,j] += 2*np.pi
        return angles

    # Cluster version of the Glauber. Note that this will just ignore changes in M and E (lazy.)
    def fast_clustglauber(self, mat, M, E, T, sweep):
        # Create a True/False array of spins (for the ones that we decide to flip.) False is unflipped.
        mark_array = np.zeros((self.lx, self.lx), dtype=types.boolean)

        # Generate a random flip for this cluster
        rand = np.array([np.random.normal(), np.random.normal()])
        randmag = np.sqrt(rand[0]**2 + rand[1]**2)
        randspin = rand/randmag

        # Pick a random point to start the cluster (the seed)
        start = np.random.randint(0, self.lx, 2)
        i, j = int(start[0]), int(start[1])

        # Take the random start point, and flip it (S -> S - 2(S*r)r...)
        mat[i, j] -= 2 * np.dot(mat[i, j], randspin) * randspin
        mark_array[i, j] = True

        # Set the tree
        current_centres = [start]
        # While the number of centres being examined is non-zero
        """
        For each centre in current_centres:
            Go to each centre i 
            Check the four surrounding points j: if they have been flipped before, do not work them. 
            Decide whether to flip them according to the Wolff probability (see papers)
            If they flip, append to new_centres as points to check on the next loop (expand the tree) 
            If they flip, keep track with mark that they've flipped 
        Once this is done, set current_centres to new_centres 
        If the length of new_centres is zero, break the loop. Else continue. 

        """
        # If len is zero, then it breaks.
        while len(current_centres) > 0:
            # Placeholder for the newly-flipped points we generate below
            new_centres = []

            # Work each recently-flipped centre in current_centres
            for centre in current_centres:
                # Set up centre index
                i, j = types.int32(centre[0]), types.int32(centre[1])

                # Grab indices for near neighbours
                leftrighttopbot_x = np.array([self.pbc(i - 1), self.pbc(i + 1), i, i])
                leftrighttopbot_y = np.array([j, j, self.pbc(j - 1), self.pbc(j + 1)])

                # Iterate over each of them
                for ii, jj in zip(leftrighttopbot_x, leftrighttopbot_y):
                    # Only work if not previously flipped
                    if mark_array[ii, jj] == False:
                        # Generate probability
                        probability = 1 - np.exp(np.min(np.array([0,(-1/T)*np.dot(-2*randspin*(np.dot(randspin, mat[ii,jj])),mat[i,j])])))
                        # Generate random number [0,1]
                        u = np.random.random()
                        # Decide whether to flip
                        if u < probability:
                            # Flip it in the array
                            mat[ii, jj] -= 2*np.dot(randspin, mat[ii,jj])*randspin
                            ## Change the mark to True
                            mark_array[ii, jj] = True
                            ## Add this newly-flipped point as a new centre
                            new_centres.append(np.array([ii, jj]))
                    else:
                        pass
            # Replace the old centres with the new ones
            current_centres = new_centres

        # Sweep complete.
        M, E = self.fast_magenergy(mat)
        sweep += 1
        return mat, M, E, T, sweep

    # Averages and errors calculation, but faster!
    def averages_errors(self, all_E, T, autocor):
        # Sample every autocor'th measurement
        all_E = all_E[0::autocor+1]
        num_samples = len(all_E)
        # Get squares
        all_EE = all_E**2
        # Get averages
        avg_E = np.mean(all_E)
        avg_EE = np.mean(all_EE)
        # Get Errors. Note sample length = Autocorrelation length, so cancels out.
        avg_E_err = np.sqrt(((avg_EE - avg_E**2)*(2))/(num_samples))
        # Estimate susceptibility and specific heat/spin
        c_true = (1/self.lx**2)*(1/T**2)*(avg_EE - avg_E**2)
        # Error estimation for chi and c via the Bootstrap method
        number_of_resamples = 100
        chi_list, c_list = np.empty(number_of_resamples), \
                           np.empty(number_of_resamples)
        for i in range(number_of_resamples):
            # Select num_samples random from num_samples
            resample = np.random.randint(0, num_samples, num_samples)
            # Grab Magnetism/Energy samples
            Qall_E = all_E[resample]
            Qall_EE = Qall_E ** 2
            # Get averages
            Qavg_E = np.mean(Qall_E)
            Qavg_EE = np.mean(Qall_EE)
            # Estimate susceptibility and specific heat/spin
            Qc = (1/self.lx**2)*(1/T**2)*(Qavg_EE - Qavg_E**2)
            # Append
            c_list[i] = Qc
        c_average = np.mean(c_list)
        cc_average = np.mean(c_list**2)
        boot_c = np.sqrt(cc_average - c_average**2)
        return avg_E/(self.lx**2), avg_E_err/(self.lx**2), c_true, boot_c

# Periodic Boundary condition.
@njit(fastmath=True)
def pbc(lx, i):
    if i >= lx:
        return 0
    if i < 0:
        return lx - 1
    else:
        return i

# Glauber Dynamics Stepping, but with FastMath instead (around 16% faster for a 128^2)
@njit(fastmath=True)
def fast_glauber(lx, mat, M, E, T, sweep):
    # Pick random point
    i, j = np.random.randint(0, lx, 2)
    # Get near neighbour vectorial sum
    nnsum = mat[pbc(lx, i - 1), j] + mat[pbc(lx, i + 1), j] + mat[i, pbc(lx, j - 1)] + mat[i, pbc(lx, j + 1)]
    # Select a new angle for the spin i,j within 0, 2pi
    new = np.array([np.random.normal(), np.random.normal()])
    new_ij = new/np.linalg.norm(new)
    # Generate difference in magnetisation
    mag_cost = new_ij - mat[i,j]
    # Generate an energy cost for doing this flip (energy new - energy old)
    energy_cost = -1 * np.dot(mag_cost, nnsum)
    # If/else. If leq 0, do the flip. Else? Decide probabilistically.
    if energy_cost <= 0:
        mat[i, j] = new_ij
        # Add change in mag/energy
        M += mag_cost
        E += energy_cost
    else:
        probability = np.exp(-1 * energy_cost / T)
        if np.random.random(1)[0] <= probability:
            mat[i, j] = new_ij
            # Add change in mag/energy
            M += mag_cost
            E += energy_cost
    # Sweep complete.
    sweep += 1
    return mat, M, E, T, sweep

# Get energy for matrix for the XY Model (uses dot products instead.)
@njit(fastmath=True)
def energy(lx, mat):
    E = 0
    for i in range(lx):
        for j in range(lx):
            # This is a vector sum
            nnsum = mat[pbc(lx, i - 1), j] \
                    + mat[pbc(lx, i + 1), j] \
                    + mat[i, pbc(lx, j - 1)] \
                    + mat[i, pbc(lx, j + 1)]
            # This is a dot product
            E += -1 * np.dot(mat[i, j],nnsum)
    E = (1 / 2) * E  # doublecounting
    return E

# Both magnetisation and energy in one for the XY Model
@njit(fastmath=True)
def fast_magenergy(lx, mat):
    E = energy(lx, mat)
    M = np.zeros(2)
    for i in range(lx):
        for j in range(lx):
            M += mat[i,j]
    return M, E

# Cluster version of the Glauber. Note that this will just ignore changes in M and E (lazy.)
@njit(fastmath=True)
def fast_clustglauber(lx, mat, T, sweep):
    # Create a True/False array of spins (for the ones that we decide to flip.) False is unflipped.
    mark_array = np.zeros((lx, lx), dtype=types.boolean)

    # Pick a random point to start the cluster (the seed)
    start = np.random.randint(0, lx, 2)
    i, j = int(start[0]), int(start[1])

    # Generate a random flip for this cluster
    rand = np.array([np.random.normal(), np.random.normal()])
    randmag = np.sqrt(rand[0]**2 + rand[1]**2)
    randspin = rand/randmag

    # Take the random start point, and flip it (S -> S - 2(S*r)r...)
    mat[i, j] -= 2 * np.dot(mat[i, j], randspin) * randspin
    mark_array[i, j] = True

    # Set the tree
    current_centres = [start]
    # While the number of centres being examined is non-zero
    """
    For each centre in current_centres:
        Go to each centre i 
        Check the four surrounding points j: if they have been flipped before, do not work them. 
        Decide whether to flip them according to the Wolff probability (see papers)
        If they flip, append to new_centres as points to check on the next loop (expand the tree) 
        If they flip, keep track with mark that they've flipped 
    Once this is done, set current_centres to new_centres 
    If the length of new_centres is zero, break the loop. Else continue. 

    """
    # If len is zero, then it breaks.
    while len(current_centres) > 0:
        # Placeholder for the newly-flipped points we generate below
        new_centres = []

        # Work each recently-flipped centre in current_centres
        for centre in current_centres:
            # Set up centre index
            i, j = types.int32(centre[0]), types.int32(centre[1])
            # Grab spin for centre
            centre_spin = mat[i, j]
            # Grab indices for periodic edges
            leftrighttopbot_x = np.array([pbc(lx, i - 1), pbc(lx, i + 1), i, i])
            leftrighttopbot_y = np.array([j, j, pbc(lx, j - 1), pbc(lx, j + 1)])

            # Iterate over each of them
            for ii, jj in zip(leftrighttopbot_x, leftrighttopbot_y):
                # Only work if not previously flipped
                if mark_array[ii, jj] == False:
                    # Generate probability
                    probability = 1 - np.exp(np.min(
                        np.array([0, (-1 / T) * np.dot(-2 * randspin * (np.dot(randspin, mat[ii, jj])), mat[i, j])])))
                    # Generate random number [0,1]
                    u = np.random.random()
                    # Decide whether to flip
                    if u < probability:
                        # Flip it in the array
                        mat[ii, jj] -= 2 * np.dot(randspin, mat[ii, jj]) * randspin
                        ## Change the mark to True
                        mark_array[ii, jj] = True
                        ## Add this newly-flipped point as a new centre
                        new_centres.append(np.array([ii, jj]))
                else:
                    pass
        # Replace the old centres with the new ones
        current_centres = new_centres

    # Sweep complete.
    M, E = fast_magenergy(lx, mat)
    sweep += 1

    return mat, M, E, T, sweep


# Compute the angles for all elements of the array
@njit(fastmath=True)
def angle(lx, mat):
    angles = np.zeros((lx, lx))
    for i in range(lx):
        for j in range(lx):
            vec = mat[i,j]
            angles[i,j] = np.arctan2(vec[1],vec[0])
    # reverse the negatives
    for i in range(lx):
        for j in range(lx):
            if angles[i,j] > 0:
                pass
            else:
                angles[i,j] += 2*np.pi
    return angles