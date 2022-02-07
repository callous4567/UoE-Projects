import numpy as np
from numba.experimental import jitclass
from numba import types

variable_specification = [
    ('good_luck', types.unicode_type),
    ('lx', types.int32),
    ('mat', types.float32[:,:,:])
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
                E += -1 * np.dot(mat[i, j],nnsum)
        E = (1 / 2) * E  # doublecounting
        return E

    # Both magnetisation and energy in one for the XY Model
    def fast_magenergy(self, mat):
        E = self.energy(mat)
        M = np.zeros(2)
        for i in range(self.lx):
            for j in range(self.lx):
                M += mat[i,j]
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