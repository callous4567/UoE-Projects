import numpy as np
from numba import njit, float64


@njit(fastmath=True)
def electrofield(pot_dist, zmid, boundary):

    # Shape of pot
    shape = np.shape(pot_dist)

    # Const
    const = -1/2

    # Return value from pot_dist accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif k >= shape[2]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        elif k < 0:
            return boundary
        else:
            return pot_dist[i,j,k]

    # Placeholder for field directions and magnitudes for this slice
    field_vecs = np.zeros((2, shape[0], shape[1]), float64)

    # Get the field vectors about zmid/etc. Only the XY components (who cares about z, am I right?)
    for i in range(shape[0]):
        for j in range(shape[1]):
            field_vector = const*np.array([pot(i,j+1,zmid) - pot(i,j-1,zmid),
                                           pot(i+1,j,zmid) - pot(i-1,j,zmid)])
            mag = np.linalg.norm(field_vector)
            field_vecs[:,i,j] = field_vector/mag

    # Return
    return field_vecs

@njit(fastmath=True)
def electrofield_nonnormal(pot_dist, zmid, boundary):
    # Shape of pot
    shape = np.shape(pot_dist)

    # Const
    const = -1/2

    # Return value from pot_dist accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif k >= shape[2]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        elif k < 0:
            return boundary
        else:
            return pot_dist[i,j,k]

    # Placeholder for field directions and magnitudes for this slice
    field_vecs = np.empty((2, shape[0], shape[1]), float64)

    # Get the field vectors about zmid/etc. Only the XY components (who cares about z, am I right?)
    for i in range(shape[0]):
        for j in range(shape[1]):
            field_vector = const*np.array([pot(i,j+1,zmid) - pot(i,j-1,zmid),
                                           pot(i+1,j,zmid) - pot(i-1,j,zmid)])
            field_vecs[:,i,j] = field_vector

    # Return
    return field_vecs