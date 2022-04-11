import numpy as np
from numba import njit, float64


@njit(fastmath=True)
def fast_gs_afield(rho, potential_distribution, shape, boundary, n):

    # Return value from potential_distribution accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        # If above or below, then return the midline of the array as the boundary (thus avoiding edge effects on wire)
        elif k >= shape[2]:
            return potential_distribution[i,j,int(shape[2]/2)]
        elif k < 0:
            return potential_distribution[i,j,int(shape[2]/2)]
        else:
            return potential_distribution[i,j,k]

    # Iterate the array forward by one step IN PLACE! (Gauss-Seidel!)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                potential_distribution[i,j,k] = (1/6)*(rho[i,j,k] +
                                         pot(i+1,j,k) + pot(i-1,j,k) +
                                         pot(i,j+1,k) + pot(i,j-1,k) +
                                         pot(i,j,k+1) + pot(i,j,k-1))

    return potential_distribution, n + 1

@njit(fastmath=True)
def fast_sor_inplace(rho, potential_distribution, shape, boundary, n, omega):

    # Const
    const = 1 - omega
    const_2 = omega/6

    # Return value from potential_distribution accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        # If above or below, then return the midline of the array as the boundary (thus avoiding edge effects on wire)
        elif k >= shape[2]:
            return potential_distribution[i,j,int(shape[2]/2)]
        elif k < 0:
            return potential_distribution[i,j,int(shape[2]/2)]
        else:
            return potential_distribution[i,j,k]

    # Iterate the array forward by one step IN PLACE! (Gauss-Seidel!)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                potential_distribution[i,j,k] = const*potential_distribution[i,j,k] + const_2*(rho[i,j,k] +
                                                                   pot(i+1,j,k) + pot(i-1,j,k) +
                                                                   pot(i,j+1,k) + pot(i,j-1,k) +
                                                                   pot(i,j,k+1) + pot(i,j,k-1))

    return potential_distribution, n + 1

@njit(fastmath=True)
def fast_sor(rho, pot_dist, shape, boundary, n, omega):

    pot_gs, new_n = fast_gs_afield(rho, pot_dist, shape, boundary, n)
    pot_fin = (1-omega)*pot_dist + omega*pot_gs

    return pot_fin, new_n

@njit(fastmath=True)
def magnetofield(potential_distribution, zmid, boundary):

    # Shape of pot
    shape = np.shape(potential_distribution)

    # Const
    const = 1/2

    # Return value from potential_distribution accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        # If above or below, then return the midline of the array as the boundary (thus avoiding edge effects on wire)
        elif k >= shape[2]:
            return potential_distribution[i,j,zmid]
        elif k < 0:
            return potential_distribution[i,j,zmid]
        else:
            return potential_distribution[i,j,k]

    # Placeholder for field directions and magnitudes for this slice
    magnetic_field_matrix = np.empty((2, shape[0], shape[1]), float64)

    # Get the field vectors about zmid/etc. Only the XY components (who cares about z, am I right?)
    for i in range(shape[0]):
        for j in range(shape[1]):
            magfieldvec = const*np.array([-1*pot(i+1,j,zmid) + pot(i-1,j,zmid),
                                          pot(i,j+1,zmid) - pot(i,j-1,zmid)])
            mag = np.linalg.norm(magfieldvec)
            magnetic_field_matrix[:,i,j] = magfieldvec/mag

    # Return
    return magnetic_field_matrix

@njit(fastmath=True)
def magnetofield_nonnormal(potential_distribution, zmid, boundary):

    # Shape of pot
    shape = np.shape(potential_distribution)

    # Const
    const = 1/2

    # Return value from potential_distribution accounting for boundary condition (ifs/elifs. Idk how else to do this.)
    def pot(i, j, k):
        if i >= shape[0]:
            return boundary
        elif j >= shape[1]:
            return boundary
        elif i < 0:
            return boundary
        elif j < 0:
            return boundary
        # If above or below, then return the midline of the array as the boundary (thus avoiding edge effects on wire)
        elif k >= shape[2]:
            return potential_distribution[i,j,zmid]
        elif k < 0:
            return potential_distribution[i,j,zmid]
        else:
            return potential_distribution[i,j,k]

    # Placeholder for field directions and magnitudes for this slice
    magnetic_field_matrix = np.empty((2, shape[0], shape[1]), float64)

    # Get the field vectors about zmid/etc. Only the XY components (who cares about z, am I right?)
    for i in range(shape[0]):
        for j in range(shape[1]):
            magfieldvec = const*np.array([-1*pot(i+1,j,zmid) + pot(i-1,j,zmid),
                                          pot(i,j+1,zmid) - pot(i,j-1,zmid)])
            magnetic_field_matrix[:,i,j] = magfieldvec

    # Return
    return magnetic_field_matrix
