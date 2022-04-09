# RELAX TWICE AS HARD
import time
import timeit

import numpy as np
from numba import njit



@njit(fastmath=True)
def converged(pot_ini, pot_fini, frac):

    # Add a bit to prevent zero errors
    pot_ini += 1e-16
    pot_fini += 1e-16

    # Difference + ratio (need average dif to be below frac)
    dif = np.abs((pot_fini - pot_ini)/pot_fini)
    avg_dif = np.mean(dif)

    if avg_dif > frac:
        return False # not converged: maximum frac is less than our limit
    else:
        return True # has converged! the maximum difference between two arrays satisfies our limit

@njit(fastmath=True)
def fast_gs(rho, pot_dist, shape, boundary, n):

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

    # Iterate the array forward by one step IN PLACE! (Gauss-Seidel!)
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                pot_dist[i,j,k] = (1/6)*(rho[i,j,k] +
                                         pot(i+1,j,k) + pot(i-1,j,k) +
                                         pot(i,j+1,k) + pot(i,j-1,k) +
                                         pot(i,j,k+1) + pot(i,j,k-1))

    return pot_dist, n + 1

@njit(fastmath=True)
def potential(rho, sqr_dist_relax):

    """
    :param rho: density matrix
    :param sqr_dist_relax: relaxation distance SQUARED for potential calculation
    :return: potential matrix

    """

    # Create a matrix
    pot_dist = np.empty_like(rho)

    # Going to evaluate it numerically straight from the density matrix
    shape = np.shape(rho)

    # Pifrac
    pifrac = 4*np.pi

    # Go over each point in the density matrix
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):

                # Sum for point i,j,k
                potsum = 0

                # Go over all other points and use the standard formula for potential
                for ii in range(shape[0]):
                    for jj in range(shape[1]):
                        for kk in range(shape[2]):

                            # Distance (skip self-contribution)
                            dist = np.sqrt((ii-i)**2 + (jj-j)**2 + (kk-k)**2 + sqr_dist_relax)

                            # Add contribution
                            potsum += (rho[ii,jj,kk]/(pifrac))/dist

                # Set the potential
                pot_dist[i,j,k] = potsum

    # Return
    return pot_dist



