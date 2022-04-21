import numpy as np
from numba import njit

# Periodic boundary conditions.
@njit()
def pbc(i, lx):
    """
    Periodic boundary along one axis for index i of axis length lx
    :param i: index
    :param lx: axis length
    :return: periodic index
    """
    if i >= lx:
        return 0
    if i < 0:
        return lx - 1
    else:
        return i

# Evaluate the chemical potential for a phimat
@njit(fastmath=True)
def mu(mat, a, k, dx):

    """

    Get the chemical potential abiding the equation
    mu = -a(phi) + a(phi**3) - k(grad**2)phi
    With centred laplacian

    :param mat: Phi-matrix
    :param a: see: ref. equation
    :param k: see: ref. equation
    :param dx: spatial discretization
    :return: chemical potential matrix

    """

    # Placeholder
    newmat = np.empty_like(mat)

    # Generate values
    lx,ly = np.shape(mat)
    for i in range(lx):
        for j in range(ly):
            newmat[i,j] = -a*(mat[i,j] - mat[i,j]**3) - k*(mat[i,pbc(j+1, ly)] + mat[pbc(i+1, lx),j]
                                                           - 4*mat[i,j] +
                                                           mat[i,pbc(j-1, ly)] + mat[pbc(i-1, lx),j])/(dx**2)

    # Return the mumat
    return newmat

# Evaluate using numpy.roll instead
def numpy_mu(mat, a, k, dx):

    """

    Get the chemical potential abiding the equation
    mu = -a(phi) + a(phi**3) - k(grad**2)phi
    With centred laplacian

    :param mat: Phi-matrix
    :param a: see: ref. equation
    :param k: see: ref. equation
    :param dx: spatial discretization
    :return: chemical potential matrix

    """

    # Define newmat
    newmat = -a*(mat - mat**3) - k*(np.roll(mat, 1, 1) + np.roll(mat, 1, 0)
                                            - 4*mat +
                                            np.roll(mat, -1, 1) + np.roll(mat, -1, 0))/(dx**2)

    # Return the mumat
    return newmat

# Iterate forward a phi matrix
@njit(fastmath=True)
def fast_cahn(mat, a, k, dx, dt, M):

    """

    Iterate the composite order matrix forward by dt according to the Cahn-Hilliard equation
    dphi/dt = M(grad**2)mu

    :param mat: comp. ord. matrix
    :param a: see: mumat
    :param k: see: mumat
    :param dx: spatial discretization
    :param dt: temporal discretization
    :param M: see: ref. equation
    :return: updated composite order matrix
    """

    # Generate the mumat for this mat
    mumat = mu(mat, a, k, dx)

    # Placeholder
    newmat = np.empty_like(mat)

    # Set up constant
    const = M*dt/(dx**2)

    # Generate values
    lx, ly = np.shape(mat)
    for i in range(lx):
        for j in range(ly):
            newmat[i, j] = mat[i,j] + const*(mumat[i, pbc(j + 1, ly)] + mumat[pbc(i + 1, lx), j]
                                             - 4 * mumat[i, j] +
                                             mumat[i, pbc(j - 1, ly)] + mumat[pbc(i - 1, lx), j])

    # Return the newmat
    return newmat

# Iterate it forward but using np.roll (NOT WITH NUMBA!) just to benchmark speed/etc
def numpy_cahn(mat, a, k, dx, dt, M):

    # Generate the mumat for this mat
    mumat = numpy_mu(mat, a, k, dx)

    # Set up constant
    const = M*dt/(dx**2)

    # Define newmat
    newmat = mat + const*(np.roll(mumat, 1, 1) + np.roll(mumat, 1, 0)
                          - 4*mumat +
                          np.roll(mumat, -1, 1) + np.roll(mumat, -1, 0))

    # Return the newmat
    return newmat

# Integrate the phi-matrix (to check for conservation of phi)
@njit(fastmath=True)
def phitegrate(mat, dx):

    """
    Get spatial integral of the composite order matrix (2D)
    :param mat: composite order matrix
    :param dx: spatial discretization
    :return: scalar value integral
    """

    # Each point has an effective area of dx**2. So just return dx**2 * array all summed up
    return (dx**2)*np.sum(mat)

# Get the free energy over the array (integrated return- no array return.)
@njit(fastmath=True)
def free_energy(mat, a, k, dx):

    # Placeholder
    newmat = np.empty_like(mat)

    # Set up constant
    const = k/2
    const /= (2*dx)**2

    # Generate values. We're using centred differencing for gradphi**2
    lx, ly = np.shape(mat)
    for i in range(lx):
        for j in range(ly):
            newmat[i,j] = -a*(0.5*(mat[i,j]**2) - 0.25*(mat[i,j]**4)) + \
                          const*((mat[i,pbc(j+1,ly)] - mat[i,pbc(j-1,ly)])**2 +
                                 (mat[pbc(i+1,lx),j] - mat[pbc(i-1,lx),j])**2)

    # Get int and return
    return np.sum(newmat)*(dx**2)