import numpy as np
from numba import njit, int32, float64


@njit(fastmath=True)
def pbc(i, lx):
    """
    :param i: grid coord oned
    :param lx: grid scale
    :return: pbc index
    """
    if i >= lx:
        return 0
    if i < 0:
        return lx - 1
    else:
        return i

@njit(fastmath=True)
def typefield(a, b, c, d):

    """
    :return: typefield [1,2,3,0] for max of [a,b,c,d]
    """

    # Get shape
    nx, ny = np.shape(a)

    # Specify the indices for each one (a,b,c,d) for whichever is maximum
    indices = np.array([1,2,3,0], int32)

    # Set up an array stack (new array is 4:nx:ny in shape)
    stack = np.zeros((4, nx, ny), float64)
    stack[0], stack[1], stack[2], stack[3] = a,b,c,d

    # Preallocate
    type = np.ones((nx,ny), int32)

    # Go over the array
    for i in range(nx):
        for j in range(ny):

            # Get the argmax for this i,j
            argmax = np.argmax(stack[:,i,j])

            # Set the type
            type[i,j] = indices[argmax]

    # Return the type array
    return type

@njit(fastmath=True)
def fast_model(a, b, c, d, D, q, p, dx, dt, sweep):

    """
    :return: a,b,c,d
    """

    # Preallocate
    new_a, new_b, new_c = np.zeros_like(a), np.zeros_like(a), np.zeros_like(a)

    # Get shape
    nx,ny = np.shape(a)

    # Set up multiplicative constants
    const_1 = (D*dt)/(dx**2)
    const_2 = (q*dt)
    const_3 = (p*dt)

    # Go over the array
    for i in range(nx):
        for j in range(ny):

            # Iterate forward for a
            new_a[i,j] = const_1*(a[pbc(i+1, nx),j] + a[pbc(i-1,nx),j]
                                  + a[i,pbc(j+1,ny)] + a[i,pbc(j-1,ny)]
                                  - 4*a[i,j]) \
                         + const_2*a[i,j]*d[i,j] \
                         - const_3*(a[i,j]*c[i,j]) \
                         + a[i,j]

            # Iterate forward for b
            new_b[i,j] = const_1*(b[pbc(i+1, nx),j] + b[pbc(i-1,nx),j]
                                  + b[i,pbc(j+1,ny)] + b[i,pbc(j-1,ny)]
                                  - 4*b[i,j]) \
                         + const_2*b[i,j]*d[i,j] \
                         - const_3*(a[i,j]*b[i,j]) \
                         + b[i,j]

            # Iterate forward for c
            new_c[i,j] = const_1*(c[pbc(i+1, nx),j] + c[pbc(i-1,nx),j]
                                  + c[i,pbc(j+1,ny)] + c[i,pbc(j-1,ny)]
                                  - 4*c[i,j]) \
                         + const_2*c[i,j]*d[i,j] \
                         - const_3*(b[i,j]*c[i,j]) \
                         + c[i,j]

    # Return the matrices a,b,c,d
    return new_a, new_b, new_c, (1 - new_a - new_b - new_c), sweep + 1

@njit(fastmath=True)
def get_typefraction(typematrix):

    """
    :param typematrix: array of types
    :return: [fraca, fracb, frac] as array
    """

    # Scale
    size = typematrix.shape[0]*typematrix.shape[1]

    # Fractions
    frac = np.zeros(3, float64)

    # Go over array and ascertain
    for i in range(typematrix.shape[0]):
        for j in range(typematrix.shape[1]):
            if typematrix[i,j] == 0:
                pass
            else:
                frac[typematrix[i,j]-1] += 1

    # Return
    return frac/size

@njit(fastmath=True)
def distance(coord_1, coord_2, nx):
    dx = np.abs(coord_2[0] - coord_1[0])
    dy = np.abs(coord_2[1] - coord_1[1])
    shape = np.shape(coord_1)
    if dx > nx/2:
        dx = nx - dx
    if dy > nx/2:
        dy = nx - dy
    return np.sqrt(dx**2 + dy**2)

@njit(fastmath=True)
def correlation_func(typetrix):
    """

    :param typetrix:
    :return: r_range and probabilities (same cardinality)
    """
    nx = np.shape(typetrix)[0]

    # Specify the r_range to do this over
    r_range = np.linspace(1, nx/2, 10)
    probs = []

    # Go over array
    for r in r_range:
        probabilities = []
        for i in range(nx):
            for j in range(nx):
                # Total number of points not paired
                nonsame = 0
                # Total same
                same = 0
                # Again
                for k in range(nx):
                    for l in range(nx):
                        if distance([i,j],[k,l],nx) < r:
                            if typetrix[i,j] == typetrix[k,l]:
                                same += 1
                            else:
                                nonsame += 1
                # Get the correlation for this point
                probability_of_same = (same/(same + nonsame))
                probabilities.append(probability_of_same)
        # average
        average_probability = np.mean(np.array(probabilities))
        probs.append(average_probability)

    # Return
    return r_range, probs


