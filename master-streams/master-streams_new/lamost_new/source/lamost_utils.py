import pickle

import numba
import numpy as np
from astropy.io import fits

# Grab / Set-up Source IDs used for matching (Gaia ones, DR2 afaik.) Returns indices, too. Astropy table.
from astropy.table import Table


# Get one big table with the source ids and indices
from numba import njit
from scipy.optimize import linear_sum_assignment
from sklearn.metrics.cluster import contingency_matrix
import astropy.units as u


def source_ids_indices(path, savepath, column='gaia_source_id', try_save=True, reset=False):

    """

    Columns of return are
    - _dr2_source_id
    - _table_index
    Both integers (signed) w/ the latter being ordered pythonically for the fits.

    :param path: of fits
    :param savepath: of save table
    :param column: -
    :param try_save: True/False bool
    :param: reset: True/False bool whether to reset the file
    :return: astropy table
    """

    # Have all the source_ids nice and saved (it's fast! :D)
    try:
        with open(savepath, "rb") as f:
            source_id_indices = pickle.load(f)

    except:

        with fits.open(path, mode='readonly', memmap=True) as hdul:

            # IDs and pythonic indices
            source_ids = hdul[1].data[column]
            source_ids = [int(d) for d in source_ids]
            source_ids = np.array(source_ids)
            source_indices = np.arange(0, len(source_ids), 1)

            # Slap them into an array together
            source_id_indices = np.ndarray(shape=(len(source_ids), 2), dtype=int)
            source_id_indices[:, 0] = source_ids
            source_id_indices[:, 1] = source_indices

            # Table up!
            columns = ["_dr2_source_id", "_table_index"]
            source_id_indices = Table(data=source_id_indices, names=columns)

        if try_save == True:
            with open(savepath, "wb") as f:
                pickle.dump(obj=source_id_indices, file=f)
                
        else:
            pass

    if reset == True:

        with fits.open(path, mode='readonly', memmap=True) as hdul:

            # IDs and pythonic indices
            source_ids = np.array(hdul[1].data[column])
            source_ids = np.array(source_ids)
            source_indices = np.arange(0, len(source_ids), 1)

            # Slap them into an array together
            source_id_indices = np.ndarray(shape=(len(source_ids), 2), dtype=np.int64)
            source_id_indices[:, 0] = source_ids
            source_id_indices[:, 1] = source_indices

            # Table up!
            columns = ["dr2_source_id_local", "_table_index"]
            source_id_indices = Table(data=source_id_indices, names=columns)

        if try_save == True:
            with open(savepath, "wb") as f:
                pickle.dump(obj=source_id_indices, file=f)

        else:
            pass

    # Return
    return source_id_indices

# Assess size/info of FITS files
def into(path):

    with fits.open(path, mode='readonly', memmap=True) as hdul:
        data = hdul[1].data

    print(("Size of {0} is {1}").format(path, len(data)))


def return_indices_of_a(a, b):
    """
    Find out which values from a (by index) can be found inside b.

    Example
        >>> a = [1,2,3]
        >>> b = [1,69,5,3,1,7]
        >>> return_indices_of_a(a,b)
        >>> (0,2)

    Works best when a 1-2-1 matching is possible, as this will not warn you of any duplicates.

    :param a:
    :param b:
    :return:
    """
    b_set = set(b)
    return [i for i, v in enumerate(a) if v in b_set]

@njit(fastmath=True)
def intmatch(a,b):

    """
    Find the mapping that takes array-of-integers a onto array-of-integers b, as tuples, i.e.

    Example 1
        >>> a = [1,4,3,6]
        >>> b = [4,6,3,1]
        >>> ret = intmatch(a,b)
        >>> [(0,3), (1,0), (2,2), (3,1)]

    Example 2
        >>> a = [1,4,3,6]
        >>> b = [1,6]
        >>> ret = intmatch(a,b)
        >>> [(0,0), (3,1)]

    Provides a set of tuples such that each tuple gives the index that maps (index_a, index_b) between (a,b),
    assuming that only a 1-2-1 mapping is possible (i.e. not multiple instances.) If a 1-2-1 ain't possible (DUPLICATES)
    then the provided "tuples" will exceed cardinality of a for example...

    Example 3
        >>> a = [1,2,3,4]
        >>> b = [1,1,1,1,1,4]
        >>> ret = intmatch(a,b)
        >>> [(0,0),(0,1),(0,2),(0,3),(0,4),(3,5)]
        >>> len(a)
        >>> 4
        >>> len(ret)
        >>> 6

    WE ASSUME THAT THE LIST a IS A SET with NO DUPLICATED VALUES
    The list b can have duplicates.

    :param a: array int64
    :param b: array int64
    :return: list of tuples
    """

    tuples = []
    for i in range(0,len(a)):
        for j in range(0,len(b)):
            if a[i] == b[j]:
                tuples.append((i,j))

    return tuples


# Identify duplicates
@njit(fastmath=True)
def catch_duplicates(default_length, list):

    """
    Takes the first array of intmatch (i.e. np.array(intmatch(a,b)).T[0]) such that a[i] = b[j]
    and we want the list of i

    default_length is the length of the non-duplicate list, i.e. len(set(a))

    Returns an array of booleans of size len(a) stating whether each element of a has duplicates in b.

    :param default_length: len(set(a))
    :param list: see above.
    :return: array size len(a) of booleans on whether duplicates exist for a[i] in b
    """

    duplicates = np.zeros(default_length, numba.core.types.bool_)
    for i in range(0, default_length):
        if list.count(i) > 1:
            duplicates[i] = True
    return duplicates

@njit(fastmath=True)
def remove_duplicates_firsttuple(tuples):
    """
    Only retain the first occurrence of a tuple and ignore the rest, where the original length of the list
    responsible for the first element of each tuple is master_len
    :param tuples: list of tuples of list matched (a,b)
    :param master_len: length of a (so as to get a 1-2-1 of a to elements of b)
    :return: list of tuples without duplicates
    """
    non_dupli = []
    non_dupli_tuples = []
    for tuple in tuples:
        if tuple[0] not in non_dupli:
            non_dupli.append(tuple[0])
            non_dupli_tuples.append(tuple)
        else:
            pass

    return non_dupli_tuples

@njit(fastmath=True)
def sigmaclip(a, low=4., high=4.):

    """

    Numba-fied Scipy Sigmaclip with the addition of index retrieval for the clipping.

    :param a: list/array
    :param low: lower-edge sigma
    :param high: upper-edge sigma
    :return: sigma-clipped a, indices of clip

    Example

        ---

        >>> print(lu.sigmaclip(np.array([1,2,3,4,10,15,50], float64), low=1, high=1))
        >>> (array([2., 3.]), 2.0, 3.0, array([0, 1], dtype=int64))

        ---

    Only the first and final argument returned are generally useful- critlower/critupper are remnants from scipy.

    """

    c = np.asarray(a).ravel()
    delta = 1
    while delta:
        c_std = c.std()
        c_mean = c.mean()
        size = c.size
        critlower = c_mean - c_std * low
        critupper = c_mean + c_std * high
        c_indices = [numba.types.int64(i) for i in range(0)]
        new_c = [numba.types.float64(i) for i in range(0)]
        for num, i in enumerate(c):
            if critlower <= i <= critupper:
                c_indices.append(num)
                new_c.append(i)
        c = np.asarray(new_c)
        delta = size - c.size

    return c, np.asarray(c_indices)

# Specify the default duplicate-removal parameters (we mainly care about kinematics, so sigma-clip by the radial vel)
rv_sigma = 1 # number of sigma to clip duplicate values by
rv_range_thresh = 15 # flat-value range of rv to accept in sigma-clip
@njit(fastmath=True)
def sigma_duplicate(rv_array, rv_err_array, feh_array, feh_err_array):
    """
    Will take the rv and feh arrays (and errors) for a given source with duplicates,
    and then decide (based on sigma-clip of the radial velocities first, and then the range
    of the resulting sigma clip) whether to retain it (in which case sigma-clipped means and
    errors are returned, with a flag "True") or refuse it (flag "False")

    Data must be provided as arrays with constrained values (i.e. pre-clean the fits files
    of interest to avoid the LAMOST -9999 flags.)

    :param rv_array:
    :param rv_err_array:
    :param feh_array:
    :param feh_err_array:
    :return: mean_rv, mean_rv_err, mean_feh, mean_feh_err, flag
    """

    # Sigma-clip the rv_array
    rv_clip, rv_clip_indices = sigmaclip(rv_array, low=rv_sigma, high=rv_sigma)

    # Now clip the feh/feh_err array
    rv_err_array_clip, \
    feh_array_clip, \
    feh_err_array_clip = rv_err_array[rv_clip_indices], \
                         feh_array[rv_clip_indices], \
                         feh_err_array[rv_clip_indices]

    # Decide if to reject the result or not based on if it exceeds our range threshold
    if np.max(rv_clip) - np.min(rv_clip) > rv_range_thresh:

        mean_rv = -9999.0
        mean_rv_err = -9999.0
        mean_feh = -9999.0
        mean_feh_err = -9999.0
        flag = False

    # Now sort out the new source information
    else:

        mean_rv = np.mean(rv_clip)
        mean_rv_err = (1/len(rv_err_array_clip))*np.sqrt(np.sum(rv_err_array_clip**2))
        mean_feh = np.mean(feh_array_clip)
        mean_feh_err = (1/len(feh_err_array_clip))*np.sqrt(np.sum(feh_err_array_clip**2))
        flag = True

    return mean_rv, mean_rv_err, mean_feh, mean_feh_err, flag
rv_dr3_difference_threshold_sigma = 5 # as multiple of dr3_rv_err
@njit(fastmath=True)
def sigma_duplicate_dr3(rv_array, rv_err_array, feh_array, feh_err_array, dr3_rv, dr3_rv_err):
    """
    Same as sigma_duplicate, except will take the dr3_rv and use that as a "base" against rv_array to assuage
    whether this source truly is the one that gaia observed (i.e. will only sigma-clip/use LAMOST rv within some
    threshold from the Gaia ID.)

    Use this in cases where a Gaia estimate for radial velocity is available.
    """

    # Get differences from dr3_rv
    difference = np.abs(rv_array - dr3_rv)
    retain = np.array([True if d < rv_dr3_difference_threshold_sigma*dr3_rv_err else False for d in difference], numba.types.bool_)
    retain_rv = rv_array[retain]

    # Now, if len(retain_rv) == 0, then none of the possible sources we have in the FITS matches this dr3_source_id
    if len(retain_rv) == 0:

        mean_rv = -9999.0
        mean_rv_err = -9999.0
        mean_feh = -9999.0
        mean_feh_err = -9999.0
        flag = False

        return mean_rv, mean_rv_err, mean_feh, mean_feh_err, flag

    # If len(retain_rv) != 0, then some of the sources provided have radial velocities that correlate to that of Gaia
    else:

        # Clip the rest of the arrays and continue
        rv_err_array, \
        feh_array, \
        feh_err_array = rv_err_array[retain],\
                        feh_array[retain],\
                        feh_err_array[retain]

        # Sigma-clip the rv_array
        rv_clip, rv_clip_indices = sigmaclip(rv_array, low=rv_sigma, high=rv_sigma)

        # Now clip the feh/feh_err array
        rv_err_array_clip, \
        feh_array_clip, \
        feh_err_array_clip = rv_err_array[rv_clip_indices], \
                             feh_array[rv_clip_indices], \
                             feh_err_array[rv_clip_indices]

        # Decide if to reject the result or not based on if it exceeds our range threshold
        if np.max(rv_clip) - np.min(rv_clip) > rv_range_thresh:

            mean_rv = -9999.0
            mean_rv_err = -9999.0
            mean_feh = -9999.0
            mean_feh_err = -9999.0
            flag = False

        # Now sort out the new source information
        else:

            mean_rv = np.mean(rv_clip)
            mean_rv_err = (1/len(rv_err_array_clip))*np.sqrt(np.sum(rv_err_array_clip**2))
            mean_feh = np.mean(feh_array_clip)
            mean_feh_err = (1/len(feh_err_array_clip))*np.sqrt(np.sum(feh_err_array_clip**2))
            flag = True

        return mean_rv, mean_rv_err, mean_feh, mean_feh_err, flag

mas = u.mas.to(u.rad)
@njit(fastmath=True)
def distance_matrix(ra1,dec1,ra2,dec2):

    """
    Take lists (as appropriate) and return a len(ra1)xlen(ra2) float64 array- assumes that the inputs are in radians.
    Element [i,j] is the distance from ra1[i] to ra2[j]

    Distances are returned in units of mas.

    :param ra1: 1D array float64
    :param dec1: 1D array float64
    :param ra2: 1D array float64
    :param dec2: 1D array float64
    :return: 2D array float64
    """

    # Empty arr
    distrix = np.empty((len(ra1), len(ra2)), numba.types.float64)

    # Generate distances
    for i in range(len(ra1)):
        for j in range(len(ra2)):

            # Set up unit vectors for 1/2
            v1 = np.array([np.cos(dec1[i])*np.cos(ra1[i]),
                      np.cos(dec1[i])*np.sin(ra1[i]),
                      np.sin(dec1[i])])
            v2 = np.array([np.cos(dec2[j])*np.cos(ra2[j]),
                      np.cos(dec2[j])*np.sin(ra2[j]),
                      np.sin(dec2[j])])
            dot_product = np.dot(v1,v2)

            # Error catch for close-unities
            if dot_product > 1:
                distrix[i,j] = 0
            elif dot_product < -1:
                distrix[i,j] = np.pi
            else:
                distrix[i,j] = np.arccos(dot_product)

                # Based on the small angle approximation - invalid
            # distrix[i,j] = np.sqrt((dec2[j]-dec1[i])**2 + ((ra2[j]-ra1[i])*np.cos(dec1[i]))**2)

    # Done
    return distrix/mas

# Resolution in mas
resolution = 1
@njit(fastmath=True)
def which_delete(distrix, matches):

    to_delete = [numba.types.int64(i) for i in range(0)]
    for num,match in enumerate(matches):
        if distrix[match[0],match[1]] >= resolution:
            to_delete.append(num)

    return to_delete

def crossmatch(ra1,dec1,ra2,dec2):
    """
    Recommend that len(ra1) < len(ra2). Returns the matches as a list of tuples. Removes matches that do not satisfy
    the resolution of which_delete. Tuple (a,b) gives the index from ra1 that matches the index in ra2.
    TODO: Deprecated- a JL function has been made for crossmatching instead. See APOGEE_STANDALONE.
    """

    # Get distrix + do linear sum assignment + set up matches tuples
    distrix = distance_matrix(ra1,dec1,ra2,dec2)
    row_ind, col_ind = linear_sum_assignment(distrix)
    matches = np.empty((len(ra1), 2), int)
    matches[:,0], matches[:,1] = row_ind, col_ind

    # Clean the matches up a notch (ones that are separated by more than the tolerable resolution of a single mas.)
    # Effectively consider a match only if it's within 1 mas- else it's an original star.
    to_delete = which_delete(distrix, matches)
    matches = np.delete(matches, to_delete, axis=0)

    # Return the matches
    return matches

def distmatch(dist1, dist2):
    """
    Recommend that len(ra1) < len(ra2). Returns the matches as a list of tuples. Removes matches that do not satisfy
    the resolution of which_delete. Tuple (a,b) gives the index from ra1 that matches the index in ra2.
    """

    # Get distrix + do linear sum assignment + set up matches tuples
    distrix = np.empty((len(dist1),len(dist2)), float)
    for i in range(len(dist1)):
        for j in range(len(dist2)):
            distrix[i,j] = np.abs(dist1[i]-dist2[j])
    row_ind, col_ind = linear_sum_assignment(distrix)
    matches = np.empty((len(dist1), 2), int)
    matches[:,0], matches[:,1] = row_ind, col_ind

    # Clean the matches up a notch (ones that are separated by more than the tolerable resolution of a single mas.)
    # Effectively consider a match only if it's within 1 mas- else it's an original star.
    to_delete = which_delete(distrix, matches)
    matches = np.delete(matches, to_delete, axis=0)

    # Return the matches
    return matches

npoints = 1000
@njit(fastmath=True)
def normals_monte_ICRSGAL(pmra, pmra_error, pmdec, pmdec_error):

    """

    All inputs floats.

    :return: Arrays [],[] with pmra/pmdec generated according to a Gaussian

    (cosfactor still exists!)

    """

    return np.random.normal(pmra, pmra_error, npoints), \
           np.random.normal(pmdec, pmdec_error, npoints)

from galpy.util import bovy_coords
def monte_ICRSGAL(ra, dec, b, pmra, pmra_error, pmdec, pmdec_error):

    """
    :param ra: deg
    :param dec: deg
    :param b: deg
    :param pmra: pmra*cos(dec)
    :param pmra_error: dpmra*cos(dec)
    :param pmdec: -
    :param pmdec_error: -
    :return: np.array([edmu_l, edmu_b]) with cos(b) factor not attached to edmu_l
    """

    # Get pmra,pmdec,vlos normals
    normals = normals_monte_ICRSGAL(pmra, pmra_error, pmdec, pmdec_error)
    ras, decs = ra*np.ones_like(normals[0]), dec*np.ones_like(normals[0])

    # Convert pmra*cos(dec), pmdec into dmu_l*cos(b), dmu_b
    dmuls, dmubs = bovy_coords.pmrapmdec_to_pmllpmbb(normals[0], normals[1], ras, decs, degree=True).T

    # Get standard deviations/etc w/ cos(b) factor removed
    edmu_l, edmu_b = np.std(dmuls)/np.cos(np.degrees(b)), np.std(dmubs)

    return np.array([edmu_l, edmu_b])

# TODO: faster.
def monte_ICRSGAL_table(table):
    table['edmu_l'], table['edmu_b'] = float(np.inf), float(np.inf)
    for row in table:
        row['edmu_l'], row['edmu_b'] = monte_ICRSGAL(*row[['ra','dec','b','pmra','pmra_error','pmdec','pmdec_error']])
    return table

