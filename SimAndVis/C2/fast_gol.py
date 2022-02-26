import numba
import numpy as np
from numba.experimental import jitclass
from numba import types, extending, njit

# <><> THIS IS NOT OUR OWN CODE! CREDIT BELOW!!!!! <><>
"""
Required for array comparisons. See https://github.com/numba/numba/pull/7067. Thanks jpivarski.
"""
@njit
def _isclose_item(x, y, rtol, atol, equal_nan):
    if np.isnan(x) and np.isnan(y):
        return equal_nan
    elif np.isinf(x) and np.isinf(y):
        return (x > 0) == (y > 0)
    elif np.isinf(x) or np.isinf(y):
        return False
    else:
        return abs(x - y) <= atol + rtol * abs(y)
@extending.overload(np.isclose)
def isclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
    if (isinstance(a, numba.types.Array) and a.ndim > 0) or (
        isinstance(b, numba.types.Array) and b.ndim > 0
    ):
        def isclose_impl(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
            # FIXME: want to broadcast_arrays(a, b) here
            x = a.reshape(-1)
            y = b.reshape(-1)
            out = np.zeros(len(y), np.bool_)
            for i in range(len(out)):
                out[i] = _isclose_item(x[i], y[i], rtol, atol, equal_nan)
            return out.reshape(b.shape)

    elif isinstance(a, numba.types.Array) or isinstance(b, numba.types.Array):
        def isclose_impl(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
            return np.asarray(
                _isclose_item(a.item(), b.item(), rtol, atol, equal_nan)
            )

    else:
        def isclose_impl(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):
            return _isclose_item(a, b, rtol, atol, equal_nan)

    return isclose_impl
# <><> THIS IS NOT OUR OWN CODE! CREDIT ABOVE!!!!! <><>

# Something specific for us to use, similar to np.all_close.
@njit
def state_compare(a,b):
    # "isclose" array. To verify if two arrays are identical or not. Useless in the face of oscillators...
    is_it_close = np.isclose(a, b)
    # Verify if all close
    if is_it_close.all() == True:
       ALL = True
    else:
       ALL = False
    return ALL

# Some notes.
"""
This is Conways Game of Life.
We're using all eight near neighbours (as typical) TOP BOT LEFT RIGHT and clockwise diagonals from TOP LEFT.
You're either alive, or dead: True or False. 
• Any live cell with less than 2 live neighbours dies.
• Any live cell with 2 or 3 live neighbours lives on to the next step.
• Any live cell with more than 3 live neighbours dies.
• Any dead cell with exactly 3 live neighbours becomes alive.

"""
variable_specification = [
    ('good_luck', types.unicode_type),
    ('lx', types.int32),
    ('sweep', types.int32),
]
@jitclass(variable_specification)
class fast_gol():
    def __init__(self, lx):
        good_luck = "Godspeed, mate. It's the second round! Go dirty with Numba- code speed doesn't matter."
        self.lx = lx

    # Periodic boundary conditions.
    def pbc(self, i):
        if i >= self.lx:
            return 0
        if i < 0:
            return self.lx - 1
        else:
            return i

    # Periodic boundary conditions, but using modulo instead. Kept here for reference. Same speed as pbc in tests!!!
    """
    Equivalency test:
    50 -> 0, 51 -> 1, 99->49, all correct pythonically for the upper region
    -1 -> 49, and so forth (so correct for lower region.) 
    """
    def pbc_mod(self, i):
        return i%self.lx

    # Iterate the Game of Life.
    def fast_gol(self, mat, sweep):
        # Empty mat.
        new_mat = np.zeros_like(mat)
        # Go through mat and exploit the rules of the Game of Life.
        for i in range(self.lx):
            for j in range(self.lx):
                # Gather sum of near neighbours
                neighbours = mat[self.pbc(i - 1), self.pbc(j - 1)] + \
                             mat[self.pbc(i - 1), self.pbc(j)] + \
                             mat[self.pbc(i - 1), self.pbc(j + 1)] + \
                             mat[self.pbc(i), self.pbc(j - 1)] + \
                             mat[self.pbc(i), self.pbc(j + 1)] + \
                             mat[self.pbc(i + 1), self.pbc(j - 1)] + \
                             mat[self.pbc(i + 1), self.pbc(j)] + \
                             mat[self.pbc(i + 1), self.pbc(j + 1)]
                # If alive
                if mat[i,j] == True:
                    if neighbours < 2:
                        new_mat[i,j] = False
                    elif 2 <= neighbours <= 3:
                        new_mat[i,j] = True
                    elif neighbours > 3:
                        new_mat[i,j] = False
                # If dead
                else:
                    if neighbours == 3:
                        new_mat[i,j] = True

        # Return (including the number of "alive"
        return new_mat, np.sum(new_mat), sweep + 1

    # Get centre of mass of a glider/object. Returns (x,y) for a graph, i.e. (row,col).
    """
    Assume that glider/object is much much smaller than the array we're dealing with
    Specifically, assume that the entire width OR height of the glider, is AT MOST (lx/2 - 1)
    Ascertain if the glider indeed crosses any boundaries
    - i.e. using difference between top and bottom-most points
    - i.e. using difference between left and rightmost points
    Roll the glider forward if it does cross a boundary 
    Determine centre of mass normally
    Roll the COM back. 
    Bingo, done.
    We assume we're not bothering with calculating active sites (and thus include this at computational cost.) 
    If you do not wish to have this cost, insert as argument and remove "active_sites += 1" 
    """
    def glide_com(self, mat):
        # Get all live coordinates, Pythonically
        active_sites = 0
        live_coords = []
        for i in range(self.lx):
            for j in range(self.lx):
                if mat[i,j] == True:
                    live_coords.append([i,j])
                    active_sites += 1
        live_coords = np.array(live_coords).T

        # Determine if a row-roll is necessary: if larger than half a grid, we can assume the glider has passed bottom
        roll_row = True if np.max(live_coords[0]) - np.min(live_coords[0]) > self.lx/2 else False
        # Roll if it is
        if roll_row == True:
            live_coords[0] += int(self.lx/2)

        # Determine if a col-roll is necessary.
        roll_col = True if np.max(live_coords[1]) - np.min(live_coords[1]) > self.lx / 2 else False
        # Roll if it is
        if roll_col == True:
            live_coords[1] += int(self.lx/2)

        # (If) we've rolled, use modulo on it all (periodic boundary conditions, effectively)
        # This will take the bit we pushed off the edge, and put it at the start
        if roll_row or roll_col == True:
            live_coords %= self.lx

        # Get the centre of mass now
        comx, comy = np.mean(live_coords[0]), np.mean(live_coords[1])

        # Revert the rolls we did (and do another % to account for negative positions of COM.)
        if roll_row == True:
            comx -= int(self.lx/2)
            comx %= self.lx
        if roll_col == True:
            comy -= int(self.lx/2)
            comy %= self.lx

        # Return is flipped, since when plotting the row is the y and column is the x.
        return np.array([comy, comx])


