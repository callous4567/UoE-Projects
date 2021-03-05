import time
import sys
import random
import numpy as np

# print() except it delays release of characters individually and makes the script look smoother on display.
# Change arguments for numbertime to 0,0 if you intend to do this vector stuff for decimal vectors instead of integer vectors... printing can take a while otherwise.
def time_delay(text):
    print()
    text = str(text)
    for c in text:
        sys.stdout.write(c)
        sys.stdout.flush()
        numbertime = random.uniform(0.005, 0.01)
        time.sleep(numbertime)
    print()
# Hard coded to allow for generation of random values in a range of -a to +a, either integers or continuous reals.
# Specify whether you want integer or continuous real value generation, and also specify the range. Hard-coded for this: no terminal input (it wasn't requested)
def rR_gen():
    integer_random = "True" # If false, decimal real is taken.
    random_range = int(10) # Specifies symmetric range of random generation
    if integer_random == "True":                                    # Generates random integers symmetrically in range
        return random.randint(-random_range, random_range)
    else:                                                           # Generates random reals symmetrically in range
        return random.uniform(-random_range, random_range)
# Produces a random vector, components satisfy range and symmetry of rR_gen(). Numpy array format instead of a list.
def rv_gen():
    v_base = np.zeros(3)
    return np.array([rR_gen() for d in v_base])
# Produce three vectors, print their magnitudes, and produce products and sums between the first two of them. Note, vectors 1,2,3 -> index 0,1,2 (normal Python).
# Returns (in list form)... [[v_1, v_2, v_3], [V_1, V_2, V_3], (v_1 + v_2), v_1*v_2, (v_1 x v_2)] where v_i is vector i and V_i is modulus of v_i.
def triple_v():
    v = np.array([rv_gen(), rv_gen(), rv_gen()]) # Matrix of vectors... makes things easier to call.
    time_delay("Three generated vectors:")
    for d in v:
        print(d)
    v_mags = [(np.inner(d,d))**0.5 for d in v] # Matrix of vector magnitudes

    time_delay("Moduli of vectors (in order):")
    time_delay(v_mags)

    v_sum = v[0] + v[1] # Sum of v_1 and v_2

    time_delay("Sum of first two vectors:")
    time_delay(v_sum)

    v_dot = np.inner(v[0], v[1]) # Dot of v_1 and v_2

    time_delay("Dot product between first two vectors:")
    time_delay(v_dot)

    v_cross = np.cross(v[0], v[1]) # v_1 cross v_2

    time_delay("Cross product of first two vectors, 1st cross 2nd:")
    time_delay(v_cross)

    return v, v_mags, v_sum, v_dot, v_cross
# is_equal adapted from vector_tester.py to work with numpy arrays (as opposed to vector class)
# see vector.py for documentation (what little there is)
def is_equal(self, other):
    index = int(0)
    for d in range(0, 2):
        str_self = ("{0:.4E}").format(self[index])
        str_othe = ("{0:.4E}").format(other[index])
        if str_self == str_othe:
            index += int(1)
        else:
            break
    if index == int(2):
        return ("True")
    else:
        return ("False")
# Identity testing for the vectors, with three identities provided. Calculates LHS and RHS of identities and compares them. Prints results of calculation.
# Also prints whether test was a success.
def tester_v():
    v = triple_v()[0]

    # First identity: v_1 x v_2 = - v_2 x v_1
    time_delay("Testing the first identity requested: v_1 x v_2 = -v_2 x v_1")
    lhs_1 = np.cross(v[0], v[1])
    rhs_1 = np.cross(v[1], v[0])*-1
    comp_1 = is_equal(lhs_1,rhs_1)
    time_delay(("LHS result of identity 1 is {}, RHS is {}").format(lhs_1, rhs_1))
    if comp_1 == "True":
        time_delay("LHS equals RHS, hence vector test 1 is a success.")
    else:
        time_delay("LHS does not equal RHS, hence vector test 2 is a failure.")

    # Second identity: v_1 x (v_2 + v_3) = (v_1 x v_2) + (v_1 x v_3)
    time_delay("Testing the second identity requested: v_1 x (v_2 + v_3) = (v_1 x v_2) + (v_1 x v_3)")
    lhs_2 = np.cross(v[0],(v[1] + v[2]))
    rhs_2 = np.cross(v[0], v[1]) + np.cross(v[0], v[2])
    comp_2 = is_equal(lhs_2, rhs_2)
    time_delay(("LHS result of identity 2 is {}, RHS is {}").format(lhs_2, rhs_2))
    if comp_2 == "True":
        time_delay("LHS equals RHS, hence vector test 2 is a success.")
    else:
        time_delay("LHS does not equal RHS, hence vector test 2 is a failure.")

    # Third identity: v_1 x (v_2 x v_3) = (v_1 * v_3)v_2 - (v_1 * v_2)v_3
    time_delay("Testing the third identity requested: v_1 x (v_2 + v_3) = (v_1*v_3)v_2 - (v_1*v_2)v_3")
    lhs_3 = np.cross(v[0], np.cross(v[1], v[2]))
    rhs_3 = (np.inner(v[0], v[2]))*v[1] - (np.inner(v[0], v[1]))*v[2]
    comp_3 = is_equal(lhs_3, rhs_3)
    time_delay(("LHS result of identity 3 is {}, RHS is {}").format(lhs_3, rhs_3))
    if comp_3 == "True":
        time_delay("LHS equals RHS, hence vector test 3 is a success.")
    else:
        time_delay("LHS does not equal RHS, hence vector test 3 is a failure.")
# Execute tester_v() function
tester_v()

# The array comparitor you're after is called "numpy.allclose(array_1, array_2, rtol=relative_tolerance, atol=absolute tolerance)
# It effectively does numerically what is_equal does... I mean, kinda?
# is_equal takes it and does it to scientific notation of four or whatever significant figures, which makes it more dynamic in a sense... not sure how well it works though.
# Worked so far though...