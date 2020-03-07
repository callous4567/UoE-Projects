import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import time
import random
import math
import sys
import scipy
# Note to self for future reference, learn to use Scipy.
# Note to self... the above imports should be on all our documents! They super useful <3

class utilities(object):
    def time_delay(text):
        for c in text:
            sys.stdout.write(c)
            sys.stdout.flush()
            numbertime = random.uniform(0.001, 0.005)
            time.sleep(numbertime)
        print()

    def value_getter(name):
        time_text = ("Please provide the value for {}").format(name)
        utilities.time_delay(time_text)
        value = float(input())
        return(value)

    def string_getter(name):
        time_text = ("Please provide a string for {}").format(name)
        utilities.time_delay(time_text)
        stringer = str(input())
        return(stringer)


    def array_former(self, rows, cols):
        row = [int(0) for d in range(int(cols))]
        array = [row for d in range(int(rows))]
        numpyarray = np.array(array)
        return numpyarray
    # Returns empty array with (n x m) n = rows and m = cols.

    def array_formerComplex(self, rows, cols):
        array = np.zeros((rows, cols), dtype=np.complex_)
        return array
    # Complex array
