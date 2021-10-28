import numpy as np
import galcentricutils
import graphutils

# You can nest dataframes with Pandas. https://stackoverflow.com/questions/51505504/pandas-nesting-dataframes
# Consider switching to this instead of Astropy tables.

def is_odd(num):
    return num & 0x1

gcc = galcentricutils.greatcount()
grapher = graphutils.threed_graph()
thetas, phis, indices = gcc.partialgen(10, 1, 45, 45)
grapher.unitsphere(thetas,phis, indices)
