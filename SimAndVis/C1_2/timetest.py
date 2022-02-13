import timeit
from SimAndVis.C1_2.xy_ising import xy_ising
import fast_xy

test_ising = xy_ising()
test_ising.M, test_ising.E = test_ising.fast.fast_magenergy(test_ising.mat)
test_ising.fast_glauber()
test_ising.fast_glauber_fastmath()
def test():
    test_ising.fast_glauber()
a = timeit.timeit(test, number=int(1e6))
print(a)