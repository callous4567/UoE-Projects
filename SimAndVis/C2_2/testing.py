import timeit
import time
import twod_sirs

# Just some generic speed tests to see whether sequential or parallel is faster.
test_sirs = twod_sirs.twod_sirs(lx=50, p1=0.1, p2=0.2, p3=0.3,
                                equilibration=int(1e6), measurements=int(25e6),
                                not_up=True)
test_sirs.I = test_sirs.fast.fast_infsum(test_sirs.mat)
test_sirs.fast_sequential()
test_sirs.fast_parallel()

# Define number of sweeps to test
nsweeps = 10000
nflips = int((test_sirs.lx**2) * nsweeps)

def seq():
    test_sirs.fast_sequential()
def par():
    test_sirs.fast_parallel()

uwu = timeit.timeit(seq, number=nflips)
owo = timeit.timeit(par, number=nsweeps)
print(("Sequential took {} while Parallel took {} per sweep.").format(uwu/nsweeps, owo/nsweeps))