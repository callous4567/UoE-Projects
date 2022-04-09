import warnings
from cahn import cahn
converge_frac = 1e-6
def convergence_cahn(zippar):
    locahn = cahn(*zippar, 100, 0.01, int(1e4), int(1e3)) # a, k, dx, dt, M, iniphi
    sweep, fe = locahn.run_multi(converge_frac)
    return sweep, fe