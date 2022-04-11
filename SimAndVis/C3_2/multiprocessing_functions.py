from SimAndVis.C3_2.relax import relaxed_poisson


def do_sor(sor_valuerho):
    sor_set, rho = sor_valuerho
    poisson = relaxed_poisson(0, rho, 1e-3, 10, 2000)
    n_value = poisson.converge_runsim(sor_set)
    return n_value