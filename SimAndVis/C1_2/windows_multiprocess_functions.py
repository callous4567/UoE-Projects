from twod_ising import twod_ising

# Params for sims
lx = 50
equilibration = int(0.25e6)
measurements = int(25e6)

def twod_run(temp_dynamic):
    temp, dynamic = temp_dynamic
    model = twod_ising(lx, temp, dynamic, equilibration, measurements)
    model.main_multi()
    return 1

def twod_regenerate_averages(temp_dynamic):
    temp, dynamic = temp_dynamic
    model = twod_ising(lx, temp, dynamic, equilibration, measurements)
    model.main_multi(run=False)