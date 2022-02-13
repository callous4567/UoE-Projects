from xy_ising import xy_ising

# Params for sims

def model_run(temp):
    model = xy_ising(temp)
    model.run_multi()
    return 1


