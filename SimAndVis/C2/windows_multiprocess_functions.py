from twod_gol import twod_gol

# Sim Params.
lx = 50
equilibration = int(0)
measurements = int(20e3)

def twod_run(identifier):
    model = twod_gol(lx=lx, equilibration=equilibration, measurements=measurements, not_up=True, identifier=identifier)
    model.run_multi()
    return 1

