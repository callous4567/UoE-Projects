from SimAndVis.C2.twod_sirs import twod_sirs

# Sim Params.
lx = 50
equilibration = int(2.5e6)
measurements = int(25e6)

def twod_run(p1p2p3):
    p1, p2, p3 = p1p2p3
    model = twod_sirs(lx=lx, p1=p1, p2=p2, p3=p3, equilibration=equilibration, measurements=measurements, not_up=True)
    model.main_multi(run=True)
    return 1

def twod_regenerate_averages(p1p2p3):
    p1, p2, p3 = p1p2p3
    model = twod_sirs(lx=lx, p1=p1, p2=p2, p3=p3, equilibration=equilibration, measurements=measurements, not_up=True)
    model.main_multi(run=False)

