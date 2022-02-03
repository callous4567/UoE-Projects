from twod_ising import twod_ising
import time
import numpy as np

def twod_run(temp_dynamic):
    temp, dynamic = temp_dynamic
    model = twod_ising()
    model.T = temp
    model.dyn = dynamic
    model.init_multiprocess()
    #model.time_delay(("Running {0:.2f} {1}").format(temp_dynamic[0], temp_dynamic[1]))
    model.main_multi()

def twod_regenerate_averages(temp_dynamic):
    temp, dynamic = temp_dynamic
    model = twod_ising()
    model.T = temp
    model.dyn = dynamic
    model.init_multiprocess()
    #model.time_delay(("Running {0:.2f} {1}").format(temp_dynamic[0], temp_dynamic[1]))
    model.main_multi(run=False)