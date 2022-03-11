import os
import pickle

import numpy as np
from twod_sirs import twod_sirs

# Sim Params.
lx = 50
equilibration = int(100*(lx**2))
measurements = int(1000*(lx**2))
number_of_runs = 5
graph = False

# Create directory for save
rootdir = os.getcwd()
savedir = "part_5"
try:
    os.mkdir(rootdir + "\\" + savedir)
except:
    pass

def twod_run(p1p2p3immu):
    p1, p2, p3, immu = p1p2p3immu
    array_dumps = []
    for ID in range(number_of_runs):
        model = twod_sirs(lx=lx, p1=p1, p2=p2, p3=p3, equilibration=equilibration, measurements=measurements, not_up=False, identifier=ID, immune_frac=immu)
        array_dump = model.main_multi()
        array_dumps.append(array_dump) # np.array([avg_I, avg_I_err, chi_true, chi_error, self.p[0], self.p[1], self.p[2]])
    # Take all the array dumps and get the data
    array_dumps = np.array(array_dumps).T
    avg_I, avg_I_err, chi_true, chi_error = array_dumps[0], array_dumps[1], array_dumps[2], array_dumps[3]
    # Get the true averages (uncorrelated so yeah- just propagate the errors.)
    avg_I, chi_true = np.mean(avg_I), np.mean(chi_true)
    # Go ahead and propagate errors
    """
    a = b + c + d / N 
    da = sqrt(db**2 + dc**2 + dd**2) / N 
    """
    avg_I_err, chi_error = avg_I_err**2, chi_error**2
    avg_I_err, chi_error = np.sqrt(np.sum(avg_I_err))/number_of_runs, np.sqrt(np.sum(chi_error))/number_of_runs
    # Save these values.
    unique_string = ("{0:.3f}").format(immu)
    with open(rootdir + "\\" + savedir + "\\" + unique_string + ".txt", 'wb') as f:
        pickle.dump(obj=np.array([avg_I, avg_I_err, chi_true, chi_error]), file=f)
    return 1
