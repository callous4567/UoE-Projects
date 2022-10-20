import numpy as np
import APOGEE_standalone

def do_strmatch(pars):

    starhorse_ID, apogee_ID = pars
    matches = np.array(
        APOGEE_standalone.strmatch(starhorse_ID, apogee_ID),
        int
    )
    return matches

