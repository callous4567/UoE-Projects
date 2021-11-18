import dill
import numpy as np
from galpy import potential
from galpy.potential import plotRotcurve, MWPotential2014
from matplotlib import pyplot as plt
import astropy.units as u

import windows_directories

# This file holds GALPY-related things- it's slow to load and thus separated from galcentricutils.

# Generate energy and other associated potential properties for data.
# Uses Astropy Tables.
from galcentricutils import angular


class energistics(object):
    def __init__(self):
        null = " null "
        self.pot = "true"

    # Evaluate potentials for Table alongside Total Energy
    """
    v^2/2 + phi = total 
    """
    def pot_eval(self, table):
        # Get azimuth (phi) and cylindrical R
        table['R'] = np.sqrt(table['x']**2 + table['y']**2)
        # Evaluate potential
        table['pot'] = potential.evaluatePotentials(self.pot,
                                                    R = table['R'],
                                                    z = table['z'],
                                                    t=0)
        # Evaluate kinetic
        table['kinetic'] = (1/2)*np.sqrt(table['vx']**2 + table['vy']**2 + table['vz']**2)

        # Set total
        table['E'] = table['pot'] + table['kinetic']

        print(table['pot'])
        # Return.
        return table
