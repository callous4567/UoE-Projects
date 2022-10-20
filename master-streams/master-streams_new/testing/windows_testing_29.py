import numpy as np
from astropy import units as u
from galpy import orbit

from energistics_new import orbigistics
orbigist = orbigistics()
integrate_time = np.linspace(0, 0.3*1e9, 500)*u.yr
element = [27,34,13,4,2,140]
element[0] *= u.deg
element[1] *= u.deg
element[2] *= u.kpc
element[3] *= u.mas / u.yr
element[4] *= u.mas / u.yr
element[5] *= u.km / u.s
orbit_forward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0] * u.kpc, vo=orbigist.rovo[1] * u.km / u.s,
                            zo=orbigist.zo * u.kpc, lb=True)
orbit_backward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0] * u.kpc,
                             vo=orbigist.rovo[1] * u.km / u.s, zo=orbigist.zo * u.kpc, lb=True)


orbit_forward.integrate(integrate_time, orbigist.pot, 'rk4_c')
orbit_backward.integrate(-integrate_time, orbigist.pot, 'rk4_c')
