import os
from juliacall import Main
import windows_directories_new
import numpy as np

Main.include(os.path.join(windows_directories_new.jl_dir, "strmatch.jl"))
def strmatch(list1, list2):
    """
    Get indices of list 2 corresponding to list1 assuming 1:1 is perfect.
    :param list1: str array
    :param list2: str array
    :return: matches
    """
    return Main.strmatch(list1, list2)

def strmatch_general(list1, list2):
    """
    Get indices of list 2 corresponding to list1 allowing non-existence, i.e.

        - matches has the indices of matches in 2 s.t list2[matches] = list1
        - found is a booltrix

    :param list1: str array
    :param list2: str array
    :return: matches, found
    """
    return Main.strmatch_general(list1, list2)

Main.include(os.path.join(windows_directories_new.jl_dir, "radecmatch.jl"))

# Note that this has been updated since use with our APOGEE code- see the JL file
# TODO: See https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu3ast/sec_cu3ast_proc/ssec_cu3ast_proc_xmatch.html
# Haversine equation is better
# Consider KDTree implementation.
def radecmatch(ra1, dec1, ra2, dec2):
    return Main.radecmatch(ra1, dec1, ra2, dec2)
def radecmatch_argmin(ra1, dec1, ra2, dec2): # new argmin_resolution argument
    return Main.radecmatch_argmin(ra1, dec1, ra2, dec2, 1)
def radecmatch_argmin_memlim(ra1, dec1, ra2, dec2, argmin_resolution):
    """
    Juliacall helper function to carry out a cross-match using angular distances on-sky.
    Computationally cheaper to have cardinality of ra1 < cardinality of ra2 / |ra1|<|ra2|.

    Returns:
        matches:-- a list len(ra1) of tuples with index(ra1),index(ra2)

        distances:-- in milliarcseconds, of the match

        truefalse:-- a truefalse 1darray for selecting matches

    :param ra1: first catalogue RADIANS
    :param dec1:
    :param ra2: second catalogue RADIANS
    :param dec2:
    :param argmin_resolution: distance threshold MILLIARCSECONDS
    :return: matches, distances, truefalse
    """
    return Main.radecmatch_argmin_memlim(ra1, dec1, ra2, dec2, argmin_resolution)


# Alright. Bad. We need to try propagating the epoch??? Not really. Example of pygaia fork.
def propagate_astrometry(phi, theta, parallax, muphistar, mutheta, vrad, t0, t1):
    """
    Same as the PyGaia documentation. Propagate the astrometric parameters of a source from the reference epoch t0 to
    the new epoch t1. Note coordinates are ra-dec/lat-long style, not theta-phi style.

    Parameters
    ----------
    phi : float
        Longitude at reference epoch (radians).
    theta : float
        Latitude at reference epoch (radians).
    parallax : float
        Parallax at the reference epoch (mas).
    muphistar : float
        Proper motion in longitude (including np.cos(latitude) term) at reference
        epoch (mas/yr).
    mutheta : float
        Proper motion in latitude at reference epoch (mas/yr).
    vrad : float
        Radial velocity at reference epoch (km/s).
    t0 : float
        Reference epoch (Julian years).
    t1 : float
        New epoch (Julian years).

    Returns
    -------
    phi1, theta1, parallax1, muphistar1, mutheta1, murad1 : float or array
        Astrometric parameters, including the "radial proper motion" (NOT the radial
        velocity), at the new epoch.
    """
    from juliacall import Main
    Main.include(os.path.join(windows_directories_new.jl_dir, "propagate_epoch.jl"))
    phi1, theta1, parallax1, \
    muphistar1, mutheta1, murad1 = Main.propagate_astrometry(phi, theta, parallax,
                                                             muphistar, mutheta, vrad,
                                                             t0, t1)
    return phi1, theta1, parallax1, muphistar1, mutheta1, murad1
"""
jkgss['parallax'] = 1/jkgss['dist']
ra1, dec1, parallax1, pmra_cosdec1, pmdec1, pmrad1 = propagate_astrometry(np.deg2rad(jkgss['ra']),
                                                                          np.deg2rad(jkgss['dec']),
                                                                          jkgss['parallax'],
                                                                          jkgss['pmra']*np.cos(np.radians(jkgss['dec'])),
                                                                          jkgss['pmdec'],
                                                                          jkgss['vlos'],
                                                                          2008,2016)
"""