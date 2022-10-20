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
def radecmatch(ra1, dec1, ra2, dec2):
    return Main.radecmatch(ra1, dec1, ra2, dec2)
def radecmatch_argmin(ra1, dec1, ra2, dec2):
    return Main.radecmatch_argmin(ra1, dec1, ra2, dec2)
def radecmatch_argmin_memlim(ra1, dec1, ra2, dec2):
    return Main.radecmatch_argmin_memlim(ra1, dec1, ra2, dec2)
