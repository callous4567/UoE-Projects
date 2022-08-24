import functools
from threading import Thread
import numpy as np


class utils(object):
    def __init__(self):
        holder = "Target"

    # Convert radians to arcminutes/sec
    def ra_arcmin(self, ra):
        return np.rad2deg(ra)*60
    def ra_arcsec(self, ra):
        return np.rad2deg(ra)*3600

    # Get distance + error from parallax + parallax error, in parsecs and MILLIARCSECONDS
    def par_to_dist(self, parallax, parallax_error):
        parallax, parallax_error = parallax*(1e-3), parallax_error*(1e-3)
        d = 1/parallax
        # 1/p^2 * dp
        d_err = parallax_error/(parallax**2)
        return d, d_err

    # Timeout caller.
    # To use, redefine a function as new_function = timeout(time_to_timeout)(old_function)
    # then use new_function as normal
    def timeout(self, seconds_before_timeout):
        def deco(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                res = [
                    Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, seconds_before_timeout))]

                def newFunc():
                    try:
                        res[0] = func(*args, **kwargs)
                    except Exception as e:
                        res[0] = e

                t = Thread(target=newFunc)
                t.daemon = True
                try:
                    t.start()
                    t.join(seconds_before_timeout)
                except Exception as e:
                    print('error starting thread')
                    raise e
                ret = res[0]
                if isinstance(ret, BaseException):
                    raise ret
                return ret

            return wrapper
        return deco

    # Converts RA:DEC string (22:22:22.333,44:43:20) to degrees and returns array [ra,dec]
    def radec_deg(self, string):
            ra, dec = string.split(",")[0], string.split(",")[1]
            ra, dec = ra.split(":"), dec.split(":")
            ra, dec = [float(d) for d in ra], [float(d) for d in dec]
            ra, dec = (15 * ra[0] + (15 / 60) * ra[1] + (15 / 3600) * ra[2]), (
                        dec[0] + (1 / 60) * dec[1] + (1 / 3600) * dec[2])
            return ([ra, dec])
            # Converts RA,DEC degree array to RA:DEC string. Returns it in the same format as radec_deg takes.

    def deg_radec(self, radecarray):
            ra, dec = radecarray[0], radecarray[1]
            # RA FIRST
            ra_hours, ra_minutes = int(ra / int(15)), ra % int(15)
            ra_minutes, ra_seconds = int(ra_minutes / (15 / 60)), ra_minutes % (15 / 60)
            ra_seconds = ra_seconds / (15 / 3600)
            # DEC SECOND
            dec_degrees, dec_minutes = int(dec), dec % 1
            dec_minutes, dec_seconds = int(dec_minutes / (1 / 60)), dec_minutes % (1 / 60)
            dec_seconds = dec_seconds / (1 / 3600)
            # Format a string
            print_format = ("{0:0.0f}:{1:0.0f}:{2:0.2f}, {3:0.0f}:{4:0.0f}:{5:0.2f}").format(ra_hours, ra_minutes,
                                                                                             ra_seconds, dec_degrees,
                                                                                             dec_minutes, dec_seconds)

            return print_format

    # Run radec_deg over array.
    def radec_degray(self, array):
        return [self.radec_deg(u) for u in array]

    # Treatise: owo = timeit.timeit(basic_func, number=100) (to time it: print owo.)
