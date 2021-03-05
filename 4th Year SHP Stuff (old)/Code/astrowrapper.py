# This file must be included with main_runetime.py for the latter to be viable, and kept in the same directory.
# Specify rootdir for the operating directory under which the file structure resides.
# Import Preamble includes packages that may/may not be used, or were used for testing purposes. Do not remove any T_T.
import nexusformat.nexus as nx
import IPython
import astropy
from MeanStars import MeanStars
import random
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import from_table, writeto
from astroquery.gaia import Gaia
import copy
from threading import Thread
import functools
import gc
import yaml
import matplotlib as mpl
from adjustText import adjust_text
from astropy.io.misc import hdf5
from astropy.io import fits
from astropy import wcs as wcs
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Table, Column
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.visualization import simple_norm
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import pandas
import h5py
import shutil
from astroquery.astrometry_net import AstrometryNet
from photutils import psf
from astropy.modeling import fitting
from photutils import background
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import SkyCircularAnnulus
from photutils import SkyCircularAperture
from photutils import datasets # Sample data to use for testing
from photutils import IRAFStarFinder
from photutils import DAOStarFinder
from photutils import aperture_photometry
import numpy as np
import matplotlib.pyplot as plt
import sys
import io
import os
import glob
import time
import scipy
import itertools
import multiprocessing
from matplotlib.patches import Patch

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
from skimage import transform
from scipy.ndimage.interpolation import shift
# Set AstrometryNet() API key.
ast = AstrometryNet()
ast.api_key = 'rkwvuedpwhdrllyp'
# Set GAIA login
gaia_username = "sstrasza"
gaia_password = "CENT@rino4657898"
# Set the root directory in WINDOWS format. This file is made for usage on PyCharm on a Windows OS.
rootdir = 'D:\Prog\pycharm\SHP2021\Code'
# This is a bunch of info on the file structure, for folders + directories, alongside the h5py file.
"""

DIRECTORIES
===========
Sci -> Cluster Folders + 'STD' Folder -> U B V -> CALIBRATED -> WCSDONE (which has the final files)
Cali -> U B V DARKS BIAS -> masterflatUBV / masterbias / masterdark


h5py/hdf5 "data.hdf5"
=====================
Groups -> Datasets 
------------------
111773/2039 -> depth_u,b,v depth_err_u,b,v alt_dif_u,b,v 

"""

# DISCLAIMER FOR REUSE.
"""
If you intend to reuse this code in part or in full, make reference to the original author and/or source.
Callicious/plasmolian@gmail.com. Proudly coded by Callicious! <3 
"""


# Holds various utilities for photometric analysis for TGP2020 Project.
# INCLUDES CALIBRATION
# INCLUDES OPTICAL DEPTH STUFF
# INCLUDES ASTROMETRY
# INCLUDES FILE UTILITIES (SO FAR)
class utils(object):
    def __init__(self):
        holder = "Target"
        self.directory = rootdir
        self.filename = "data.hdf5"

    # Convert apparent to absolute given a distance (parsecs)
    def abs_to_app(self, m, m_err, d, d_err):
        try:
            M = m - 5 * np.log10(d / 10)
            M_err = np.sqrt((m_err)**2 + (((np.log10(np.exp(1)))**2)*((5*d_err/d)**2)))
            return M, M_err
        except:
            return np.NaN, np.NaN

    # Get distance + error from parallax + parallax error, in parsecs and MILLIARCSECONDS
    def par_to_dist(self, parallax, parallax_error):
        try:
            parallax, parallax_error = parallax*(1e-3), parallax_error*(1e-3)
            d = 1/parallax
            # 1/p^2 * dp
            d_err = parallax_error/(parallax**2)
            return d, d_err
        except:
            return np.NaN, np.NaN


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

    # Extracts variables [VARIS] from 'name' FITS file in subdirectory DIRECTORY.
    def test_function(self):
        print("Success.")
    # Misc utility to identify all the IMMEDIATE subdirectories for a directory (with no files).
    def dir_lister(self, directory):
        os.chdir(directory)
        subdirs = os.listdir()
        subdirs = [rootdir + '\\' + directory + '\\' + d for d in subdirs]
        os.chdir(rootdir)
        return subdirs
    # Misc utility to identify all 'endpoint' subdirectories inside a directory. Directory. <3
    # Ex: run on 'Sci' directory, returns 'Sci\...other_roots. Endpoints have no subdirectories.
    def dir_sub_lister(self, directory):
        end_point_directories = []
        for directories, subdirectories, files in os.walk(directory):
            if subdirectories == []:
                to_add = directories
                is_true = int(0)
                for u in end_point_directories:
                    if to_add == u:
                        is_true += int(1)
                if is_true == int(0):
                    end_point_directories.append(to_add)
        return end_point_directories





    # Gets various statistics/histogram/etc for a FITS file background sigma clip.
    def stat_reduce(self, name):
        # Load in the file and reduce to data + header.
        fitso = fits.open(name)
        fitsodata, fitsoheader = fitso[0].data, fitso[0].header

        # Sigma-clipped Stats.
        # mean, median, sigma = sigma_clipped_stats(fitsodata)

        # Non-sigma Stats.
        # mean, median = np.mean(fitsodata), np.median(fitsodata)


        # Generate histograms about the median.
        sigma_data = sigma_clip(fitsodata, sigma=3, cenfunc='median', stdfunc='std')
        plt.hist(sigma_data, bins=30)
        #plt.show()
        plt.savefig("backogram.png",dpi=300)
        # Generate 1D y(x) for the 2D array.
        sigma_data_int = np.int32(sigma_data)
        max, min = np.max(sigma_data_int), np.min(sigma_data_int)
        # Generate y(x) values.
        x_vals = np.arange(min, max + 1, 1)
        y_vals = [0 for d in x_vals]
        shape = np.shape(sigma_data_int)
        for i in range(shape[0]):
            for j in range(shape[1]):
                value = sigma_data_int[i,j]
                if type(value) == np.int32:
                    value_corrected = value - min
                    y_vals[value_corrected] += 1
                else:
                    pass

        plt.plot(x_vals, y_vals, label='Background Median Histogram')
        #plt.show()
    # Varis array from filename name (fits)
    def fithead(self, name, varis):
        header = np.array([0])
        with fits.open(name) as f:
            header = f[0]

        # Variable Array
        vararray = [header.header[d] for d in varis]
        vararray = np.array(vararray)
        return vararray
    # Utility to get all current files in directory and fuck right off the non-fits files.
    def fits_identifier(self):
        files = os.listdir()
        new_files = []
        for file in files:
            filesplit = file.split(".")
            fitso = filesplit[len(filesplit) - 1]
            if fitso == "fits":
                new_files.append(file)
        return new_files

    # zeros all coordinates except those within the cuboid formed by the two square coordinates afforded (in the image)
    # Returns the data (explicitly) and saves the file for debug optionally.
    # Coords should be quoted pythonically in "XY FITS" format from Aladin
    def fittrim_coords(self, file, coord1, coord2):
        # Load data and shape.
        data, header = fits.open(file)[0].data, fits.open(file)[0].header
        datashape = np.shape(data)
        datashape = np.array(datashape)

        # Take coord1 and coord2, split into x,y coords. Rudimentary check to see if coordinates are messed up.
        xy1,xy2 = coord1.split(","),coord2.split(",")
        x1,x2,y1,y2 = xy1[0],xy2[0],xy1[1],xy2[1]
        x1,x2,y1,y2 = np.int(x1),np.int(x2),np.int(y1),np.int(y2)
        if x1 >= x2:
            while True:
                print("Coordinates are pogged.")
                time.sleep(1)
        if y1 >= y2:
            while True:
                print("Coordinates are pogged.")
                time.sleep(1)

        # Trim the data by rows (y) and columns (x), shape is (y,x), pythonically.
        for x in range(datashape[1]):
            for y in range(datashape[0]):
                if x >= x2:
                    data[y,x] = 0
                elif x <= x1:
                    data[y,x] = 0
                elif y <= y1:
                    data[y,x] = 0
                elif y >= y2:
                    data[y,x] = 0

        # Save data for debug.
        #newdata = np.int32(data)
        #newfitso = fits.PrimaryHDU(newdata)
        #newfitso.header = header
        #newfitsodata = fits.HDUList(newfitso)
        #newfitsodata.writeto("trimtest.fits", overwrite="True")

        # Return data
        return data

    # Same but just for data
    def fittrim_data(self, data, coord1, coord2):
        # Load data and shape.
        data = data
        datashape = np.shape(data)
        datashape = np.array(datashape)

        # Take coord1 and coord2, split into x,y coords. Rudimentary check to see if coordinates are messed up.
        xy1,xy2 = coord1, coord2
        x1,x2,y1,y2 = xy1[0],xy2[0],xy1[1],xy2[1]
        x1,x2,y1,y2 = np.int(x1),np.int(x2),np.int(y1),np.int(y2)
        if x1 >= x2:
            while True:
                print("Coordinates are pogged.")
                time.sleep(1)
        if y1 >= y2:
            while True:
                print("Coordinates are pogged.")
                time.sleep(1)

        # Trim the data by rows (y) and columns (x), shape is (y,x), pythonically.
        for x in range(datashape[1]):
            for y in range(datashape[0]):
                if x >= x2:
                    data[y,x] = 0
                elif x <= x1:
                    data[y,x] = 0
                elif y <= y1:
                    data[y,x] = 0
                elif y >= y2:
                    data[y,x] = 0

        # Save data for debug.
        #newdata = np.int32(data)
        #newfitso = fits.PrimaryHDU(newdata)
        #newfitso.header = header
        #newfitsodata = fits.HDUList(newfitso)
        #newfitsodata.writeto("trimtest.fits", overwrite="True")

        # Return data
        return data


    # Trims down an array around coordinate (x,y) by square radius (region_clip)
    def array_trim(self, data, xy, region_clip):
        xy = [int(d) for d in xy]
        data = np.array(data)
        newdata = data[xy[1] - region_clip:xy[1] + region_clip, xy[0] - region_clip: xy[0] + region_clip]
        return newdata


    # Take a small data clip, and use a local IRAF starfind to improve accuracy of positioning the centroid (compared to manual apertures w/ blanket inaccurate wcs.)
    # Sig-clipped stats method is also included but not used due to inaccuracy (ostensible)
    def array_gauss_phot(self, data_clip, aprad, apanin, apanout, debug_folder, debug_identifier, save_figs, default_fwhm):
        #norm = simple_norm(data_clip, 'sqrt', percent=99)
        #plt.imshow(data_clip, norm=norm, origin='lower', interpolation='nearest')
        #plt.show()
        # Run a preliminary starfind on the data (trimming the edge coords by 1/3th the image size)
        shape = np.shape(data_clip)[0]
        trimshape = int(shape/3)
        trimcoords = [trimshape, trimshape], [shape - trimshape, shape - trimshape]
        trimdata = self.fittrim_data(copy.deepcopy(data_clip), trimcoords[0], trimcoords[1])
        starfinder = IRAFStarFinder(threshold=5, fwhm=3, minsep_fwhm=2, sigma_radius=3, sharplo=0, roundlo=0,
                                    sharphi=999, roundhi=999, exclude_border=True)
        sourced = starfinder(trimdata)  # Astropy of sources.
        #sourced = starfinder(data_clip)

        # Sort the list in descending order + also an "if" incase none are found
        try:
            sourced.sort('flux', reverse='true')
            x, y, iraf_fwhm = sourced['xcentroid'][0], sourced['ycentroid'][0], default_fwhm #sourced['fwhm'][0]
        except:
            x, y, iraf_fwhm = int(shape / 2), int(shape / 2), default_fwhm


        # Deprecated.
        """
        # Trim down the data clip for the actual gauss fit.
        gauss_clip = int(0.6 * 0.5 * shape)

        # Try and fit the data clip
        timed_fit = self.timeout(15)(self.gauss_fit)
        try:
            gauss_params = timed_fit(self.array_trim(copy.deepcopy(data_clip), [x,y], gauss_clip))
            # Get FWHM
            stdev = gauss_params[3]
            fwhm = np.abs(stdev) * 2.355  # Approximately
        except:
            if save_figs == True:
                os.chdir(self.directory + "\\Debug\\" + debug_folder)
                norm = simple_norm(copy.deepcopy(data_clip), 'sqrt', percent=99)
                plt.imshow(copy.deepcopy(data_clip), norm=norm, origin='lower', interpolation='nearest')
                plt.savefig(str(debug_identifier) + "_ERROR.png")
            if iraf_fwhm != 0:
                fwhm = iraf_fwhm
            if iraf_fwhm == 0:
                fwhm = 4 # Safe "reasonable" value.
        """

        # DEBUG STUFF
        # norm = simple_norm(data_clip_gauss, 'sqrt', percent=99)
        # plt.imshow(data_clip_gauss, norm=norm, origin='lower', interpolation='nearest')
        # plt.scatter([gauss_params[0]], [gauss_params[1]], s=4, lw=2, color="red")
        # plt.show()
        # plt.close()

        # Quick check to see if aprad*fwhm >= size of clip. If it does, roughly scale it down to the "safe bet" (enable if using iraf starfind here)
        #if shape <= 2.2 * apanout * iraf_fwhm:
        #    iraf_fwhm = default_fwhm





        """
        # (SIGMA-CLIPPED MEAN BACKGROUND ESTIMATION (possibly? deprecated)
        # Do photometry
        aperture = CircularAperture((x, y), r=aprad * iraf_fwhm)
        annulus = CircularAnnulus((x, y), r_in=iraf_fwhm * apanin, r_out=iraf_fwhm * apanout)
        aperturinos = [aperture]
        phottable = aperture_photometry(copy.deepcopy(data_clip), aperturinos)
    
        # (SIGMA-CLIPPED MEAN BACKGROUND ESTIMATION
    
        # Create an annular mask to do some sigma clipping for the background.
        annular_mask = annulus.to_mask(method='center')
        annulus_data = annular_mask.multiply(copy.deepcopy(data_clip))
        mean, med, std = sigma_clipped_stats(annulus_data)  # Per pixel/unit area. Take std as the error in this.

        # Get fluxes + errors 
        apflux = np.abs(phottable['aperture_sum_0'][0])
        aparea = np.pi * ((aprad * iraf_fwhm) ** 2)
        anfluxinaparea = aparea * mean
        apflux_corrected = apflux - anfluxinaparea

        # Calculate errors
        apflux_err, anflux_err = np.sqrt(apflux), std * aparea
        apflux_corrected_error = np.sqrt(apflux_err ** 2 + anflux_err ** 2)
        """




        # FOR SIMPLE APERTURE METHOD (no sigma clip bs) possibly (?) deprecated
        # Do photometry
        aperture = CircularAperture((x, y), r=aprad * iraf_fwhm)
        annulus = CircularAnnulus((x, y), r_in=iraf_fwhm * apanin, r_out=iraf_fwhm * apanout)
        aperturinos = [aperture, annulus]
        phottable = aperture_photometry(copy.deepcopy(data_clip), aperturinos)

        # Calculate fluxes
        apflux, anflux = np.abs(phottable['aperture_sum_0'][0]), phottable['aperture_sum_1']
        aparea, anarea = np.pi * ((aprad * iraf_fwhm) ** 2), np.pi * ((apanout * iraf_fwhm) ** 2 - (apanin * iraf_fwhm) ** 2)
        apareadivanarea = aparea / anarea
        anfluxinaparea = anflux * apareadivanarea
        apflux_corrected = apflux - anfluxinaparea

        # Calculate errors
        apflux_err, anflux_err = np.sqrt(apflux), np.sqrt(np.abs(anflux))
        anfluxinaparea_err = anflux_err * apareadivanarea
        apflux_corrected_error = np.sqrt(apflux_err ** 2 + anfluxinaparea_err ** 2)

        # Create a plot and save it (if save_figs is true)
        if save_figs == True:
            plt.clf()
            os.chdir(self.directory + "\\Debug\\" + debug_folder)
            norm = simple_norm(copy.deepcopy(data_clip), 'sqrt', percent=99)
            plt.imshow(copy.deepcopy(data_clip), norm=norm, origin='lower', interpolation='nearest')
            #plt.show()
            plt.scatter([x], [y], s=4, color="green", lw=2)
            plt.xlim([0, shape])
            plt.ylim([0, shape])
            aperture.plot(color='red', lw=1, alpha=0.5)
            annulus.plot(color='pink', lw=1, alpha=0.5)
            plt.savefig(str(debug_identifier) + ".png")
            plt.close()

        # save_array = [apflux_corrected, apflux_corrected_error, fwhm]
        # filer = hdf5_writer(rootdir, "data.hdf5")
        # filer.write(datagroup, dataset, save_array)
        return apflux_corrected, apflux_corrected_error, iraf_fwhm




    # Trims single FIT image. (40 pix all sides, hard coded.)
    # Produces various statistics for the FITS file provided, including means/medians and sigma clipped histograms.
    # You'll need to edit this for what you want, specifically. It does not return anything.
    def fittrim(self, file):
        with fits.open(file, mode='update') as fitopened:
                newdata = fitopened[0].data[40:4056, 40:4056]
                newfitso = fits.PrimaryHDU(newdata)
                newfitso.header = fitopened[0].header
                newfitsodata = fits.HDUList(newfitso)
                newfitsodata.writeto(file + 'trimmed.fits', overwrite="True")
    # Trims down all FITS files within (parent) directory DIRECTORY (recursive). Verifies FITS validity.
    def fittrim_recursive(self, directory):
        oswalk = os.walk(directory, topdown=True)
        for root, directories, files in oswalk:
            for u in files:
                h = u.split(".")
                h = h[len(h) - 1]
                if str(h) == "fits":
                    os.chdir(root)
                    self.fittrim(u)
                    os.remove(u)
                    os.chdir(rootdir)
    # Generate flats/bias/darks. Run Standalone.
    # masterbias.fits, masterdark.fits in DARKS, BIAS
    # masterflat[u,b,v].fits in U B V
    # DARKS are for 60 seconds. Masterdark header has this information. Future scale for exposure.
    def calib_gen(self, bands):
        os.chdir("Cali")
        # Biasgen
        os.chdir('BIAS')
        filenames = os.listdir()
        totalstack = np.zeros(shape=(4016, 4016))
        for f in filenames:
            with fits.open(f) as file:
                indifit = np.float64(file[0].data)
                totalstack += indifit
        biasdata = totalstack / len(filenames)
        biasfits = fits.PrimaryHDU(np.float64(biasdata))
        biasfitshduldata = fits.HDUList([biasfits])
        biasfitshduldata.writeto('masterbias' + '.fits', overwrite='true')
        os.chdir("..")
        # Flatgen
        for d in bands:
            os.chdir(d)
            filenames = os.listdir()
            stacker = []
            for f in filenames:
                with fits.open(f) as file:
                    indifit = np.float64(file[0].data)
                    indifit -= biasdata
                    indifit = indifit/np.mean(indifit)
                    stacker.append(indifit)
            totalstack = np.median(stacker, axis=0)
            flatfits = fits.PrimaryHDU(np.float64(totalstack))
            flatfitshduldata = fits.HDUList([flatfits])
            flatfitshduldata.writeto('masterflat' + d + '.fits', overwrite='true')
            os.chdir("..")
        # Darkgen
        os.chdir('DARKS')
        filenames = os.listdir()
        totalstack = np.zeros(shape=(4016, 4016))
        for f in filenames:
            with fits.open(f) as file:
                indifit = np.float64(file[0].data)
                indifit -= biasdata
                totalstack += indifit
        totalstack = totalstack/len(filenames)
        darkfits = fits.PrimaryHDU(np.float64(totalstack))
        darkfits.header['EXPOSURE'] = '60'# Set exposure to that of the dark frame for calibration purposes.
        darkfitshduldata = fits.HDUList([darkfits])
        darkfitshduldata.writeto('masterdark' + '.fits', overwrite='true')
        os.chdir(rootdir)
    # Calibrates an image and saves it in a subdirectory: "CALIBRATED"
    def cali_img(self, name):
        fitso = fits.open(name) # Image Data
        fitso[0].data = np.float64(fitso[0].data) # Convert to float64

        darko = fits.open(rootdir + '\Cali\DARKS\masterdark.fits')
        biaso = fits.open(rootdir + '\Cali\BIAS\masterbias.fits')
        filter_band = fitso[0].header['FILTER']

        # Adjust dark frame for the exposure and apply DARK + BIAS
        darkojusted = darko[0].data * (float(fitso[0].header['EXPOSURE'])/float(darko[0].header['EXPOSURE']))
        newfitsodata = fitso[0].data - (darkojusted + biaso[0].data)

        # Locate and apply flatfield.
        flatto = fits.open((rootdir + '\Cali\{0}\masterflat{0}.fits').format(filter_band))
        newfitsodata = np.int32(np.divide(newfitsodata,flatto[0].data))

        # Make new file in "CALIBRATED" subdirectory
        newfitso = fits.PrimaryHDU(newfitsodata)
        newfitso.header = fitso[0].header
        newfitsodata = fits.HDUList(newfitso)
        # Delete old file.
        fitso.close()
        os.remove(name)
        try:
            os.chdir("CALIBRATED")
            newfitsodata.writeto(name + 'fixed.fits', overwrite="True")
        except:
            os.mkdir("CALIBRATED")
            os.chdir("CALIBRATED")
            newfitsodata.writeto(name + 'fixed.fits', overwrite="True")
    # Recursively runs cali_img for all FITS within directory DIRECTORY
    def recur_cali(self, directory):
        oswalk = os.walk(directory, topdown=True)
        for root, directories, files in oswalk:
            for u in files:
                h = u.split(".")
                h = h[len(h) - 1]
                if str(h) == "fits":
                    os.chdir(root)
                    self.cali_img(u)
                    os.chdir(rootdir)
    # Average-stacks all FITS images in CWD. Hard-coded 4016,4016 array.
    # Assumes already aligned images (hence same WCS)
    def average_stacker(self, identifier):
        # Identify files and CWD
        all_files = os.listdir()
        files = []
        filename = 'average_stack_' + identifier + '.fits'

        # Loop to avoid stacking the old image stack
        for d in all_files:
            if d == filename:
                pass
            else:
                # Also check if its actually a fits file.
                splitter = d.split(".")
                if splitter[len(splitter) - 1] != 'fits':
                    pass
                else:
                    files.append(d)


        # Stack images and write file.
        # GET THE SHAPE!

        shape = 0
        with fits.open(files[0]) as f:
            shape = np.shape(f[0].data)
        average_stack = np.zeros(shape, dtype=np.float64)

        for d in files:
            with fits.open(d) as f:
                average_stack += f[0].data


        # Write file and all that stuff.
        new_header = fits.open(files[0])[0].header
        average_stack = average_stack/len(files)
        avgstack = fits.PrimaryHDU(np.int32(average_stack))
        avgstack.header = new_header  # Set exposure to that of the dark frame for calibration purposes.
        avgstack_data = fits.HDUList([avgstack])
        avgstack_data.writeto(filename, overwrite='true')
    # Recursively runs average_stacker for all cluster subdirectories provided. Must be in array.
    # Alter for the average_stacker identifier argument (needs to be done/tdb)
    def average_stacker_recursive(self, cluster_identifiers):
        for d in cluster_identifiers:
            endpoints = self.dir_sub_lister(rootdir + '\\Sci\\' + d)
            for endpoint in endpoints:
                os.chdir(endpoint)
                self.average_stacker("stacked_up")


    # Visualize data and aperture sources. Data visual fails on calibrated frames due to negative/darkpix/etc.
    def source_visualizer(self, data, sources, fwhm, save_index):
        # [x,y] coordinates, all in a long array, i.e. [r1, r2, r3...]
        fig = plt.figure(figsize=(5,5), dpi=300)
        # Quick error catch for numpy ndarray fwhm.
        try:
            fwhm = np.sum(fwhm)
        except:
            pass
        xymarkers = np.transpose((sources['xcentroid'], sources['ycentroid']))
        apertures = CircularAnnulus(xymarkers, r_in=fwhm, r_out=2*fwhm)
        norm = simple_norm(data, 'sqrt', percent=99)
        plt.imshow(data, norm=norm, origin='lower', interpolation='nearest')
        apertures.plot(color='red', lw=2, alpha=0.6)
        #plt.show()
        if save_index != False:
            try:
                cwd = os.getcwd()
                os.chdir(rootdir + "\\Debug")
                fig.savefig(save_index + ".png")
                os.chdir(cwd)
            except:
                pass
    # Takes FITS data , not file, and returns astropy table of sources using IRAFSTARFIND. This is useful exclusively for astrometry and the like.
    # Do not use this for something that requires more precision, like the source extraction for data use.
    def source_extractor(self, data, threshold):
        # You might need to decrease the threshold dependent on how faint your image is.
        #print(np.min(data), np.max(data))
        mean,med,dev = sigma_clipped_stats(data, sigma=2.0) # For background/etc; mean:median:st_dev
        # FOR STANDARD STARS starfinder = IRAFStarFinder(threshold=20, fwhm=4, minsep_fwhm=10, sigma_radius=4)
        starfinder = IRAFStarFinder(threshold=threshold, fwhm=3, minsep_fwhm=2, sigma_radius=1.5, sharplo=0, roundlo=0, sharphi=999, roundhi=999) # FOR CLUSTERS
        sourced = starfinder(data) # Astropy of sources.

        # Sort the list in descending order.
        sourced.sort('flux', reverse='true')

        # Unhash if you want to see the sources you're handling while sextracting.
        # self.source_visualizer(data, sourced, 5, False) # fwhm by default set to 5 for this part.
        #print("Sources extracted successfully...")
        return sourced



    # Returns WCS_HEADER given source list. plasmolian@gmail.com API key.
    def wcsgetter(self, astrotable):
        # RA,DEC are the old RA and DEC, taken from the pointing data.

        print("Running WCS solve.")

        wcs_header = "null"
        try_again = True
        submission_id = None
        hw = 4016
        varidict = {
            'scale_upper':120,
            'scale_lower':2,
            'scale_units':'arcminwidth',
        }

        while try_again:
            try:
                if not submission_id:
                    print("Uploading...")
                    wcs_header = ast.solve_from_source_list(astrotable['xcentroid'], astrotable['ycentroid'], hw, hw, solve_timeout=1200, submission_id=submission_id, **varidict)
                else:
                    wcs_header = ast.monitor_submission(submission_id,
                                                        solve_timeout=1200)
            except:
                continue
            else:
                # got a result, so terminate
                try_again = False

        return wcs_header


        #wcs_header = ast.solve_from_source_list(astrotable['xcentroid'], astrotable['ycentroid'], hw, hw, solve_timeout=300)
        #return wcs_header
    # Returns WCS_HEADER given the entire FITS file to upload.
    def wcsgetter_fits(self, entire_file):
        # RA,DEC are the old RA and DEC, taken from the pointing data.
        hw = 4016
        wcs_header = "null"
        try_again = True
        submission_id = None
        varidict = {
            'scale_upper': 120,
            'scale_lower': 2,
            'scale_units': 'arcminwidth',
        }

        print("Running WCS solve. Uploading.")

        while try_again:
            try:
                if not submission_id:
                    print("Uploading initiating.")
                    wcs_header = ast.solve_from_image(entire_file, hw, hw,
                                                            solve_timeout=1200, submission_id=submission_id, **varidict)
                    print("Waiting.")
                else:
                    print("Waiting.")
                    wcs_header = ast.monitor_submission(submission_id,
                                                        solve_timeout=1200)
            except:
                print("Failed.")
                continue
            else:
                # got a result, so terminate
                try_again = False

        return wcs_header
    # Astrometry-calibrates a FITS file. New FITS retains EXPOSURE, FILTER, ALTITUDE.
    # Checks to see if old one exists. If it does, then skip.
    def astro_metrifier(self, name):
        print(name)
        # Check for existence of already-metrified file
        cwd = os.getcwd()
        tokeno = 0
        try:
            os.chdir("WCSDONE")
            file_list = self.fits_identifier()
            tokeno = 0
            for file in file_list:
                if file == name + "wcs.fits":
                    tokeno += 1
        except:
            pass


        if tokeno != 0:
            print("Exists. Passing " + name)
            pass
        else:
            os.chdir(cwd)
            old_header, old_data, wcs_header = "null", "null", "null"
            print(name)
            print()
            with fits.open(name) as f:
                old_data = f[0].data
                old_header = f[0].header
                print("Opened " + name)

            band = old_header['FILTER']
            threshold = "dolt"
            if band == 'U':
                threshold = 150
            elif band == 'B':
                threshold = 1000
            elif band == 'V':
                threshold = 1000
            elif band == 'Luminance':
                threshold = 150
            else:
                print("Band fucked up mate.")

            if old_header['OWNERID'] != "WCS_FITTED":
                with fits.open(name) as f:
                    # Use this to get the WCS via extracting sources and then uploading the astropy catalogue. Not 100% accurate, sadly.
                    data = f[0].data
                    print("Got data.")
                    sourcelist = self.source_extractor(data, threshold)
                    wcs_header = self.wcsgetter(sourcelist)
                    print("Gotten WCS. Plugging it in.")

                    # Use this to upload the entire file. Most accurate result, takes an obscene amount of time.
                    # cwd = os.getcwd()
                    # wcs_header = self.wcsgetter_fits(cwd + "\\" + name)
                new_data = fits.PrimaryHDU(old_data)
                new_data.header = wcs_header
                new_data.header['OWNERID'] = "WCS_FITTED"
                new_data.header['EXPOSURE'], new_data.header['ALTITUDE'], new_data.header['FILTER'] = old_header[
                                                                                                          'EXPOSURE'], \
                                                                                                      old_header[
                                                                                                          'ALTITUDE'], \
                                                                                                      old_header[
                                                                                                          'FILTER']
                new_data = fits.HDUList([new_data])
                try:
                    os.chdir("WCSDONE")
                    new_data.writeto(name + "wcs.fits", overwrite='true')
                    print("Done")
                except:
                    os.mkdir("WCSDONE")
                    os.chdir("WCSDONE")
                    new_data.writeto(name + "wcs.fits", overwrite='true')
                    print("Done")
            elif old_header['OWNERID'] == "WCS_FITTED":
                print("Already Fitted. Passing.")
            else:
                print("Header Error with" + name)

    # Uploads entire image instead of just doing with source extraction.
    def astro_metrifier_entire_file(self, name):
        print(name)
        # Check for existence of already-metrified file
        cwd = os.getcwd()
        tokeno = 0
        try:
            os.chdir("WCSDONE")
            file_list = self.fits_identifier()
            tokeno = 0
            for file in file_list:
                if file == name + "wcs.fits":
                    tokeno += 1
        except:
            pass


        if tokeno != 0:
            print("Exists. Passing " + name)
            pass
        else:
            os.chdir(cwd)
            old_header, old_data, wcs_header = "null", "null", "null"
            print(name)
            print()
            with fits.open(name) as f:
                old_data = f[0].data
                old_header = f[0].header
                print("Opened " + name)

            band = old_header['FILTER']
            threshold = "dolt"
            if band == 'U':
                threshold = 150
            elif band == 'B':
                threshold = 1000
            elif band == 'V':
                threshold = 1000
            else:
                print("Band fucked up mate.")

            if old_header['OWNERID'] != "WCS_FITTED":
                with fits.open(name) as f:
                    # Use this to get the WCS via extracting sources and then uploading the astropy catalogue. Not 100% accurate, sadly.
                    #data = f[0].data
                    #print("Got data.")
                    #sourcelist = self.source_extractor(data, threshold)
                    #wcs_header = self.wcsgetter(sourcelist)
                    #print("Gotten WCS. Plugging it in.")

                    # Use this to upload the entire file. Most accurate result, takes an obscene amount of time.
                    cwd = os.getcwd()
                    wcs_header = self.wcsgetter_fits(cwd + "\\" + name)
                new_data = fits.PrimaryHDU(old_data)
                new_data.header = wcs_header
                new_data.header['OWNERID'] = "WCS_FITTED"
                new_data.header['EXPOSURE'], new_data.header['ALTITUDE'], new_data.header['FILTER'] = old_header[
                                                                                                          'EXPOSURE'], \
                                                                                                      old_header[
                                                                                                          'ALTITUDE'], \
                                                                                                      old_header[
                                                                                                          'FILTER']
                new_data = fits.HDUList([new_data])
                try:
                    os.chdir("WCSDONE")
                    new_data.writeto(name + "wcs.fits", overwrite='true')
                    print("Done")
                except:
                    os.mkdir("WCSDONE")
                    os.chdir("WCSDONE")
                    new_data.writeto(name + "wcs.fits", overwrite='true')
                    print("Done")
            elif old_header['OWNERID'] == "WCS_FITTED":
                print("Already Fitted. Passing.")
            else:
                print("Header Error with" + name)

    # Recursively runs astro_metrifier for all images in subdirectories.
    def astro_recursifier(self, directory):
            oswalk = os.walk(directory, topdown=True)
            for root, directories, files in oswalk:
                if files != []:
                    for u in files:
                        h = u.split(".")
                        h = h[len(h) - 2]
                        if str(h) == "fitsfixed":
                            os.chdir(root)
                            self.astro_metrifier(u)
                            os.chdir(rootdir)
                            os.chdir(root)
                            os.remove(u)
                            os.chdir(rootdir)

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


    # Does photometry for a ["RA,DEC"]. Must be encased by a list.
    # RA:DEC should be given as a string, i.e. ["15:27:16.5,60:23:05.6"]
    # Returns AstroPy table with 'aperture_sum_0_(corrected)' + 'backapsu' + 'aperture_sum_1'
    # Inst count (corrected), uncorrected, background count, annular background count respectively
    # Also returns error (aperture_sum_0_corrected_error in the astropy table) based on poisson errors, i.e.
    # FILTER and ALTITUDE arguments are also returned inside the astropy table.
    # THIS CODE IS SPECIFICALLY MADE FOR A SINGLE STAR TARGET TO BE USED, NOT MULTIPLE (although it should *technically* work.)
    def wcs_photom(self, name, coordinates, clip, threshold, aprad, anradin, anradout):
        # Load data + headers.
        fitso = fits.open(name)
        fitsodata, fitsoheader = fitso[0].data, fitso[0].header
        directory = os.getcwd()

        # The scale is about 0.6302739312067 for the PIRATE sensor (DATA)
        pix_scale = 1 / 0.6302739312067

        # Set up WCS
        fitsowcs = wcs.WCS(fitsoheader)

        try:
            fwhm = self.star_profiler(name, coordinates[0], clip, threshold, aprad, anradin, anradout)

            # Tab Code 6958
            aprad_radius, anradin_radius, anradout_radius = aprad*fwhm, anradin*fwhm, anradout*fwhm
        except RuntimeWarning:
            print("FWHM solution failed. See wcs_photom to rectify error.")

        # Convert coords + generate aperture centers
        coordinates = self.radec_degray(coordinates)

        # Generate apertures and annuli
        positions = SkyCoord(coordinates, unit='deg')
        apertures = SkyCircularAperture(positions, r=aprad_radius * u.arcsec)
        annuli = SkyCircularAnnulus(positions, r_in=anradin_radius*u.arcsec, r_out=anradout_radius*u.arcsec)
        combined = [apertures, annuli]

        # Photometry using simple mean of local background. [aperture_sum_0,1] are the aperture and annular sums respec.
        photdata = aperture_photometry(fitsodata, combined, wcs=fitsowcs)

        # Area calculation for the annulus
        anarea = np.pi*(anradout_radius**2 - anradin_radius**2)/(pix_scale**2)
        backmean = photdata['aperture_sum_1']/anarea
        aparea = np.pi*(aprad_radius**2)/(pix_scale**2)
        photdata['backapsu'] = backmean * aparea

        # Subtract the background from the main apertures for the final sum. Attach ALTITUDE and FILTER data.
        photdata['aperture_sum_0_corrected'] = photdata['aperture_sum_0'] - photdata['backapsu']
        photdata['ALTITUDE'],photdata['FILTER'], photdata['EXPOSURE'] = fitsoheader['ALTITUDE'], fitsoheader['FILTER'], fitsoheader['EXPOSURE']
        photdata['ZENITH'] = 90 - photdata['ALTITUDE']
        # Error calculation
        # Error over any aperture count = SQRT(count over aperture)
        # Propagation is given here.
        # CORRECTED = APERTURE - ANNULAR*AREA_SCALE_FACTOR
        # CORRECTED_ERROR = SQRT(APERTURE - ANNULAR*AREA_SCALE_FACTOR) = SQRT(CORRECTED)
        # Easy enough since error in count is the sqrt of that count, via poisson statistics.
        photdata['aperture_sum_0_corrected_error'] = np.sqrt(photdata['aperture_sum_0_corrected'])
        photdata['fwhm'] = fwhm
        os.chdir(directory)
        return photdata



    # Returns [depths, errors, combinatoric_alt_difference, band, fluxes] for a directory of FITS files.
    # The directory should lead all the way from 'Sci'. Fluxes is all the fluxes used for airmass solve.
    # Give it a (parent) directory and aperture/annulus radius in arcseconds.
    # The RADEC should be in standard format inside a list, i.e. ['88:88:88,44:44:44', -....] where -.... is absent.
    # # Tab Code 6958 TO SPECIFY THE ANNULI/APRAD FOR THE AIRMASS EXTRACTION/ETC
    def airmass_solver(self, directory, RADEC, clip, threshold, aprad, anradin, anradout):
        # Navigate to the directory from the main folder and extract files.
        os.chdir(directory)

        # Gather the STD star identifier for the target.
        dirsplit = directory.split("\\")
        band = dirsplit[3]

        # Identify Files
        files = os.listdir()

        # Only use FITS.
        corrected_files = []
        for d in files:
            d_cor = d.split('.')
            if d_cor[len(d_cor) - 1] == 'fits':
                corrected_files.append(d)
        files = corrected_files

        # Extract Astropy tables for all stars.
        astrodata = [self.wcs_photom(d, RADEC, clip, threshold, aprad, anradin, anradout) for d in files]

        # Gather up the FWHM's, band + zeniths for statistics.
        fwhms = [d['fwhm'][0] for d in astrodata]
        zeniths = [d['ZENITH'][0] for d in astrodata]


        # Generate flux (corrected) list for "return"
        astrofluxes = [d['aperture_sum_0_corrected'][0] for d in astrodata]

        # Generate combinatorics for all possible combinations of airmass.
        base_list = np.arange(0,len(astrodata))
        base_combis = list(itertools.combinations(base_list, 2))

        # Define lambda funcs for optical depth and optical error.
        optidepth = lambda I_1,I_2,AM_1,AM_2: (np.log(I_1/I_2))/(AM_2 - AM_1)
        optierror = lambda I_1,I_2,AM_1,AM_2,ERR_1,ERR_2: (((ERR_1/I_1)**2 + (ERR_2/I_2)**2)/((AM_2 - AM_1)**2))**0.5

        # Just quickly add some empty lists for the variables.
        depths, errors, alt_difference = [],[],[]
        for d in base_combis:
            phottables = [self.wcs_photom(files[u], RADEC, clip, threshold, aprad, anradin, anradout) for u in d]
            fluxes, flux_errors, altitudes, filters = [d['aperture_sum_0_corrected'][0] for d in phottables], [d['aperture_sum_0_corrected_error'][0] for d in phottables], [d['ALTITUDE'][0] for d in phottables], [d['FILTER'][0] for d in phottables]
            zenith_angles = [np.deg2rad(90 - d) for d in altitudes]
            secant_airmass = [1/np.cos(d) for d in zenith_angles]
            #true_airmass = lambda airmass: airmass - (0.0018167)*(airmass-1) - 0.002875*((airmass - 1)**2) - 0.0008083*((airmass - 1)**3)
            #alt_true_airmass = lambda airmass: airmass*(1 - 0.0012*(airmass**2 - 1)
            # Photometry Cookbook = Sauce.
            depth = optidepth(fluxes[0], fluxes[1], secant_airmass[0], secant_airmass[1])
            depth_error = optierror(fluxes[0], fluxes[1], secant_airmass[0], secant_airmass[1], flux_errors[0], flux_errors[1])
            depths.append(depth)
            errors.append(depth_error)
            # Calculates the altitude difference for the combinatoric, in degrees.
            altidif = altitudes[1] - altitudes[0]
            alt_difference.append(np.absolute(altidif))
        # Astropy LaTeX table generation.
        """
        # Generate an Astropy Q-Table with all the variables/depths/etc.
        table = Table([depths, errors, combinatoric, band],
                       names=('optical_depth', 'depth_error', 'combi_tuple', 'FILTER'),
                       meta={"title": directory + 'depths'})
        try:
            os.remove('depth_data.dat')
            table.write('depth_data.csv')
            os.chdir(rootdir)
        except:
            table.write('depth_data.csv', format='ascii')
            os.chdir(rootdir)
        """
        return [depths, errors, alt_difference, band, astrofluxes, fwhms, zeniths]
    # Quick utility to set the optical depth values from the standard "IDENTIFIER" to be the "true" optical depth
    def true_depth(self, identifier):
        band_indices = ['U','B','V']
        writer = hdf5_writer(rootdir, 'data.hdf5')
        for band in band_indices:
            value = writer.read(identifier, 'optical_depth_' + band)
            writer.write('true_depth', 'optical_depth_' + band, value)

    # quick utility to set the z_p values from the standard "IDENTIFIER" to be the "true" zero point
    def true_zpers(self, identifier):
        band_indices = ['U','B','V']
        writer = hdf5_writer(rootdir, 'data.hdf5')
        for band in band_indices:
            value = writer.read(identifier, 'zero_point_' + band)
            writer.write('true_zero_point', 'zero_point_' + band, value)

    # Quick statistics for the zeropoints, instmags, etc, for fits_zp_solver_recursive. If squaresum_true is true, error is done via sqrt(sum_of_errs) and not stdev on the values.
    def zp_rec_stats(self, zeros, zeroerrs, insts, insterrs, squaresum_true):
        # Averages
        zero, inst = np.average(zeros), np.average(insts)

        # Errors
        if squaresum_true == False:
            zero_err, inst_err = np.std(zeros), np.std(insts)

        else:
            zero_err, inst_err = np.sqrt(np.sum([d**2 for d in zeroerrs])), np.sqrt(np.sum([d**2 for d in insterrs]))

        return zero, zero_err, inst, inst_err


    # Statistics on sets of bands, fwhms, zeniths (for a given standard star) with save directory specified inside "Debug"
    # Ripped from https://stackoverflow.com/questions/34280444/python-scatter-plot-with-multiple-y-values-for-each-x
    def standard_fwhmstats(self, bands, fwhms, zeniths, standard):
        x_values = 0,1,2

        texts = []
        # Generate texts
        for x_val in x_values:
            zenith_vals = zeniths[x_val]
            fwhm_vals = fwhms[x_val]
            for num, d in enumerate(zenith_vals):
                texts.append(plt.text(x=x_val, y=fwhm_vals[num], s=("{0:.1f}").format(d)))

        plt.title(("FWHM vs. Band for {}, zenith in deg.").format(standard))

        for xv, yv in zip(x_values, fwhms):
            plt.scatter([xv] * len(yv), yv, s=4, lw=2)

        # Sort out figure
        plt.xticks([0, 1, 2])
        plt.axes().set_xticklabels(bands)

        plt.grid(True, which='major', color="blue", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        plt.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='green', alpha=0.8, lw=0.5))


        try:
            os.mkdir(rootdir + "\\" + "Debug" + "\\" + standard)
        except:
            pass

        try:
            plt.savefig(rootdir + "\\" + "Debug" + "\\" + standard + "\\" + standard + "_fwhmfit.png")
        except:
            pass

        plt.clf()
        plt.close()


    # Statistical reduction for the airmass values (not very mathematical... just preliminary...)
    """
    - clip negative values
    - clip values >= 1 (just for safe measure) and also <= 0.2 (based on a "visual guesstimate" of what constitutes a good value.)
    - clip values >= some_error_fraction 
    """
    # Some thoughts on reduction...
    """
    - new values are a,b,c with errors da,db,dc
    - take average: AVG = (a + b + c)/3
    - take error: D_AVG = SQRT((da/3)^2 + (db/3)^2 + (dc/3)^2) 
    """
    def airmass_stats(self, arr, err, clip_range, err_frac, maxdev):
        arr_corr = []
        err_corr = []
        for num, u in enumerate(err):
            if arr[num] <= clip_range[0]:
                pass
            elif arr[num] >= clip_range[1]:
                pass
            elif err[num]/arr[num] >= err_frac:
                pass
            else:
                arr_corr.append(arr[num]), err_corr.append(err[num])
        length, mean, std_dev = 0,0,0
        try:
            length, mean, std_dev = len(arr_corr), np.mean(arr_corr), np.std(arr_corr)
        except RuntimeWarning:
            print("RuntimeWarning #1 for airmass_stats")
        err_mean = np.sqrt(np.sum([(d/3)**2 for d in err_corr]))
        new_arr_corr, new_err_corr = [],[]
        for num, u in enumerate(arr_corr):
            if np.abs(u - mean) >= maxdev:
                pass
            elif np.abs(u - mean) < maxdev:
                new_arr_corr.append(u), new_err_corr.append(err_corr[num])
        try:
            new_mean, new_std_dev, new_length = np.mean(new_arr_corr), np.std(new_arr_corr), len(new_arr_corr)
            new_err_mean = np.sqrt(np.sum([(d/3)**2 for d in new_err_corr]))
            return new_arr_corr, new_err_corr, new_length, new_mean, new_std_dev, new_err_mean
        except RuntimeWarning:
            print("airmass_stats has encountered a RuntimeWarning, #2, likely due to no suitable depth data. Passing.")
            new_length, new_mean, new_std_dev, new_err_mean = len(new_arr_corr), 0,0,0
            return new_arr_corr, new_err_corr, new_length, new_mean, new_std_dev, new_err_mean



    # Returns the atmos corrected FLUX, FLUX_ERR, INST_MAGNITUDE, INST_MAGNITUDE_ERR for file target (singular).
    # ERROR PROPAGATION:
    # I_corr = I * exp(optical_depth * airmass)
    # I_corr_I = (dI) * exp(optical_depth * airmass)
    # I_corr_DEPTH = (dD * airmass) * (I * exp(optical_depth * airmass))
    # I_corr_ERR = SQRT(sum of above^2)
    # INST_MAG = -2.5log10(I_corr)
    # INST_MAG_ERR = -2.5 * 1/I_corr * dI_corr * log_10(e)
    # Also returns the filter band.
    # Note that threshold specifies threshold for source extraction.
    def star_parametrifier(self, name, coordinates, clip, threshold, STD_identifier, aprad, anradin, anradout):
        # Get photometric data
        photomdata = self.wcs_photom(name, coordinates, clip, threshold, aprad, anradin, anradout)
        flux, flux_err, alti = photomdata['aperture_sum_0_corrected'][0], photomdata['aperture_sum_0_corrected_error'][0], photomdata['ALTITUDE'][0]


        # Obtain exposure time and relevant band header.
        exposure, band = "null", "null"
        with fits.open(name) as fitso:
            exposure = fitso[0].header['EXPOSURE']
            band = fitso[0].header['FILTER']

        # Get the optical depth and associated error.
        hdf = hdf5_writer(rootdir, 'data.hdf5')
        optical_depth, optical_depth_error = hdf.read('true_depth', 'optical_depth_' + band)

        # Return fluxes and magnitudes, fluxes normalized to the exposure period.
        airmass = 1 / np.cos(np.deg2rad(90 - alti))

        #print(exposure, airmass, flux, optical_depth, band)

        flux_atmospherically_done = flux * np.exp(optical_depth * airmass) * 1/exposure
        flux_atmospherically_done_error = (1/exposure) * np.sqrt((flux_err * np.exp(optical_depth * airmass))**2 + (optical_depth_error*airmass*flux_atmospherically_done)**2)

        inst_mag = -1*(5/2)*np.log10(flux_atmospherically_done)
        inst_mag_error = -1*(5/2)*np.log10(np.exp(1))*flux_atmospherically_done_error/flux_atmospherically_done
        return [flux_atmospherically_done, flux_atmospherically_done_error, inst_mag, inst_mag_error, photomdata['FILTER'][0]]
    # RADEC IS THE INDIVIDUAL STRING HERE! NOT AN ARRAY! Breaking with convention, baby!
    # You should make sure the star is not in a crowded region. Alter "Threshold" and ensure it's the brightest.
    # Verify in the produced visualization that it is indeed what you want.
    # Clip is the +- clip around the object. I.e. clip = 20 gives 40x40 rect aperture.
    # Threshold defines threshold for source extraction (above which the source is considered)
    def star_profiler(self, name, radec, clip, threshold, aprad, anradin, anradout, wcs=wcs):
        # Load file.
        fitso = fits.open(name)
        fitsodata, fitsoheader = fitso[0].data, fitso[0].header

        # Extract X,Y coordinate
        wcs = wcs.WCS(fitsoheader)
        radec_m = self.radec_deg(radec)
        radec_sky = astropy.coordinates.SkyCoord(ra=radec_m[0],dec=radec_m[1],unit='deg')
        pixy = astropy.wcs.utils.skycoord_to_pixel(radec_sky, wcs=wcs)
        pixy = [int(d) for d in pixy] # 0,1 is X,Y

        # Initial clip for stars.
        clip = fitsodata[pixy[1] - clip:pixy[1] + clip, pixy[0] - clip:pixy[0] + clip] # Numpy slicing: first Y, then X. Matrically! (row, then column)


        # Run the IRAF star finder on the image to estimate the centroid.
        sources = self.source_extractor(clip, threshold)

        # Crop the clip even more to increase the accuracy of the centering.
        pixynew = np.array([sources['xcentroid'], sources['ycentroid']])
        pixynew = [int(d) for d in pixynew]
        manuclip = 25
        clip = clip[pixynew[1] - manuclip: pixynew[1] + manuclip, pixynew[0] - manuclip:pixynew[0] + manuclip]
        plt.imshow(clip)
        #plt.show()

        # Pass the data and various boundary conditions to the solver and await results.
        try:
            debug_limiter = name.replace("_", "")
            solved = self.gauss_fit(clip)
            # Calculate and return the FWHM for the data (to be used for photometry/etc)
            fwhm = solved[3] * 2.355 # (approximation)

            # More debug
            # Try to make directory for debug files.
            try:
                os.mkdir(rootdir + "\\" + "Debug")
            except:
                pass

            # Bit of debug
            cwd = os.getcwd()
            os.chdir(rootdir + "\\Debug")
            filer = hdf5_writer(rootdir, "data.hdf5")
            filolial = name.split("OSL_ROE")[0]
            filolial = filolial.replace("_", "")
            filer.write("wcsphotom_fwhm", filolial, fwhm)

            plt.imshow(clip)
            aperture = CircularAperture((solved[0], solved[1]), r=aprad*fwhm)
            annulus = CircularAnnulus((solved[0], solved[1]), r_in=anradin*fwhm, r_out=anradout*fwhm)
            aperture.plot(color="w", lw=1)
            annulus.plot(color="red", lw=1)
            plt.text(1, 1, fwhm)
            cwd = os.getcwd()
            os.chdir(rootdir + "\\Debug")
            plt.savefig(name + ".png")
            plt.close()

            os.chdir(cwd)
            return fwhm
        except RuntimeWarning:
            print("FWHM solve encountered a RuntimeWarning. See star_profiler.")
    # Generate 1D gaussian slice y(x) from [-max_x,+max_x]
    def oned_gauss(self, peak, sigma, max_x):
        gausser = lambda x: peak*np.exp(-1*((x/sigma)**2))
        x_vals = np.arange(-1*max_x, max_x, 0.001)
        y_vals = [gausser(x) for x in x_vals]

        # The half-point of a gaussian occurs where (x/sigma) = log(2)
        # hence where x = sigma*log(2)

        x_half = sigma*np.sqrt(np.log(2))

        self.grapher(x_vals,y_vals, [-x_half,x_half])
    # Generate z(x,y), assumes sigma_x = sigma_y. Single value or run for tupled meshgrid.
    def twod_gauss(self, gridcoord, x_0, y_0, peak, sigma):
        x,y = gridcoord
        uwu = peak * np.exp(-1 * ((x - x_0)**2 + (y - y_0)**2) / (2*(sigma**2)))
        return uwu
    # Create 2D gaussian with shape SHAPE and other *params. Returns [X,Y,Z]
    def threed_gauss(self, x_0, y_0, peak, sigma, shape):
        x,y = np.arange(0, shape[1], 1), np.arange(0, shape[0], 1)
        X,Y = np.meshgrid(x,y)
        Z = self.twod_gauss((X,Y), x_0, y_0, peak, sigma)
        plt.imshow(Z)
        #plt.show()
        return X,Y,Z
    # Fancy Grapher for Gaussians (alter variables or make it generic later. tbd.)
    def grapher(self, x, y, v_lines):
        fig = plt.figure(1)  # Define a figure as fig.
        ax = fig.add_subplot(1, 1, 1)  # Defines our axis within this figure.
        ax.set(title=("Graph for {}").format("FWHM of the target star fitted gaussian."), ylabel="(Fitted) Intensity)",
               xlabel="Radius in Pix")  # Defines variabels of our axis. r'$\mu$'
        # We used plt.xlabel and plt.ylabel, or plt.title, and plt.axis([xstart, xfinish]) and etc.... this just puts all that into ax.set, setting the variables of ax, hence our subplot.
        # We can also plot on this axis. So, example. ax.plot([x], [y], label, etc) and ax.scatter...
        ax.plot(x, y, label="Scatter", lw=1)
        # We can also set xticks and yticks. ax.set(xticks=[], yticks=[]) sets ticks at each point defined, i.e. on the sides, but doesn't set the grid. So, a if xticks=[1] you get one singular tick on the entire x axis, only at one, LABELLED as 1.
        # So... to get custom minor ticks on our axis...
        ax.minorticks_on()  # Enables the minor ticks.
        ax.grid(True, which='minor', color='r', alpha=0.2)  # Enables minor ticks on the grid
        ax.grid(True, which='major', color='g', alpha=0.2)  # Enables major ticks on the grid
        # To set custom minor and major locators... First define the locator(s) using pltick, matplotlib.ticker!

        # Add vertical lines.
        for d in v_lines:
            plt.axvline(d)
        #plt.show()
    # Fit a gaussian to a 2D numpy array.
    # Bounds should be nested list, with each sublist having the lower/upper limit. peak and sigma bound should be singular floats.
    # Returns [X_0, Y_0, PEAK_VALUE, STDEV] fitted.
    def gauss_fit(self, data):
        # Generate the requisite meshgrid/coordinate-tuples.
        shapez = np.shape(data)
        x,y = np.arange(0, shapez[1], 1), np.arange(0, shapez[0], 1)
        X,Y = np.meshgrid(x,y)
        Z = data

        # Sap in an initial guess for the actual values of the boundaries.
        guess = [shapez[1]/2, shapez[0]/2, np.max(Z), 3]
        # Infinite Boundaries baby.
        boundaries = [-np.inf, np.inf]

        popt, pcov = scipy.optimize.curve_fit(self.twod_gauss, (X.ravel(),Y.ravel()), Z.ravel(), p0=guess, bounds=boundaries, maxfev=100)

        # X_0,Y_0,PEAK,STDEV
        return popt


    # Generates the zero-point (and associated error) for a given flux value.
    # Returns: [z_p, err, inst_mag, err]
    def zp_solver(self, corrected_flux, corrected_flux_error, true_mag):
        # Get inst_mag + error on inst_mag
        inst_mag = -(5/2)*np.log10(corrected_flux)
        inst_mag_err = -2.5*corrected_flux_error*np.log10(np.exp(1))/corrected_flux
        """
        inst_mag = -2.5log_10(FLUX) 
        log_e(flux) = log_10(flux)/log_10(e) 
        log_10(flux) = log_10(e) * log_e(flux) 
        d(log_10(flux))/d(flux) = (1/flux) * log_10(e)
        hence error in inst_mag is:
        -2.5 * (dFLUX/FLUX) * log_10(e) 
        """
        # true = inst + z_p
        z_p = true_mag - inst_mag
        # error = error in inst_mag (since true_mag error is assumed very small, ~0)
        return z_p, inst_mag_err
    # Takes a FITS, RADEC, CLIp/THRESHOLD FROM LEAST SQUARES GAUSS; generates the instrumental magnitude and zeropoint + errors.
    # Gets z_p and error, inst_mag and error, filter_band, for a given FITS file + target. Requires airmass calculations to be done.
    def fits_zp_solver(self, name, coordinates, clip, threshold, std_identifier, aprad, anradin, anradout):
        # Generate photometry data.
        data = self.star_parametrifier(name, coordinates, clip, threshold, std_identifier, aprad, anradin, anradout)
        # Extract the band.
        band = data[4]
        # Get the true_magnitude from data.
        file_utility = hdf5_writer(rootdir, 'data.hdf5')
        true_mag = file_utility.read(std_identifier, 'true_mag_' + band)
        # Generate zero_point
        z_p = self.zp_solver(data[0], data[1], true_mag)
        z_p = [d for d in z_p]

        return [z_p[0], z_p[1], data[2], data[3], data[4], data[0], data[1]]
    # Runs fits_zp_solver for all files in directory.
    # Writes all inst_mags, errors, z_p, errors, to file, under:
    # inst_mags_u,b,v, inst_mags_errs_u,b,v, z_ps_u,b,v, z_ps_errs_u,b,v
    # You don't need to specify the band, but the magnitude should be band specific for this particular identifier.
    # VERY BASIC STDDEV/AVERAGING FOR THE ERRORS ON THE ZEROPOINT FINAL VALUE!!!
    def fits_zp_solver_recursive(self, STD_identifier, gamma, aprad, anradin, anradout):
        # Instantiate reader and utils.
        utility, hdf = utils(), hdf5_writer(rootdir, 'data.hdf5')

        # Get files and also the coordinates/clip/threshold
        file_list_all = os.listdir()
        file_list = []
        for u in file_list_all:
            uwu = u.split('.')
            if uwu[len(uwu) - 1] == 'fits':
                file_list.append(u)
        radec, clip, threshold = hdf.read(STD_identifier, 'radec'),hdf.read(STD_identifier, 'clip'),hdf.read(STD_identifier, 'threshold')

        # Convert radec back to RA,DEC format, as wcs_photom takes it.
        radec = [utility.deg_radec(radec)]

        # Quick fix to the fact that file reader exports an ndarray. Clean this up later.
        threshold = np.mean(threshold)

        # Get back to directory (sorry, side-product of my file editing method. We're stuck with it now :/
        os.chdir(rootdir + '\\' + gamma)

        # Generate z_p for all folder members.
        zperrinstmagerr = []
        for d in file_list:
            array = self.fits_zp_solver(d, radec, clip, threshold, STD_identifier, aprad, anradin, anradout)
            zperrinstmagerr.append(array)
            os.chdir(rootdir + '\\' + gamma)

        # Transpose the zperrinstmagerr array.
        zperrinstmagerr = np.array(zperrinstmagerr)
        zperrinstmagerr_trans = zperrinstmagerr.T

        # Note the z_ps, z_perrs, inst_mags, inst_mag_errs and another quick fix for transpose making strings.
        zeros, zero_errors, insts, inst_errors, bands, fluxvals, fluxerrs = zperrinstmagerr_trans
        zeros, zero_errors, insts, inst_errors, fluxvals, fluxerrs = [np.float(d) for d in zeros], [np.float(d) for d in zero_errors], [np.float(d) for d in insts], [np.float(d) for d in inst_errors], [np.float(d) for d in fluxvals], [np.float(d) for d in fluxerrs]


        # Write down all the above information.
        hdf.write(STD_identifier, 'zeropoints_' + bands[0], zeros)
        hdf.write(STD_identifier, 'zeropoints_' + bands[0] + '_errors', zero_errors)
        hdf.write(STD_identifier, 'inst_mags_' + bands[0], insts)
        hdf.write(STD_identifier, 'inst_mags_' + bands[0] + '_errors', inst_errors)

        # Also run some basic statistics on the z_p's. CURRENTLY USING STANDARD DEVIATION. You can alter for square-sum errors.
        statzero, statzeroerr, statinst, statinsterr = self.zp_rec_stats(zeros, zero_errors, insts, inst_errors, False)
        zeropoint, inst = [statzero, statzeroerr], [statinst, statinsterr]

        # Write statistically-reduced zero and instmag
        hdf.write(STD_identifier, 'zero_point_' + bands[0], zeropoint)
        hdf.write(STD_identifier, 'inst_mag_' + bands[0], inst)

        # Also for debug purposes, write down the fluxes/errors for the star that was z_p solved. The flux/error is normalized.
        #print(bands[0])
        hdf.write(STD_identifier, 'zp_flux_' + bands[0], fluxvals)
        hdf.write(STD_identifier, 'zp_flux_err_' + bands[0], fluxerrs)

# Because I'm sick of dealing with a huge "utils" class, I'm writing a new class that draws from utils.
# This is exclusively for alignment of FITS files. Takes "fwhm" and "threshold" as initialization parameter, for DAOPHOT starfinder.
# Threshold is now deprecated (uses static values based on filter band)
# INCLUDES FITS ALIGNMENT PROCESSES
class fits_alignment(object):
    def __init__(self, fwhm):
        self.utilities = utils()
        self.filer = hdf5_writer(rootdir, 'data.hdf')
        self.fwhm = fwhm
        self.identifier = 'fitsfitsfitsaligned21'
    # Rotates a numpy array by 180 degrees
    def rotate_2d(self,array):
        shape = np.shape(array)
        n,m = shape
        retray = np.zeros(shape=shape)
        for i in range(n):
            for j in range(m):
                retray[i,j] = array[n - i - 1, m - j - 1]
        return np.int32(retray)

    # Utility explicitly for stacking aligned files in a given directory with name + identifier, here "'aligned21.fits'"
    def identifier_stacker(self):
        # Get files that are going to be stacked
        files = self.utilities.fits_identifier()
        new_files = []
        for file in files:
            filesplit = file.split(".")
            if filesplit[len(filesplit) - 2] == self.identifier:
                new_files.append(file)

        # Get an identifier to uniquely specify the stack
        directory = os.getcwd()
        dirsplit = directory.split("\\")
        length = len(dirsplit)
        cluster,offset,band = dirsplit[length - 4:length - 1]

        # Stack files new_files
        full_stack = np.zeros(shape=(4016,4016), dtype=np.int32)
        for file in new_files:
            data = fits.open(file)[0].data
            full_stack += data
        full_stack = full_stack/len(new_files)
        full_stack = np.int32(full_stack)

        # Grab an image header (the first one) for airmass reference
        header = fits.open(new_files[0])[0].header
        # Package Data. Overwrites the old FITS2, remember this.
        newfitso = fits.PrimaryHDU(full_stack)
        newfitso.header = header
        newfitsodata = fits.HDUList(newfitso)
        newfitsodata.writeto("fully_stacked" + "_" + cluster + "_" + offset + "_" + band + ".fits", overwrite=True)
    # Align image 2 to image 1, i.e. 2 will now be parallel to 1. Overwrites old fits..
    # fits1/fits2 same shape.
    def fits_aligner(self, fits1,fits2):
        # Print debug
        print("Processing " + fits1 + fits2)

        # Note the name for what the ideal "saved" file should be.
        filename = fits2 + "." + self.identifier + '.fits'
        # Check for the required identifier to run an alignment.
        fits2split = fits2.split(".fitstrimmed.")
        fits2_alignedidentifier = fits2split[1]
        # The identifier should be 'fitsfixed' if it goes in line with the data processing so far.
        if fits2_alignedidentifier == 'fitsfixed.fits':
            # Import the fits
            fitso1, fitso2 = fits.open(fits1), fits.open(fits2)
            # Get header of #2 and thus the band, too.
            header2 = fitso2[0].header
            band = header2['FILTER']
            threshold = "dolt"
            if band == 'U':
                threshold = 80
            elif band == 'B':
                threshold = 200
            elif band == 'V':
                threshold = 200
            else:
                print("Band fucked up mate.")

            # Extract data
            data1, data2 = fitso1[0].data, fitso2[0].data

            # Close fitso2
            fitso2.close()


        # Data is fine until here.


            data1_source,data2_source = data1,data2

            #plt.imshow(data2)
            #plt.show()

            # Record shape
            new_data_shape = np.shape(data2)
            datey, datex = new_data_shape

            # Adulterate the data such that no sources can be discovered within 10 pixels of the boundary of the image (a tolerable error).
            # Adulteration should use old shape.
            # The shape of the new array. (for new data down the line. We need the variable here, though.)
            # Both arrays should be same shape.

            # Check for if FITS2 is from the rotaligner. If it is, then we must trim by at least 700 from each side.
            clip = 10
            try:
                rot_value = header2['ROTER']
                if rot_value == 'TRUE':
                    print("File is rotated. Clipping for central stars and aligning.")
                    clip = 700
                else:
                    pass
            except:
                pass

            # Trim for source extraction within the bounds. x-bounds not equal to y-bounds.


            # First trim off the left/right
            x_low,x_high = np.arange(0,clip,1), np.arange(datex - clip, datex, 1)
            xdomain = np.append(x_low, x_high)
            for j in xdomain:
                for i in range(datey):
                    data1_source[i, j] = 0
                    data2_source[i, j] = 0

            # Then trim off the top/bottom
            y_low,y_high = np.arange(0,clip,1), np.arange(datey - clip, datey, 1)
            ydomain = np.append(y_low, y_high)
            for i in ydomain:
                for j in range(datex):
                    data1_source[i, j] = 0
                    data2_source[i, j] = 0


            # Get sources from both images using DAOPHOT (optional, just use IRAF starfind instead, anyway)
            src1,src2 = "null","null"
            try:
                starfinder = DAOStarFinder(threshold, fwhm=2*self.fwhm, sigma_radius=3, sharplo=0, roundlo=0, sharphi=999, roundhi=999, exclude_border=True)
                src1, src2 = starfinder(data1_source), starfinder(data2_source)
            except:
                print("Daophot starfinder does not have enough sources. Lower the threshold mate.")


            # src_xy will be a list of lists of vectors to the identified sources.
            sources = [src1, src2]
            src_xy, src_phots = [],[]

            # Extract x,y lists.
            for source in sources:
                fluxes = np.array(source['flux'])
                x, y = np.array(source['xcentroid']), np.array(source['ycentroid'])

                # Debug
                """
                if len(x) != len(y):
                    print("Length Error. See src_xy in fits_alignment. Stopping process for 20 seconds before error code.")
                    time.sleep(20)
                if len(x) != len(fluxes):
                    print("Length Error in src_xy vs. src_phots, as the flux list is not the same length as x,y list.")
                """

                vector_array = np.array([x, y]).T
                flux_array = fluxes

                if len(flux_array) != len(vector_array):
                    print("LENGTH ERROR IN THE SOURCE EXTRACTOR BABY! See fits_alignment")
                elif len(flux_array) == len(vector_array):
                    src_phots.append(flux_array)
                    src_xy.append(vector_array)
                else:
                    print("Fucked up.")
                    time.sleep(999)







            # src_xy1 is dif size to src_phots1.
            # Pair sources: smallest vector distance. We want vector (2 - 1), i.e. the vector that maps 1 onto 2.
            #pairwise_combinatorics = []
            pairwise_vectors = []
            list_pairwise = []
            magnitude = lambda x, y: np.sqrt(x ** 2 + y ** 2)


            # Get a translation to try alignment.
            for num, u in enumerate(src_xy[0]):

                # Get vec seps and their lengths
                separation_vectors = [d - u for d in src_xy[1]]
                magnitudes = []
                for d in separation_vectors:
                    dx,dy = d
                    magnitudes.append(magnitude(dx,dy))

                # Get min argument and attach this vector and other information.
                min_arg = np.argmin(magnitudes)
                #combinatoric = [num, min_arg]
                combivector = src_xy[1][min_arg] - src_xy[0][num]
                list_pairwise.append(combivector)
                if np.sqrt((np.linalg.norm(combivector))**2) <= 5:
                    try:
                        if np.abs(src_phots[0][num] - src_phots[1][min_arg]) <= (2*threshold):
                            pairwise_vectors.append(combivector)#, pairwise_combinatorics.append(combinatoric)
                        else:
                            #print("Failed vector test of length and separation.")
                            null = "null"
                    except IndexError:
                        print("WE DONE FUCKED UP NOW. See fits_aligner.") # Somehow src_xy[1] longer than src_phots[1].
                        print(len(src_xy[1]),len(separation_vectors), len(magnitudes), min_arg, len(src_phots[1]))
                        time.sleep(500)
                else:
                    #print("failed first test of length and separation")
                    null = "null"

            # Generic statistics list for the pairwise vectors for debug purposes.
            #sourcesxy = src_xy[0].T
            #pairwise = np.array(list_pairwise).T
            #x,y = sourcesxy
            #plt.clf()
            #plt.close()
            #plt.quiver(x,y,pairwise[0],pairwise[1])
            #plt.savefig("quiverplot" + fits2split[0] + ".png",dpi=300)
            #plt.show()

            #plt.clf()
            #plt.close()
            #lengths = [np.linalg.norm(d) for d in pairwise_vectors]
            #plt.hist(lengths, bins=100)
            #plt.savefig("histplot" + fits2split[0] + ".png",dpi=300)
            #plt.show()

            # Sigma-clip stats for the X,Y and Z components of the translation, getting the sigma-clipped mean/etc/STDEV. Sigmaclip = 2 sigma.
            xy = np.array(pairwise_vectors).T

            try:
                xm, xmed, xdev = sigma_clipped_stats(xy[0], sigma=1)
                ym, ymed, ydev = sigma_clipped_stats(xy[1], sigma=1)
            except Exception as e:
                print("Error has occurred in finding xm/ym/etc for fits_aligner. Sleeping.")
                print(e)
                print(xy)
                print(fits1, fits2)
                print("Error has occurred in finding xm/ym/etc for fits_aligner. Sleeping.")
                time.sleep(500)
            # Get the integer (rounded) values for xm,ym
            xm, ym = int(np.round(xm)), int(np.round(ym))






            # Next move the data back... we specify a "miniclip" to trim the data edge values slightly.
            # We don't expect a shift by more than ~5 pixels, so clipping by 5 keeps it nice and clean.
            # We also specify the "new_data_shape", a TUPLE!!!. Shape X should be Shape Y: whole tuple is a formality.

            # The amount by which to "clip" the data, i.e. to set the edge elements to a neutral value that can be ignored (the mean of the original data.)
            miniclip = 5
            # New shape was specified above before "adulteration" of the data.
            # Create new data array for aligned fits2.
            new_data2 = np.zeros(shape=new_data_shape, dtype=np.int32)
            mean, median, dev = sigma_clipped_stats(data2)

            # Restate data2.
            fitso2 = fits.open(fits2)
            oridata2 = fitso2[0].data

            #newfitso = fits.PrimaryHDU(oridata2)
            #newfitso.header = header2
            #newfitsodata = fits.HDUList(newfitso)
            #newfitsodata.writeto("testsaveddata.fits", overwrite=True)

            for i in range(0 + miniclip, datey - miniclip):
                for j in range(0 + miniclip, datex - miniclip):
                    try:
                        new_data2[i - ym, j - xm] = oridata2[i, j]
                    except:
                        print("We done fucked up now bro.")
                        continue

            #plt.imshow(new_data2)
            #plt.show()

            # Close original files
            fitso1.close(), fitso2.close()

            # Package Data. Overwrites the old FITS2, remember this.
            newfitso = fits.PrimaryHDU(new_data2)
            newfitso.header = header2
            newfitsodata = fits.HDUList(newfitso)
            newfitsodata.writeto(filename, overwrite=True)
        else:
            #print("Already done this one.", fits2)
            null = "null"
    # Runs fits_aligner for the provided directory. Uses zeroth image that matches parameters of fitsfixed as a file.
    def fits_dir_aligner(self, directory):
        # Navigate, get files, find out if the chosen reference image is suitable (i.e. an unaltered one.)
        os.chdir(directory)
        print(directory)
        files = self.utilities.fits_identifier()
        fits_choice_1 = files[0].split(".")
        fits_choice_1 = fits_choice_1[len(fits_choice_1) - 2]
        # Based on our choice of file arrangement, this should work. If not, come back and adapt to select an unregistered/unmoved image. ### BIT OF A RED AREA ###
        if fits_choice_1 == 'fitsfixed':
            for d in files:
                self.fits_aligner(files[0], d)
        else:
            print("See 'fits_dir_aligner' for error. Reference image is unsuitable.")
    # Gathers up the offsets accordingly.
    # NOTES ON A TRUE OFFSET ALIGNER BASED ON WCS TO MIRROR WREGISTER
    """
    Import two FITS files
    Get WCS for both of them (assuming they already have, if not then download/get it/etc)
    Work out the polar-coordinate rotation necessary on the ICRS to move from Fits 2 to Fits 1
    Apply that transformation to Fits 2
    Align those mother fucking offsets.
    """
    def offset_stackwcs(self, clusters):
        # Generates the new filestructure.
        bands = ['U', 'B', 'V']
        offsets = [str(d) for d in np.arange(0, 3, 1)]

        # Generate file structure.
        try:
            os.mkdir(rootdir + "\\Sci\\ClusterstacksDone")
        except:
            pass

        for cluster in clusters:
            try:
                os.mkdir(rootdir + "\\Sci\\ClusterstacksDone" + "\\" + cluster)
            except:
                continue
            for band in bands:
                try:
                    os.mkdir(rootdir + "\\Sci\\ClusterstacksDone" + "\\" + cluster + "\\" + band)
                except:
                    continue

        # Gather together all the various offset files.
        # Generate offset fixes
        for cluster in clusters:
            for band in bands:
                # Navigate to working dir.
                os.chdir(rootdir + "\\Sci\\ClusterstacksDone\\" + cluster + "\\" + band)
                # Specify the "original" dir
                oridir = rootdir + "\\Sci\\" + cluster
                # Specify CWD in case we need it.
                cwd = os.getcwd()

                # Generate list of files to copy over for this particular cluster and band...
                # Copy over the offset files we are going to stack in full directory format for this particular band.
                for offset in offsets:
                    # You'll want to alter this as the crow flies. This is all just provisional for the current file setup. Sorry not sorry for the fuckiness.
                    # aligned_0_CLUSTER_OFFSET_BAND.fits.fits is the base image.
                    # aligned_1_CLUSTER_OFFSET_BAND.fits.fits is the "aligned" image.
                    # Aligned images are for offsets [1:]/all non-initial elements of the offset list, which starts at offset 0.
                    file = ("{}_{}_{}_{}.fits").format("fully_stacked", cluster, offset, band)
                    shutil.copyfile(oridir + "\\" + offset + "\\" + band + "\\CALIBRATED\\" + file, cwd + "\\" + file)
                    # Run astrometry on the file that was just plonked in.
                    os.chdir(cwd)

                    # General debug stuff
                    directory = os.getcwd()
                    files = os.listdir()
                    print("Metrifying ", directory, " ", files)
                    self.utilities.astro_metrifier(file)
    # This will attempt to rotate one fits image to align with the other, using the WCS as a basis for this.
    # Assumes that the images are already aligned to be fairly close (i.e. using WREGISTER or otherwise)
    # Direct output from offset_aligner is doable but you'd need to twist the files to be the same size.
    # METHOD HERE!!!
    # The newly rotated file is called "fits2.rotated.fits" where fits2 is the full name of the old file.
    """
    This assumes that the WCS has already been generated. 
    Use two coordinate points somewhat near the centre from the 2nd on the first
    Get the vectors corresponding to them
    Work out the rotation used on them 
    Rotate 2 back onto 1
    Save the files in a rot-fixed subdirectory.
    """
    def fit_rot_aligner(self, fits1, fits2):
        # Load files, data, and headers.
        fitso1, fitso2 = fits.open(fits1), fits.open(fits2)
        data1, data2 = fitso1[0].data, fitso2[0].data
        heado1, heado2 = fitso1[0].header, fitso2[0].header
        # Get wcs
        wcs1, wcs2 = wcs.WCS(heado1), wcs.WCS(heado2)
        # Get dimensions. Assumed the same shape.
        shape = np.shape(data1)
        # Pick midpoint (parallel to x) and a point to the right by ~100pix.
        pytho_shape = [d - 1 for d in shape] # Pythonic shape (4016 ~ 4015, index wise.) *ROWS:COLUMNS*
        pytho_mid = np.array([int((d - 1)/2) for d in pytho_shape]) # Pythonic centre (y_mid, x_mid)
        pytho_midl = np.array([pytho_mid[1], pytho_mid[0]]) # Vectorial centre (x_mid, y_mid)
        right = pytho_midl + np.array([1000,0]) # Right after pixel coords to the left and right by 100 respectively.

        # Get WCS coords for centre,right in first frame (radec)
        centre, right_1 = wcs1.wcs_pix2world(pytho_midl[0], pytho_midl[1], 0), wcs1.wcs_pix2world(right[0], right[1], 0)
        # relative_vector_11 = right - pytho_midl
        # relative_vector_1 = relative_vector_11/np.linalg.norm(relative_vector_11) Normalized

        # Get the pixel coordinates in the 2nd frame for first
        centre2pix, right2pix = wcs2.wcs_world2pix(centre[0], centre[1], 0), wcs2.wcs_world2pix(right_1[0], right_1[1], 0)
        centre2pix, right2pix = np.array([float(d) for d in centre2pix]), np.array([float(d) for d in right2pix])

        relative_vector_22 = right2pix - centre2pix  # Vector from left to right (hopefully)
        # relative_vector_2 = relative_vector_22/np.linalg.norm(relative_vector_22) # Normalized


        # Just some notation to make things easier to understand.
        x_1,y_1 = relative_vector_22[0], relative_vector_22[1]

        ## Note, here theta is positive counter clockwise, negative clockwise. We need to account for negative theta.
        #sintheta = -1*(y_0*x_1 - y_1*x_0)/(y_0**2 + x_0**2) #Negative the theta! We're doing the active transform, not passive.
        #theta = np.arcsin(sintheta)
        # ,x_1, ,y_1
        """
        Need a smarter way to deduce/solve for theta, which should be between 0 and 360 degrees (or higher. Just, the full 360.)
        Theta is the (positive) counter clockwise angle from relvect 1 onto relvect 2. 
        There are four areas we have: A S T C.
        If x>=0 y>0 then theta is just arctan(y/x) (it'll be 0>pi/2)
        If x>=0 y<0 then theta is again just arctan(y/x) (it'll be 0>-pi/2)
        If x<0 y>0 then theta is going to be pi/2 + arctan(y/(-x)) 
        If x<0 y<0 then theta is going to be pi + arctan(y/x) 
        """
        theta = "null"
        if y_1 >= 0:
            if x_1 >= 0:
                theta = np.arctan(y_1/x_1)
            if x_1 < 0:
                theta = np.pi - np.arctan(y_1/(-1*x_1))
        if y_1 < 0:
            if x_1 >= 0:
                theta = np.arctan(y_1/x_1)
            if x_1 < 0:
                theta = np.pi + np.arctan(y_1/x_1)
        thetadeg = np.rad2deg(theta)
        thetarot = -thetadeg
        thetarot = -thetarot # Try the negiflip.

        # Save thetarot to a temp file
        with open(fits2 + "rotated.txt", 'w') as f:
            f.write(("Rotation is {}").format(thetarot))

        # Rotate
        data2_new = transform.rotate(data2, thetarot, preserve_range=True, resize=False, )

        # DEPRECATED CODE (n^2 manual rotation of array by angle theta).
        # data2_new is the new array name you should use for alternative methods.
        """
        # Sine and cos for future use.
        siner = np.sin(theta)
        coser = np.cos(theta)

        # Next create a new empty array with the same shape as original.
        data2_new = np.zeros(shape=np.shape(data2))

        # Mapping guesses...

        # Generate the relative vector that points from our chosen reference point to the given pixel that we are mapping
        # Rotate the vector
        # Map the old coordinate to the new one in the array data2_new
        for i in range(shape[0]):
            for j in range(shape[1]):
                # We want to translate the array.
                blcorneroriginvector = np.array([j, shape[0] - 1 - i]) # Pythonically speaking.
                relvector = blcorneroriginvector - pytho_midl

                # Map it (rotate it!)
                # new_x = x_1*cos(theta) - y_1*sin(theta)
                # new_y = x_1*sin(theta) + y_1*cos(theta)

                # The new (derotated) relative vector from the midpoint.
                new_relvector = np.array([relvector[0]*np.cos(theta) + relvector[1]*np.sin(theta), (-1)*relvector[0]*np.sin(theta) + relvector[1]*np.cos(theta)])
                # Convert back to the pythonic relative to the origin
                new_blcorneroriginvector = new_relvector + pytho_midl

                # Convert back to i,j notation
                new_index = np.array([shape[0] - new_blcorneroriginvector[1] - 1, new_blcorneroriginvector[0]])
                new_index = np.array([np.round(d) for d in new_index]) # This is the derotated index . [i,j] is the original index.
                new_index = [d + 0.1 for d in new_index]
                new_index = np.array([int(d) for d in new_index])

                # Map values. Take original [i,j] and dump it at the new [i,j]
                try:
                    data2_new[new_index[0],new_index[1]] = data2[i,j]
                except:
                    continue
        """

        # Note the new filename
        rotated_name = fits2 + ".rotated.fits"


        # For debug purposes. Write down teh newly rotated array.
        # Write file and all that stuff.
        data2new = fits.PrimaryHDU(np.int32(data2_new))
        data2new.header = heado2
        data2new.header['OWNERID'] = "NOT_FITTED"
        data2new = fits.HDUList([data2new])
        data2new.writeto(rotated_name, overwrite='true')

        # Save CWD
        cwd = os.getcwd()

        # Run the astrometrifier for the file
        self.utilities.astro_metrifier(rotated_name)

        # Return to CWD (in case we have lost our original place.)

        os.chdir(cwd)

        # Magically bring back the rotated fits2 file.
        os.replace(cwd + "\\WCSDONE\\" + rotated_name + "wcs.fits", cwd + "\\" + rotated_name + "wcs.fits")

        # Attempt to align the now-rotated fits file to the original fits file.
        # Load files, data, and headers.
        rotated_name = fits2 + ".rotated.fits"
        fitso1, fitso2 = fits.open(fits1), fits.open(rotated_name + "wcs.fits")
        data1, data2 = fitso1[0].data, fitso2[0].data
        heado1, heado2 = fitso1[0].header, fitso2[0].header
        # Get wcs
        wcs1, wcs2 = wcs.WCS(heado1), wcs.WCS(heado2)
        # Get dimensions. Assumed the same shape.
        shape = np.shape(data1)
        # Pick midpoint (parallel to x)
        pytho_shape = [d - 1 for d in shape]  # Pythonic shape (4016 ~ 4015, index wise.) *ROWS:COLUMNS*
        pytho_mid = np.array([int((d - 1) / 2) for d in pytho_shape])  # Pythonic centre (y_mid, x_mid)
        pytho_midl = np.array([pytho_mid[1], pytho_mid[0]])  # Vectorial centre (x_mid, y_mid)

        # Get WCS coords for centre in first frame (radec)
        centre1wcs = wcs1.wcs_pix2world(pytho_midl[0], pytho_midl[1], 0)

        # Get the pixel coordinates in the 2nd frame
        centre2pix = wcs2.wcs_world2pix(centre1wcs[0], centre1wcs[1], 0)

        # Get the translation required to move 1 onto 2 (x,y)
        relative_vector_21 = np.array([int(d) for d in centre2pix]) - np.array([int(d) for d in pytho_midl])

        with open(fits2 + "translated.txt", 'w') as f:
            f.write(("Translation is {}").format(relative_vector_21))

        # Shift using scipy.ndimage.shift
        data2_new = shift(data2, shift=np.array([-1*relative_vector_21[1], -1*relative_vector_21[0]]), cval=0.0)

        """
        # Deprecated code to manually translate. Replaced with skimage translate.  data2_new is new array.

        data2_new = np.zeros(shape=np.shape(data2))

        # Translate the fits2 by -1*relative_vector_21 (to align it.)
        # A few pixels worth of error is acceptable.
        for i in range(shape[0]):
            for j in range(shape[1]):
                # i is the y coordinate, j is the x coordinate.
                # Let X and Y be the translations needed to take 2 onto 1
                # new_i will be i - y
                # new_j will be j + x
                try:
                    data2_new[i, j] = data2[i + relative_vector_21[1], j + relative_vector_21[0]] # Remember that we're translating by NEGATIVE the relative_vector_21, hence the fact that we have opposite signs.
                except:
                    continue

        """

        # Write file and all that stuff.
        rot_aligned_name = fits2 + "rot_trans_aligned.fits"
        data2new = fits.PrimaryHDU(np.int32(data2_new))
        data2new.header = heado1 # USE THE HEADER FROM THE FIRST FILE! This is for the sake of having equal WCS/Aladin ease.
        data2new.header['OWNERID'] = "NOPE"
        data2new.header['ROTER'] = 'TRUE'
        data2new = fits.HDUList([data2new])
        data2new.writeto(rot_aligned_name, overwrite='true')

    # Runs fit_rot_aligner for all the FITS in the given directory.
    def fit_rot_aligner_all(self, dir_file_concats):
        directory, files = dir_file_concats[0], dir_file_concats[1]
        os.chdir(directory)
        cwd = os.getcwd()
        print(files)
        for file in files:
            # Run aligner
            self.fit_rot_aligner(files[0],file)
            print("Aligning ",file)
            # New name for rotated/aligned file
            new_file = file + "rot_trans_aligned.fits"
            # Rename new_file
            os.chdir(cwd)
            even_newer = new_file + ".fitstrimmed.fitsfixed.fits" # fully_stacked_NGC7789_2_U.fitswcs.fitsrot_trans_aligned.fits.fitstrimmed.fitsfixed EXAMPLE NAME
            try:
                os.replace(new_file, even_newer)
            except:
                print("Error has occurred in file replacement. See fits_rot_aligner_all")
                continue
            # Attempt alignment XY
            self.fits_aligner(files[0], even_newer)


            # New file name for this.
            new_name = even_newer + "." + self.identifier + '.fits' # fully_stacked_NGC7789_0_U.fitswcs.fitsrot_trans_aligned.fits.fitstrimmed.fitsfixed.fits.fitsfitsfitsaligned21 EXAMPLE NAME
            os.chdir(cwd)
            try:
                os.mkdir(cwd + "\\FinalAligns")
                os.replace(cwd + "\\" + new_name, cwd + "\\FinalAligns\\" + new_name)
            except:
                os.replace(cwd + "\\" + new_name, cwd + "\\FinalAligns\\" + new_name)

        # Once all aligned, navigate to the final subdir and do a quick preliminary average stack.
        # An identifier is achieved based on the file directory.
        dirsplit = directory.split("\\")
        idi1,idi2 = dirsplit[len(dirsplit) - (1 + 1)],dirsplit[len(dirsplit) - (1 + 2)] # BAND THEN CLUSTER IN CURRENT SETUP
        stack_identifier = ("fully_stacked_{}_{}").format(idi2,idi1)
        os.chdir(directory + "\\FinalAligns")
        self.utilities.average_stacker(stack_identifier)


    # Old Fits Offset Aligner (deprecated)
    # Takes two fits, attempts WCS alignment, accounts for rotation, gets them within spurting distance, trims down the axes instead of keeping fixed axis size.
    def fits_offset_aligner(self, fits1, fits2):
        # New filename c:
        fitsos = [fits1, fits2]

        # Get data
        data = []
        for d in fitsos:
            with fits.open(d) as f:
                data.append(f[0].data)

        # Get original headers
        heados = []
        for u in fitsos:
            with fits.open(u) as f:
                heados.append(f[0].header)

        # Get threshold
        band = heados[0]['FILTER']
        threshold = "dolt"
        if band == 'U':
            threshold = 150
        elif band == 'B':
            threshold = 800
        elif band == 'V':
            threshold = 1000
        else:
            print("Band fucked up mate.")

        # Specify CWD
        cwd = os.getcwd()

        # Note the true-false noter
        not_succeeded = True
        infinitoken = 0

        # Try-except loop (A BIG ONE!) to account for the possibility of field rotation. Will attempt to correct using 180-degree rotator! :)
        while not_succeeded:
            try:
                # Currently hashed out for debugging.
                # Get WCS header and repackage data and headers.

                src_1, src_2 = self.utilities.source_extractor(data[0], threshold), self.utilities.source_extractor(
                    data[1], threshold)

                # HASH OUT HERE FOR DEBUG

                wcs1, wcs2 = self.utilities.wcsgetter(src_1), self.utilities.wcsgetter(src_2)
                os.chdir(cwd)

                """
                print("RUNNING SEX")
                time.sleep(5)

                # Write WCS to file for time saving.
                newfitso = fits.PrimaryHDU(np.int32(data[0]))
                newfitso.header = wcs1
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("wcs1test.fits", overwrite=True)

                newfitso = fits.PrimaryHDU(np.int32(data[1]))
                newfitso.header = wcs2
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("wcs2test.fits", overwrite=True)

                # HASH IN HERE TO SKIP WRITING
                """
                """
                # Fileread wcs/etc for debug. Clean this up later.
                wcss1,wcss2 = "null","null"
                with fits.open("wcs1test.fits") as f:
                    wcss1 = f[0].header
                with fits.open("wcs2test.fits") as f:
                    wcss2 = f[0].header


                wcs1,wcs2 = wcss1, wcss2"""

                # HASH IN HERE FOR NO DEBUG

                """
                # Repackage the WCS with original FITS (for debug purposes. You can mop this up/not do this. Just do it if you feel the WCS fit needs to be checked.)
                wcss = [wcs1,wcs2]
                for num, dataxi in enumerate(data):
                    newfitso = fits.PrimaryHDU(dataxi)
                    newfitso.header = wcss[num]
                    newfitsodata = fits.HDUList(newfitso)
                    # New WCS fit file (for debug purposes)
                    newfitsodata.writeto("temp_" + str(num) + ".fitsbbtrimmedWCSfit." + "fitsfixed.fits", overwrite=True)"""

                # Get the WCS for further processing.
                wcs11, wcs22 = wcs.WCS(wcs1), wcs.WCS(wcs2)

                # We need to work out the optimal shift to preserve the most data, if that makes sense.
                """
                Get WCS of centre of each field (approx)
                Get vector that goes from centre 1 to centre 2
                divide by two to get a vector that goes from 1 to the midway of 1-2
                Translate matrices to match at midpoint (to attempt a close overlap)
                Once at the overlap point, attempt to run standard aligner (assuming sub-threshold overlap, currently at 5 pixels as of writing.
                """

                # The midpoint (in pixels) is given by this argument. 4016 - 1 = array argument length.
                # We assume fits are in the range of having a width of 4016 pixels. Whatever, fix this if it gets to be a problem.
                midx, midy = int((4015 - 1) / 2), int(
                    (4015 - 1) / 2)  # Long, but incase we add an argument to represent the clipped length.
                midpoint = np.array([midx, midy])
                # Get [RA,DEC] in degrees for FITS2 centre
                radec2 = wcs22.wcs_pix2world(midpoint[0], midpoint[1], 0)
                radec2 = np.array([float(d) for d in radec2])
                # Find [PIX_X,PIX_Y] centre of 2 in image of 1
                radec2on1 = wcs11.wcs_world2pix(radec2[0], radec2[1], 0)
                radec2on1 = np.array([float(d) for d in radec2on1])
                # Get vector that points from coord in 1 to coord in 2
                translation = midpoint - radec2on1

                # The method will be to clip from the second offset and the first offset to their overlap area. The X,Y components of this clip are:
                x, y = int(np.round(translation[0])), int(np.round(translation[1]))

                # Both x and y may be positive or negative. In any case...
                # Clip arrays.
                shape = np.shape(data[0])
                shape = [d for d in shape]
                shapex, shapey = shape[1], shape[0]

                # For instantaneous debug.
                print(translation)
                print(type(translation[0]))
                print(x, y)
                print(type(x), type(y))
                print(shape)
                print(type(shapex), type(shapey))

                # x,y can be positive or negative. Account for this. We're going it using NP h and v splits.
                # Notes on numpy slicing.
                """
                Here is example rotation code.

                # Attempt a rotation
                print(np.shape(data))
                rotata = fitso.rotate_2d(data)
                rotata = rotata[500:,:]
                rotata = fitso.rotate_2d(rotata)

                The above trims off 500 from the top of data. 

                If you don't want to trim off 500 from the top but instead want to trim off 500 from the bottom... 
                rotata = rotata[500:,:] 

                Why? Fucking hell if I know. That 500 shouldn't even exist in the end, but hey, whatever, fucking fuck.
                X-Slicing is justttt fine. <3 

                IN THE CASE OF X BEING POSITIVE:
                Trim off the back of 2 and front of 1 (x_wise)
                X BEING NEGATIVE:
                Trim off the back of 1 and front of 2 (x_wise)

                IN THE CASE OF Y BEING POSITIVE:
                Trim off the top of 1 and bot of 2 (y_wise)
                Y BEING NEGATIVE:
                Trim off the top of 2 and bot of 1 (y_wise) 
                """

                # SLICING FOR X. This works as normal and isn't a problem.
                data1, data2 = "null", "null"  # debug identify.
                if x < 0:  # x<0
                    data1, data2 = data[0][:, -x:], data[1][:, :x + shapex]
                elif x >= 0:
                    data1, data2 = data[0][:, :shapex - x], data[1][:, x:]
                else:
                    print("posneg error in offset alignment.")
                newdat1, newdat2 = "null", "null"
                data1, data2 = np.array(data1), np.array(data2)

                """
                # Save files for various debug purposes (in the case of a crash occurring)
                newfitso = fits.PrimaryHDU(data1)
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("data1debug.fits", overwrite=True)

                newfitso = fits.PrimaryHDU(data2)
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("data2debug.fits", overwrite=True)  # DEBUG"""

                # Go for Y. We're going to have to use rotations here to handle this shit.
                if y < 0:  # y<0
                    # Trim off the top of 2 and bottom of 1
                    # FIRST DO 1
                    newdat1 = data1[-y:, :]
                    # THEN DO 2
                    data2 = self.rotate_2d(data2)
                    data2 = data2[-y:, :]
                    newdat2 = self.rotate_2d(data2)
                elif y >= 0:
                    # Trim off the top from 1, bottom from 2.
                    # FIRST DO 2
                    newdat2 = data2[y:, :]  # Trim Y from the bottom of data2
                    # THEN DO 1
                    data1 = self.rotate_2d(data1)  # Rotate data
                    data1 = data1[y:, :]  # Trim from the bottom (hence the top when you rotate it back)
                    newdat1 = self.rotate_2d(data1)  # Rotate back
                else:
                    print("posneg error in offset alignment.")

                # Bit of an artefact from debugging.
                data1, data2 = np.array(newdat1), np.array(newdat2)

                # Debug Crap
                # SOMETHING IS REAL GOD DAMN FUCKY WHAT THE HELL IS WRONG WITH MY SLICING DAMN IT ALL TO GOD DAMN HELL FUCK YOU NUMPY
                newfitso = fits.PrimaryHDU(newdat1)
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("data1debug2.fits", overwrite=True)

                newfitso = fits.PrimaryHDU(newdat2)
                newfitsodata = fits.HDUList(newfitso)
                # New WCS fit file (for debug purposes)
                newfitsodata.writeto("data2debug2.fits", overwrite=True)  # DEBUG

                # Deprecated code that is not suitable for this. Keep for future method reference. This is for a "half clip" (naivity at its finest.)
                """
                transvector = 0.5*translation
                transvector_int, = np.array([round(d) for d in transvector])
                x,y = transvector_int[0], transvector_int[1] 

                x,y are the translation vectors.
                we're moving array 1 by +x and +y, 2 by -x and -y 
                thus the relevant data for array 1 will be array_1[0:4016-y, x:4016]
                for array two the relevant data will be array_2[y:4016, 0:4016-x] 
                # Let's try it. 

                # Clip arrays.
                shape = np.shape(data[0])
                shape = [d for d in shape]
                shapex,shapey = shape[1],shape[0]

                # x,y can be positive or negative. Account for this.
                #---------------------------------------------------
                if x >= 0:
                    data1,data2 = data[0][0:shapey,x:shapex], data[1][0:shapey, 0:shapex - x]
                else: # x<0
                    data1,data2 = data[0][0:shapey,0:shapex + x], data[1][0:shapey, -x:shapex]
                # retrieve new shapes
                shape1,shape2 = np.shape(data1),np.shape(data2)
                oldx1,oldx2 = shape1[1],shape2[1]
                # go for y
                if y >= 0:
                    data1,data2 = data1[0:shapey - y, 0:oldx1], data2[y:shapey, 0:oldx2]
                else: # y<0
                    data1,data2 = data1[-y:shapey, 0:oldx1], data2[0:shapey + y, 0:oldx2]
                """

                # For list comprehension
                datas = [data1, data2]

                # Debug.
                """
                print(np.shape(data1), np.shape(data2), np.shape(data[0]), np.shape(data[1]))
                for datax in datas:
                print(np.shape(datax))"""

                # Save CWD for purposes
                os.chdir(cwd)
                for num, daton in enumerate(datas):
                    newfitso = fits.PrimaryHDU(daton)
                    newfitso.header = heados[num]
                    newfitsodata = fits.HDUList(newfitso)
                    # Tempfile is temp_1,2.fits"
                    newfitsodata.writeto("temp_" + str(num) + ".fitstrimmed." + "fitsfixed.fits", overwrite=True)

                # tempfiles
                os.chdir(cwd)
                tempfiles = ['temp_0.fitstrimmed.fitsfixed.fits', 'temp_1.fitstrimmed.fitsfixed.fits']
                self.fits_aligner(tempfiles[0], tempfiles[1])

                # Remove temps and rename original file for uniqueness
                cwd = os.getcwd()
                # Get filename
                fits2splittifier = fits2.split("fully_stacked_")[1]

                os.replace(tempfiles[0], "aligned_0_" + fits2splittifier + ".fits")
                os.remove(tempfiles[1])
                os.replace(tempfiles[1] + "." + self.identifier + '.fits', "aligned_1_" + fits2splittifier + ".fits")
                # All done.
                not_succeeded = False
            except Exception as e:
                print(e)
                print("Index error: field may be rotated 180 degrees. Attempting to rotate 180 and then resolve.")
                print(
                    "Note, if this takes too long, an infinite loop has occurred and no solutions have been found. Go back and recode to avoid this.")
                if infinitoken >= 1:
                    print("Reached max limit of rotation attempts.")
                    not_succeeded = False
                data[1] = self.rotate_2d(data[1])
                infinitoken += 1
    # This assumes you've tried to correct rotation using the the IRAF WREGISTER.
    # Same file structure as offset_stackwcs output files except under "ClusterstacksDoneIRAF"
    # Not due for use since IRAF will no longer be utilized.
    def offset_wcsfinish(self, clusters):
            # Fitser
            fitser = fits_alignment(3)
            utilities = utils()
            # Generates the new filestructure.
            bands = ['U', 'B', 'V']
            offsets = [str(d) for d in np.arange(0, 3, 1)]

            # Align all images to eachother and then stack.
            for cluster in clusters:
                for band in bands:
                    # Navigate to working dir.
                    os.chdir(rootdir + "\\Sci\\ClusterstacksDoneIRAF\\" + cluster + "\\" + band + "\\WCSDONE")

                    # Specify CWD in case we need it.
                    cwd = os.getcwd()

                    # Get the file list. This is a bit of a point that will need to be returned to.
                    # The "try-except" loop is for renaming the files to match the format fits_align requires.
                    file_list = []
                    for offset in offsets:
                        try:

                            # This is the file identifier we currently have. Totally random, just for the test data.
                            file = ("{}_{}_{}_{}.fitswcs.fits").format("fully_stacked", cluster, offset, band)
                            # Attempt a rename to the format that fits_aligner will accept.

                            """
                            # Note the name for what the ideal "saved" file should be.
                            filename = fits2 + "." + self.identifier + '.fits'

                            # Check for the required identifier to run an alignment.
                            fits2split = fits2.split(".fitstrimmed.")
                            fits2_alignedidentifier = fits2split[1]

                            # The identifier should be 'fitsfixed' if it goes in line with the data processing so far.
                            if fits2_alignedidentifier == 'fitsfixed.fits':
                            """

                            # Note a new filename
                            new_file = "stack_" + offset + ".fitstrimmed." + "fitsfixed.fits"
                            file_list.append(new_file)
                            # Try to rename.
                            os.rename(file, new_file)
                        except:
                            # The file already has the required name, so, we can just go ahead and append the name to the list.
                            new_file = "stack_" + offset + ".fitstrimmed." + "fitsfixed.fits"
                            file_list.append(new_file)

                    # This is a bit of a quick fix edit we're doing to make sure all the FITS have edge artefacts trimmed off.
                    for file in file_list:
                        data, header = 0, 0
                        f = fits.open(file)
                        data = f[0].data
                        header = f[0].header
                        x, y = np.shape(data)
                        data = data[200:y - 200, 200:x - 200]
                        newfitso = fits.PrimaryHDU(data)
                        newfitso.header = header
                        newfitsodata = fits.HDUList(newfitso)
                        newfitsodata.writeto(file + ".holder", overwrite=True)
                    print("done")

                    try:
                        os.mkdir("Holders")
                    except:
                        pass

                    cwd = os.getcwd()
                    for file in file_list:
                        filer = file + ".holder"
                        try:
                            os.replace(cwd + "\\" + filer, cwd + "\\Holders\\" + file)
                        except:
                            continue

                    # Align fits.
                    # for file in file_list:
                    #    fitser.fits_aligner(file_list[0], file)

                    # Generate the new fits file names. This is just based on what we already have, sorry.
                    file_fits_list = [d + ".fitsfitsfitsaligned21.fits" for d in file_list]
                    # We're going to move to a subdirectory so that we can make use of the average_stacker utility we already have.
                    cwd = os.getcwd()
                    try:
                        os.mkdir("ALIGNED")
                    except:
                        pass
                    for file in file_fits_list:
                        try:
                            os.replace(cwd + "\\" + file, cwd + "\\ALIGNED\\" + file)
                        except:
                            continue
                    utilities.average_stacker("ALIGNED")
    # Recursively runs fits_offset_aligner for all individual cluster bands under our file structure. Offset N0 and Band Identifiers are hard coded as usual.
    def offset_recursifier(self, clusters):
        # Example: CLUSTER + OFFSET + BAND + CALIBRATED + (FILES)
        # Generate endpoints, manually this time.
        # Also hard coded is the identifier: "fully_stacked" is the final stacked image for this offset directory.
        bands = ['U', 'B', 'V']
        offsets = [str(d) for d in np.arange(0, 3, 1)]

        # Generate file structure.
        try:
            os.mkdir(rootdir + "\\Sci\\Clusterstacks")
        except:
            pass

        for cluster in clusters:
            try:
                os.mkdir(rootdir + "\\Sci\\Clusterstacks" + "\\" + cluster)
            except:
                continue
            for band in bands:
                try:
                    os.mkdir(rootdir + "\\Sci\\Clusterstacks" + "\\" + cluster + "\\" + band)
                except:
                    continue

        # Generate offset fixes
        for cluster in clusters:
            for band in bands:
                # Navigate to working dir.
                os.chdir(rootdir + "\\Sci\\Clusterstacks\\" + cluster + "\\" + band)
                cwd = os.getcwd()

                # Generate the offset files we are going to stack in full directory format for this particular band.
                file_lister = []
                for num, offset in enumerate(*offsets):
                    file_lister.append(
                        ("{}\\Sci\\{}\\{}\\{}\\CALIBRATED\\{}" + ".fits").format(rootdir, cluster, offset, band,
                                                                                 ("fully_stacked_{}_{}_{}").format(
                                                                                     cluster, offset, band)))
                    shutil.copyfile(file_lister[num],
                                    cwd + "\\" + ("fully_stacked_{}_{}_{}.fits").format(cluster, offset, band))

                # Files are of the form "fully_stacked_M52_0_V" in the current directory.
                # Use Offet #0 as the reference. Offset #1,2...etc against that.
                print(file_lister)
                for num, file in enumerate(file_lister):
                    # Check for alignment already...
                    file_list = self.utilities.fits_identifier()
                    splitter = file.split("fully_stacked_")[1]
                    split_concat = ("{0}{1}.fits").format("aligned_1_", splitter)
                    tokeno = 0

                    if file == file_lister[0]:
                        print("Passing " + file)
                        pass
                    else:
                        for filo in file_list:
                            if filo == split_concat:
                                tokeno += 1
                        if tokeno != 0:
                            print("Passing")
                            pass
                        else:
                            print("Solving " + file_lister[num])
                            print(file)
                            os.chdir(cwd)
                            self.fits_offset_aligner(file_lister[0], file_lister[num])


# h5py/hdf5 file-writing utilities for the 'data.hdf5' file in the specified directory.
# Directory is explicit (i.e. D:\TGP2020\etc). Only singular arrays to each.
# INCLUDES HDF5 FILE HANDLING UTILITIES
class hdf5_writer(object):
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename
    # Creates the file if it doesn't exist, does nothing otherwise.
    def create(self):
        os.chdir(self.directory)
        file_created = h5py.File(self.filename, 'a')
        file_created.close()
        os.chdir(self.directory)
    # Returns the list of groups for the file.
    def info(self):
        os.chdir(self.directory)
        f = nx.nxload(self.filename)
        print(f.tree)
    # Writes ARRAY to the DATASET inside GROUP of FILENAME. Can also write SINGLE STRINGS.
    def write(self, group, dataset, array):
        os.chdir(self.directory)
        # First try to create group and add to dataset.
        # Except: The group exists. Access the group instead.
        # Further Except: The group exists, AND the dataset exists. Delete old dataset and retry.

        # In each case, consider whether the data is in the form of a single string or not.
        try:
            with h5py.File(self.filename, 'a') as f:
                grouper = f.create_group(group, 'a')
                if type(array) == str:
                    stringer = array.encode("ascii", "ignore")
                    datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                else:
                    datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)
        except:
            try:
                with h5py.File(self.filename, 'a') as f:
                    grouper = f[group]
                    if type(array) == str:
                        stringer = array.encode("ascii", "ignore")
                        datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                    else:
                        datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)
            except:
                with h5py.File(self.filename, 'a') as f:
                    del f[group][dataset]
                    grouper = f[group]
                    if type(array) == str:
                        stringer = array.encode("ascii", "ignore")
                        datasetter = grouper.create_dataset(dataset, shape=(1, 1), dtype='S10', data=stringer)
                    else:
                        datasetter = grouper.create_dataset(dataset, shape=np.shape(array), data=array)


    # Write an astropy table
    def write_table(self, group, dataset, astropytable):
        os.chdir(self.directory)
        astropy.io.misc.hdf5.write_table_hdf5(astropytable, self.filename, path=(group + "/" + dataset), append=True,overwrite=True, serialize_meta=True)

    # Read an astropy table
    def read_table(self, group, tablename):
        os.chdir(self.directory)
        tab = astropy.io.misc.hdf5.read_table_hdf5(input=self.filename, path=(group + "/" + tablename), character_as_bytes=False)
        return tab

    # Reads ARRAY from DATASET inside GROUP of FILENAME
    def read(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a np array.
            data = np.array(f[group][dataset])
        return data

    # Read a string instead, singular, not tested on arrays.
    def read_string(self, group, dataset):
        os.chdir(self.directory)
        data = "null"
        with h5py.File(self.filename, 'r') as f:
            # As a string (well, bytes.)
            data = str(np.array(f[group][dataset]), "ascii")
        return data

    # Converts a table to votable xml format + exports. Exportname automatically gets .xml export key.
    def write_votable(self, group, dataset, exportname):
        table = self.read_table(group, dataset)
        votable = from_table(table)
        writeto(votable,exportname + ".xml")

# Class to handle miscellaneous file import, i.e. (now deprecated, use query method) GAIA votables, and the Hipparcos votable.
# Also has a method to read in the Pleiades WEBDA catalogue CSV file (which has been processed for this purpose).
class misc_file_handler(object):
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename
        self.filters = ['g', 'bp', 'rp'] # GAIA filter keys

    # Import Shazzles text files (# R;A[0]; ;D;E;C[1]; ;V[2]; ;V;e;r;r[3]; ;U;B[4]; ;U;B;e;r;r[5]; ;B;V[6]; ;B;V;e;r;r[7] is the order/format)
    def syazza_porter(self, directory, filename, group, set):
        file = open(directory + "\\" + filename)
        arr = []
        for num,i in enumerate(file):
            if num == 0: # skip first line which has info on contents
                pass
            else:
                split = i.split(";") # ; is delimiter
                split = [float(d) for d in split] # Make float
                arr.append(split)
        arr = np.array(arr)
        table = Table()
        table['ra'], table['dec'], table['V'], table['V_err'], table['UB'], table['UB_err'], table['BV'], table['BV_err'], table['par'], table['par_err'] = arr.T
        filer = hdf5_writer(self.directory, self.filename)
        filer.write_table(group, set, table)





        file.close()

    # Import any old generic VOTABLE file.
    # ColFlags is optional (if it throws an error, ColFlags lets you specifically select columns to use)
    def votable_import(self, directory, filename, group, set, ColFlags):
        # Set up filer
        filer = hdf5_writer(self.directory, self.filename)
        # Get votable, also in the rootdir/etc.
        votable = parse_single_table(directory + "\\" + filename).to_table(use_names_over_ids=True)

        # Write the table
        filer.write_table(group, set, votable)

    # Makes use of https://github.com/dsavransky/MeanStars to obtain a typical colour vs. M/Spectral_Type/etc table (+ saves it to file)
    def meanstar_dwarf(self):
        MSdata = MeanStars().data
        filer = hdf5_writer(self.directory, self.filename)
        filer.write_table("MEANSTARS", "spec_table", MSdata)

    # Reads in the Fitzgerald CSV, ; delimiter
    # Sp. Type	U-B	(B-V)	(V-R)C	(V-I)C	(V-J)	(V-H)	(V-K)	(V-L)	(V-M)	(V-N)	(V-R)J	(V-I)J
    def fitz_read(self, csvdir, csvname):
        # Change to dir, get file, import with numpy, get a nice table of SP, UB, BV
        filer = hdf5_writer(self.directory, self.filename)
        os.chdir(csvdir)
        csv = open(csvname)
        table = Table()
        SP,UB,BV = [],[],[]
        for i in csv:
            i = i.split(";")
            sp, ub, bv = i[0], float(i[1]), float(i[2])
            SP.append(sp), UB.append(ub), BV.append(bv)
        csv.close()
        table['SP'], table['UB'], table['BV'] = SP, UB, BV
        # Map a numerical scale to the spectral type. It's roughly in order from HOT to COLD (B0 to M4).
        # As long as we retain the mapping throughout we should be fine (since we're mapping the observed BV on the graph to SP Type, subsequently to the MS age and other parameters.
        index = np.arange(0, len(SP), 1)
        table['ID'] = index

        # BAFGKM: B0.0 is example.
        spec_indices = ["B", "A", "F", "G", "K", "M"]

        # Save
        filer.write_table("FITZGERALD", "spec_table", table)

    # Reads in the GLIESE text catalogue and writes the U,B,V to file. SPECIFIC TO THE TEXT FILE OBTAINED.
    def gliese_read(self, textdir, textname):
        # Get filer
        filer = hdf5_writer(self.directory, self.filename)

        # Move to dir + read + clip off random crap at start
        os.chdir(textdir)
        gliese = open(textname, "r")
        gliese_list = []
        for x in gliese:
            gliese_list.append(x)
        gliese.close()
        # Then we need to get the data
        for num, i in enumerate(gliese_list):
            isplit = i.split(";")
            isplit_nought = isplit[0].replace(" ", "")
            if isplit_nought == str(1):
                gliese_list = gliese_list[num:len(gliese_list)]
                break

        # 6,7,8 are pythonic indices for U,B,V in the extracted table.
        gliese_index = []
        u, b, v = [], [], []
        for num, i in enumerate(gliese_list):
            istrip = i.replace(" ", "")
            isplit = istrip.split(";")
            try:
                if isplit[6] != "~":
                    if isplit[7] != "~":
                        if isplit[8] != "~":
                            u.append(float(isplit[6])), b.append(float(isplit[7])), v.append(
                                float(isplit[8])), gliese_index.append(int(isplit[0]))
            except:
                pass

        # Create a phottable and save all the indices.
        phottable = Table()
        phottable["id"], phottable["U"], phottable["B"], phottable["V"] = gliese_index, u, b, v
        filer.write_table("GLIESE", "phottable_table", phottable)

    # Should have 5 lists. RA, DEC, and bands 1,2,3, labelled in capitals for what they are. Tricolour will force
    def import_votable_gaia(self, filename, datagroup, dataset):
        # Set up filer
        filer = hdf5_writer(self.directory, self.filename)
        # Get votable, also in the rootdir/etc.
        votable = parse_single_table(filename).to_table(use_names_over_ids=True)
        # Save votable as Astropy Table. Bit of dancing, here...
        names, tabled = votable.colnames, Table()
        for name in names:
            tabled[name] = votable[name]

        filer.write_table(datagroup, dataset, tabled)
        filer.write(datagroup, dataset + "group_names", names) # Weird bug with this... ".mask"

        # Will also go ahead and write in a table using the band names individually, for the sake of it.
        phottable = Table()
        g, rp, bp = [],[],[]
        for i in tabled:
            if np.isnan(i["phot_g_mean_mag"]) == False:
                if np.isnan(i["phot_bp_mean_mag"]) == False:
                    if np.isnan(i["phot_rp_mean_mag"]) == False:
                        g.append(i["phot_g_mean_mag"]), rp.append(i["phot_rp_mean_mag"]), bp.append(i["phot_bp_mean_mag"])

        phottable['g'], phottable['rp'], phottable['bp'] = g,rp,bp
        phottable['id'] = np.arange(0, len(phottable))
        filer.write(datagroup, "mag_table", phottable)

    # Note that it rejects any stars with d_parallax/parallax >= 0.1  and error on B-V >= 0.025 (# 0.1/0.025 are the default from release paper)
    def import_votable_hipviz(self, filename, datagroup, dataset, par_fraction, bv_error):
        # Set up filer
        filer = hdf5_writer(self.directory, self.filename)
        # Get votable, also in the rootdir/etc.
        votable = parse_single_table(filename).to_table(use_names_over_ids=True)
        # Save votable as Astropy Table. Bit of dancing, here...
        names, tabled = votable.colnames, Table()
        for name in names:
            tabled[name] = votable[name]

        print(tabled.info)

        filer.write_table(datagroup, dataset, tabled)
        filer.write(datagroup, dataset + "group_names", names)  # Weird bug with this... ".mask"

        # Will also go ahead and write in a table using the band names individually, for the sake of it.
        # The V-Mag here is ABSOLUTE V MAGNITUDE using the parallax.
        # Parallax in mas.
        # This section trims down so that the parallax error is less than 10% and the bv error is less than 0.025
        phottable = Table()
        vt, bt, vtabs, plx, dist, distsec = [], [], [], [], [], []
        for i in tabled:
            if np.isnan(i["BTmag"]) == False:
                if np.isnan(i["VTmag"]) == False:
                    parratio = i['e_Plx']/i['Plx']
                    if np.abs(parratio) <= par_fraction:
                        bv_err = np.sqrt((i['e_BTmag'])**2 + (i['e_VTmag'])**2)
                        if bv_err <= bv_error:
                            v_mag = i["VTmag"]
                            parallax = (1*10**-3)*(np.pi/180) * (1/3600) * i['Plx']
                            distance = 1.496*(10**11)/parallax
                            parsec = 3.26*(3*(10**8)*(365*86400))
                            distsec_ratio = distance/(10*parsec)
                            v_mag_absolute = v_mag - 5*np.log10(distsec_ratio)
                            vt.append(i["VTmag"]), bt.append(i["BTmag"]), vtabs.append(v_mag_absolute), plx.append(i['Plx']), dist.append(distance), distsec.append(distsec_ratio)

        phottable['vt'], phottable['bt'], phottable['vtabs'], phottable['Plx'], phottable['Dist'], phottable['Distsec'] = vt, bt, vtabs, plx, dist, distsec
        phottable['id'] = np.arange(0, len(phottable))
        filer.write(datagroup, "mag_table_errored", phottable)

        # This section will do the above, but will NOT consider error/etc.
        phottable = Table()
        vt, bt, vtabs, plx, dist, distsec = [], [], [], [], [], []
        for i in tabled:
            if np.isnan(i["BTmag"]) == False:
                if np.isnan(i["VTmag"]) == False:
                    v_mag = i["VTmag"]
                    parallax = (1 * 10 ** -3) * (np.pi / 180) * (1 / 3600) * i['Plx']
                    distance = 1.496 * (10 ** 11) / parallax
                    parsec = 3.26 * (3 * (10 ** 8) * (365 * 86400))
                    distsec_ratio = distance / (10 * parsec)
                    v_mag_absolute = v_mag - 5 * np.log10(distsec_ratio)
                    vt.append(i["VTmag"]), bt.append(i["BTmag"]), vtabs.append(v_mag_absolute), plx.append(i['Plx']), dist.append(distance), distsec.append(distsec_ratio)

        phottable['vt'], phottable['bt'], phottable['vtabs'], phottable['Plx'], phottable['Dist'], phottable['Distsec'] = vt, bt, vtabs, plx, dist, distsec
        phottable['id'] = np.arange(0, len(phottable))
        filer.write(datagroup, "mag_table", phottable)

    # Import the M31 stuff. We only want the colour-colour relationship (U-B, B-V) values. V is inconsequential.
    # Format is N REF V B-V U-B
    def plei_import(self, filename, group, dataset):
        # Import table
        table = np.genfromtxt(filename, delimiter=";")
        V,B,U = [],[],[]
        for i in table:
            if np.isnan(i[3]) == False:
                if np.isnan(i[4]) == False:
                    V.append(i[2])
                    b = i[3] + i[2]
                    B.append(b)
                    u = b + i[4]
                    U.append(u)
        # Phottable the result
        table = Table()
        table['V'], table['B'], table['U'] = V,B,U
        # Save
        filer = hdf5_writer(rootdir, "data.hdf5")
        filer.write_table(group, dataset, table)

# GAIA interface class
# G, GBP, GRP ~ R,V,I FILTER BANDS (approximately...)
class GAIA(object):
    def __init__(self, directory, filename):
        self.filename = filename
        self.directory = directory
        self.group = "GAIA"
        self.username = gaia_username
        self.password = gaia_password

    # Brings back ALL targets within RADIUS (arcmin) of target ra,dec (degrees) and saves to appropriate folder under "GAIA" group as a query table. Will proceed to try to "clean up" the data + save it under dataset.
    # Also gets proper motions/etc in mas, per year
    def obtain_data(self, radius, ra, dec, dataset):
        # Get data + instantiate the filer.
        filer = hdf5_writer(self.directory, self.filename)
        Gaia.ROW_LIMIT = -1
        coordinate = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
        w, h = u.Quantity(radius, u.arcmin), u.Quantity(radius, u.arcmin)
        query = Gaia.query_object_async(coordinate=coordinate, width=w, height=h)
        table = Table(names=['ra','dec','bp','g','rp','par','par_err', 'pmra', 'pmdec', 'pmra_error', 'pmdec_error'])
        # Go through and filter out the bad apples.
        for i in query:
            if np.isnan(i['phot_g_mean_mag']) != True:
                if np.isnan(i['phot_bp_mean_mag']) != True:
                    if np.isnan(i['phot_rp_mean_mag']) != True:
                        # G, BP, RP equivalent (kinda like) B U V. Anyway.
                        ra,dec,bp,g,rp,par,par_err, pmra, pmdec, pmra_error, pmdec_error = np.array(i['ra']), np.array(i['dec']), np.array(i['phot_bp_mean_mag']), np.array(i['phot_g_mean_mag']), np.array(i['phot_rp_mean_mag']), np.array(i['parallax']), np.array(i['parallax_error']), np.array(i['pmra']), np.array(i['pmdec']), np.array(i['pmra_error']), np.array(i['pmdec_error'])
                        table.add_row([ra,dec,bp,g,rp,par,par_err, pmra, pmdec, pmra_error, pmdec_error])
        # Write table
        filer.write_table(self.group, dataset, table)
    """
    # Crossmatch service. Much code curtosy of https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
    def cross_match(self, input_table, radii, savename):
        # Log in to the GAIA service to authenticate
        radii_new = radii #u.Quantity(copy.copy(radii), u.arcsec)
        Gaia.login(user=gaia_username, password=gaia_password)


        # Quick check
        #tables = Gaia.load_tables(only_names=True, include_shared_tables=True)
        #print(tables)

        try:
            Gaia.delete_user_table(table_name='my_sources', force_removal=True)
        except:
            pass

        table_named = 'my_sources'
        Gaia.upload_table(upload_resource=input_table, table_name=table_named)
        time.sleep(5)

        Gaia.update_user_table(table_name=table_named,
                               list_of_changes = [["ra", "flags", "Ra"],
                                                  ["dec", "flags", "Dec"]])

        # Other misc stuff that I don't get but works
        full_qualified_table_name = "user_" + gaia_username + "." + table_named
        crossmatched_table_name = savename
        crossmatched_table = "user_" + gaia_username + "." + crossmatched_table_name

        # Set up the job
        Gaia.cross_match(full_qualified_table_name_a = full_qualified_table_name,
                         full_qualified_table_name_b = 'gaiadr2.gaia_source',
                         results_table_name = crossmatched_table_name, radius=radii_new) # Radius defines the radius to search for a match

        # Obtain the results
        query = ('SELECT c."dist"*3600 as dist, a.*, b.* FROM gaiadr2.gaia_source AS a, '
                 'full_qualified_table_name AS b, '
                 'crossmatched_table_name AS c '
                 'WHERE (c.gaia_source_source_id = a.source_id AND '
                 'c.my_sources_my_sources_oid = b.my_sources_oid)')
        job = Gaia.launch_job(query=query)
        results = job.get_results()
        print(f"results = {results}")
        return results"""

# Class for extraction and atmospheric correction of sources from FITS files.
# Also has bicolour/tricolour graphing functions for extracted sources.
class source_analysis(object):
    # Directory should be the root directory for the project, in which the hdf5 file of name "filename" resides.
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename

    # Tool to quickly generate file tree for the source analyser.
    # Note that the initial files must manually be positioned appropriately, alone, in the start directories.
    # For automatic positioning, add argument in main_runtime.py to copy over astrometrified average_stack_fully_stacked_cluster_band.fitswcs.fits files.
    def anatree_maker(self, clusters, bands):
        # Generate file structure.
        try:
            os.mkdir(self.directory + "\\Sci\\SEX")
        except:
            pass
        try:
            os.mkdir(self.directory + "\\Sci\\SEX\\Images")
        except:
            pass

        for cluster in clusters:
            try:
                os.mkdir(self.directory + "\\Sci\\SEX" + "\\" + cluster)
            except:
                continue
            for band in bands:
                try:
                    os.mkdir(self.directory + "\\Sci\\SEX" + "\\" + cluster + "\\" + band)
                except:
                    continue





    # Returns a list of sources from an image NAME given a menagerie of parameters, functions using IRAF starfinder.
    # Ignores sources past the linearity threshold of 55,000
    # "Identifier" specifies the directory in which to save the extracted data, including the astropy table from the starfinder.
    # squarecoords define a square region to use in the image for source extraction as a list: ["x1,y1","x2,y2"] as FITS IMAGE FORMAT from Aladin, bottom left and top right of image rectangle (pythonically)
    """
    The roundness statistic is defined as follows roundness = (hx - hy) / (hx + hy) where and hx and hy are the amplitudes of the best fitting 1-D gaussian
    """
    def anasource_extraction(self, name, identifier, squarecoords, threshold, fwhm, minsep_fwhm, sigma_radius, sharplo, roundlo, sharphi, roundhi, minflux):
        # Define the starfinder operation alongside other miscellanies
        starfinder = IRAFStarFinder(peakmax=55000, threshold=threshold, sigma_radius=sigma_radius, fwhm=fwhm, minsep_fwhm=minsep_fwhm, sharplo=sharplo, sharphi=sharphi, roundlo=roundlo, roundhi=roundhi)
        #starfinder = DAOStarFinder(peakmax=55000, threshold=threshold, sigma_radius=sigma_radius, fwhm=fwhm, sharplo=sharplo, sharphi=sharphi, roundlo=roundlo, roundhi=roundhi, ratio=1)
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()

        # Load data and find stars in image
        file = fits.open(name)
        data, header = file[0].data, file[0].header

        # Trim data down to the appropriate square area for analysis, retaining original shape of array.
        data = utilities.fittrim_coords(name, squarecoords[0], squarecoords[1])

        # Extract sources using the provided FWHM (blanket extraction for all sources using same FWHM.
        startable = starfinder(data)

        # Sort the sources in terms of the instrumental flux, in descending order.
        startable.sort('flux', reverse='true')

        # Find the index beneath which flux drops under minflux and clip the table
        for num,i in enumerate(startable):
            flux = i['flux']
            if flux <= minflux:
                startable = startable[0:num]
                break

        # Alongside that, we go ahead and retrieve the RA,DEC coordinates for all the targets...
        w = astropy.wcs.WCS(header)
        x,y = startable['xcentroid'],startable['ycentroid']
        xy = np.array([x,y])
        xy = xy.T # Coordinate list of all centroids.
        for i in range(len(xy)):
            xy[i] = w.wcs_pix2world(xy[i][0], xy[i][1], 0)

        # We'll save those too! Gotta split 'em up first, though.
        id = np.arange(0, len(x), 1)
        startable['id'] = id
        xy = xy.T
        ra,dec = xy[0],xy[1]
        filer.write(identifier, "starfind_ra",ra)
        filer.write(identifier, "starfind_dec", dec)

        # Attach the ra,dec to the starfind table
        startable['ra'] = ra
        startable['dec'] = dec
        startable['id'] = np.arange(1, len(dec) + 1, 1)

        # Iraf starfinder retrieves the FWHM's and other parameters. We can write entire table to hdf5 ^.^
        astropy.io.misc.hdf5.write_table_hdf5(startable, self.directory + "\\" + self.filename, path=identifier + "/" + "starfind_table", append=True, overwrite=True)

        # We also attach the determined sigma-clipped mean FWHM to the table
        mean,median,std = sigma_clipped_stats(startable['fwhm'])
        filer.write(identifier, "sigclip_fwhm", mean)

        # For completeness, we also return the table, for any usage that may be needed
        return startable

    # Generate a smoothed 2D background map for a fits image NAME
    def anaback_2d(self, name, segment_size):
        # Load file and data
        data, header = fits.open(name)[0].data, fits.open(name)[0].header
        # Generate bkimg
        bkg = background.Background2D(data, box_size=segment_size)
        background_image = bkg.background

        # Package data background (for debug purposes)
        newfitso = fits.PrimaryHDU(background_image)
        newfitso.header = header
        newfitsodata = fits.HDUList(newfitso)
        newfitsodata.writeto("debug_background_map.fits", overwrite=True)

    # This takes a list of RA and DECs for targets, including their FWHM's, and runs simple aperture photometry for them, ONE AT A TIME.
    # No area trim here, the skycoords for targets should lie in the stack overlap area already.
    # Very slow process, but photutils lacks the required tools to give a list of aprads/anrads, sadly.
    # aprad and anrad are in units of FWHM, recommended to be equal to (3,4) respectively by Ross
    # All results are saved to the data.hdf5 file (or whatever file was used in instantiating this class)
    def anasource_simpaertures(self, name, identifier, radec_identifier, aprad, anradin, anradout):
        # Load requisites
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()


        # Get data, header, and starfind table
        data, header = fits.open(name)[0].data, fits.open(name)[0].header
        startable = filer.read_table(radec_identifier, "starfind_table")
        startable.sort('flux', reverse='true')

        # Record the altitude/filter/exposure for the file under the IDENTIFIER group in the data file
        # Random sleep to prevent parallel read-write by processes
        random_value = random.randint(0, 20)
        time.sleep(random_value)

        filer.write(identifier, "altitude", header['ALTITUDE'])
        filer.write(identifier, "filter", header['FILTER'])
        filer.write(identifier, "exposure", header['EXPOSURE'])

        # Get sigma-clipped-stats for the data for uniform background consideration
        sig_mean, sig_median, std = sigma_clipped_stats(data)

        # Load the radec_identifier starfind_ra,starfind_dec information alongside the sig_clip FWHM
        ra,dec = filer.read(radec_identifier, "starfind_ra"),filer.read(radec_identifier, "starfind_dec")
        radec = np.array([ra,dec])
        radec = radec.T
        fwhm = filer.read(radec_identifier, "sigclip_fwhm")


        # Generate WCS and convert ra,dec into xy coordinates for this
        w = wcs.WCS(header)
        xy = []
        for rd in radec:
            xy_wcs = w.wcs_world2pix(rd[0],rd[1],0)
            xy.append(xy_wcs)
        xy = np.array(xy)

        # Generate the apertures and annulars
        apertures, annulars = CircularAperture(xy, r=aprad*fwhm), CircularAnnulus(xy, r_in=fwhm*anradin, r_out=anradout*fwhm)
        aperturinos = [apertures, annulars]

        # Go ahead and generate photometry data and attach various quantitites to phottable
        phottable = aperture_photometry(data, aperturinos)
        phottable['ra'], phottable['dec'] = ra, dec

        # Next go through and subtract the background using the annuli.
        anu_corrected_count = []
        anu_corrected_count_err = []
        sig_corrected_count = []
        sig_corrected_count_err = []

        for i in phottable:
            # Calculate fluxes
            apflux,anflux = np.abs(i['aperture_sum_0']),i['aperture_sum_1']
            aparea,anarea = np.pi*((aprad*fwhm)**2), np.pi*((anradout*fwhm)**2 - (anradin*fwhm)**2)
            apareadivanarea = aparea/anarea
            anfluxinaparea = anflux*apareadivanarea
            apflux_corrected = apflux - anfluxinaparea
            anu_corrected_count.append(apflux_corrected)

            # Calculate errors
            apflux_err, anflux_err = np.sqrt(apflux), np.sqrt(np.abs(anflux))
            anfluxinaparea_err = anflux_err*apareadivanarea
            apflux_corrected_error = np.sqrt(apflux_err**2 + anfluxinaparea_err**2)
            anu_corrected_count_err.append(apflux_corrected_error)

            # For the sake of considering static background, we'll also use SCM (entire field) for background corrections.
            bkg_pua = sig_median
            anfluxinaparea = aparea*bkg_pua
            apflux_corrected = apflux - anfluxinaparea
            apflux_err, anfluxinaparea_err = np.sqrt(apflux), np.sqrt(np.abs(anfluxinaparea))
            apflux_corrected_error = np.sqrt((apflux_err**2) + (anfluxinaparea_err**2))
            sig_corrected_count.append(apflux_corrected), sig_corrected_count_err.append(apflux_corrected_error)





        # Add the corrected counts/errors to the phottable
        phottable['apflux_annuli_corr'], phottable['apflux_annuli_corr_err'] = anu_corrected_count, anu_corrected_count_err


        # And for the sigma background
        phottable['apflux_sig_corr'], phottable['apflux_sig_corr_err'] = sig_corrected_count, sig_corrected_count_err

        # Write table
        astropy.io.misc.hdf5.write_table_hdf5(phottable, self.directory + "\\" + self.filename, path=identifier + "/" + "phottable_table", append=True, overwrite=True, serialize_meta=True)

        # Voluntary stat reduction (just a quick check to see if we can use a flat-out background instead of anrads
        #utilities.stat_reduce(name)

        # Voluntary source visualization/etc
        utilities.source_visualizer(data, startable, fwhm, identifier)


    # Re-written
    def anasource_manuapertures_redo(self, name, identifier, radec_identifier, aprad, anradin, anradout, region_clip, save_figs):

        # Set up requisites
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()

        # Load in the data
        data, header = fits.open(name)[0].data, fits.open(name)[0].header
        phottable = copy.deepcopy(filer.read_table(identifier, "phottable_table"))

        # Also grab a "default fwhm" to use (the one simpaurtures or whatever used)
        default_fwhm = filer.read(radec_identifier, "sigclip_fwhm")
        # Grab the iraf fwhms
        fwhm_list = np.array(filer.read_table(radec_identifier, "starfind_table")['fwhm'])

        # Set up WCS
        w = wcs.WCS(header)

        # Attempt to make the necessary debug directory.
        print("Attempting to write Debug file directory")
        try:
            os.mkdir(self.directory + "\\Debug")
        except:
            pass
        try:
            os.mkdir(self.directory + "\\Debug" + "\\" + identifier + "_manu")
        except:
            pass

        alist, aelist, fwlist = [],[],[]

        # Go through list, do photometry, etc.
        for num,i in enumerate(phottable):
            ra, dec = i['ra'], i['dec']
            x, y = w.wcs_world2pix(ra, dec, 0)
            xy = [x,y]

            # Redefine default_fwhm according to the iraf starfind fwhm
            default_fwhm = fwhm_list[num]

            clip = utilities.array_trim(data, xy, int(region_clip))
            a, ae, fw = utilities.array_gauss_phot(clip, aprad, anradin, anradout, identifier + "_manu", str(num + 1), save_figs, default_fwhm)
            alist.append(a), aelist.append(ae), fwlist.append(fw)

        print(identifier, "...clip photometry success.")

        # The resultant ap/apfl/fwhm are written to group_id under "num_clipdata" in that order. We need to re-read from the file + collate into a table.

        print("Attaching Data Step, boom. ", identifier)
        phottable['apflux_annuli_manu_corr'], phottable['apflux_annuli_manu_corr_err'], phottable['manu_fwhm'] = alist, aelist, fwlist

        print("Writing Table!")
        # Write table
        random_sleep = random.randint(0, 20)
        time.sleep(random_sleep)
        filer.write_table(identifier, "phottable_table", phottable)







    # Given the group identifier, will go and correct for ATMOSPHERIC ABSORPTION, altering the phottable "phottable_table" under the "identifier" group
    def ana_atmos_correction(self, identifier):
        # Get required classes going
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()

        # Load in the phottable
        print("Atmospherically correcting ", identifier)
        phottable = astropy.io.misc.hdf5.read_table_hdf5(input=self.directory + "\\" + self.filename, path=identifier + "/" + "phottable_table")

        # Identify the filter, load in altitude/exposure, and get optical depth, and also get the zenith.
        altitude, filter, exposure = filer.read(identifier, "altitude"), filer.read_string(identifier, "filter"), filer.read(identifier, "exposure")
        zenith = np.deg2rad(90 - altitude)
        airmass = 1/(np.cos(zenith))
        optical_depth, optical_depth_error = filer.read("true_depth", ("optical_depth_{}").format(filter)) # [depth error]

        # Quick Debug for B-Band M52
        #if identifier == "M52_B":
        #    print(exposure, filter, altitude, zenith, airmass, optical_depth, optical_depth_error)
        #    time.sleep(50)



        # Also grab the zeropoint for a preliminary test of getting scientific magnitudes
        zero_point = filer.read("true_zero_point", ("zero_point_{}").format(filter)) # Get z_p and err


        # Indices for types of fluxes, so far just the sigma-clipped-bkg and annular, and errors.
        # index has "_corr" for the flux, and "_corr_err" for the flux error, in the dataset. See data.hdf5 for samples.
        indices = ["apflux_annuli", "apflux_sig", "apflux_annuli_manu"]

        # Flux and Flux_Err lists
        fluxes, fluxerrs, instmags, instmagerrs = [[] for d in indices],[[] for d in indices],[[] for d in indices],[[] for d in indices]

        # Also zeropoint/sci_mag lists
        scimags, scimagerrs = [[] for d in indices], [[] for d in indices]

        # For all the stars inside the phottable, calculate atmospherically corrected flux/etc.
        # Remember to consider exposure (normalize to exposure), and do this for both types of background subtraction.
        for num, index in enumerate(indices):
            for i in phottable:
                # Get fluxes/errors
                fluxdex, fluxerrdex = index + "_corr", index + "_corr_err"
                flux, flux_err = i[fluxdex], i[fluxerrdex]

                # Correct and get error
                flux_atmospherically_done = flux * np.exp(optical_depth * airmass) * 1 / exposure
                flux_atmospherically_done_error = (1 / exposure) * np.sqrt((flux_err * np.exp(optical_depth * airmass)) ** 2 + (optical_depth_error * airmass * flux_atmospherically_done) ** 2)
                inst_mag = -1 * (5 / 2) * np.log10(flux_atmospherically_done)
                inst_mag_error = -1 * (5 / 2) * np.log10(np.exp(1)) * flux_atmospherically_done_error / flux_atmospherically_done

                inst_mag = inst_mag
                flux_atmospherically_done_error = np.abs(flux_atmospherically_done_error)
                inst_mag_error = np.abs(inst_mag_error)

                # Get the scientific magnitude and associated error
                sci_mag = inst_mag + zero_point[0]
                sci_mag_err = np.sqrt(zero_point[1]**2 + inst_mag_error**2)

                # Attach to the phottable
                fluxes[num].append(flux_atmospherically_done), fluxerrs[num].append(flux_atmospherically_done_error), instmags[num].append(inst_mag), instmagerrs[num].append(inst_mag_error)

                # And the sci_mag
                scimags[num].append(sci_mag), scimagerrs[num].append(sci_mag_err)

        # Attach the atmofluxes/errs to the list. THESE ARE NORMALIZED TO THE EXPOSURE REMEMBER THAT YOU CUCK
        for num, index in enumerate(indices):
            fluxdex, fluxerrdex, instdex, insterrdex = index + "_norm_atmos_corr", index + "_norm_atmos_corr_err", index + "_norm_instmag_corr", index + "_norm_instmag_corr_err"
            phottable[fluxdex], phottable[fluxerrdex], phottable[instdex], phottable[insterrdex] = fluxes[num], fluxerrs[num], instmags[num], instmagerrs[num]

            # We'll also do the scientific magnitudes.
            scimagdex, scimagerrdex = index + "atmos_scimag", index + "atmos_scimag_err"
            phottable[scimagdex], phottable[scimagerrdex] = scimags[num], scimagerrs[num]

        # Write table
        astropy.io.misc.hdf5.write_table_hdf5(phottable, self.directory + "\\" + self.filename,path=identifier + "/" + "phottable_table_pre_reddenning", append=True,overwrite=True, serialize_meta=True)

    # Collect up all the catalogues for CLUSTER, BANDS. NOT CLUSTER(S), but CLUSTER.
    # Condense the data, keeping the RA, DEC, ATMOS_FLUX, ATMOS_FLUX_ERR, SCI_MAG, SCI_MAG_ERR
    def ana_mag_cleanup_1(self, cluster, bands):
        # Get file reader
        filer = hdf5_writer(self.directory, self.filename)
        # Create table
        phottable_new = Table(names=['id'])
        # Import bands, and attach the details.
        for band in bands:
            table_path = cluster + "_" + band
            table_name = "phottable_table_pre_reddenning" # Pre-reddenned data
            phottable = astropy.io.misc.hdf5.read_table_hdf5(input=self.directory + "\\" + self.filename, path=table_path + "/" + table_name)
            # Snatch up the ra,dec,etc
            phottable_new['id'], phottable_new['ra'], phottable_new['dec'] = phottable['id'], phottable['ra'], phottable['dec']
            # Also snatch up fluxes/mags/errs
            indices = ["apflux_annuli", "apflux_sig", "apflux_annuli_manu"]
            for index in indices:
                phottable_new[("{}_{}_{}").format(band, index, "norm_atmos_corr")]  = phottable[("{}_{}").format(index, "norm_atmos_corr")]
                phottable_new[("{}_{}_{}").format(band, index, "norm_atmos_corr_err")] = phottable[("{}_{}").format(index, "norm_atmos_corr_err")]
                phottable_new[("{}_{}_{}").format(band, index, "atmos_scimag")] = phottable[("{}{}").format(index, "atmos_scimag")]
                phottable_new[("{}_{}_{}").format(band, index, "atmos_scimag_err")] = phottable[("{}{}").format(index, "atmos_scimag_err")]

        # Save the file
        astropy.io.misc.hdf5.write_table_hdf5(phottable_new, self.directory + "\\" + self.filename,path=cluster + "_" "PREREDCATALOGUE" + "/" + "phottable_table_prered", append=True,overwrite=True, serialize_meta=True)

    # Create bi-colour diagram (B-V, V, etc.)
    # group/dataset is table location/name in file. band1/2 are the band prefixes (capitalized, i.e. "U", "B").
    # band_suffix describes the suffix to the band for the magnitudes, i.e. "B_apflux_annuli_atmos_scimag"... error is auto taken as "_err" on the end of this.
    # By default, colour is put along the X axis, and the magnitude along the Y axis. "flip" = True will flip this.
    # Colour is calculated via BAND[1] - BAND[0], with BAND[0] put along the brightness axis.
    # Limits should be set for x, y. If autolims is False, they will be used. Else, pyplot sets lims automatically.
    # If absolute_V is anything but "False" then this provided band (no suffix, just from the file straight-up) will be used for the Y axis.
    def ana_bicolour(self, group, dataset, band_suffix, band_1_2, flip, autolims, xlims, ylims, titler, absolute_V, dotwidth):
        # Define classes
        filer = hdf5_writer(self.directory, self.filename)

        # Obtain table
        phottable = filer.read_table(group, dataset)

        # Mag_lists
        mag_1_2 = [[] for d in band_1_2]

        # Obtain all the indices, magnitudes, etc, for band1/band2 in order.
        for num,band in enumerate(band_1_2):
            mag = phottable[band + band_suffix]
            mag_1_2[num] = mag

        # Trim the list down a bit...
        #for num, i in enumerate(mag_1_2):
        #    mag_1_2[num] = i[0:50]

        # Array-up
        mag_1_2 = [np.array(d) for d in mag_1_2]

        # Calculate colour
        colour_2_1 = mag_1_2[1] - mag_1_2[0] # Band 2 - Band 1

        # Set axes (just to avoid confusion later)
        x_values, y_values = colour_2_1, mag_1_2[0]

        # Set the labels
        xlabel, ylabel = ("{} - {}").format(band_1_2[1], band_1_2[0]), ("{} (Apparent)").format(band_1_2[0])

        # Check for absolute_V
        if absolute_V != False:
            ylabel = ("{} (Absolute)").format(absolute_V)
            y_values = phottable[absolute_V]

        # Check for flip...
        if flip == True:
            x_values, xlabel, y_values, ylabel, xlims, ylims = y_values, ylabel, x_values, xlabel, ylims, xlims


        # Create the plot.
        fig = plt.figure(figsize=(3, 4), dpi=300)
        ax1 = fig.subplots()
        ax1.scatter(x=x_values, y=y_values, marker=".", lw=0.05, color="black", s=dotwidth)
        ax1.set(title=titler, # Insert a title, dummy is here. Baaaaaaaka nya!
                ylabel = ylabel,
                xlabel = xlabel)
        if autolims == False:
            ax1.set_ylim(ylims)
            ax1.set_xlim(xlims)

        ax1.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax1.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        ax1.plot()
        fig.show()
        os.chdir(self.directory + "\\Sci\\SEX\\Images")
        fig.savefig(group + "_" + dataset + ("_{}-{}").format(band_1_2[1], band_1_2[0])+"_BICOLOUR.png")

    # Create a tri-colour diagram (U - B, B - V, etc.). Suffix should include a preliminary _
    # By default will automatically save to SEX/Images
    def ana_tricolour(self, group, dataset, band_suffix, band_1_2_3, flip, autolims, point_label, xlims, ylims, listclip, titler, dotwidth):
        # Define classes
        filer = hdf5_writer(self.directory, self.filename)

        # Obtain table
        phottable = filer.read_table(group, dataset)


        # Mag_lists
        mag_1_2_3 = [[] for d in band_1_2_3]

        # Obtain all the indices, magnitudes, etc, for band1/band2 in order.
        for num, band in enumerate(band_1_2_3):
            mag = phottable[band + band_suffix]
            mag_1_2_3[num] = mag

        #Trim the list down a bit...
        if listclip != False:
            for num, i in enumerate(mag_1_2_3):
                mag_1_2_3[num] = i[0:listclip]

        # Array-up
        mag_1_2_3 = [np.array(d) for d in mag_1_2_3]

        # Calculate colour
        colour_2_1 = mag_1_2_3[1] - mag_1_2_3[0]  # Band 2 - Band 1
        colour_3_2 = mag_1_2_3[2] - mag_1_2_3[1] # Band 3 - Band 1

        # For debug, we save.
        phottable_new = Table()
        phottable_new["colour_21"] = colour_2_1
        phottable_new["colour_32"] = colour_3_2
        filer.write_table(group, dataset + band_suffix +"_.", phottable_new)

        # Set axes (just to avoid confusion later)
        x_values, y_values = colour_2_1, colour_3_2

        # Set the labels
        xlabel, ylabel = ("{} - {}").format(band_1_2_3[1], band_1_2_3[0]), ("{} - {}").format(band_1_2_3[2], band_1_2_3[1])



        # Check for flip...
        if flip == True:
            x_values, xlabel, y_values, ylabel, xlims, ylims = y_values, ylabel, x_values, xlabel, ylims, xlims


        # Create the plot.
        fig = plt.figure(dpi=300)
        fig.set_size_inches(5,7)
        ax1 = fig.subplots()

        # Set labels for each point, if wished.
        if point_label == True:
            # Also set ID against all the individual tags, if you wish.
            ids = phottable['id']
            if listclip != False:
                ids = ids[0:listclip]
            texts = [ax1.text(x_values[i], y_values[i], ids[i], ha='center', va='center') for i in range(len(x_values))]
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

        ax1.scatter(x=x_values, y=y_values, marker=".", lw=0.05, color="black", s=dotwidth)
        ax1.set(title=titler + "..." + band_suffix,  # Insert a title, dummy is here. Baaaaaaaka nya!
                ylabel=ylabel,
                xlabel=xlabel)
        if autolims == False:
            ax1.set_ylim(ylims)
            ax1.set_xlim(xlims)

        ax1.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax1.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        ax1.figure.set_size_inches(5,7)
        ax1.plot()
        fig.show()
        os.chdir(self.directory + "\\Sci\\SEX\\Images")
        fig.savefig(group + "_" + dataset + "_" + band_suffix + "_TRICOLOUR.png")
        plt.clf()
        plt.close()



    # Takes cluster(s), band(s), band suffix, and estimated distance(s), and quickly runs a correction on them.
    # Spits out a new band suffix/etc (attaches testabs behind the end, i.e. _apflux_annuli_manu_atmos_scimag goes to _apflux_annuli_manu_atmos_scimag_testabs)
    def ana_disttest_corr(self, datagroups, datasets, clusters, bands, band_suffix, distance_over_10pc):
        # Instantiate filer
        filer = hdf5_writer(self.directory, self.filename)
        # Correct for distance
        for num, cluster in enumerate(clusters):
            phottable = copy.deepcopy(filer.read_table(datagroups[num], datasets[num]))
            for band in bands:
                uncor_label = band + band_suffix
                uncor_mag = phottable[uncor_label]
                cor_value = -5*np.log10(distance_over_10pc[num])
                cor_mag = uncor_mag + cor_value
                cor_label = uncor_label + "_testabs"
                phottable[cor_label] = cor_mag
            filer.write_table(datagroups[num], datasets[num], phottable)

    """
    # Dataclip writer. DEPRECATED!
    def anasource_clipwriter(self, identifier, data, xy, region_clip):
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()
        # Generate data clips for each of the xy. The clips will have to be saved to disk, way too many to keep.
        # Clips are used for photometry, and even smaller clips (of the clips) are used for the gauss fit.
        save_identifier = identifier + "_clips"
        print("Writing dataclips for ", identifier)
        for num, coord in enumerate(xy):
            filer.write(save_identifier, str(num), utilities.array_trim(data, coord, region_clip))
    # Long-winded simpaertures. DEPRECATED AND REPLACED!!! (This code does not work.)
    # Will clip data by region_clip (the rectangular radius) and run a gauss fit on this region, and then run photometry on it, and return sigma-clipped background for it.
    def anasource_manuapertures(self, name, identifier, radec_identifier, aprad, anradin, anradout, region_clip, save_figs):
        # PARALLELISATION TEST
        #name, identifier, radec_identifier, aprad, anradin, anradout, region_clip, save_figs = splitstring.split(";")
        #aprad, anradin, anradout, region_clip = float(aprad), float(anradin), float(anradout), int(region_clip)
        # , name, identifier, radec_identifier, aprad, anradin, anradout, region_clip, save_figs):
        # Set up requisites
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()

        # Load in the data
        data, header = fits.open(name)[0].data, fits.open(name)[0].header
        phottable = filer.read_table(radec_identifier, "starfind_table")
        phottable.sort('flux', reverse='true')
        phot_new = Table()
        phot_new['id'], phot_new['ra'], phot_new['dec'] = phottable['id'], phottable['ra'], phottable['dec']

        # Record the altitude/filter/exposure for the file under the IDENTIFIER group in the data file
        filer.write(identifier, "altitude", header['ALTITUDE'])
        filer.write(identifier, "filter", header['FILTER'])
        filer.write(identifier, "exposure", header['EXPOSURE'])

        # Get ra,dec
        ra,dec = phottable['ra'], phottable['dec']
        radec = np.array([ra, dec])
        radec = radec.T # in the form of [[ra1,dec1], [ra2,dec2]...etc]
        fwhm = filer.read(radec_identifier, )

        # Generate WCS and convert ra,dec into xy coordinates for this
        w = wcs.WCS(header)
        xy = []
        for rd in radec:
            xy_wcs = w.wcs_world2pix(rd[0], rd[1], 0)
            xy.append(xy_wcs)
        xy = np.array(xy) # in the form of [[x1,y1], [x2,y2],...]

        # Write dataclips
        #self.anasource_clipwriter(identifier, data, xy, region_clip)

        # Also make the necessary debug directory.
        print("Attempting to write Debug file directory")
        try:
            os.mkdir(self.directory + "\\Debug")
        except:
            pass
        try:
            os.mkdir(self.directory + "\\Debug" + "\\" + identifier + "_manu")
        except:
            pass

        # Run for each data clip, the solver.
        cwd = os.getcwd()
        # Specify the directory for the "manu" data...
        group_id = identifier + "_manu"
        save_identifier = identifier + "_clips" # Group for num_clipdata

        ap, apfl, fwhm = [],[],[]
        # Get photometry for clips
        print(identifier, "...Running clip photometry")
        for num in range(len(xy)):
            clip = filer.read(save_identifier, str(num))
            print(num)
            # Try-except loop
            try:
                timed_phot = utilities.timeout(15)(utilities.array_gauss_phot)
                a,ae,fw = timed_phot(clip, aprad, anradin, anradout, identifier + "_manu", num, save_figs, group_id, str(num) + "_clipdata")
                a, ae, fw = utilities.array_gauss_phot(clip, aprad, anradin, anradout, identifier + "_manu", num, save_figs)
                ap.append(a), apfl.append(ae), fwhm.append(fw)
            except:
                # In the case of error, try a stronger clip.
                clip = utilities.array_trim(data, xy[num], int(0.5*region_clip))
                try:
                    timed_phot = utilities.timeout(15)(utilities.array_gauss_phot)
                    a,ae,fw = timed_phot(clip, aprad, anradin, anradout, identifier + "_manu", num, save_figs, group_id, str(num) + "_clipdata")
                    a, ae, fw = utilities.array_gauss_phot(clip, aprad, anradin, anradout, identifier + "_manu", num, save_figs)
                    ap.append(a), apfl.append(ae), fwhm.append(fw)
                except:
                    # The stronger clip still failed. Use nan values.
                    print("Extended clip failed in manu, appending NaN.")
                    a,ae,fw = [np.nan, np.nan, np.nan]
                    ap.append(a), apfl.append(ae), fwhm.append(fw)
                    print("Failed clip...", num)
                    time.sleep(20)
            os.chdir(cwd)

        print(identifier, "...clip photometry success.")

        # The resultant ap/apfl/fwhm are written to group_id under "num_clipdata" in that order. We need to re-read from the file + collate into a table.

        print("Attaching Data Step, boom. ", identifier)
        phot_new['apflux_annuli_manu_corr'], phot_new['apflux_annuli_manu_corr_err'], phot_new['manu_fwhm'] = ap, apfl, fwhm

        print("Writing Table!")
        # Write table
        filer.write_table(identifier + "_TEST","phottable_table", phot_new)

        # Voluntary source visualization/etc
        #utilities.source_visualizer(data, phottable, fwhm, identifier)

        return True
    """

# Class containing functions related to determination of cluster ages via turnoff point deduction.
class turnoff_analysis(object):
    # Directory should be the root directory for the project, in which the hdf5 file of name "filename" resides.
    def __init__(self, directory, filename):
        self.directory = directory
        self.filename = filename

    # Generate file hierarchy for resources associated with class. Once generated, add in any other requisites.
    def turn_treemaker(self):
        # Rootdir + Sci + TURNOFF + Images (only directory as of yet. More may follow.)
        try:
            os.mkdir(rootdir + "\\Sci\\TURNOFF")
        except:
            pass
        try:
            os.mkdir(rootdir + "\\Sci\\TURNOFF\\Images")
        except:
            pass

    # Generate absolute V-band magnitudes given the table
    # Megan/Olivia has parallax/etc + Syazza has V-band mag, so you need to do this for two dif tables.
    def turn_absmag_gen(self, MAGPARgroup, MAGPARset):
        # Get filer + table + utilities
        filer = hdf5_writer(self.directory, self.filename)
        utilities = utils()

        # Get table that has pars/mags
        par_table = filer.read_table(MAGPARgroup, MAGPARset)
        parname, parerrname, Vname, Verrname = 'par', 'par_err', 'V', 'V_err'

        # Abs mag list
        abs_mags, abs_mags_errs, dists, disterrs = [], [], [], []

        # Get absmags/etc
        for num,i in enumerate(par_table):
            par, par_err = i[parname], i[parerrname] # In arcsec/''
            distance, distance_error = utilities.par_to_dist(par, par_err)
            v_app, v_app_err = i[Vname], i[Verrname]
            v_abs, V_abs_err = utilities.abs_to_app(v_app, v_app_err, distance, distance_error)
            abs_mags.append(v_abs), abs_mags_errs.append(V_abs_err), dists.append(distance), disterrs.append(distance_error)

        # Attach abs_mags
        par_table['V_abs'], par_table['V_abs_err'], par_table['dist'], par_table['dist_err'] = abs_mags, abs_mags_errs, dists, disterrs

        # Resave the table (with the abs-mags now included)
        # clust[num] + "_MEMBERS", "corrected_ubv"
        filer.write_table(MAGPARgroup, MAGPARset, par_table)


    # Bicolour plotter for provided V_abs/BV colour
    def turn_bicolour(self, group, dataset, V_ABS, BV_COLOUR, autolims, xlims, ylims, titler, dotwidth):
        # Define classes
        filer = hdf5_writer(self.directory, self.filename)

        # Obtain table
        phottable = filer.read_table(group, dataset)

        # Mag_lists
        BV, ABS_V = phottable[BV_COLOUR], phottable[V_ABS]

        # Set axes (just to avoid confusion later)
        x_values, y_values = BV, ABS_V

        # Set the labels
        xlabel, ylabel = "B - V", "V"

        # Invert the Y axis
        #ylims = [min(y_values - 0.1), max(y_values + 0.1)]

        # Create the plot.
        fig = plt.figure(figsize=(5, 4), dpi=300)
        ax1 = fig.subplots()
        ax1.scatter(x=x_values, y=y_values, marker=".", lw=0.05, color="black", s=dotwidth)
        ax1.set(title=titler, # Insert a title, dummy is here. Baaaaaaaka nya!
                ylabel = ylabel,
                xlabel = xlabel)

        if autolims == False:
            ax1.set_ylim(ylims)
            ax1.set_xlim(xlims)

        ax1.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax1.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        plt.gca().invert_yaxis()

        ax1.plot()
        fig.show()
        os.chdir(self.directory + "\\Sci\\TURNOFF\\Images")
        fig.savefig(group + "_" + dataset + "_" + "VBV_BICOLOUR.png")

    # Bicolour plotter for provided V_abs/BV colour, that overplots with hipparcos
    def turn_bicolour_hipparcos(self, group, dataset, V_ABS, BV_COLOUR, autolims, xlims, ylims, titler, dotwidth):
        # Define classes
        filer = hdf5_writer(self.directory, self.filename)

        # Obtain table
        phottable = filer.read_table(group, dataset)

        # Mag_lists
        BV, ABS_V = phottable[BV_COLOUR], phottable[V_ABS]

        # Set axes (just to avoid confusion later)
        x_values, y_values = BV, ABS_V

        # Set the labels
        xlabel, ylabel = "B - V", "V"

        # Grab the Hipparcos stuff (hardcoded)
        hip_table = filer.read_table("HIPPARCOS", "mag_table_errored")
        vabs, bv = hip_table['vtabs'], hip_table['bt'] - hip_table['vt']

        # Create the plot.
        fig = plt.figure(figsize=(5, 4), dpi=300)
        ax1 = fig.subplots()
        ax1.scatter(x=bv, y=vabs, marker=".", lw=0.1, color="black", s=dotwidth)
        ax1.scatter(x=x_values, y=y_values, marker=".", lw=0.05, color="red", s=dotwidth)

        legend_elements = [Patch(facecolor='black', edgecolor='black', label="Hipparcos"),
                           Patch(facecolor='red', edgecolor='black', label="Obtained")]
        ax1.legend(handles=legend_elements, loc='upper right')



        ax1.set(title=titler, # Insert a title, dummy is here. Baaaaaaaka nya!
                ylabel = ylabel,
                xlabel = xlabel)

        if autolims == False:
            ax1.set_ylim(ylims)
            ax1.set_xlim(xlims)

        ax1.grid(True, which='major', color="pink", alpha=1, linestyle='dotted', lw=0.5)  # Enable grids on subplot
        ax1.grid(True, which='minor', color="pink", alpha=1, linestyle='dotted', lw=0.5)

        plt.gca().invert_yaxis()

        ax1.plot()
        fig.show()
        os.chdir(self.directory + "\\Sci\\TURNOFF\\Images")
        fig.savefig(group + "_" + dataset + "_" + "VBV_BICOLOUR.png")

    # Misc plotting
    def turn_meanstarplotter(self):
        filer = hdf5_writer(rootdir, "data.hdf5")
        data = filer.read_table("MEANSTARS", "spec_table")
        plt.plot(data['Mv'], data['logL'])
        plt.show()


    # Estiamtes Teff/Luminosity/Mass/other_parameters for all stars inside a table
    def turn_mspectype_estimator(self, group, set):
        # Set all indicators
        BV = "BV"

        # Grab table of data for processing (i.e. our data)
        filer = hdf5_writer(rootdir, "data.hdf5")
        table = filer.read_table(group, set)
        BV = table['BV']

        # Get sample table + generate T(BV) function
        sample_table = filer.read_table("MEANSTARS", "spec_table")

        # Teff(BV)
        teff_bv = interp1d(sample_table['B-V'], sample_table['Teff'], fill_value="extrapolate")
        # M(Mv)
        M_Mv = interp1d(sample_table['Mv'], sample_table['Msun'])
        # M(BV)
        M_bv = interp1d(sample_table['B-V'], sample_table['Msun'])


        # For our data, generate T_eff
        TEFF = []
        for bv in BV:
            TEFF.append(teff_bv(bv))
        table['Teff'] = TEFF

# Class containing function to estimate cluster membership. Written hastily and dirtily.
class member_finder(object):
    def __init__(self, root, file):
        self.directory = root
        self.filename = file
        self.filer = hdf5_writer(root, file)

    # Take table (group, set) and create proper motion scatter plots.
    def pmscat(self, group, set, pmra, pmdec):
        table = self.filer.read_table(group, set)
        pmr, pmd = table[pmra], table[pmdec]
        fig, ax = plt.subplots(1)
        ax.scatter(pmr, pmd, s=0.5, color="black", xlabel="pmra", ylabel="pmdec", title="PM Scatter for " + group)







#owo = turnoff_analysis(rootdir, "data.hdf5")
#owo.turn_meanstarplotter()
#owo.turn_mspectype_estimator("M52_MEMBERS","corrected_ubv")

#heyo = misc_file_handler(rootdir, "data.hdf5")

#heyo.syazza_porter(rootdir + "\\Sci\\TURNOFF", "corrected_ubv_M52.txt", "CLUST_TRIM", "M52_TEST")
#heyo.meanstar_dwarf()
#owo = source_analysis(rootdir, "data.hdf5")
#owo.ana_tricolour("GLIESE", "phottable_table", "", ["V","B","U"], False, False, False, [-0.5, 2], [2.5, -1], False, "Gliese Example CC", 4)
#owo.ana_tricolour("NGC7789_GAIA_TEST", "mag_table", "", ["rp", "g", "bp"], False, False, False, [0.4,1.8], [1.5, -0.5], False, "Near 7789", 4)

#votableha = misc_file_handler(rootdir, "data.hdf5")
#votableha.import_votable_hipviz("vizier_votable.vot", "HIPPARCOS", "hip_votable", 0.1, 0.025) # 0.1/0.025 are the default from release paper

# Finito Bonito Baby.
print("Loaded utilities task complete.")
