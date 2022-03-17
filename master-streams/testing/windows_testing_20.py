import copy
import os
import pickle
import time
import astropy.units as u
from matplotlib.patches import Patch
import galpy
from galpy.util import bovy_coords
import numpy as np
from galpy import orbit
from matplotlib import pyplot as plt
import ascii_info
import hdfutils
import windows_directories
from energistics import orbigistics, orbifitter


"""
Least Squares Preliminary Fitting.
- iteration defines number of orbtis to generate. Orphan is good for 4,000
- time_to_integrate is number of years forward or backward for the orbit. Orphan is good for 0.1e9
- number_of_steps defines the number of points to generate for the integration. Orphan is good for 1000
- clust_to_fit defines which cluster you want to fit, from the data. 
- try_load specifies whether we should try reloading old orbits to save on generation time: saves a lot of time for MC.
"""
def least_squares_preliminary_orbit_fitting_thing_with_a_long_name(table,
                                                                   membership_table,
                                                                   clust_to_fit,
                                                                   iterations,
                                                                   time_to_integrate,
                                                                   number_of_steps,
                                                                   try_load=True,
                                                                   graph=False):
    """
    Method for orbit fitting:
    - Establish the range of orbits that fit the data
    - Set up parameter space (assume gaussians, use covariance matrix and means, generate set.)
    - Monte-carlo generate a set of orbits within this space (n-iterations to test.)
    - Get least-squares estimate of best fit
    This will get you the orbital parameters of the best-fitting orbit for that set.
    Galpy orbit generation is the bottleneck: make sure to save the monte-carlo'd orbits for re-use.
    Once a good orbit has been generated, a markov-chain can be used to search near it to see if an improved orbit exists.
    The same basic monte-carlo'd orbital elements should realistically be usable for all sets of the data, saving comp time.
    """


    # Galactocentric Coordinate Frame Least-squares Orbit Fitting for the data.
    """
    TODO:
    Change fit to fit inside galactic coordinates, instead of ra/dec space. 
    """

    # Set up integration time
    integrate_time = np.linspace(0, time_to_integrate, number_of_steps)*u.yr  # in years.

    # Grab the clustering you're fitting from the membership table and probabilities for weighting.
    clustering = membership_table['probable_clust']
    clustering_probability = membership_table['probability']
    clustselec = [True if d == clust_to_fit else False for d in clustering]

    # Try to craft specific directory
    try:
        os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables_guess")
    except:
        pass
    savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables_guess" + "\\" + str(clust_to_fit)
    try:
        os.mkdir(savedir)
    except:
        pass

    # Set up orbigistics (for potential) and cylindrification, and also the orbit fitter
    orbigist = orbigistics()
    # Set up the orbit fitter
    orbifitt = orbifitter()

    # Clip data
    data_to_fit = table[clustselec]

    # Get data as array in left-handed galpy, for the fitting.
    R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
    # Renormalize our data to Galpy ROVO, since that's what orbits are outputted in.
    R /= orbigist.rovo[0]
    vR /= orbigist.rovo[1]
    vT /= orbigist.rovo[1]
    z /= orbigist.rovo[0]
    vz /= orbigist.rovo[1]
    # Set up table values of these re-normalized values.
    data_to_fit['R_galpy'], \
    data_to_fit['vR_galpy'], \
    data_to_fit['vT_galpy'], \
    data_to_fit['z_galpy'], \
    data_to_fit['vz_galpy'], \
    data_to_fit['phi_galpy'] = R, \
                               vR, \
                               vT, \
                               z, \
                               vz, \
                               phi
    data_array = np.empty((6, len(R)))
    data_array[0,:], \
    data_array[1,:], \
    data_array[2,:], \
    data_array[3,:], \
    data_array[4,:], \
    data_array[5,:] = R, vR, vT, z, vz, phi  # replace xyzvxvyvz with R,vR,vT,z,vz,phi (in radians)
    data_array = data_array.T # convert to vectors. note that we're now in a left-handed system.

    # TODO: Important Note for the Report!!! Defines local solar.
    """
    Use the solarmotion='schoenrich' local solar velocities. 
    (as noted by Drimmel and Poggio as being decent- these are the defaults.)
    solarmotion='schoenrich'
    """

    # Try to load orbits if try_load is true.
    try:
        if try_load == True:
            with open(
                    windows_directories.orbitsdir + "\\" + ("preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_{2:.0f}_timesteps").format(clust_to_fit,iterations,number_of_steps) + ".txt", 'rb') as f:
                orbit_list = pickle.load(file=f)
            with open(
                    windows_directories.orbitsdir + "\\" + ("preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_elements_{2:.0f}_timesteps").format(clust_to_fit,iterations,number_of_steps)
                            .format(clust_to_fit, iterations) + ".txt", 'rb') as f:
                orbit_elements = pickle.load(file=f)
            print("Orbits successfully loaded for " + str(clust_to_fit) + " with an iteration number of " + str(iterations))
        else:
            raise ValueError("Not loading orbits, Generating them.")
    # Failed to load them. Generate them.
    except Exception as e:
        if try_load != False:
            print("Orbits failed to load. Generating preliminary Monte orbits.")
        else:
            print(e)

        # Grab the data in galactic coords
        l, b, dist, dmu_l_cosdec, dmu_b, vlos = data_to_fit['l'],\
                                                data_to_fit['b'],\
                                                data_to_fit['dist'],\
                                                data_to_fit['dmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                data_to_fit['dmu_b'],\
                                                data_to_fit['vlos'] # deg deg kpc mas/yr mas/yr km/s

        # Set up the means and deviations
        lm, bm, distm, dmulcosdecm, dmubm, vlosm = np.mean(l), np.mean(b), np.mean(dist), \
                                                   np.mean(dmu_l_cosdec), np.mean(dmu_b), np.mean(vlos)
        ldev, bdev, distdev, dmulcosdecdev, dmubdev, vlosdev = l - lm, b - bm, dist - distm, \
                                                               dmu_l_cosdec - dmulcosdecm, dmu_b - dmubm, vlos - vlosm
        devlist = [ldev, bdev, distdev, dmulcosdecdev, dmubdev, vlosdev]

        # Set up covariance matrix and mean vector, then generate points
        mean_vector = np.array([lm, bm, distm, dmulcosdecm, dmubm, vlosm])
        covtrix = np.empty(shape=(6,6))
        for i in range(6):
            for j in range(6):
                covtrix[i,j] = np.mean(devlist[i]*devlist[j])
        orbit_elements = np.random.default_rng().multivariate_normal(mean=mean_vector, cov=covtrix, size=iterations)

        # For all this orbit elements, generate orbit data (list of arrays of orbit vectors: each array is [vec1, vec2...]
        start = time.time()
        orbit_list = []
        for element in orbit_elements:
            element = list(element)
            element[0]*=u.deg
            element[1]*=u.deg
            element[2]*=u.kpc
            element[3]*=u.mas/u.yr
            element[4]*=u.mas/u.yr
            element[5]*=u.km/u.s
            orbit_forward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0]*u.kpc, vo=orbigist.rovo[1]*u.km/u.s, zo=orbigist.zo*u.kpc, lb=True)
            orbit_backward = orbit.Orbit(vxvv=element, ro=orbigist.rovo[0]*u.kpc, vo=orbigist.rovo[1]*u.km/u.s, zo=orbigist.zo*u.kpc, lb=True)
            orbit_forward.integrate(integrate_time, orbigist.pot, 'rk4_c')
            orbit_backward.integrate(-integrate_time, orbigist.pot, 'rk4_c')
            orbits = np.concatenate([orbit_forward.getOrbit(), orbit_backward.getOrbit()], axis=0)
            orbit_list.append(orbits)
            # Galpy returns are in units of RO and VO from orbigist, here at least.
            # Phi is also [-pi,pi] instead of [0,2pi.] Be aware of this.
            # To rectify this, re-normalize our data to rovo: we took care of this earlier for this reason.

            """
            # Try Plotting (for test.) DEBUG DEBUG DEBUG 
            R, vR, vT, z, vz, phi = orbits.T  # the fitted parameters, which now need to be converter.
            phi = [np.pi - d for d in phi]  # convert galpy left to astropy right
            phi = np.array(phi)
            x, y = R * np.cos(phi), R * np.sin(phi)
            table = Table()
            table['x'], table['y'], table['z'] = x, y, z
            table = galcentricutils.angular().get_polar(table)  # get the fit in galactocentric polars
            plt.scatter(table['theta'], table['phi'], color='green', marker='x', s=0.05)  # plot the fit. """
        end = time.time()
        print("Integration of initial orbits too " + str(end - start) + " seconds, for cluster " + str(clust_to_fit))
        orbit_list = np.array(orbit_list)

        # Try to save the orbits that have been generated
        try:
            with open(windows_directories.orbitsdir + "\\" + ("preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_{2:.0f}_timesteps").format(clust_to_fit,iterations,number_of_steps) + ".txt", 'wb') as f:
                pickle.dump(obj=orbit_list, file=f)
        except:
            pass

        # Also try to save the orbital elements we've generated.
        try:
            with open(windows_directories.orbitsdir + "\\" + ("preliminary_fit_cluster_{0:.0f}_with_{1:.0f}_MCdraws_elements_{2:.0f}_timesteps").format(clust_to_fit,iterations,number_of_steps) + ".txt", 'wb') as f:
                pickle.dump(obj=orbit_elements, file=f)
        except:
            pass

    # Generate preliminary MC least squares fits
    argmin, least_list = orbifitt.least_multisquares(data_array, orbit_list, clustering_probability)

    # Grab the correct orbit element and integrate/etc, get the fit, plot it, etc, in galactocentric polars.
    argmin_element = orbit_elements[argmin] # l, b, dist, dmu_l_cosdec, dmu_b, vlos
    argmin_element = list(argmin_element)
    argmin_element[0]*=u.deg
    argmin_element[1]*=u.deg
    argmin_element[2]*=u.kpc
    argmin_element[3]*=u.mas/u.yr
    argmin_element[4]*=u.mas/u.yr
    argmin_element[5]*=u.km/u.s
    best_fit = orbit.Orbit(vxvv=argmin_element, ro=orbigist.rovo[0]*u.kpc, vo=orbigist.rovo[1]*u.km/u.s, zo=orbigist.zo*u.kpc, lb=True)

    # Produce graphs if true.
    if graph == True:
        best_fit_reverse = orbit.Orbit(vxvv=argmin_element, ro=orbigist.rovo[0]*u.kpc, vo=orbigist.rovo[1]*u.km/u.s, zo=orbigist.zo*u.kpc, lb=True)
        best_fit.integrate(integrate_time, orbigist.pot, method='rk4_c')
        best_fit_reverse.integrate(-integrate_time, orbigist.pot, method='rk4_c')
        best_fits_full = np.concatenate([np.flipud(best_fit_reverse.getOrbit()), best_fit.getOrbit()], axis=0)
        R, vR, vT, z, vz, phi = best_fits_full.T  # the fitted parameters, which now need to be converted.
        R*=orbigist.rovo[0]
        vR*=orbigist.rovo[1]
        vT*=orbigist.rovo[1]
        z*=orbigist.rovo[0]
        vz*=orbigist.rovo[1]
        modelfits = [R, vR, vT, z, vz, phi]
        X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
        l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T

        # Main example sky plot in galactic longitude/latitude
        data_array = data_array.T
        Rdata, vRdata, vTdata, zdata, vzdata, phidata = data_array[0,:], \
                                                        data_array[1,:], \
                                                        data_array[2,:], \
                                                        data_array[3,:], \
                                                        data_array[4,:], \
                                                        data_array[5,:]
        Rdata*=orbigist.rovo[0]
        vRdata*=orbigist.rovo[1]
        vTdata*=orbigist.rovo[1]
        zdata*=orbigist.rovo[0]
        vzdata*=orbigist.rovo[1]
        datafits = [Rdata, vRdata, vTdata, zdata, vzdata, phidata]
        Xd, Yd, Zd = bovy_coords.galcencyl_to_XYZ(Rdata, phidata, zdata, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
        ld, bd, dd = bovy_coords.XYZ_to_lbd(Xd, Yd, Zd, degree=True).T
        plt.scatter(ld, bd, color='black', marker='x')
        plt.scatter(l, b, color='red', marker='x', s=0.5) # plot the fit.
        plt.xlabel('l')
        plt.ylabel('b')
        plt.savefig(savedir + "\\" + str(clust_to_fit) + "_all_orbits_testgalactic.png", dpi=300)
        plt.close()

        # Set up permutations for comparative plots
        permutations = [[0,1],[0,2],[0,3],[0,4],[0,5],[1,2],[1,3],[1,4],[1,5],[2,3],[2,4],[2,5],[3,4],[3,5],[4,5]]
        titles = ["R vs. vR",
                  "R vs. vT",
                  "R vs. z",
                  "R vs. vz",
                  "R vs. phi",
                  "vR vs. vT",
                  "vR vs. z",
                  "vR vs, vz",
                  "vR vs. phi",
                  "vT vs. z",
                  "vT vs. vz",
                  "vT vs. phi",
                  "z vs. vz",
                  "z vs. phi",
                  "vz vs. phi"]

        # Create 1-1 plots to demonstrate error directions
        for permutation, title in zip(permutations, titles):
            x,y = permutation
            plt.scatter(datafits[x], datafits[y], color='black', marker='x', s=0.1)
            plt.scatter(modelfits[x], modelfits[y], color='red', marker='x', s=0.1)
            plt.title(title)
            plt.savefig(savedir + "\\" + str(clust_to_fit) + "_" + title.replace(" ", "_") + ".png", dpi=300)
            plt.close()

        # Return the best-fitting orbit object (to use as a "base" for the galpy orbit fitting suite.)
        return best_fit
    # Else, just return the best_fit.
    else:
        return best_fit

"""
Galpy Final Fitting
- Use least_squares as initial fit/init_vxvv
- Generate lots of plots
- profit.
- Needs membership table.
"""
def galpy_final_fitting(clust_to_fit,
                        iterations,
                        time_to_integrate,
                        number_of_steps,
                        try_load=True,
                        graph=False):
    # Load in the "mean" percenttable and map
    writer = hdfutils.hdf5_writer(windows_directories.datadir, ascii_info.asciiname)
    membership_table = writer.read_table(ascii_info.fullgroup, "percent_table_greatfitted")

    # Grab the data table that you are fitting
    table = hdfutils.hdf5_writer(windows_directories.datadir,
                                 ascii_info.asciiname).read_table(ascii_info.fullgroup,
                                                                  ascii_info.fullset)

    # Get an initial guess for the orbit with out code (will also produce plots for the orbit- see folder.)
    guess = least_squares_preliminary_orbit_fitting_thing_with_a_long_name(table,
                                                                           membership_table,
                                                                           clust_to_fit,
                                                                           iterations,
                                                                           time_to_integrate,
                                                                           number_of_steps,
                                                                           try_load)

    # Get clustering data and select for relevant data.
    clustering = membership_table['probable_clust']
    clustselec = [True if d == clust_to_fit else False for d in clustering]
    data_to_fit = table[clustselec]

    # Set up orbigistics (for potential and various parameters related to it/rovozo.)
    orbigist = orbigistics()

    # Use the solarmotion='schoenrich' local solar velocities
    # (as noted by Drimmel and Poggio as being decent- these are the defaults in Galpy anyway but we're specifying.)
    solarmotion='schoenrich'

    # Set up the position vectors for galpy, which will use galactic coordinates straight-up, and also use errors
    vxvv = np.empty(shape=(len(data_to_fit), 6))
    vxvv[:,0], vxvv[:,1], vxvv[:,2], vxvv[:,3],vxvv[:,4], vxvv[:,5] = data_to_fit['l'],\
                                                                      data_to_fit['b'],\
                                                                      data_to_fit['dist'],\
                                                                      data_to_fit['dmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                                                                      data_to_fit['dmu_b'],\
                                                                      data_to_fit['vlos']


    # Also set up the errors. Don't forget: edmu is multiplied by cos(dec) in Galpy.
    evxvv = np.empty(shape=(len(data_to_fit), 6))
    evxvv[:,0], \
    evxvv[:,1], \
    evxvv[:,2], \
    evxvv[:,3], \
    evxvv[:,4], \
    evxvv[:,5] = 0.05*np.ones_like(data_to_fit['l']),\
                 0.05*np.ones_like(data_to_fit['l']),\
                 data_to_fit['edist'],\
                 data_to_fit['edmu_l']*np.cos(np.radians(data_to_fit['b'])),\
                 data_to_fit['edmu_b'],\
                 data_to_fit['evlost'] # Jo recommended using a few arcmins of error for l/b

    # Set up the initial orbit guess. Per Bovy, I'm just guessing that pmll here already has the cosdec factor... guessing.
    init_vxvv = [guess.ll(ro=orbigist.rovo[0]),
                 guess.bb(ro=orbigist.rovo[0]),
                 guess.dist(ro=orbigist.rovo[0]),
                 guess.pmll(ro=orbigist.rovo[0],
                            vo=orbigist.rovo[1]),
                 guess.pmbb(ro=orbigist.rovo[0],
                            vo=orbigist.rovo[1]),
                 guess.vlos(ro=orbigist.rovo[0],
                            vo=orbigist.rovo[1])]

    load_fit = False
    if load_fit == False:
        # Run and save the fit
        fit = galpy.orbit.Orbit.from_fit(init_vxvv,
                                         vxvv,
                                         evxvv,
                                         orbigist.pot,
                                         ro=orbigist.rovo[0]*u.kpc,
                                         vo=orbigist.rovo[1]*u.km/u.s,
                                         zo=orbigist.zo*u.kpc,
                                         solarmotion=solarmotion,
                                         lb=True)
        with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt", 'wb') as f:
            pickle.dump(obj=fit, file=f)
    else:
        # Load the fit
        with open(windows_directories.datadir + "\\" + "clustered" + "\\" + str(clust_to_fit) + ".orbit_fit.txt", 'rb') as f:
            fit = pickle.load(file=f)

    # If Graph is true, generate graphs
    if graph == True:
        # Get a "flip" for backward orbit
        fit_backward = copy.copy(fit)

        # Get integrals
        fit.integrate(np.linspace(0, time_to_integrate, number_of_steps)*u.yr, orbigist.pot, 'rk4_c')
        fit_backward.integrate(-np.linspace(0, time_to_integrate, number_of_steps)*u.yr, orbigist.pot, 'rk4_c')
        print(fit.e())
        # Get orbits
        orbits = fit.getOrbit()
        orbits_backward = fit_backward.getOrbit()
        orbits_full = np.concatenate((orbits, orbits_backward), axis=0)
        R, vR, vT, z, vz, phi = orbits_full.T  # the fitted parameters, in units of rovo/etc.
        R*=orbigist.rovo[0]
        z*=orbigist.rovo[0]
        x, y = R*np.cos(phi), R*np.sin(phi)
        X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
        l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
        #r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
        fig = plt.figure()
        plt.grid(which="major", color="pink")
        plt.scatter(l, b, color='red', marker='x', s=0.1)

        # Clip data and get galpy cylindrical coordinates, alongside getting data as an array.
        R, vR, vT, z, vz, phi = orbigist.get_leftgalpy(data_to_fit)
        X, Y, Z = bovy_coords.galcencyl_to_XYZ(R, phi, z, Xsun=orbigist.rovo[0], Zsun=orbigist.zo).T
        l, b, d = bovy_coords.XYZ_to_lbd(X, Y, Z, degree=True).T
        #r, theta, phi = np.array(galcentricutils.angular().right_numba_polar(R, z, phi)).T
        plt.scatter(data_to_fit['l'], data_to_fit['b'], color='black', marker='x')

        # Generate the range for the plotting.
        xlength = np.max(data_to_fit['l']) - np.min(data_to_fit['l'])
        ylength = np.max(data_to_fit['b']) - np.min(data_to_fit['b'])
        xlim = [np.min(data_to_fit['l']) - 0.05*xlength, np.max(data_to_fit['l']) + 0.05*xlength]
        ylim = [np.min(data_to_fit['b']) - 0.05*ylength, np.max(data_to_fit['b']) + 0.05*ylength]
        plt.xlim(xlim)
        plt.ylim(ylim)
        legend_elements = [Patch(edgecolor='gray',facecolor='red',label='Fit'),
                           Patch(edgecolor='gray',facecolor='black',label='Data')]
        plt.legend(handles=legend_elements, loc='upper right')

        try:
            os.mkdir(windows_directories.imgdir + "\\" + "orbit_fitting_variables")
        except:
            pass
        savedir = windows_directories.imgdir + "\\" + "orbit_fitting_variables" + "\\" + str(clust_to_fit)
        try:
            os.mkdir(savedir)
        except:
            pass
        plt.savefig(savedir + "\\" + str(clust_to_fit) + "_final_galpy_fit.png", dpi=300)
        plt.xlabel('l')
        plt.ylabel("b")
        plt.show()


galpy_final_fitting()