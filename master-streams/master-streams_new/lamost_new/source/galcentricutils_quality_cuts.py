from astropy.stats import sigma_clipped_stats
from numba import njit, types

import windows_directories_new
from galcentricutils_new import galconversion
from energistics_new import fast_energistics_new, orbigistics
import numpy as np
import os
import astropy.units as u
import pickle

class quality_cuts(object):

    def __init__(self):

        self.converter = galconversion()
        self.converter.solinfo_grab(windows_directories_new.sourcedir, "solar_info.dat")
        self.orbigist = orbigistics()

        # Path for data cache
        self.savedir = os.path.join(windows_directories_new.lamodir, "qc_cache")
        try:
            os.mkdir(self.savedir)
        except:
            pass

    def LSR_cut(self, table, cut, extra_text=None):

        """
        Do an LSR cut based on the provided LSR velocity (given to the converter).
        - Takes into account peculiar motions based on Schoenrich 2010.

        :param table:
        :param cut:
        :return:
        """

        # https://arxiv.org/abs/0912.3693
        LSR_VELOCITY = np.array([11.1, 12.24, 7.25], float)
        LSR_VELOCITY_ERR = np.array([0.74,0.47,0.37], float)

        absdifs = np.array([table['vx'],table['vy'],table['vz']]).T - \
                  (self.converter.sol_params[2]-LSR_VELOCITY)
        absdifs = np.linalg.norm(absdifs, axis=1)
        indices_to_keep = np.where(absdifs >= cut)[0]
        indices_to_remove = np.where(absdifs < cut)[0]

        # Create a plot if required
        if extra_text != None:

            R, vR, vT, z, vz, phi = self.orbigist.get_leftgalpy(table)  # Get R in units of kpc
            feh = table['feh']  # get fehs
            import matplotlib.pyplot as plt

            # Set up figure
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.set(xlim=[-550,550],
                    ylim=[-5/2,1/2],
                    xlabel=r'$\textrm{v}_\textrm{los}$',
                    ylabel=r'$[\textrm{FeH}]$')
            axs.grid(True, which='major', color='pink')

            # Plot data
            axs.scatter(vT[indices_to_remove], feh[indices_to_remove], color='red', marker='x', s=0.2)
            axs.scatter(vT[indices_to_keep], feh[indices_to_keep], color='green', marker='x', s=0.2)
            savepath = os.path.join(self.savedir, extra_text + "_LSR.png")
            plt.savefig(savepath, dpi=300)
            plt.close()
            plt.clf()

        return indices_to_keep

    def feh_cut(self, table, cut):

        """

        Cut below the provided feh.

        :param table: or pandas
        :param cut:
        :return: indices
        """

        indices = np.where(table['feh'] < cut)

        return indices

    def fancy_feh_cut(self, table, extra_text=None):

        R, vR, vT, z, vz, phi = self.orbigist.get_leftgalpy(table) # Get R in units of kpc
        feh = table['feh'] # get fehs

        # Rescale VT to feh (it is pretty damn different you know)
        mean_vT = np.mean(np.abs(vT))
        mean_feh = np.mean(np.abs(feh))
        scale = mean_vT/mean_feh
        vT /= scale

        # Cluster
        twod_vectors = np.array([feh, vT]).T
        from hdbscan import flat
        clusterer = flat.HDBSCAN_flat(X=twod_vectors,
                                      n_clusters=2,
                                      min_cluster_size=int(len(feh)/10),
                                      min_samples=7,
                                      metric='l2',
                                      algorithm='best',
                                      prediction_data=True)
        labels = clusterer.labels_
        available_clusters = np.array(list(set(labels)), int)
        available_clusters = available_clusters[np.where(available_clusters != -1)[0]]
        size = [len(np.where(labels==label)[0]) for label in available_clusters]
        disk_star_index = available_clusters[np.argmin(size)]
        indices_to_remove = np.where(labels==disk_star_index)[0]
        indices_to_keep = np.where(labels!=disk_star_index)[0]

        # Create a plot if required
        if extra_text != None:

            import matplotlib.pyplot as plt
            vT *= scale

            # Set up figure
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.set(xlim=[-550,550],
                    ylim=[-5/2,1/2],
                    xlabel=r'$\textrm{v}_\textrm{los}$',
                    ylabel=r'$[\textrm{FeH}]$')
            axs.grid(True, which='major', color='pink')

            # Plot data
            axs.scatter(vT[indices_to_remove], feh[indices_to_remove], color='red', marker='x', s=0.2)
            axs.scatter(vT[indices_to_keep], feh[indices_to_keep], color='green', marker='x', s=0.2)
            savepath = os.path.join(self.savedir, extra_text + "_fancy_feh.png")
            plt.savefig(savepath, dpi=300)
            plt.close()
            plt.clf()

        return indices_to_keep


    def bound_cut(self, table):

        """

        Get rid of unbound stars (E>0) according to the model given to fast_energistics.

        :param table: or pandas, really.
        :return: indices
        """

        E = fast_energistics_new().default_E_c_vals(table)
        indices = np.where(E<=0)[0]

        return indices

    def zmax_cuts(self, table, zmax,
                  time_to_integrate,
                  number_of_steps,
                  try_load=False,
                  try_save=False,
                  filename=None,
                  extra_text=None):

        """

        Cut based on ZMAX obtained through an appropriate orbit integration.

            - Take table
            - Integrate all stars backward in time by time_to_integrate
            - For number_of_steps
            - Get np.max() for galactocentric Z
            - If greater than ZMAX, accept

        Use load/save if you are running the same table in succession and want to avoid reloading
        the orbit integrations in the future (good for debug.)

        :param table: or pandas
        :param ZMAX: float
        :param time_to_integrate: int
        :param number_of_steps: int
        :param try_load: bool
        :param try_save: bool
        :param filename: the filename for savedir
        :return: indices
        """

        savepath = os.path.join(self.savedir, filename)

        try:

            if try_load == True:

                with open(savepath, "rb") as f:
                    orbits = pickle.load(file=f)

            else:

                raise RuntimeError("Trying to generate orbits for quality cuts...")

        except:

            orbits = self.orbigist.orbits(table)
            times = np.linspace(0, -1 * time_to_integrate, number_of_steps) * u.yr
            orbits.integrate(times, self.orbigist.pot, 'rk4_c')

            if try_save == True:

                with open(savepath, "wb") as f:
                    pickle.dump(file=f, obj=orbits)


        R, vR, vT, z, vz, phi = orbits.getOrbit().T
        R *= self.orbigist.rovo[0]
        vR *= self.orbigist.rovo[1]
        vT *= self.orbigist.rovo[1]
        z *= self.orbigist.rovo[0]
        vz *= self.orbigist.rovo[1]

        """
        
        A note on the format of the integration results.
        The result of the above (including the tranpose) is a vector of the shape...
        
            (number_of_steps, len(table))
            
        I.e. each column is a given stars orbit integrated over. So, to get the max
        of R for example,
        
            max_R_of_zero = np.max(R[:,0])
            
        Gives the max of all rows (all steps) for star 0.
        
        """

        # Get the zmaxes
        zmax_list = np.max(z, axis=0)
        indices_to_keep = np.where(zmax_list > zmax)
        indices_to_remove = np.where(zmax_list <= zmax)

        # Create a plot if required
        if extra_text != None:

            import matplotlib.pyplot as plt

            feh=table['feh']
            R, vR, vT, z, vz, phi = self.orbigist.get_leftgalpy(table)  # Get R in units of kpc

            # Set up figure
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.set(xlim=[-550,550],
                    ylim=[-5/2,1/2],
                    xlabel=r'$\textrm{v}_\textrm{los}$',
                    ylabel=r'$[\textrm{FeH}]$')
            axs.grid(True, which='major', color='pink')

            # Plot data
            axs.scatter(vT[indices_to_remove], feh[indices_to_remove], color='red', marker='x', s=0.2)
            axs.scatter(vT[indices_to_keep], feh[indices_to_keep], color='green', marker='x', s=0.2)
            savepath = os.path.join(self.savedir, extra_text + "_ZMAX_CUT.png")
            plt.savefig(savepath, dpi=300)
            plt.close()
            plt.clf()

        # Return
        return indices_to_keep

    @staticmethod
    @njit(fastmath=True)
    def _fancy_LSR_cut(vR, vT, vz, phi, vcirc):

        # Set up empties for the circular vectors and the actual data vectors
        vec_vcirc = np.empty((len(vcirc),3), types.float64)
        vec_datas = np.empty_like(vec_vcirc)

        # Set up the vcirc vectors
        vec_vcirc[:, 0] = -vcirc*np.sin(phi)
        vec_vcirc[:, 1] = vcirc*np.cos(phi)
        vec_vcirc[:, 2] = 0

        # Set up the vec_true vectors
        vec_datas[:, 0] = -vT*np.sin(phi) + vR*np.cos(phi)
        vec_datas[:, 1] = vT*np.cos(phi) + vR*np.sin(phi)
        vec_datas[:, 2] = vz

        # Difs
        vel_difs = vec_vcirc - vec_datas
        return np.sqrt(vel_difs[:,0]**2 + vel_difs[:,1]**2 + vel_difs[:,2]**2)

    def fancy_LSR_cut(self, table, cut, extra_text=None):

        """

        Do a fancy-cut. This will...

            - For star
            - Elucidate local disk velocity vector (i.e. LSR equivalent, basically)
            - Cut as usual, but relative to this LSR value (rotational velocity curve vector basically.)

        Will hopefully eliminate stars that belong to the disk, but be sensitive to different areas of the disk
        around the galaxy that may otherwise have different LSR to our local LSR.

        Note: Does not work that well- not very conservative, it seems.

        :param table: or pandas
        :param cut: float
        :return: indices

        """

        R, vR, vT, z, vz, phi = self.orbigist.get_leftgalpy(table) # Get R in units of kpc
        from galpy.potential.Potential import vcirc
        vcirc = np.array(vcirc(self.orbigist.pot, R*u.kpc), float)
        norms = self._fancy_LSR_cut(vR, vT, vz, phi, vcirc)
        feh = table['feh']
        indices_to_keep = np.where(norms>=cut)[0]
        indices_to_remove = np.where(norms<cut)[0]

        # Create a plot if required
        if extra_text != None:

            import matplotlib.pyplot as plt

            # Set up figure
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.set(xlim=[-550,550],
                    ylim=[-5/2,1/2],
                    xlabel=r'$\textrm{v}_\textrm{los}$',
                    ylabel=r'$[\textrm{FeH}]$')
            axs.grid(True, which='major', color='pink')

            # Plot data
            axs.scatter(vT[indices_to_remove], feh[indices_to_remove], color='red', marker='x', s=0.2)
            axs.scatter(vT[indices_to_keep], feh[indices_to_keep], color='green', marker='x', s=0.2)
            savepath = os.path.join(self.savedir, extra_text + "_fancy_LSR.png")
            plt.savefig(savepath, dpi=300)
            plt.close()
            plt.clf()

        return indices_to_keep

    def GSE_cut(self, table, extra_text=None):

        """
        Do a cut using L-space to within |L| of cut about 0,0,0 (target the GSE, basically)
        only on stars near the FeH of the GSE

        https://academic.oup.com/mnras/article-abstract/508/1/1489/6370605
        About -1.15 (we will use a width of about 0.25.)

        :param table:
        :param cut:
        :return:

        """

        GSE_FEH = -1.15
        GSE_DISP = 0.25
        L_MAG_CUT = np.sqrt(3) * 1500
        L_MAG_DEFINITE_CUT = 2000

        L = np.sqrt(table['Lx']**2 + table['Ly']**2 + table['Lz']**2)
        feh = table['feh']
        feh_disp = np.abs(feh-GSE_FEH)
        truefalse = [False if A < GSE_DISP and B < L_MAG_CUT else True
                     for (A,B) in zip(feh_disp, L)]
        truefalse_secondary = [False if A < L_MAG_DEFINITE_CUT else True for A in L]
        truefalse = np.array(truefalse, bool)
        truefalse_secondary = np.array(truefalse_secondary, bool)
        truefalse = truefalse*truefalse_secondary
        indices_to_keep = np.arange(0, len(L), 1)[truefalse]
        truefalse = -1*np.array(truefalse, bool)
        indices_to_remove = np.arange(0, len(L), 1)[truefalse]

        # Create a plot if required
        if extra_text != None:

            import matplotlib.pyplot as plt
            R, vR, vT, z, vz, phi = self.orbigist.get_leftgalpy(table)  # Get R in units of kpc

            # Set up figure
            fig, axs = plt.subplots(nrows=1, ncols=1)
            axs.set(xlim=[-550,550],
                    ylim=[-5/2,1/2],
                    xlabel=r'$\textrm{v}_\textrm{los}$',
                    ylabel=r'$[\textrm{FeH}]$')
            axs.grid(True, which='major', color='pink')

            # Plot data
            axs.scatter(vT[indices_to_remove], feh[indices_to_remove], color='red', marker='x', s=0.2)
            axs.scatter(vT[indices_to_keep], feh[indices_to_keep], color='green', marker='x', s=0.2)
            savepath = os.path.join(self.savedir, extra_text + "_GSE_REMOVAL.png")
            plt.savefig(savepath, dpi=300)
            plt.close()
            plt.clf()

        return indices_to_keep