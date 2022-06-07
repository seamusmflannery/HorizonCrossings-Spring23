# Author: Nathaniel Ruhl

import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d

# This should be derined in crossing_time_range in the future

class OrbitModel:

    # INPUTS: obs_dict (dict), orbit_model = "mkf", "aster", or "circle"
    # OUTPUTS: position and time at intervals of time step
    @staticmethod
    def define_orbit_model(obs_dict, orbit_model, time_step):
        t0 = obs_dict["crossing_time_range"][0]
        tf = obs_dict["crossing_time_range"][1]
        t_array = np.arange(t0, tf, time_step)
        if orbit_model == "mkf":
            r_mkf, t_mkf = OrbitModel.readMKF(obs_dict["mkf_path"])
            r_array = OrbitModel.interpolate_mkf_position(r_mkf, t_mkf, t_array)
        # temporary solution until we implement read_aster_orbit()
        elif orbit_model == "aster":
            r_mkf, t_mkf = OrbitModel.readMKF(obs_dict["mkf_path"])
            r_array = OrbitModel.interpolate_mkf_position(r_mkf, t_mkf, t_array)
        return r_array, t_array

    # Reads the orbital state from NICER's mkf file
    @staticmethod
    def readMKF(fn_string):
        tab_mkf = Table.read(fn_string, hdu=1)
        r = np.array(tab_mkf['POSITION'])
        t = np.array(tab_mkf['TIME'])
        # v = np.array(tab_mkf['VELOCITY'])
        return r, t

    # time_array is the list of times that corresponds to r_array
    @staticmethod
    def interpolate_mkf_position(r_mkf, t_mkf, t_array):
        x_interpolator = interp1d(t_mkf, r_mkf[:, 0])
        y_interpolator = interp1d(t_mkf, r_mkf[:, 1])
        z_interpolator = interp1d(t_mkf, r_mkf[:, 2])

        rx = x_interpolator(t_array).reshape((len(t_array), 1))
        ry = y_interpolator(t_array).reshape((len(t_array), 1))
        rz = z_interpolator(t_array).reshape((len(t_array), 1))
        r_array = np.hstack((rx, ry, rz))

        return r_array
