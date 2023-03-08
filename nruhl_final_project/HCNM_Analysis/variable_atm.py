# Author: Seamus Flannery
# Compare the precision vs altitude plots for two planets
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import random
# set path for local modules
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")
# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide_double_gaus import CurveComparison, generate_crossings
from HC_precision_vs_altitude import read_data, median_zero_remover, curve_comp_fit, plot_inverse_log_fit, write_data, \
    plot_data

N = 5378 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99]  # range of transmittance in which to compare the curves
E_kev = 1.5 # keV
hc_type = "rising"
cb_str = ""
cwd = os.getcwd()  # Get the current working directory (cwd)
wd = cwd.replace("HCNM_Analysis", "")
print(wd)


def variable_write(cb_str_list, min_alt, max_alt, alt_interval, how_many):
    global cb_str
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros((len(altitude_list), how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list), how_many), dtype="float")
    fail_counter = 0
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            cb_str = cb_str_list[random.randrange(len(cb_str_list))]  # rotates through the planet ephem files
            sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
            try:
                comp_obj = CurveComparison(sat, hc_type, N)
                dt_list[i][j] = comp_obj.dt_e
                dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega
            except RuntimeError:
                print(str(sat.H) + " failed to fit")
                fail_counter += 1
                print("fail counter: " + str(fail_counter))
            except ValueError:
                print(str(sat.H) + " failed on a value error")
                fail_counter += 1
                print("fail counter: " + str(fail_counter))
            print("iteration: " + str(j) + "/" + str(how_many) + " altitude: " + str(alt) +
                  ", " + str(
                round((j * len(altitude_list) + i + 1) * 100 / (how_many * len(altitude_list)), 2)) + "% complete")
    total_runs = how_many * len(altitude_list)
    print("percent failure: " + str(fail_counter / total_runs * 100) + "%")
    dt_path = "sample_data/variable_" + cb_str_list[0] + "_dt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print(dt_path)
    dr_path = "sample_data/variable_" + cb_str_list[0] + "_dr_int_" + str(alt_interval) + "_iter_" + str(how_many)
    alt_path = "sample_data/variable_" + cb_str_list[0] + "_alt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    np.save(dt_path, dt_list)
    np.save(dr_path, dr_list)
    np.save(alt_path, altitude_list)


def plot_variable_data(planet, interval, iter):
    suffix = "_int_" + str(interval) + "_iter_" + str(iter) + ".npy"
    dt_name = "sample_data/variable_" + planet + "_dt" + suffix
    dr_name = "sample_data/variable_" + planet + "_dr" + suffix
    alt_name = "sample_data/variable_" + planet + "_alt" + suffix
    dt_list = read_data(dt_name)
    dr_list = read_data(dr_name)
    altitude_list = read_data(alt_name)
    while np.median(dt_list[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list = np.delete(dt_list, 0, axis=0)
        dr_list = np.delete(dr_list, 0, axis=0)
        altitude_list = np.delete(altitude_list, 0)
    plot_data(dt_list, dr_list, altitude_list)


cb_str_list = ["Jupiter", "Jupiter1", "Jupiter2", "Jupiter3"]
variable_write(cb_str_list, 300, 1000, 10, 10)
plot_variable_data("Jupiter", 10, 10)
