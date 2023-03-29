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
from tcc_double_slide_fitless import CurveComparison
from HC_precision_vs_altitude import read_data, median_zero_and_outlier_remover

N = 5378  # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99]  # range of transmittance in which to compare the curves
E_kev = 1.5  # keV
hc_type = "rising"
cb_str = ""
cwd = os.getcwd()  # Get the current working directory (cwd)
wd = cwd.replace("HCNM_Analysis", "")
if wd[-1] != "/":
    wd += "/"
print(wd)

def random_iterate(list):
    shuffled_list = list
    random.shuffle(shuffled_list)
    return shuffled_list


def variable_write(cb_str_list, min_alt, max_alt, alt_interval, how_many):
    global cb_str
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    t_e_list = np.zeros((len(altitude_list), how_many), dtype="float")
    fail_counter = 0
    total_runs = how_many * len(altitude_list) * len(cb_str_list)
    current_runs = 0
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            rand_list = random_iterate(np.copy(cb_str_list))  # rotates through the planet ephem files, once each, random order
            for k, cb in enumerate(rand_list):
                cb_str = rand_list[k]
                rand_alt = np.random.choice(np.arange(alt, alt+alt_interval, 1))
                sat = AnalyzeCrossing(cb=cb_str, H=rand_alt, E_kev=E_kev)
                current_runs += 1
                try:
                    comp_obj = CurveComparison(sat, hc_type, N)
                    t_e_list[i][j] = comp_obj.t0_e
                except RuntimeError:
                    print(str(sat.H) + " failed to fit")
                    fail_counter += 1
                    print("fail counter: " + str(fail_counter))
                except ValueError:
                    print(str(sat.H) + " failed on a value error")
                    fail_counter += 1
                    print("fail counter: " + str(fail_counter))
                print("iteration: " + str(j) + "/" + str(how_many) + " altitude: " + str(rand_alt) +
                      ", " + str(
                        round(current_runs * 100 / total_runs, 2)) + "% complete")
    print("percent failure: " + str(fail_counter / total_runs * 100) + "%")
    t_e_path = wd + "sample_data/variable_" + cb_str_list[0] + "_t_e_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print(t_e_path)
    alt_path = wd + "sample_data/variable_" + cb_str_list[0] + "_alt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    np.save(t_e_path, t_e_list)
    np.save(alt_path, altitude_list)
    plt.show()


def plot_variable_data(planet, interval, iter):
    global cb_str_list
    suffix = "_int_" + str(interval) + "_iter_" + str(iter) + ".npy"
    t_e_name = wd + "sample_data/variable_" + planet + "_t_e" + suffix
    alt_name = wd + "sample_data/variable_" + planet + "_alt" + suffix
    t_e_list = np.subtract(read_data(t_e_name), 2000)
    altitude_list = read_data(alt_name)
    t_e_sort = abs(np.sort(t_e_list))

    t_e_median = median_zero_and_outlier_remover(t_e_sort)

    plt.title(r"$\delta t$ uncertainty as a function of orbital altitude, " + str(len(cb_str_list)) + " atmospheric composition(s)")
    # plt.fill_between(altitude_list, dt_66_invlog_fit, dt_33_invlog_fit)
    plt.plot(altitude_list, t_e_median, label=fr"median, " + str(iter*4) + " iterations", color="red")
    # plt.plot(altitude_list, dt_med_comp_fit, label=fr"median invlog fit", color="red")
    plt.ylabel(
        fr"$\delta t$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.show()



cb_str_list = ["Jupiter", "Jupiter1", "Jupiter2", "Jupiter3"]
# variable_write(cb_str_list, 600, 10000, 200, 11)
plot_variable_data("Jupiter", 200, 11)
cb_str_list = ["Jupiter"]
# variable_write(cb_str_list, 600, 10000, 200, 400)
# plot_variable_data("Jupiter", 200, 400)
# variable_write(cb_str_list, 600, 10000, 200, 12)
# plot_variable_data("Jupiter", 200, 12)
