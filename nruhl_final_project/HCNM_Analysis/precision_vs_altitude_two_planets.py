# Author: Seamus Flannery
# Compare the precision vs altitude plots for two planets
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
# set path for local modules
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")
# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide_double_gaus import CurveComparison, generate_crossings
from HC_precision_vs_altitude import read_data, median_zero_remover, curve_comp_fit, plot_inverse_log_fit, write_data
N = 5378 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99]  # range of transmittance in which to compare the curves
E_kev = 1.5 # keV
hc_type = "rising"
cb_str = ""
cwd = os.getcwd()  # Get the current working directory (cwd)
wd = cwd.replace("HCNM_Analysis", "")
print(wd)


def plot_compare_planets(planet1, planet2, min_alt, max_alt, alt_interval, iterations, read=False):
    global cb_str
    cb_str = planet1
    if not read:
        write_data(min_alt, max_alt, alt_interval, iterations)
    cb_str = planet2
    if not read:
        write_data(min_alt, max_alt, alt_interval, iterations)
    suffix = "_int_" + str(alt_interval) + "_iter_" + str(iterations) + ".npy"
    # wd = "/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project/"
    dt_name1 = wd + "sample_data/" + planet1 + "_dt" + suffix
    dr_name1 = wd + "sample_data/" + planet1 + "_dr" + suffix
    alt_name1 = wd + "sample_data/" + planet1 + "_alt" + suffix
    dt_list1 = read_data(dt_name1)
    dr_list1 = read_data(dr_name1)
    altitude_list1 = read_data(alt_name1)
    while np.median(
            dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list1 = np.delete(dt_list1, 0, axis=0)
        dr_list1 = np.delete(dr_list1, 0, axis=0)
        altitude_list1 = np.delete(altitude_list1, 0)
    dt_name2 = wd + "sample_data/" + planet2 + "_dt" + suffix
    dr_name2 = wd + "sample_data/" + planet2 + "_dr" + suffix
    alt_name2 = wd + "sample_data/" + planet2 + "_alt" + suffix
    dt_list2 = read_data(dt_name2)
    dr_list2 = read_data(dr_name2)
    altitude_list2 = read_data(alt_name2)
    while np.median(
            dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list2 = np.delete(dt_list2, 0, axis=0)
        dr_list2 = np.delete(dr_list2, 0, axis=0)
        altitude_list2 = np.delete(altitude_list2, 0)
    dt_sort1 = np.sort(dt_list1)
    dr_sort1 = np.sort(dr_list1)
    dt_sort2 = np.sort(dt_list2)
    dr_sort2 = np.sort(dr_list2)
    dt_median1 = median_zero_remover(dt_sort1)
    dr_median1 = median_zero_remover(dr_sort1)
    dt_median2 = median_zero_remover(dt_sort2)
    dr_median2 = median_zero_remover(dr_sort2)
    dt_med_comp_fit1 = curve_comp_fit(altitude_list1, dt_median1, True)  # TODO clarify whats happening with this fit
    dt_med_comp_fit2 = curve_comp_fit(altitude_list2, dt_median2, True)  # TODO clarify whats happening with this fit
    dr_median_invlog_fit1 = curve_comp_fit(altitude_list1, dr_median1, True)  # TODO clarify whats happening with this fit
    dr_median_invlog_fit2 = curve_comp_fit(altitude_list2, dr_median2, True)  # TODO clarify whats happening with this fit

    plt.subplot(221)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude" + planet1)
    plt.plot(altitude_list1, dt_median1, label=fr"median", color="orange")
    plt.plot(altitude_list1, dt_med_comp_fit1, label=fr"median invlog fit", color="red")
    plt.ylabel(
        fr"Temporal uncertaintainty in HCNM measurement, $\delta t_e$ (sec), {E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.subplot(222)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude" + planet2)
    plt.plot(altitude_list2, dt_median2, label=fr"median", color="orange")
    plt.plot(altitude_list2, dt_med_comp_fit2, label=fr"median invlog fit", color="red")
    plt.ylabel(
        fr"Temporal uncertaintainty in HCNM measurement, $\delta t_e$ (sec), {E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.subplot(223)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude" + planet1)
    plt.plot(altitude_list1, dr_median1, label=fr"median", color="orange")
    plt.plot(altitude_list1, dr_median_invlog_fit1, label=fr"median invlog fit", color="red")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.subplot(224)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude" + planet2)
    plt.plot(altitude_list2, dr_median2, label=fr"median", color="orange")
    plt.plot(altitude_list2, dr_median_invlog_fit2, label=fr"median invlog fit", color="red")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.show()


def write_data(min_alt, max_alt, alt_interval, how_many):
    print(wd)
    dt_path = wd + "sample_data/" + cb_str + "_dt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print(dt_path)
    dr_path = wd + "sample_data/" + cb_str + "_dr_int_" + str(alt_interval) + "_iter_" + str(how_many)
    alt_path = wd + "sample_data/" + cb_str + "_alt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print("save path:: " + dt_path)
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros((len(altitude_list),how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list),how_many), dtype="float")
    fail_counter = 0
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
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
                  ", " + str(round((j*len(altitude_list)+i+1)*100/(how_many*len(altitude_list)), 2)) + "% complete")
    total_runs = how_many*len(altitude_list)
    print("percent failure: " + str(fail_counter/total_runs*100) + "%")
    np.save(dt_path, dt_list)
    np.save(dr_path, dr_list)
    np.save(alt_path, altitude_list)


# plot_compare_planets("Earth", "Jupiter", 600, 10000, 100, 100, read=True)  # stable
plot_compare_planets("Jupiter", "Earth", 600, 10000, 200, 300)  # TODO stabilize on Maria
print(wd)