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
from tcc_double_slide_minfit import CurveComparison
from HC_precision_vs_altitude import read_data, median_zero_and_outlier_remover, write_data

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
    rand_alt_name1 = wd + "sample_data/" + planet1 + "_randalt" + suffix
    print(rand_alt_name1)
    dt_list1 = read_data(dt_name1)
    dr_list1 = read_data(dr_name1)
    dt_scatter1 = np.copy(dt_list1)
    dr_scatter1 = np.copy(dr_list1)
    altitude_list1 = read_data(alt_name1)
    rand_alt_list1 = read_data(rand_alt_name1)
    while np.median(
            dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list1 = np.delete(dt_list1, 0, axis=0)
        dr_list1 = np.delete(dr_list1, 0, axis=0)
        altitude_list1 = np.delete(altitude_list1, 0)
    dt_name2 = wd + "sample_data/" + planet2 + "_dt" + suffix
    dr_name2 = wd + "sample_data/" + planet2 + "_dr" + suffix
    alt_name2 = wd + "sample_data/" + planet2 + "_alt" + suffix
    rand_alt_name2 = wd + "sample_data/" + planet1 + "_randalt" + suffix
    dt_list2 = read_data(dt_name2)
    dr_list2 = read_data(dr_name2)
    dt_scatter2 = np.copy(dt_list2)
    dr_scatter2 = np.copy(dr_list2)
    altitude_list2 = read_data(alt_name2)
    rand_alt_list2 = read_data(rand_alt_name2)
    while np.median(
            dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list2 = np.delete(dt_list2, 0, axis=0)
        dr_list2 = np.delete(dr_list2, 0, axis=0)
        altitude_list2 = np.delete(altitude_list2, 0)
    dt_sort1 = np.sort(dt_list1)
    dr_sort1 = np.sort(dr_list1)
    dt_sort2 = np.sort(dt_list2)
    dr_sort2 = np.sort(dr_list2)
    dt_median1 = median_zero_and_outlier_remover(dt_sort1)
    dr_median1 = median_zero_and_outlier_remover(dr_sort1)
    dt_median2 = median_zero_and_outlier_remover(dt_sort2)
    dr_median2 = median_zero_and_outlier_remover(dr_sort2)

    plt.subplot(221)
    plt.title(r"$\delta t$ vs. altitude - " + planet1)
    plt.plot(altitude_list1, dt_median1, label=fr"median", color="blue")
    plt.scatter(rand_alt_list1, dt_scatter1, c="green")
    plt.ylabel(
        fr"$\delta t$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.xlim(min_alt, max_alt - alt_interval)
    plt.legend()

    plt.subplot(222)
    plt.title(r"$\delta t$ vs. altitude - " + planet2)
    plt.plot(altitude_list2, dt_median2, label=fr"median", color="red")
    plt.scatter(rand_alt_list2, dt_scatter2, c="orange")
    plt.ylabel(
        fr"$\delta t$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.xlim(min_alt, max_alt - alt_interval)
    plt.legend()

    plt.subplot(223)
    plt.title(r"$\delta r$ vs. altitude - " + planet1)
    plt.plot(altitude_list1, dr_median1, label=fr"median", color="blue")
    plt.scatter(rand_alt_list1, dr_scatter1, c="green")
    plt.ylabel(r"$\delta r$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.xlim(min_alt, max_alt-alt_interval)
    plt.legend()

    plt.subplot(224)
    plt.title(r"$\delta r$ vs. altitude - " + planet2)
    plt.plot(altitude_list2, dr_median2, label=fr"median", color="red")
    plt.scatter(rand_alt_list2, dr_scatter2, c="orange")
    plt.ylabel(r"$\delta r$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.xlim(min_alt, max_alt - alt_interval)
    plt.legend()

    plt.tight_layout()
    plt.show()


def write_data(min_alt, max_alt, alt_interval, how_many):
    print(wd)
    dt_path = wd + "sample_data/" + cb_str + "_dt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print(dt_path)
    rand_alt_path = wd + "sample_data/" + cb_str + "_randalt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    dr_path = wd + "sample_data/" + cb_str + "_dr_int_" + str(alt_interval) + "_iter_" + str(how_many)
    alt_path = wd + "sample_data/" + cb_str + "_alt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print("save path:: " + dt_path)
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    rand_alt_list = np.zeros((len(altitude_list), how_many), dtype="float")
    dt_list = np.zeros((len(altitude_list), how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list), how_many), dtype="float")
    fail_counter = 0
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            rand_alt = np.random.choice(np.arange(alt, alt + alt_interval, 1))
            rand_alt_list[i][j] = rand_alt
            sat = AnalyzeCrossing(cb=cb_str, H=rand_alt, E_kev=E_kev)
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
    np.save(rand_alt_path, rand_alt_list)
    np.save(dt_path, dt_list)
    np.save(dr_path, dr_list)
    np.save(alt_path, altitude_list)


plot_compare_planets("Earth", "Jupiter", 600, 10200, 200, 300, read=False)  # stable
