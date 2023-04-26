# Author: Seamus Flannery
# Compare the precision vs altitude plots for two planets, at a variety of photon energies
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
from tcc_double_slide_minfit import CurveComparison
from HC_precision_vs_altitude import read_data, median_zero_and_outlier_remover

N = 5378  # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99]  # range of transmittance in which to compare the curves
# E_kev as a list for this file
hc_type = "rising"
cb_str = ""
cwd = os.getcwd()  # Get the current working directory (cwd)
wd = cwd.replace("HCNM_Analysis", "")
if wd[-1] != "/":
    wd += "/"
print(wd)


def plot_compare_planets_band_test(planet1, planet2, min_alt, max_alt, alt_interval, iterations, E_kev_list, write=False):
    global cb_str
    cb_str = planet1
    if write:
        write_data_band_test(min_alt, max_alt, alt_interval, iterations, E_kev_list)
    cb_str = planet2
    if write:
        write_data_band_test(min_alt, max_alt, alt_interval, iterations, E_kev_list)
    for E_kev in E_kev_list:
        suffix = "_int_" + str(alt_interval) + "_iter_" + str(iterations) + "_kev_" + str(E_kev) + ".npy"
        # wd = "/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project/"
        dt_name1 = wd + "sample_data/" + planet1 + "_dt" + suffix
        dr_name1 = wd + "sample_data/" + planet1 + "_dr" + suffix
        alt_name1 = wd + "sample_data/" + planet1 + "_alt" + suffix
        rand_altname1 = wd + "sample_data/" + planet1 + "_rand_alt" + suffix
        dt_list1 = read_data(dt_name1)
        dr_list1 = read_data(dr_name1)
        altitude_list1 = read_data(alt_name1)
        rand_alt_list1 = read_data(rand_altname1)
        while np.median(
                dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
            dt_list1 = np.delete(dt_list1, 0, axis=0)
            dr_list1 = np.delete(dr_list1, 0, axis=0)
            altitude_list1 = np.delete(altitude_list1, 0)
        dt_name2 = wd + "sample_data/" + planet2 + "_dt" + suffix
        dr_name2 = wd + "sample_data/" + planet2 + "_dr" + suffix
        alt_name2 = wd + "sample_data/" + planet2 + "_alt" + suffix
        rand_altname2 = wd + "sample_data/" + planet1 + "_rand_alt" + suffix
        dt_list2 = read_data(dt_name2)
        dr_list2 = read_data(dr_name2)
        altitude_list2 = read_data(alt_name2)
        rand_alt_list2 = read_data(rand_altname2)
        while np.median(
                dt_list1[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
            dt_list2 = np.delete(dt_list2, 0, axis=0)
            dr_list2 = np.delete(dr_list2, 0, axis=0)
            altitude_list2 = np.delete(altitude_list2, 0)
        # plot each dr and dt random plot for each planet for this KeV
        plt.subplot(221)
        plt.scatter(rand_alt_list1, dt_list1, label=fr"median " + str(E_kev) + " KeV")

        plt.subplot(222)
        plt.scatter(rand_alt_list2, dt_list2, label=fr"median " + str(E_kev)  + " KeV")

        plt.subplot(223)
        plt.scatter(rand_alt_list1, dr_list1, label=fr"median " + str(E_kev)  + " KeV")

        plt.subplot(224)
        plt.scatter(rand_alt_list2, dr_list2, label=fr"median " + str(E_kev)  + " KeV")
    # after all KeV are plotted, add titles and such
    plt.subplot(221)
    plt.title(r"$\delta t$ vs. altitude - " + planet1)
    plt.ylabel(fr"$\delta t$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.subplot(222)
    plt.title(r"$\delta t$ vs. altitude - " + planet2)
    plt.ylabel(fr"$\delta t$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.subplot(223)
    plt.title(r"$\delta r$ vs. altitude - " + planet1)
    plt.ylabel(r"$\delta r$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.subplot(224)
    plt.title(r"$\delta r$ vs. altitude - " + planet2)
    plt.ylabel(r"$\delta r$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.tight_layout()
    plt.show()


def write_data_band_test(min_alt, max_alt, alt_interval, how_many, E_kev_list):
    for E_kev in E_kev_list:
        print(wd)
        dt_path = wd + "sample_data/" + cb_str + "_dt_int_" + str(alt_interval) + "_iter_" \
                  + str(how_many) + "_kev_" + str(E_kev)
        print(dt_path)
        dr_path = wd + "sample_data/" + cb_str + "_dr_int_" + str(alt_interval) + "_iter_" + \
                  str(how_many) + "_kev_" + str(E_kev)
        alt_path = wd + "sample_data/" + cb_str + "_alt_int_" + str(alt_interval) + \
                   "_iter_" + str(how_many) + "_kev_" + str(E_kev)
        rand_alt_path = wd + "sample_data/" + cb_str + "_rand_alt_int_" + str(alt_interval) + \
                   "_iter_" + str(how_many) + "_kev_" + str(E_kev)
        print("save path:: " + dt_path)
        altitude_list = np.arange(min_alt, max_alt, alt_interval)
        dt_list = np.zeros((len(altitude_list), how_many), dtype="float")
        dr_list = np.zeros((len(altitude_list), how_many), dtype="float")
        rand_alt_list = np.zeros((len(altitude_list), how_many), dtype="float")
        fail_counter = 0
        for j in range(how_many):
            for i, alt in enumerate(altitude_list):
                rand_alt = np.random.choice(np.arange(alt, alt+alt_interval, 1))
                rand_alt_list[i,j]=rand_alt
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
                      ", " + str(E_kev) + " keV " + str(round((j * len(altitude_list) + i + 1) *
                                                              100 / (how_many * len(altitude_list)), 2)) + "% complete")
        total_runs = how_many * len(altitude_list)
        print("percent failure: " + str(fail_counter / total_runs * 100) + "%")
        np.save(rand_alt_path, rand_alt_list)
        np.save(dt_path, dt_list)
        np.save(dr_path, dr_list)
        np.save(alt_path, altitude_list)


plot_compare_planets_band_test("Earth", "Jupiter", 600, 10000, 200, 100, [1.0, 1.5, 2, 5, 8, 10], write=True)  # stable
# plot_compare_planets("Jupiter", "Earth", 600, 10000, 200, 300)  # TODO stabilize on Maria
print(wd)
