# Authors: Nathaniel Ruhl + Seamus Flannery
# This script uses Nate's toy model to determine the effects of orbital altitude on the precision
# of a horizon crossing, and takes different planets, and a variety of plotting parameters
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import datetime
import sys
import math
# set path for local modules
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")
# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide_double_gaus import CurveComparison, generate_crossings


# Global variables
N = 5378 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99] # range of transmittance in which to compare the curves
E_kev = 1.5 # keV
hc_type = "rising"
cb_str = "Earth" # planet being plotted
# allows for command-line arguments to determine planet
if len(sys.argv) == 2:
    cb_str = str(sys.argv[1])


def main():
    np.random.seed(3)
    write_data(10000, 20000, 101, 101)
    plot_read_data(cb_str, 101, 101)
    return 0


def test():
    # do_one()
    # do_a_bunch(400, 2100, 25, 100)
    # do_a_bunch_median(400, 2100, 25, 100)
    # plot_a_bunch(400,2100, 25, 100)
    # print(plot_inverse_root_fit([1, 2, 3, 4], [1, 0.7, 0.57, 0.5]))
    # do_a_bunch_max_min(400, 2100, 100, 1)
    # print(poly_fit([1, 2, 3], [1, 4, 9], 2))
    # write_data(300, 10000, 25, 1000)
    # array = [[0, 0, 3, 5, 7], [1, 2, 3, 4, 5], [0, 0, 3, 4, 6], [0, 0, 0, 0, 1]]
    # median_zero_remover(array)
    plot_read_data(cb_str, 25, 1000)
    return 0


def poly_fit(x, y, degree, printout=False):  #for polynomial fitting, returns fit data in array
    fit_params = np.polyfit(x, y, degree)
    if printout and degree == 3:
        print(str(fit_params[0]) + " x^3 + " + str(fit_params[1]) + " x^2 + "
              + str(fit_params[2]) + " x + " + str(fit_params[3]))
    fit = np.poly1d(fit_params)
    out_array = fit(x)
    return out_array


def sqrt_inverse(x, a):
    return a/np.sqrt(x)


def plot_inverse_root_fit(x, y, printout=False):
    fit_params = curve_fit(sqrt_inverse, x, y)
    if printout:
        print("a = " + str(fit_params[0][0]) + "\ncovariance: " + str(fit_params[1]))
    out_array = sqrt_inverse(x, fit_params[0][0])
    return out_array


def sqrt_func(x, a, k):
    return (a*np.sqrt(x))+k


def plot_root_fit(x, y, printout=False):
    fit_params = curve_fit(sqrt_func, x, y)
    if printout:
        print("a = " + str(fit_params[0][0]) + "\nk = " + str(fit_params[0][1]) + "\ncovariance: " + str(fit_params[1]))
    out_array = sqrt_func(x, fit_params[0][0], fit_params[0][1])
    return out_array


def exponential(x, a, b, k):
    return (a*(math.e**(b*x)))+k


def plot_exponential_fit(x, y, printout=False):
    fit_params = curve_fit(exponential, x, y, p0=[0.01, -200, 0])
    if printout:
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
        print("error: " +str(fit_params[1]))
    out_array = exponential(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
    return out_array


def inverse_log(x, a, b, k):
    return b/(np.log10(a*x))+k


def plot_inverse_log_fit(x, y, printout=False):
    fit_params = curve_fit(inverse_log, x, y, bounds=([0, -np.inf, -np.inf], np.inf))
    if printout:
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
    out_array = inverse_log(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
    return out_array


def curve_compaction_form(altitude, fit_factor):
    grav_const = 6.6743 * (10 ** (-20))  # kilometers cubed, inverse kilograms, inverse seconds
    if cb_str == "Earth":
        r_planet, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt = [6378, 5.972 * (10 ** 24), 85, 160]
    elif cb_str == "Mars":
        r_planet, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt = [3390, 6.39 * (10 ** 23), 100, 200]
    elif cb_str == "Venus":
        r_planet, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt = [6051, 4.867 * (10 ** 24), 150, 300]
    elif cb_str == "Moon":
        r_planet, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt = [1737, 7.35 * (10 ** 22), 0.0001, 0.00011]
    r_orbit = r_planet+altitude  # km
    b = zero_transmit_tang_alt  # km
    a = one_transmit_tang_alt  # km
    term1 = ((altitude+r_planet)**(3/2))
    term2 = np.arccos((r_planet+b) / r_orbit)-np.arccos((r_planet+a) / r_orbit)
    term3 = np.sqrt(grav_const*m_planet)
    duration = term1 * term2 / term3 * fit_factor
    return duration


def curve_comp_fit(x, y, printout=False):
    popt, pcov = curve_fit(curve_compaction_form, x, y)
    if printout:
        print("fit factor= " + str(popt[0]))
    out_array = curve_compaction_form(x, *popt)
    return out_array


def convert_to_log_scales(x, y):
    logx = np.zeros(len(x))
    logy = np.zeros(len(y))
    for i in range(len(x)):
        logx[i]=np.log10(x[i])
        logy[i]=np.log10(y[i])
    return [logx, logy]


def do_a_bunch_max_min(min_alt, max_alt, alt_interval, how_many):
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros((len(altitude_list),how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list),how_many), dtype="float")
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
            comp_obj = CurveComparison(sat, hc_type, N)
            dt_list[i][j] = comp_obj.dt_e
            dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega
            print("iteration: " + str(j) + "/" + str(how_many) + " altitude: " + str(alt) +
                  ", " + str(
                round((j * len(altitude_list) + i + 1) * 100 / (how_many * len(altitude_list)), 2)) + "% complete")
    plot_data(dt_list, dr_list, altitude_list)
    return 0


def read_data(pathname):
    return np.load(pathname)


def write_data(min_alt, max_alt, alt_interval, how_many):
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
    dt_path = "sample_data/" + cb_str + "_dt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    print(dt_path)
    dr_path = "sample_data/" + cb_str + "_dr_int_" + str(alt_interval) + "_iter_" + str(how_many)
    alt_path = "sample_data/" + cb_str + "_alt_int_" + str(alt_interval) + "_iter_" + str(how_many)
    np.save(dt_path, dt_list)
    np.save(dr_path, dr_list)
    np.save(alt_path, altitude_list)

def plot_read_data(planet, interval, iter):
    suffix = "_int_" + str(interval) + "_iter_" + str(iter) + ".npy"
    dt_name = "sample_data/" + planet + "_dt" + suffix
    dr_name = "sample_data/" + planet + "_dr" + suffix
    alt_name = "sample_data/" + planet + "_alt" + suffix
    dt_list = read_data(dt_name)
    dr_list = read_data(dr_name)
    altitude_list = read_data(alt_name)
    while np.median(dt_list[0]) == 0:  # handles situations where you generated data too low and need to edit out zeroes
        dt_list = np.delete(dt_list, 0, axis=0)
        dr_list = np.delete(dr_list, 0, axis=0)
        altitude_list = np.delete(altitude_list, 0)
    plot_data(dt_list, dr_list, altitude_list)


def median_zero_remover(sorted_array):  # removes zeroes from a 2-D array, then generates a list of the median of the remaining data
    out_array = np.zeros(len(sorted_array))
    for i in range(len(sorted_array)):
        row_zeroes = 0
        for j in range(len(sorted_array[1])):
            if sorted_array[i][j] == 0:
                row_zeroes += 1
        out_array[i] = np.median(sorted_array[i][row_zeroes:])
        if row_zeroes == len(sorted_array[1]):
            out_array[i] = 0
    return out_array


def plot_data(dt_list, dr_list, altitude_list, save=False):
    dt_sort = np.sort(dt_list)
    dr_sort = np.sort(dr_list)
    # dt_median = np.median(dt_sort, axis=1)
    # dr_median = np.median(dr_sort, axis=1)
    dt_median = median_zero_remover(dt_sort)
    dr_median = median_zero_remover(dr_sort)
    dt_66 = np.percentile(dt_list, 66, axis=1)
    dt_33 = np.percentile(dt_list, 33, axis=1)
    dr_66 = np.percentile(dr_list, 66, axis=1)
    dr_33 = np.percentile(dr_list, 33, axis=1)
    # dt fits
    # print("dt median inverse log fit: ")
    # dt_median_invlog_fit = plot_inverse_log_fit(altitude_list, dt_median, True)
    # print("dt 66 inverse log fit: ")
    # dt_66_invlog_fit = plot_inverse_log_fit(altitude_list, dt_66, True)
    # plt.plot(altitude_list, dt_66_invlog_fit)
    # print("dt 33 inverse log fit: ")
    # dt_33_invlog_fit = plot_inverse_log_fit(altitude_list, dt_33, True)
    dt_med_comp_fit = curve_comp_fit(altitude_list, dt_median, True)
    # dr fits
    print("dr median inverse log fit: ")
    dr_median_invlog_fit = plot_inverse_log_fit(altitude_list, dr_median, True)
    print("dr 66 inverse log fit: ")
    dr_66_invlog_fit = plot_inverse_log_fit(altitude_list, dr_66, True)
    print("dr 33 inverse log fit: ")
    # dr_33_invlog_fit = plot_inverse_log_fit(altitude_list, dr_33, True)
    # log scaling plots/fits
    dt_log_data = convert_to_log_scales(altitude_list, dt_median)
    fit_y_dt = poly_fit(dt_log_data[0], dt_log_data[1], 1)
    dr_log_data = convert_to_log_scales(altitude_list, dr_median)
    fit_y_dr = poly_fit(dr_log_data[0], dr_log_data[1], 1)

    plt.figure(1)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    # plt.fill_between(altitude_list, dt_66_invlog_fit, dt_33_invlog_fit)
    plt.plot(altitude_list, dt_median, label=fr"median", color="orange")
    plt.plot(altitude_list, dt_med_comp_fit, label=fr"median invlog fit", color="red")
    plt.ylabel(
        fr"Temporal uncertaintainty in HCNM measurement, $\delta t_e$ (sec), {E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    # plt.xscale("log")
    # plt.yscale("log")
    if save:
        plt.savefig("plots/dt_v_alt_" + str(datetime.datetime.now()) + ".png")

    plt.figure(2)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    # plt.fill_between(altitude_list, dr_33_invlog_fit, dr_66_invlog_fit)
    plt.plot(altitude_list, dr_median, label=fr"median", color="orange")
    plt.plot(altitude_list, dr_median_invlog_fit, label=fr"median invlog fit", color="red")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    # plt.xscale("log")
    # plt.yscale("log")
    if save:
        plt.savefig("plots/dr_v_alt_" + str(datetime.datetime.now()) + ".png")

    plt.figure(3)
    plt.title("log log dt")
    plt.plot(dt_log_data[0], dt_log_data[1])
    plt.plot(dt_log_data[0], fit_y_dt)

    plt.figure(4)
    plt.title("log log dr")
    plt.plot(dr_log_data[0], dr_log_data[1])
    plt.plot(dr_log_data[0], fit_y_dr)

    plt.show()


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
