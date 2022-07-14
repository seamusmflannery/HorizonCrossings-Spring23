# Author: Nathaniel Ruhl
# This script uses the toy model to determine the effects of orbital altitude on the precision of a horizon crossing
from scipy.optimize import curve_fit
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import datetime
from math import e
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")

# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide import CurveComparison, generate_crossings

# Global variables
N = 5378 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99] # range of transmittance in which to compare the curves
cb_str = "Mars"
E_kev = 1.5 # keV
hc_type = "rising"


def main():
    write_data(500, 10000, 25, 1000)
    plot_read_data(cb_str, 25, 1000)
    #write_data(200, 10000, 25, 1000)
    #test()
    return 0


def test():
    # do_one()
    # do_a_bunch(400, 2100, 25, 100)
    # do_a_bunch_median(400, 2100, 25, 100)
    # plot_a_bunch(400,2100, 25, 100)
    # print(plot_inverse_root_fit([1, 2, 3, 4], [1, 0.7, 0.57, 0.5]))
    # do_a_bunch_max_min(400, 2100, 100, 1)
    # print(poly_fit([1, 2, 3], [1, 4, 9], 2))
    write_data(500, 10000, 1000, 5)
    plot_read_data(1000, 5)
    return 0


def poly_fit(x, y, degree, printout="false"):  #for polynomial fitting, returns fit data in array
    fit_params = np.polyfit(x, y, degree)
    if printout == "true" and degree == 3:
        print(str(fit_params[0]) + " x^3 + " + str(fit_params[1]) + " x^2 + "
              + str(fit_params[2]) + " x + " + str(fit_params[3]))
    fit = np.poly1d(fit_params)
    out_array = fit(x)
    return out_array


def sqrt_inverse(x, a):
    return (a/np.sqrt(x))

def plot_inverse_root_fit(x, y, printout="false"):
    fit_params = curve_fit(sqrt_inverse, x, y)
    if printout == "true":
        print("a = " + str(fit_params[0][0]) + "\ncovariance: " + str(fit_params[1]))
    out_array = sqrt_inverse(x, fit_params[0][0])
    return out_array


def sqrt_func(x, a, k):
    return (a*np.sqrt(x))+k

def plot_root_fit(x, y, printout="false"):
    fit_params = curve_fit(sqrt_func, x, y)
    if printout == "true":
        print("a = " + str(fit_params[0][0]) + "\nk = " + str(fit_params[0][1]) + "\ncovariance: " + str(fit_params[1]))
    out_array = sqrt_func(x, fit_params[0][0], fit_params[0][1])
    return out_array


def exponential(x, a, b, k):
    return (a*(e**(b*x)))+k

def plot_exponential_fit(x, y, printout="false"):
    fit_params = curve_fit(exponential, x, y, p0=[0.01, -200, 0])
    if printout == "true":
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
        print("error: " +str(fit_params[1]))
    out_array = exponential(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
    return out_array


def inverse_log(x, a, b, k):
    return b/(np.log10(a*x))+k

def plot_inverse_log_fit(x, y, printout="false"):
    fit_params = curve_fit(inverse_log, x, y, bounds=([0, -np.inf, -np.inf], np.inf))
    if printout == "true":
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
    out_array = inverse_log(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
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


def read_data(pathname): #TODO finish this?
    return np.load(pathname)


def write_data(min_alt, max_alt, alt_interval, how_many):
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
                  ", " + str(round((j*len(altitude_list)+i+1)*100/(how_many*len(altitude_list)), 2)) + "% complete")
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
    plot_data(dt_list, dr_list, altitude_list)


def plot_data(dt_list, dr_list, altitude_list, save="false"):
    dt_sort = np.sort(dt_list)
    dr_sort = np.sort(dr_list)
    dt_median = np.median(dt_sort, axis=1)
    dr_median = np.median(dr_sort, axis=1)
    dt_66 = np.percentile(dt_list, 66, axis=1)
    dt_33 = np.percentile(dt_list, 33, axis=1)
    dr_66 = np.percentile(dr_list, 66, axis=1)
    dr_33 = np.percentile(dr_list, 33, axis=1)
    # dt fits
    print("dt median inverse log fit: ")
    dt_median_invlog_fit = plot_inverse_log_fit(altitude_list, dt_median, "true")
    print("dt 66 inverse log fit: ")
    dt_66_invlog_fit = plot_inverse_log_fit(altitude_list, dt_66, "true")
    plt.plot(altitude_list, dt_66_invlog_fit)
    print("dt 33 inverse log fit: ")
    dt_33_invlog_fit = plot_inverse_log_fit(altitude_list, dt_33, "true")
    # dr fits
    print("dr median inverse log fit: ")
    dr_median_invlog_fit = plot_inverse_log_fit(altitude_list, dr_median, "true")
    print("dr 66 inverse log fit: ")
    dr_66_invlog_fit = plot_inverse_log_fit(altitude_list, dr_66, "true")
    print("dr 33 inverse log fit: ")
    dr_33_invlog_fit = plot_inverse_log_fit(altitude_list, dr_33, "true")

    dt_log_data = convert_to_log_scales(altitude_list, dt_median)
    fit_y_dt = poly_fit(dt_log_data[0], dt_log_data[1], 1)
    dr_log_data = convert_to_log_scales(altitude_list, dr_median)
    fit_y_dr = poly_fit(dr_log_data[0], dr_log_data[1], 1)

    plt.figure(1)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.fill_between(altitude_list, dt_66_invlog_fit, dt_33_invlog_fit)
    plt.plot(altitude_list, dt_median, label=fr"median", color="orange")
    plt.plot(altitude_list, dt_median_invlog_fit, label=fr"median invlog fit", color="red")
    plt.ylabel(
        fr"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec), {E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    if save == "true":
        plt.savefig("plots/dt_v_alt_" + str(datetime.datetime.now()) + ".png")
    # plt.xscale("log")
    # plt.yscale("log")

    plt.figure(2)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.fill_between(altitude_list, dr_33_invlog_fit, dr_66_invlog_fit)
    plt.plot(altitude_list, dr_median, label=fr"median", color="orange")
    plt.plot(altitude_list, dr_median_invlog_fit, label=fr"median invlog fit", color="red")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    if save == "true":
        plt.savefig("plots/dr_v_alt_" + str(datetime.datetime.now()) + ".png")
    plt.legend()
    # plt.xscale("log")
    # plt.yscale("log")

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
