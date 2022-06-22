# Author: Nathaniel Ruhl
# This script uses the toy model to determine the effects of orbital altitude on the precision of a horizon crossing
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")

# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide import CurveComparison, generate_crossings

# Global variables
N = 244 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99] # range of transmittance in which to compare the curves
cb_str = "Earth"
E_kev = 1.5 # keV
hc_type = "rising"


def main():
    test()
    # do_a_bunch_max_min(400, 2500, 25, 1000)
    return 0


def test():
    # do_one()
    # do_a_bunch(400, 2100, 25, 100)
    # do_a_bunch_median(400, 2100, 25, 100)
    # plot_a_bunch(400,2100, 25, 100)
    # print(plot_inverse_root_fit([1, 2, 3, 4], [1, 0.7, 0.57, 0.5]))
    # do_a_bunch_max_min(400, 2100, 100, 1)
    # print(plot_fit([1, 2, 3], [1, 4, 9], 2))
    write_data(400, 2100, 25, 10)
    return 0


def plot_fit(x, y, degree, printout="false"):  #for polynomial fitting, returns fit data in array
    fit_params = np.polyfit(x, y, degree)
    if printout == "true":
        print(str(fit_params[0]) + " x^3 + " + str(fit_params[1]) + " x^2 + "
              + str(fit_params[2]) + " x + " + str(fit_params[3]))
    fit = np.poly1d(fit_params)
    out_array = fit(x)
    return out_array


def sqrt_inverse(x, a):
    return a/np.sqrt(x)


def plot_inverse_root_fit(x, y, printout="false"):
    fit_params = curve_fit(sqrt_inverse, x, y)
    if printout == "true":
        print("a= " + str(fit_params[0]) + " covariance: " + str(fit_params[1]))
    out_array = sqrt_inverse(x, fit_params[0])
    return out_array


def exponential(x, a, b, k):
    return (a*(b**x))+k


def plot_exponential_fit(x, y, printout="false"):
    fit_params = curve_fit(exponential, x, y)
    if printout == "true":
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
    out_array = exponential(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
    return out_array


def hyperbolic(x, a, b, k):
    return b/((a*x)+k)


def plot_hyperbolic_fit(x, y, printout="false"):
    fit_params = curve_fit(hyperbolic, x, y)
    if printout == "true":
        print("a= " + str(fit_params[0][0]) + "\nb= " + str(fit_params[0][1]) + "\nk= " +str(fit_params[0][2]))
    out_array = hyperbolic(x, fit_params[0][0], fit_params[0][1], fit_params[0][2])
    return out_array


def do_one():
    altitude_list = np.arange(400, 2100, 100)
    dt_list = np.zeros_like(altitude_list, float)
    dr_list = np.zeros_like(dt_list)  # lists containing uncertainties corresponding to altitude_list

    for i, alt in enumerate(altitude_list):
        sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
        comp_obj = CurveComparison(sat, hc_type, N)
        dt_list[i] = comp_obj.dt_e
        dr_list[i] = comp_obj.dt_e * sat.R_orbit * sat.omega

    # Plot results
    plt.figure(1)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_list, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()

    plt.figure(2)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_list,
             label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.show()


def do_a_bunch(min_alt, max_alt, alt_interval, how_many):
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros_like(altitude_list, float)
    dr_list = np.zeros_like(dt_list)
    for j in range(how_many):
        #np.random.seed(j)
        for i, alt in enumerate(altitude_list):
            sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
            comp_obj = CurveComparison(sat, hc_type, N)
            dt_list[i] += comp_obj.dt_e / how_many
            dr_list[i] += comp_obj.dt_e * sat.R_orbit * sat.omega / how_many
            print("seed: " + str(j) + " altitude: " + str(alt))
    plt.figure(1)
    # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_list, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.savefig("bunches/dt_v_alt.png")

    plt.figure(2)
    # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_list,
             label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.savefig("bunches/dr_v_alt.png")
    plt.show()

    dt = open("bunches/dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dt.write(str(dt_list))
    dt.close()
    dr = open("bunches/dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dr.write(str(dr_list))
    dr.close()
    altlist = open("bunches/alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    altlist.write(str(altitude_list))
    altlist.close()
    return 0


def do_a_bunch_median(min_alt, max_alt, alt_interval, how_many):
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros((len(altitude_list),how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list),how_many), dtype="float")
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
            comp_obj = CurveComparison(sat, hc_type, N)
            dt_list[i][j] = comp_obj.dt_e / how_many
            dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega / how_many
            print("iteration: " + str(j) + " altitude: " + str(alt))

    dt_median = np.median(dt_list, axis=1)
    dr_median = np.median(dr_list, axis=1)
    plt.figure(1)
    # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_median, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.savefig("bunches/med_dt_v_alt.png")

    plt.figure(2)
    # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_median, label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.savefig("bunches/med_dr_v_alt.png")
    plt.show()

    dt_data = open("bunches/med_dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dt_data.write(str(dt_list))
    dt_data.close()
    dr_data = open("bunches/med_dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dr_data.write(str(dr_list))
    dr_data.close()
    altlist = open("bunches/med_alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    altlist.write(str(altitude_list))
    altlist.close()
    return 0


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
            print("iteration: " + str(j) + " altitude: " + str(alt))
    dt_sort = np.sort(dt_list)
    dr_sort = np.sort(dr_list)
    dt_median = np.median(dt_sort, axis=1)
    dr_median = np.median(dr_sort, axis=1)
    dt_66 = np.percentile(dt_list, 66, axis=1)
    dt_33 = np.percentile(dt_list, 33, axis=1)
    dr_66 = np.percentile(dr_list, 66, axis=1)
    dr_33 = np.percentile(dr_list, 33, axis=1)
    print("dt median poly fit: ")
    dt_medianfit = plot_fit(altitude_list, dt_median, 3, "true")
    dt_33fit = plot_fit(altitude_list, dt_33, 3)
    dt_66fit = plot_fit(altitude_list, dt_66, 3)
    print("dr median poly fit: ")
    dr_medianfit = plot_fit(altitude_list, dr_median, 3, "true")
    dr_33fit = plot_fit(altitude_list, dr_33, 3)
    dr_66fit = plot_fit(altitude_list, dr_66, 3)
    print("dt median sq-1 fit: ")
    dt_mediansqfit = plot_inverse_root_fit(altitude_list, dt_median, "true")
    dt_33sqfit = plot_inverse_root_fit(altitude_list, dt_33)
    dt_66sqfit = plot_inverse_root_fit(altitude_list, dt_66)
    print("dr median sq-1 fit: ")
    dr_mediansqfit = plot_inverse_root_fit(altitude_list, dr_median, "true")
    dr_33sqfit = plot_inverse_root_fit(altitude_list, dr_33)
    dr_66sqfit = plot_inverse_root_fit(altitude_list, dr_66)
    print("dt median hyperbolic fit: ")
    dt_median_hyper_fit = plot_hyperbolic_fit(altitude_list, dt_median, "true")
    #dt_median_exponential_fit = plot_exponential_fit(altitude_list, dt_median, "true")

    plt.figure(1)
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_66, label=fr"upper 66 data")
    plt.plot(altitude_list, dt_33, label=fr"lower 33 data")
    plt.plot(altitude_list, dt_median, label=fr"median data")
    plt.plot(altitude_list, dt_66fit, label=fr"upper 66 poly fit")
    plt.plot(altitude_list, dt_33fit, label=fr"lower 33 poly fit")
    plt.plot(altitude_list, dt_medianfit, label=fr"median poly fit")
    plt.plot(altitude_list, dt_66sqfit, label=fr"upper 66  sq-1 fit")
    plt.plot(altitude_list, dt_33sqfit, label=fr"lower 33 sq-1 fit")
    plt.plot(altitude_list, dt_mediansqfit, label=fr"median sq-1 fit")
    plt.plot(altitude_list, dt_median_hyper_fit, label=fr"median hyper fit")
    #plt.plot(altitude_list, dt_median_exponential_fit, label=fr"median exponential fit")
    plt.ylabel(fr"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec), {E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    # plt.savefig("bunches/curves_dt_v_alt.png")

    plt.figure(2)
    plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_66, label=fr"upper 66 data")
    plt.plot(altitude_list, dr_33, label=fr"lower 33 data")
    plt.plot(altitude_list, dr_median, label=fr"median data")
    plt.plot(altitude_list, dr_66fit, label=fr"upper 66 poly fit")
    plt.plot(altitude_list, dr_33fit, label=fr"lower 33 poly fit")
    plt.plot(altitude_list, dr_medianfit, label=fr"median poly fit")
    plt.plot(altitude_list, dr_66sqfit, label=fr"upper 66  sq-1 fit")
    plt.plot(altitude_list, dr_33sqfit, label=fr"lower 33 sq-1 fit")
    plt.plot(altitude_list, dr_mediansqfit, label=fr"median sq-1 fit")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    # plt.savefig("bunches/curves_dr_v_alt.png")
    plt.legend()
    plt.show()

    dt_data = open("bunches/curves_dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dt_data.write(str(dt_list))
    dt_data.close()
    dr_data = open("bunches/curves_dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dr_data.write(str(dr_list))
    dr_data.close()
    altlist = open("bunches/curves_alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    altlist.write(str(altitude_list))
    altlist.close()
    return 0


def plot_a_bunch(min_alt, max_alt, alt_interval, how_many):
    altitude_list = np.arange(min_alt, max_alt, alt_interval)
    dt_list = np.zeros((len(altitude_list), how_many), dtype="float")
    dr_list = np.zeros((len(altitude_list), how_many), dtype="float")
    for j in range(how_many):
        for i, alt in enumerate(altitude_list):
            sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
            comp_obj = CurveComparison(sat, hc_type, N)
            dt_list[i][j] = comp_obj.dt_e
            dr_list[i][j] = comp_obj.dt_e * sat.R_orbit * sat.omega
            print("iteration: " + str(j) + " altitude: " + str(alt))

    plt.figure(1)
    # plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_list, '.', label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Temporal uncertaintainty in HCNM meauremental, $\delta t_e$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.savefig("bunches/points_dt_v_alt.png")

    plt.figure(2)
    # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_list, '.', label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")
    plt.savefig("bunches/points_dr_v_alt.png")
    plt.show()

    dt_data = open("bunches/points_dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dt_data.write(str(dt_list))
    dt_data.close()
    dr_data = open("bunches/points_dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dr_data.write(str(dr_list))
    dr_data.close()
    altlist = open("bunches/points_alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    altlist.write(str(altitude_list))
    altlist.close()
    return 0

#TODO finish this
def read_data(pathname):
    with open(str(pathname)) as dt_data:
        lines = dt_data.readlines()


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
            print("iteration: " + str(j) + " altitude: " + str(alt))
    dt_data = open("sample_data/dt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dt_data.write(str(dt_list))
    dt_data.close()
    dr_data = open("sample_data/dr_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    dr_data.write(str(dr_list))
    dr_data.close()
    altlist = open("sample_data/alt_int_" + str(alt_interval) + "_iter_" + str(how_many), "w")
    altlist.write(str(altitude_list))
    altlist.close()


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
