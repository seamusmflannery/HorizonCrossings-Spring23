import sys
import datetime
from math import e
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
sys.path.append("/homes/smflannery/HorizonCrossings-Summer22/nruhl_final_project")
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project")

# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide import CurveComparison, generate_crossings

# Global variables
N = 5378 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99] # range of transmittance in which to compare the curves
cb_str = "P2" #P2 for no atm, Earth for Full atmosphere
E_kev = 1.5 # keV
hc_type = "rising"
H=420


def main():
    #count_vs_uncertainty_alts(60, 5410,[420, 1000, 2500, 10000] , interval=200, axes_type="linear") #Earth Test
    count_vs_uncertainty(50, 2000) #P2 test
def sqrt_inverse(x, a):
    return a/np.sqrt(x)

def nicer_uncertainty():
    counts = [19, 40, 89, 244, 275, 584, 1406, 5378]
    stdErr = [0.38, 0.27, 0.20, 0.14, 0.09, 0.07, 0.05, 0.03]
    xdata = np.linspace(counts[0], counts[-1], 300)
    b, b_err = curve_fit(sqrt_inverse, counts, stdErr)
    nicer = [xdata, b]
    return nicer

def count_vs_uncertainty(N0_min, N0_max, interval=50):
    count_list = np.arange(N0_min, N0_max, interval)
    dr_list = np.zeros(len(count_list), dtype="float")
    for i in range(len(count_list)):
        sat = AnalyzeCrossing(cb=cb_str, H=H, E_kev=E_kev)
        comp_obj = CurveComparison(sat, hc_type, count_list[i], "linear")
        dr_list[i] = comp_obj.dt_e * sat.R_orbit * sat.omega
        print("N=" + str(count_list[i]))
    plt.plot(count_list, dr_list)
    plt.show()


def count_vs_uncertainty_alts(N0_min, N0_max, alt_list, interval=50, axes_type="linear"):
    count_list = np.arange(N0_min, N0_max, interval)
    dr_list = np.zeros((6, len(count_list)), dtype="float")
    local_h = H
    plt.figure(1)
    plt.title(r"$\delta r_e$ uncertainty as a function of count-rate, per altitude")
    for j in range(len(alt_list)):
        for i in range(len(count_list)):
            sat = AnalyzeCrossing(cb=cb_str, H=alt_list[j], E_kev=E_kev)
            comp_obj = CurveComparison(sat, hc_type, count_list[i], "linear")
            dr_list[j][i] = comp_obj.dt_e * sat.R_orbit * sat.omega
            print("N=" + str(count_list[i]) + " H=" + str(alt_list[j]))
        plt.plot(count_list, dr_list[j], label=str(alt_list[j]) + " km")
    # nicer = nicer_uncertainty()
    # b=nicer[1]
    # plt.plot(nicer[0], sqrt_inverse(nicer[0], b), label=r"$b/ \sqrt{N_0}$ trendline, b= %.2f" % b)
    plt.legend()
    if axes_type == "log":
        plt.yscale("log")
        plt.xscale("log")
    plt.xlabel("Count Rate")
    plt.ylabel(r"$\delta r_e$ uncertainty")
    plt.show()

if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))