# Author: Nathaniel Ruhl
# This script uses the toy model to determine the effects of orbital altitude on the precision of a horizon crossing

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/seamusflannery/Documents/HorizonCrossings-Summer22/nruhl_final_project/")

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
    do_a_bunch(400, 2100, 100, 10)
    return 0


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

    plt.figure(2)
    # plt.title(r"$\delta r_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dr_list,
             label=fr"{E_kev} keV {cb_str} {hc_type} crossing, $N_0$ = {N}")
    plt.ylabel(r"Positional uncertainty in HCNM measurement, $\delta r_e$ (km)")
    plt.xlabel("Orbital altitude (km)")

    plt.show()
    f = open("dt_bunches.csv", "w")
    f.write(dt_list)
    f.close()
    g = open("dr_bunches.csv", "w")
    g.write(dr_list)
    g.close()
    return 0


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
