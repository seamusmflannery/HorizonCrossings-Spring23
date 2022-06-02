# Author: Nathaniel Ruhl
# This script uses the toy model to determine the effects of orbital altitude on the precision of a horizon crossing

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/nathanielruhl/PycharmProjects/HorizonCrossings-local/nruhl_final_project")

# import local modules
from AnalyzeCrossing import AnalyzeCrossing
from tcc_slide import CurveComparison, generate_crossings

# Global variables
N = 500 # average number of unattenuated counts in data
bin_size = 1
comp_range = [0.01, 0.99] # range of transmittance in which to compare the curves
cb_str = "Earth"
E_kev = 2 # keV
hc_type = "setting"

def main():
    altitude_list = np.arange(300, 2100, 100)
    dt_list = np.zeros_like(altitude_list, np.float)  # list containing uncertainties corresponding to altitude_list

    for i, alt in enumerate(altitude_list):
        sat = AnalyzeCrossing(cb=cb_str, H=alt, E_kev=E_kev)
        comp_obj = CurveComparison(sat, hc_type)
        dt_list[i] = comp_obj.dt_e

    plt.figure()
    plt.title(r"$\delta t_e$ uncertainty as a function of orbital altitude")
    plt.plot(altitude_list, dt_list, label=f"{E_kev} keV {cb_str} {hc_type} crossing, N = {N}")
    plt.ylabel(r"$\delta t_e$ (sec)")
    plt.xlabel("Orbital altitude (km)")
    plt.legend()
    plt.show()
    return 0

if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
