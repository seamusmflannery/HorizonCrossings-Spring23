# Author Nathaniel Ruhl
# This script times my adaptive integration method and compares it to simpson's rule and Gaussian quadrature

import numpy as np
import time
import matplotlib.pyplot as plt

from AnalyzeCrossing import AnalyzeCrossing

def main():

    E_kev = 4.0

    ES = AnalyzeCrossing(cb="Earth", H=420)

    time_array = np.arange(0, ES.time_final+1, 1)
    transmit_adaptive = np.zeros_like(time_array)
    transmit_gauss = np.zeros_like(time_array)
    transmit_simpson = np.zeros_like(time_array)

    # Time the full horizon crossing (over 300 evaluations of the integral)
    start_time = time.time()
    tol_adapt = 1e-5
    for i, t in enumerate(time_array):
        tau_adaptive, dx_list, x_midpoints = ES.tau_adaptive_simpson(t, tol=tol_adapt)
        transmit_adaptive[i] = np.exp(-tau_adaptive)
    print(f"Adaptive quadrature with tol = {tol_adapt} takes {time.time() - start_time} seconds")

    start_time = time.time()
    N_gauss = 10
    for i, t in enumerate(time_array):
        transmit_gauss[i] = np.exp(-ES.tau_gauss(t, N_gauss))
    print(f"Gaussian quadrature with N={N_gauss} takes {time.time() - start_time}")

    start_time = time.time()
    N_simp = 1000
    for i, t in enumerate(time_array):
        transmit_simpson[i] = np.exp(-ES.tau_simpson(t, N=N_simp))
    print(f"Simpson with N = {N_simp} takes {time.time() - start_time}")

    # Calculate a gaussian quadrate with N = 500 for a "best estimate" to machine precision
    transmit_gauss500 = np.zeros_like(time_array)
    for i, t in enumerate(time_array):
        transmit_gauss500[i] = np.exp(-ES.tau_gauss(t, N=500))


    plt.plot(time_array, transmit_adaptive, label=f"adaptive (tol={tol_adapt})")
    plt.plot(time_array, transmit_gauss, label=f"gauss (N={N_gauss})")
    plt.plot(time_array, transmit_simpson, label=f"simpson (N={N_simp})")
    plt.plot(time_array, transmit_gauss500, label=f"gauss (N=500)")
    plt.legend()
    plt.show()

    # I could find when there is convergence between these methods and pick step size based on that. I could also take a practical perspective... after that, I can look at the total numer of iterations in the adaptive method

    return 0

if __name__ == '__main__':
    main()
