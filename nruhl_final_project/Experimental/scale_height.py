# Author: Nathaniel Ruhl
# This script considers the effects of changing scale height on Earth and Mars, and looks for a generalizeable relationship between the two

from AnalyzeCrossing import AnalyzeCrossing

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

# Log function for hsar vs cale height, L
def hstar_vs_scale(L, A, C):
    return A*np.log(L) + C

def analyze_scale_height(SAT):
    # List altitude at 50% transmittance
    hstar50_list = []
    h50_list = []

    L_list = [4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5]   # km, scale heights
    E_kev = 4.0

    for L in L_list:
        SAT.scale_height = L
        time_array = np.arange(0, SAT.time_final+1, 1)
        transmit_array = np.zeros_like(time_array)
        tan_alt_array = np.zeros_like(time_array)
        for i, t in enumerate(time_array):
            transmit_array[i] = np.exp(-SAT.tau_gauss(t, N=100))
            tan_alt_array[i] = SAT.tan_alt(t)
        hstar_array = tan_alt_array/SAT.scale_height
        lambdastar_array = SAT.sigma * SAT.rho0 * tan_alt_array

        # Determine hstar at 50% transmission (I could solve this analytically)
        # Transmittance array must be reduced
        a = np.where(transmit_array >= 0.01)[0][0]
        b = np.where(transmit_array <= 0.99)[0][-1]
        hstar_vs_transmit = interp1d(transmit_array[a:b], hstar_array[a:b], kind="cubic")
        hstar50_list.append(hstar_vs_transmit(0.5))
        h_vs_transmit = interp1d(transmit_array[a:b], tan_alt_array[a:b], kind="cubic")
        h50_list.append(h_vs_transmit(0.5))

        plt.figure(1)
        plt.plot(tan_alt_array, transmit_array, label=f"{SAT.cb}: scale height = {SAT.scale_height} km")
        plt.figure(2)
        plt.plot(hstar_array, transmit_array, label=f"{SAT.cb}: scale height = {SAT.scale_height} km")
        plt.figure(3)
        plt.plot(lambdastar_array, transmit_array,
                 label=f"{SAT.cb}: scale height = {SAT.scale_height} km")


    plt.figure(1)
    plt.title("Transmittance vs tangent altitude")
    plt.ylabel("Transmittance")
    plt.xlabel(r"$h$ (km)")
    plt.legend()
    plt.figure(2)
    plt.title("Transmittance vs non-dimensional altitude")
    plt.ylabel("Transmittance")
    plt.xlabel(r"$h^*$ (km)")
    plt.legend()

    plt.figure(3)
    plt.title("Transmittance vs non-dimensional mean-free path")
    plt.ylabel("Transmittance")
    plt.xlabel(r"$\lambda^*$ (km)")
    plt.legend()


    # fit a log curve to the h^* plot
    popt, pcov = curve_fit(hstar_vs_scale, L_list, hstar50_list)

    plt.figure(4)
    plt.plot(L_list, hstar50_list, label=fr'{SAT.cb} satellite, $h^*_{50}$={popt[0]:.3f}ln($L$)+{popt[1]:.3f}')
    plt.ylabel(rf"$h^*$ of 50\% transmission point")
    plt.xlabel(r"Atmospheric Scale Height (km), $L$")
    plt.legend()

    plt.figure(5)
    plt.plot(L_list, h50_list, label=f"{SAT.cb} satellite")
    plt.ylabel(rf"$h$ of 50\% transmission point")
    plt.xlabel("Atmospheric Scale Height (km)")
    plt.legend()

    return popt

if __name__ == '__main__':
    plt.rc("text", usetex=True)

    sat_alt = 420   # km, satellite altitude
    ES = AnalyzeCrossing(cb="Earth", H=sat_alt)
    MS = AnalyzeCrossing(cb="Mars", H=sat_alt)
    VS = AnalyzeCrossing(cb="Venus", H=sat_alt)
    PS1 = AnalyzeCrossing(cb="P1", H=sat_alt)

    # The shape of the log looks like they all agree very well... there might be a slight correction for atmospheric make-up or energy
    A_earth, C_earth = analyze_scale_height(ES)
    A_mars, C_mars = analyze_scale_height(MS)
    A_venus, C_venus = analyze_scale_height(VS)

    print(f"Earth: A = {A_earth}, C = {C_earth}")
    print(f"Mars: A = {A_mars}, C = {C_mars}")
    print(f"Mars: A = {A_venus}, C = {C_venus}")
    # analyze_scale_height(PS1)

    plt.show()
