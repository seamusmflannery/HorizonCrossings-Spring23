# Author: Nathaniel Ruhl
# This script uses the 2d geometry of the horizon crossing and makes preliminary plots of the curves that I want to integrate with variable step sizes -- gamma = optical depth per km vs distance on the LOS

# Note that the default in this project is to use units of km for lengths, and I will specify when I use different units

import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)

# import local libraries
from AnalyzeCrossing import AnalyzeCrossing

def make_plots():
    ISS = AnalyzeCrossing(cb="Earth", H=420)
    print(f"Period of the orbit = {ISS.T} sec")

    # define the array of times in a horizon crossing
    time_array = np.arange(0, ISS.time_final + 1, 1)
    h_array = ISS.tan_alt(time_array) # km, tangent altitudes during the horizon crossing
    transmit_array = np.zeros_like(time_array) # Transmission at each time during the horizon crossing

    for i, t in enumerate(time_array):
        transmit_array[i] = np.exp(-ISS.tau_gauss(t, N=10))

    plt.figure(1)
    plt.title("Tangent Altitude (km) vs Time (sec)")
    plt.plot(time_array, h_array)

    # Calculate the densities across the LOS at time t = 40 sec, h~40 km
    t = 50

    dtot = ISS.d_tot(t)
    x_array_km = np.arange(0, dtot+1, 1)
    z_array_km = ISS.x_to_z(x_array_km, t)
    densities = ISS.rho_vs_x(x_array_km, t)
    gamma_array = ISS.gamma_vs_x(x_array_km, t)

    plt.figure(2)
    plt.title(fr"Density ($g/cm^3$) vs x coord on LOS (km) at t={t} sec")
    plt.plot(x_array_km, densities)

    plt.figure(3)
    plt.title(
        fr"Optical Depth per km vs x coord on LOS (km) at t={t} sec, scale height = {ISS.scale_height}")
    plt.plot(x_array_km, gamma_array)

    plt.figure(4)
    plt.title(fr"Optical Depth per km vs z coord on LOS (km) at t={t} sec, scale height = {ISS.scale_height} km")
    plt.plot(z_array_km, gamma_array)

    plt.figure(6)
    plt.title(fr"Optical Depth per km vs x coord on LOS (km) at different times in the crossing")
    t_range = np.arange(50, 60, 1)
    for t_i in t_range:
        dtot = ISS.d_tot(t_i)
        h_i = ISS.tan_alt(t_i)
        x_array_i = np.arange(0, dtot+1, 1)
        gamma_array_i = ISS.gamma_vs_x(x_array_i, t_i)
        plt.plot(x_array_i, gamma_array_i, label=fr"t={t_i} sec, h = {h_i:.2f} km")

    plt.legend()

    plt.figure(7)
    # This is the plot that characterizes the necessary step size!
    plt.title(fr"Optical Depth per km vs density on LOS (km) at different times in the crossing")
    for t_i in t_range:
        dtot = ISS.d_tot(t_i)
        h_i = ISS.tan_alt(t_i)
        x_array_i = np.arange(0, dtot+1, 1)
        z_array_i = ISS.x_to_z(x_array_i, t_i)
        gamma_array_i = ISS.gamma_vs_x(x_array_i, t_i)
        density_i = ISS.rho_vs_z(z_array_i, t_i)
        plt.plot(density_i, gamma_array_i,
                 label=fr"t={t_i} sec, h = {h_i:.2f} km")
    plt.legend()

    plt.figure(8)
    plt.title("Transmission vs Time")
    plt.plot(time_array, transmit_array)

    plt.show()

    return 0

# Calculate area under gamma_array vs x with variable step sizes

def calc_area():
    ISS = AnalyzeCrossing(cb="Earth", H=420)
    E_kev = 4.0
    t = 50
    dtot = ISS.d_tot(t)
    N = 4000
    a = 0.0
    b = dtot/2
    h = (b-a)/N
    x_array = np.arange(a, b+h, h)
    gamma_list = ISS.gamma_vs_x(x_array, t)   # Optical depth per km

    plt.plot(x_array, gamma_list)
    plt.show()

    s_odd = 0
    s_even = 0
    for i in range(len(x_array)):
        if i % 2 == 0:
            s_even += gamma_list[i]
        else:
            s_odd += gamma_list[i]

    simp_sum = 2*(h/3)*(gamma_list[0] + gamma_list[-1] +
                        4*s_odd + 2*s_even)   # value of integral

    print(fr'simp_sum = {simp_sum}')

    return 0

if __name__=="__main__":
    make_plots()
    # calc_area()
