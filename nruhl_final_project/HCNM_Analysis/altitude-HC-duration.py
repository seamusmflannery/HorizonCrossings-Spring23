# Author: Nathaniel Ruhl
# This script plots the relationship between orbital altitude and the duration of a HC curve

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
sys.path.append("/Users/nathanielruhl/Documents/HorizonCrossings-Summer22/nruhl_final_project/")  # add working directory, str(Path(__file__).parents[1])

from AnalyzeCrossing import AnalyzeCrossing

altitude_list = np.arange(400, 1300, 100)
dt_list = []   # list containing duration of HC's
m50_list = []   # list containing the slope at the 50% point

# Calculate and plot transmittance curves at each orbital altitude
for altitude in altitude_list:
    sat = AnalyzeCrossing(cb="Earth", H=altitude, E_kev=1.5)
    time_array = np.arange(0, sat.time_final, 1)
    transmit_array = sat.calc_transmit_curve(time_array)

    # Calculate dt = time between 99% transmit and 1% transmittance
    # First reduce the range to ensure we have a 1:1 function
    frange = np.where((transmit_array > 0.005) & (transmit_array < 0.995))[0]

    time_vs_transmit = interp1d(x=transmit_array[frange], y=time_array[frange])

    dt = time_vs_transmit(0.99) - time_vs_transmit(0.1)
    dt_list.append(dt)

    # Calculate the slope at the 50% points
    m50 = (51 - 49)/(time_vs_transmit(0.51)-time_vs_transmit(0.49))
    m50_list.append(m50)

    # Plot the tranmsittance curves
    plt.figure(1)
    plt.plot(time_array, transmit_array,
        label=fr"H = {altitude} km, $\Delta t_{{99}}$ = {dt:.1f} sec")

plt.figure(1)
# plt.title("1.5 keV Earth Horizon Crossing Curves at Different Orbital Altitudes")
plt.ylabel("Transmittance")
plt.xlabel("Time")
plt.legend()

plt.figure(2)
plt.plot(altitude_list, dt_list)
# plt.title("Duration of 1.5 keV Earth Horizon Crossing Curves at Different Orbital Altitudes")
plt.xlabel("Orbital Altitude (km)")
plt.ylabel("Time from 1% to 99% Transmittance (sec)")
plt.legend()

plt.figure(3)
# plt.title("Slope of 1.5 keV Earth Horizon Crossing Curves at Different Orbital Altitudes")
plt.plot(altitude_list, m50_list)
plt.xlabel("Orbital Altitude (km)")
plt.ylabel(r"Slope of transmittance vs time curve at $t_{50}$")

plt.show()
