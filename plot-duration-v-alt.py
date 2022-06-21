# Author Seamus Flannery
import matplotlib.pyplot as plt
import numpy as np
import math

r_earth = 6378 # kilometers
Grav_Const = 6.6743 * (10 ** (-20))  # kilometers cubed, inverse kilograms, inverse seconds
m_earth = 5.972 * (10 ** 24)  # kilograms


def geometric_function(altitude, r_planet, grav_const, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt):
    r_orbit = r_planet+altitude  # km
    b = zero_transmit_tang_alt  # km
    a = one_transmit_tang_alt  # km
    term1 = ((altitude+r_planet)**(3/2))
    term2 = math.acos((r_planet+b) / r_orbit)-math.acos((r_planet+a) / r_orbit)
    term3 = math.sqrt(grav_const*m_planet)
    duration = term1 * term2 / term3
    return duration


def plotter_earth(min_alt, max_alt, alt_interval):
    alt_list = np.arange(min_alt, max_alt, alt_interval)
    duration_list = np.zeros_like(alt_list, float)
    for i, alt in enumerate(alt_list):
        duration_list[i] = geometric_function(alt, r_earth, Grav_Const, m_earth, 85, 160)
    print("Min Duration: " + str(np.min(duration_list)) + " at " + str(np.argmin(duration_list)*alt_interval) + " km altitude.")
    plt.figure(1)
    plt.title("Duration of Transmittivity Curve vs. Orbit Altitude")
    plt.plot(alt_list, duration_list)
    plt.ylabel("Time (s)")
    plt.xlabel("Altitude (km)")
    plt.show()



plotter_earth(180, 20000, 5)
