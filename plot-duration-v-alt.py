# Author Seamus Flannery
import matplotlib.pyplot as plt
import numpy as np
import math
Grav_Const = 6.6743 * (10 ** (-20))  # kilometers cubed, inverse kilograms, inverse seconds
Earth = [6378, 5.972*(10**24), 85, 160, "Earth"]
# approx 6378km radius, 5.972*10**24kg mass, 85km lower lim x-ray penetration, 160km for upper lim for Earth
Mars = [3390, 6.39*(10**23), 10, 20, "Mars"]
# approx 3390km radius, 6.39*(10**23)kg mass, 10km lower lim x-ray penetration, 20km for upper lim for Mars
Venus = [6051, 4.867*(10**24), 150, 300, "Venus"]
# approx 6051km radius, 4.867*(10**24)kg mass, 150km lower lim x-ray penetration, 300km for upper lim for Venus
Moon = [1737, 7.35*(10**22), 0.0001, 0.00011, "Moon"]
def main():
    duration_v_alt(Earth, 300, 20000)
    duration_v_alt(Mars, 150, 20000)
    duration_v_alt(Venus, 300, 20000)
    duration_v_alt(Moon, 50, 20000)

def geometric_function(altitude, r_planet, grav_const, m_planet, zero_transmit_tang_alt, one_transmit_tang_alt):
    r_orbit = r_planet+altitude  # km
    b = zero_transmit_tang_alt  # km
    a = one_transmit_tang_alt  # km
    term1 = ((altitude+r_planet)**(3/2))
    term2 = math.acos((r_planet+b) / r_orbit)-math.acos((r_planet+a) / r_orbit)
    term3 = math.sqrt(grav_const*m_planet)
    duration = term1 * term2 / term3
    return duration


def duration_v_alt(planet, min_alt, max_alt, alt_interval=5, save = "false"):
    alt_list = np.arange(min_alt, max_alt, alt_interval)
    duration_list = np.zeros_like(alt_list, float)
    for i, alt in enumerate(alt_list):
        duration_list[i] = geometric_function(alt, planet[0], Grav_Const, planet[1], planet[2], planet[3])
    print("Min Duration: " + str(np.min(duration_list)) + " seconds at " + str(np.argmin(duration_list)*
                    alt_interval) + " +/- " + str(alt_interval/2) + "km altitude.")
    plt.figure(1)
    plt.title("Duration of " + planet[4] + " Transmittivity Curve vs. Orbit Altitude")
    plt.plot(alt_list, duration_list)
    plt.ylabel("Time (s)")
    plt.xlabel("Altitude (km)")
    if save == "true":
        plt.savefig("duration_v_alt_plots/" + str(planet[4]) + "_large.png")
    plt.show()


if __name__ == '__main__':
    import time
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))

